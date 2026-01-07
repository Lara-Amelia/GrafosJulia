using Metaheuristics
using LinearAlgebra
using Statistics
using Printf
using DataFrames
using CSV
using ProgressMeter
using Random

import Metaheuristics: initialize!, update_state!, final_stage!, gen_initial_state

include("../colorGul.jl")

const TIME_FITNESS    = Ref(0.0)
const TIME_CROSSOVER  = Ref(0.0)
const TIME_MUTATION   = Ref(0.0)
const TIME_SELECTION = Ref(0.0)
const TIME_INITIALIZE = Ref(0.0)
const TIME_UPDATE_STATE = Ref(0.0)

function reset_timers!()
    for t in (TIME_FITNESS, TIME_CROSSOVER, TIME_MUTATION,
              TIME_SELECTION, TIME_INITIALIZE, TIME_UPDATE_STATE)
        t[] = 0.0
    end
end

global ADJ = Vector{Vector{Int}}()
global V = 0

const FITNESS_CACHE = Dict{UInt64, Float64}()

function hash_perm(p::Vector{Int})
    h = UInt64(1469598103934665603)
    for x in p
        h ⊻= UInt64(x)
        h *= 1099511628211
    end
    return h
end

function fitness_harmonious_coloring(ind::Vector{Float64})
    t = @elapsed begin
        perm = sortperm(ind, rev = true)
        h = hash_perm(perm)

        if haskey(FITNESS_CACHE, h)
            res = FITNESS_CACHE[h]
        else
            cores = coloracaoHarmonicaAdjVetAux!(ADJ, perm)
            res = maximum(cores)
            FITNESS_CACHE[h] = res
        end
    end
    TIME_FITNESS[] += t
    return Float64(res)
end

function tournament_select(pop)
    t = @elapsed begin
        a, b = rand(pop, 2)
        res = (a.f <= b.f) ? a.x : b.x
    end
    TIME_SELECTION[] += t
    return res
end

function crossover_media_simples(p1, p2)
    t = @elapsed child = (p1 .+ p2) ./ 2
    TIME_CROSSOVER[] += t
    return child
end

function graph_swap_mutation!(x)
    t = @elapsed begin
        while true
            v1 = rand(1:length(x))
            if !isempty(ADJ[v1])
                v2 = rand(ADJ[v1])
                x[v1], x[v2] = x[v2], x[v1]
                return
            end
        end
    end
    TIME_MUTATION[] += t
end

mutable struct CustomGAParams <: Metaheuristics.AbstractParameters
    N::Int
    p_mutation::Float64
    stag_limit::Int
    last_best::Float64
    stag_iters::Int
end

CustomGAParams(; N, p_mutation=0.5, stag_limit) =
    CustomGAParams(N, p_mutation, stag_limit, Inf, 0)

function Metaheuristics.initialize!(
    status,
    params::CustomGAParams,
    problem,
    information,
    options
)
    t = @elapsed begin
        params.last_best = Inf
        params.stag_iters = 0
        state = gen_initial_state(problem, params, information, options, status)
    end
    TIME_INITIALIZE[] += t
    return state
end

function Metaheuristics.update_state!(
    state,
    params::CustomGAParams,
    problem,
    information,
    options
)
    t = @elapsed begin
        pop = state.population
        offspring = Metaheuristics.xf_solution[]

        for _ in 1:div(params.N, 2)
            p1 = tournament_select(pop)
            p2 = tournament_select(pop)

            child = crossover_media_simples(p1, p2)

            if rand() < params.p_mutation
                graph_swap_mutation!(child)
            end

            f = problem.f(child)
            push!(offspring, Metaheuristics.xf_solution(child, f))
        end

        # elitist replacement
        combined = vcat(pop, offspring)
        sort!(combined, by = s -> s.f)
        state.population = combined[1:params.N]

        gen_best = state.population[1]

        # GLOBAL BEST 
        if state.best_sol === nothing || gen_best.f < state.best_sol.f
            state.best_sol = deepcopy(gen_best)
        end

        # stagnation
        if gen_best.f < params.last_best
            params.last_best = gen_best.f
            params.stag_iters = 0
        else
            params.stag_iters += 1
        end
    end
    TIME_UPDATE_STATE[] += t

    return params.stag_iters < params.stag_limit
end

Metaheuristics.final_stage!(state, params, problem, information, options) = state

function run_ga_experiment(file_name::String, k_limit::Int, N_pop::Int)
    reset_timers!()
    empty!(FITNESS_CACHE)

    global ADJ, V

    num_v, num_a = leInfo!(file_name)
    ADJ = [Int[] for _ in 1:num_v]
    leArestasLista!(file_name, ADJ)
    V = num_v

    bounds = [zeros(V) ones(V)]'

    params = CustomGAParams(N = N_pop, stag_limit = k_limit)
    algo = Metaheuristics.Algorithm(params)

    t = @elapsed result =
        Metaheuristics.optimize(fitness_harmonious_coloring, bounds, algo)

    return t,
           Int(Metaheuristics.minimum(result)),
           num_v,
           num_a,
           TIME_FITNESS[],
           TIME_CROSSOVER[],
           TIME_MUTATION[],
           TIME_SELECTION[],
           TIME_INITIALIZE[],
           TIME_UPDATE_STATE[]
end

function main()
    all_files = filter(f -> startswith(f, "galaxy_") && endswith(f, ".col"), readdir())
    sort!(all_files)

    results_main, results_prof = [], []
    n_rep = 5
    K_STAG = 50
    N_POP = 200

    @showprogress 1 "Processando: " for (idx, file) in enumerate(all_files)
        println("\n[$idx/$(length(all_files))] Processing: $file")
        
        t_tot, chi = Float64[], Int[]
        t_fit, t_cross, t_mut, t_sel, t_init, t_upd = [Float64[] for _ in 1:6]
        
        # Regex para extrair parâmetros do padrão galaxy_aX_bX_cX_vX
        m_a = match(r"a(\d+)", file); a = m_a !== nothing ? parse(Int, m_a.captures[1]) : 0
        m_b = match(r"b(\d+)", file); b = m_b !== nothing ? parse(Int, m_b.captures[1]) : 0
        m_c = match(r"c(\d+)", file); c = m_c !== nothing ? parse(Int, m_c.captures[1]) : 0
        m_v = match(r"v(\d+)", file); v = m_v !== nothing ? parse(Int, m_v.captures[1]) : 0

        local nv, na
        for i in 1:n_rep
            tt, ch, nv, na, tf, tx, tm, ts, ti, tu = run_ga_experiment(file, K_STAG, N_POP)
            push!(t_tot, tt); push!(chi, ch)
            push!(t_fit, tf); push!(t_cross, tx); push!(t_mut, tm)
            push!(t_sel, ts); push!(t_init, ti); push!(t_upd, tu)
            @printf("   Run %d Done (Chi: %d)\n", i, ch)
        end

        se(v) = round(std(v)/sqrt(n_rep), digits=4)
        mn(v) = round(mean(v), digits=5)

        push!(results_main, Dict(
            :instancia => file, :a => a, :b => b, :c => c, :v => v, :N => nv, :M => na,
            :mean_chi => mean(chi), :se_chi => se(chi),
            :mean_time => mn(t_tot), :se_time => se(t_tot)
        ))

        push!(results_prof, Dict(
            :instancia => file, :a => a, :b => b, :c => c, :v => v, :N => nv, :M => na,
            :mean_init_t => mn(t_init), :se_init_t => se(t_init),
            :mean_upd_t  => mn(t_upd),  :se_upd_t  => se(t_upd),
            :mean_fit_t  => mn(t_fit),  :se_fit_t  => se(t_fit),
            :mean_sel_t  => mn(t_sel),  :se_sel_t  => se(t_sel),
            :mean_cross_t => mn(t_cross), :se_cross_t => se(t_cross),
            :mean_mut_t  => mn(t_mut),  :se_mut_t  => se(t_mut)
        ))
    end

    CSV.write("results_GA_Final.csv", DataFrame(results_main))
    CSV.write("results_GA_Profiling_Summary.csv", DataFrame(results_prof))

    println("\nCSVs Finalizados.")
end

main()