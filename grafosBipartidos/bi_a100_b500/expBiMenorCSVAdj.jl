using Metaheuristics
using LinearAlgebra
using Statistics
using Printf
using DataFrames
using CSV
using ProgressMeter
using Random

import Metaheuristics: initialize!, update_state!, final_stage!, gen_initial_state

include("../../colorGul.jl")

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

# método para extração de parâmetros do nome do arquivo de entrada
function extract_bipartite_params(file_name)
    match_result = match(r"bi_a(\d+)_b(\d+)_p(\d+)%_v(\d+)\.col", file_name)
    if match_result !== nothing
        return Dict(
            :a_param => parse(Int, match_result.captures[1]),
            :b_param => parse(Int, match_result.captures[2]),
            :p_param => parse(Int, match_result.captures[3]),
            :v_param => parse(Int, match_result.captures[4])
        )
    end
    return nothing
end

function main()
    N_REPETITIONS = 5
    K_LIMIT = 50
    N_POP = 200
    
    TARGET_A = 100
    TARGET_B = 1000

    # filtragem de arquivos de entrada de acordo com os parâmetros desejados
    all_files = filter(f -> endswith(f, ".col"), readdir())
    filtered_files = String[]

    for file in all_files
        p = extract_bipartite_params(file)
        if p !== nothing && p[:a_param] == TARGET_A && p[:b_param] == TARGET_B
            push!(filtered_files, file)
        end
    end

    if isempty(filtered_files)
        println("AVISO: Nenhum arquivo corresponde aos critérios (a=$TARGET_A, b=$TARGET_B).")
        return
    end

    sort!(filtered_files)
    results_main = []
    results_prof = []

    @showprogress 1 "Processando: " for (idx, file) in enumerate(filtered_files)
        println("\n[$idx/$(length(filtered_files))] Processing: $file")

        params = extract_bipartite_params(file)
        t_tot, chi = Float64[], Int[]
        t_fit, t_cross, t_mut, t_sel, t_init, t_upd = 
            Float64[], Float64[], Float64[], Float64[], Float64[], Float64[]

        local nv, na
        for i in 1:N_REPETITIONS
            tt, ch, nv, na, tf, tx, tm, ts, ti, tu = run_ga_experiment(file, K_LIMIT, N_POP)

            push!(t_tot, tt); push!(chi, ch)
            push!(t_fit, tf); push!(t_cross, tx)
            push!(t_mut, tm); push!(t_sel, ts)
            push!(t_init, ti); push!(t_upd, tu)

            @printf("   Run %d Done (Chi: %d)\n", i, ch)
        end

        se(v) = std(v) / sqrt(length(v))
        meanr(v) = mean(v)

        # dicionário na ordem desejada para a saída do .csv final
        push!(results_main, Dict(
            :a => params[:a_param],
            :b => params[:b_param],
            :N => nv,
            :p => params[:p_param],
            :M => na,
            :v => params[:v_param],
            :mean_time => meanr(t_tot),
            :se_time => se(t_tot),
            :mean_chi => meanr(chi),
            :se_chi => se(chi),
            :instancia => file
        ))

        push!(results_prof, Dict(
            :instancia => file,
            :a => params[:a_param],
            :b => params[:b_param],
            :mean_fit_t => meanr(t_fit),
            :mean_sel_t => meanr(t_sel),
            :mean_cross_t => meanr(t_cross),
            :mean_mut_t => meanr(t_mut),
            :mean_init_t => meanr(t_init),
            :mean_upd_t => meanr(t_upd)
        ))
    end

    df_main = DataFrame(results_main)
    col_order = [:a, :b, :N, :p, :M, :v, :mean_time, :se_time, :mean_chi, :se_chi, :instancia]
    select!(df_main, col_order)

    CSV.write("results_GA_Final_a100_b1000.csv", df_main)
    CSV.write("results_GA_Profiling_a100_b1000.csv", DataFrame(results_prof))

    println("\n--- Experimentos concluídos ---")
end

main()