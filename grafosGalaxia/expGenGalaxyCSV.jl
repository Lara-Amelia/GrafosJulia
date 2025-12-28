using Metaheuristics
using LinearAlgebra
using Statistics
using Printf
using DataFrames
using CSV
using ProgressMeter

# cronômetros para medição de tempo em funções do GA
const TIME_FITNESS = Ref(0.0)
const TIME_CROSSOVER = Ref(0.0)
const TIME_MUTATION = Ref(0.0)
const TIME_SELECTION = Ref(0.0)
const TIME_INITIALIZE = Ref(0.0)
const TIME_UPDATE_STATE = Ref(0.0)

function reset_timers!()
    TIME_FITNESS[] = 0.0
    TIME_CROSSOVER[] = 0.0
    TIME_MUTATION[] = 0.0
    TIME_SELECTION[] = 0.0
    TIME_INITIALIZE[] = 0.0
    TIME_UPDATE_STATE[] = 0.0
end

import Metaheuristics: initialize!, update_state!, final_stage!, gen_initial_state
include("../colorGul.jl")

mutable struct CustomGAParams <: Metaheuristics.AbstractParameters
    N::Int
    p_crossover::Float64
    p_mutation::Float64
    stag_limit::Int        
    last_best::Int          
    stag_iters::Int 
end

CustomGAParams(; N = 1000, p_crossover = 0.5, p_mutation = 0.5, stag_limit = 50, 
                  last_best = -1, stag_iters = 0) =
    CustomGAParams(N, p_crossover, p_mutation, stag_limit, last_best, stag_iters)

# métodos personalizados com medição de tempo
function tournament_select_timed(pop)
    t = @elapsed begin
        candidates = rand(pop, 2)
        res = candidates[argmin([c.f for c in candidates])].x
    end
    TIME_SELECTION[] += t
    return res
end

function crossover_media_simples(parents::Matrix{Float64})
    t = @elapsed begin
        p1 = parents[1, :]
        p2 = parents[2, :]
        child = (p1 .+ p2) ./ 2.0
        res = [child']
    end
    TIME_CROSSOVER[] += t
    return res
end

function graph_swap_mutation!(x::Vector{Float64}, matriz_adj_local)
    t = @elapsed begin
        n = length(x)
        if n >= 2
            v1 = rand(1:n)
            vizinhos = [i for i in 1:n if matriz_adj_local[v1, i] == 1]
            if !isempty(vizinhos)
                v2 = rand(vizinhos)
                x[v1], x[v2] = x[v2], x[v1]
            end
        end
    end
    TIME_MUTATION[] += t
    return x
end

function Metaheuristics.update_state!(state, parameters::CustomGAParams, problem, information, options)
    t_upd = @elapsed begin
        # Seleção
        p1 = tournament_select_timed(state.population)
        p2 = tournament_select_timed(state.population)
        parents = vcat(p1', p2')  
        
        # Crossover
        offsprings = crossover_media_simples(parents)

        # Mutação e Fitness
        for i in axes(offsprings, 1)  
            vec = offsprings[i, :]
            fit = Float64(problem.f(vec)) 
            push!(state.population, Metaheuristics.xf_solution(vec, fit))
        end

        # Ordenação e Elite
        sort!(state.population, by = s -> s.f)
        state.population = state.population[1:min(parameters.N, length(state.population))]

        # Estagnação
        melhor_atual = state.population[1].f
        if melhor_atual < parameters.last_best || parameters.last_best == -1
            parameters.last_best = Int(melhor_atual)
            parameters.stag_iters = 0
        else
            parameters.stag_iters += 1
        end
        cont_exec = parameters.stag_iters < parameters.stag_limit
    end
    TIME_UPDATE_STATE[] += t_upd
    return cont_exec
end

function Metaheuristics.initialize!(status, parameters::CustomGAParams, problem, information, options)
    t_init = @elapsed begin
        state = Metaheuristics.gen_initial_state(problem, parameters, information, options, status)
        parameters.last_best = Int(minimum([s.f for s in state.population]))
        parameters.stag_iters = 0
    end
    TIME_INITIALIZE[] += t_init
    return state
end

function Metaheuristics.final_stage!(state, p, prob, info, opts) return state end

# lógica dos experimentos 
function run_ga_experiment(file_name, k_limit, n_pop)
    reset_timers!()
    nv, na = leInfo!(file_name)
    matriz_adj = zeros(Int, nv, nv)
    leArestas!(file_name, matriz_adj)

    fitness_func(ind) = begin
        t_fit = @elapsed begin
            res = maximum(NOVOcoloracaoHarmonicaGuloso!(matriz_adj, sortperm(ind, rev=true)))
        end
        TIME_FITNESS[] += t_fit
        return res
    end

    bounds = [zeros(nv) ones(nv)]'
    params = CustomGAParams(N = n_pop, stag_limit = k_limit)
    algo = Metaheuristics.Algorithm(params, options = Metaheuristics.Options(iterations=3000, verbose=false))
    
    total_time = @elapsed begin
        result = Metaheuristics.optimize(fitness_func, bounds, algo)
    end
    
    return total_time, Int(Metaheuristics.minimum(result)), 
           TIME_FITNESS[], TIME_CROSSOVER[], TIME_MUTATION[], 
           TIME_SELECTION[], TIME_INITIALIZE[], TIME_UPDATE_STATE[]
end

function main()
    all_files = filter(f -> startswith(f, "galaxy_") && endswith(f, ".col"), readdir())
    sort!(all_files)

    results_main, results_prof = [], []
    n_rep = 5

    @showprogress 1 "Processando: " for (idx, file) in enumerate(all_files)
        println("\n[$idx/$(length(all_files))] Processing: $file")
        
        t_tot, chi = Float64[], Int[]
        t_fit, t_cross, t_mut, t_sel, t_init, t_upd = Float64[], Float64[], Float64[], Float64[], Float64[], Float64[]
        
        nv, na = leInfo!(file)
        m_a = match(r"a(\d+)", file); a = m_a !== nothing ? parse(Int, m_a.captures[1]) : 0
        m_b = match(r"b(\d+)", file); b = m_b !== nothing ? parse(Int, m_b.captures[1]) : 0
        m_c = match(r"c(\d+)", file); c = m_c !== nothing ? parse(Int, m_c.captures[1]) : 0
        m_v = match(r"v(\d+)", file); v = m_v !== nothing ? parse(Int, m_v.captures[1]) : 0

        for i in 1:n_rep
            tt, ch, tf, tx, tm, ts, ti, tu = run_ga_experiment(file, 50, 1000)
            push!(t_tot, tt); push!(chi, ch)
            push!(t_fit, tf); push!(t_cross, tx); push!(t_mut, tm)
            push!(t_sel, ts); push!(t_init, ti); push!(t_upd, tu)
            @printf("  Run %d Done\n", i)
        end

        se(v) = round(std(v)/sqrt(n_rep), digits=4)

        push!(results_main, Dict(
            :instancia => file, :a => a, :b => b, :c => c, :v => v, :N => nv, :M => na,
            :mean_chi => round(mean(chi), digits=2), :se_chi => round(std(chi)/sqrt(n_rep), digits=2),
            :mean_time => round(mean(t_tot), digits=5), :se_time => se(t_tot)
        ))

        push!(results_prof, Dict(
            :instancia => file, :a => a, :b => b, :c => c, :v => v, :N => nv, :M => na,
            :mean_init_t => round(mean(t_init), digits=5), :se_init_t => se(t_init),
            :mean_upd_t  => round(mean(t_upd), digits=5),  :se_upd_t  => se(t_upd),
            :mean_fit_t  => round(mean(t_fit), digits=5),  :se_fit_t  => se(t_fit),
            :mean_sel_t  => round(mean(t_sel), digits=5),  :se_sel_t  => se(t_sel),
            :mean_cross_t => round(mean(t_cross), digits=5), :se_cross_t => se(t_cross),
            :mean_mut_t  => round(mean(t_mut), digits=5),  :se_mut_t  => se(t_mut)
        ))
    end

    # Salva CSV principal (resultados do genético para cada instância)
    df_m = DataFrame(results_main)
    select!(df_m, [:instancia, :a, :b, :c, :v, :N, :M, :mean_chi, :se_chi, :mean_time, :se_time])
    sort!(df_m, [:a, :b, :c, :v])
    CSV.write("results_GA_Final.csv", df_m)

    # Salva CSV de profiling (tempos para cada função do ga)
    df_p = DataFrame(results_prof)
    select!(df_p, [:instancia, :a, :b, :c, :v, :N, :M, 
                   :mean_init_t, :se_init_t, :mean_upd_t, :se_upd_t, 
                   :mean_fit_t, :se_fit_t, :mean_sel_t, :se_sel_t,
                   :mean_cross_t, :se_cross_t, :mean_mut_t, :se_mut_t])
    sort!(df_p, [:a, :b, :c, :v])
    CSV.write("results_GA_Profiling_Summary.csv", df_p)

    println("\nCSVs Finalizados.")
end

main()