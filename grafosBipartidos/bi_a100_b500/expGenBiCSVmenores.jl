using Metaheuristics
using LinearAlgebra
using Statistics
using Printf
using DataFrames 
using CSV         
using ProgressMeter

include("../../colorGul.jl") 

# --- Timers para Profiling ---
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

global matriz_adj = Matrix{Int}(undef, 0, 0)
global num_vertices = 0

import Metaheuristics: initialize!, update_state!, final_stage!, gen_initial_state

mutable struct CustomGAParams <: Metaheuristics.AbstractParameters
    N::Int
    p_crossover::Float64
    p_mutation::Float64
    stag_limit::Int
    last_best::Float64
    stag_iters::Int
end 

# Construtor
CustomGAParams(; N = 1000, p_crossover = 0.5, p_mutation = 0.5,
                 stag_limit = 50, last_best = Inf, stag_iters = 0) =
    CustomGAParams(N, p_crossover, p_mutation, stag_limit, last_best, stag_iters)

# --- Métodos com Cronometragem (Profiling) ---

function tournament_select_timed(pop)
    t = @elapsed begin
        candidates = rand(pop, 2)
        res = candidates[argmin([c.f for c in candidates])].x
    end
    TIME_SELECTION[] += t
    return res
end

function crossover_media_simples_timed(parents::Matrix{Float64})
    t = @elapsed begin
        p1 = parents[1, :]
        p2 = parents[2, :]
        child = (p1 .+ p2) ./ 2.0
        res = [child'] 
    end
    TIME_CROSSOVER[] += t
    return res
end

function graph_swap_mutation_timed!(x::Vector{Float64}, m_adj)
    t = @elapsed begin
        n = length(x)
        if n >= 2
            v1 = rand(1:n)
            vizinhos = [i for i in 1:n if m_adj[v1, i] == 1]
            if !isempty(vizinhos)
                v2 = rand(vizinhos)
                x[v1], x[v2] = x[v2], x[v1]
            end
        end
    end
    TIME_MUTATION[] += t
    return x
end 

function Metaheuristics.initialize!(status, parameters::CustomGAParams, problem, information, options)
    t_init = @elapsed begin
        initial_state = Metaheuristics.gen_initial_state(problem, parameters, information, options, status)
        parameters.last_best = minimum([s.f for s in initial_state.population])
        parameters.stag_iters = 0
    end
    TIME_INITIALIZE[] += t_init
    return initial_state
end 

function Metaheuristics.update_state!(state::Metaheuristics.State, parameters::CustomGAParams, problem, information, options)
    global matriz_adj
    t_upd = @elapsed begin
        p1 = tournament_select_timed(state.population)
        p2 = tournament_select_timed(state.population)
        parents = vcat(p1', p2')

        offsprings = crossover_media_simples_timed(parents)

        for i in axes(offsprings, 1)
            vec = offsprings[i, :]
            graph_swap_mutation_timed!(vec, matriz_adj)
            fit = float(problem.f(vec))
            push!(state.population, Metaheuristics.xf_solution(vec, fit))
        end 

        sort!(state.population, by = s -> s.f)
        state.population = state.population[1:min(parameters.N, length(state.population))]

        melhor_atual = state.population[1].f
        if melhor_atual < parameters.last_best
            parameters.last_best = melhor_atual
            parameters.stag_iters = 0
        else
            parameters.stag_iters += 1
        end 
        continuar = parameters.stag_iters < parameters.stag_limit
    end
    TIME_UPDATE_STATE[] += t_upd
    return continuar
end 

function Metaheuristics.final_stage!(state, parameters::CustomGAParams, problem, information, options)
    return state
end

# --- Lógica do experimento ---

function run_ga_experiment(file_name::String, k_limit::Int, N_pop::Int)
    reset_timers!()
    global matriz_adj
    num_v, num_a = leInfo!(file_name)
    matriz_adj = zeros(Int, num_v, num_v)
    leArestas!(file_name, matriz_adj) 

    fitness_func(ind) = begin
        t_fit = @elapsed begin
            res = float(maximum(NOVOcoloracaoHarmonicaGuloso!(matriz_adj, sortperm(ind, rev=true))))
        end
        TIME_FITNESS[] += t_fit
        return res
    end

    bounds = [zeros(num_v) ones(num_v)]'
    algo = Metaheuristics.Algorithm(CustomGAParams(N = N_pop, stag_limit = k_limit))
    
    time_taken = @elapsed begin
        result = Metaheuristics.optimize(fitness_func, bounds, algo)
    end
    
    return time_taken, Int(Metaheuristics.minimum(result)), num_v, num_a,
           TIME_FITNESS[], TIME_CROSSOVER[], TIME_MUTATION[], 
           TIME_SELECTION[], TIME_INITIALIZE[], TIME_UPDATE_STATE[]
end

function main()
    N_REPETITIONS = 5
    K_LIMIT = 50
    N_POP = 1000 

    # 1. Coleta todos os arquivos bipartidos na pasta
    all_files = filter(f -> startswith(f, "bi_") && endswith(f, ".col"), readdir())
    filtered_file_names = String[]

    limite_a = 500
    limite_b = 1000

    for file_name in all_files
        # Extrai os parâmetros do nome do arquivo via Regex
        m = match(r"bi_a(\d+)_b(\d+)_p(\d+)%_v(\d+)\.col", file_name)
        if m === nothing continue end

        # Converte as capturas para inteiros
        a_param = parse(Int, m.captures[1])
        b_param = parse(Int, m.captures[2])

        # Lógica de decisão para inclusão no experimento
        if a_param <= limite_a && b_param <= limite_b
            push!(filtered_file_names, file_name)
        end
    end

    sort!(filtered_file_names)

    results_main = []
    results_prof = [] # Lista para o segundo CSV de profiling

    # 2. Progress bar e logs de processamento
    @showprogress 1 "Processando: " for (idx, file) in enumerate(filtered_file_names)
        println("\n[$idx/$(length(filtered_file_names))] Processing: $file")
        
        t_tot, chi = Float64[], Int[]
        t_fit, t_cross, t_mut, t_sel, t_init, t_upd = Float64[], Float64[], Float64[], Float64[], Float64[], Float64[]
        
        # 3. Extração de parâmetros corrigida para o padrão bi_a100_b100_p10%_v1.col
        m = match(r"bi_a(\d+)_b(\d+)_p(\d+)%_v(\d+)\.col", file)
        if m === nothing continue end
        a, b, p, v = parse(Int, m[1]), parse(Int, m[2]), parse(Int, m[3]), parse(Int, m[4])

        local nv, na
        for i in 1:N_REPETITIONS
            tt, ch, nv, na, tf, tx, tm, ts, ti, tu = run_ga_experiment(file, K_LIMIT, N_POP)
            
            push!(t_tot, tt); push!(chi, ch)
            push!(t_fit, tf); push!(t_cross, tx); push!(t_mut, tm)
            push!(t_sel, ts); push!(t_init, ti); push!(t_upd, tu)
            
            @printf("   Run %d Done (Chi: %d)\n", i, ch)
        end 

        se(vec) = round(std(vec)/sqrt(N_REPETITIONS), digits=5)
        rmean(vec) = round(mean(vec), digits=5)

        # 4. Resultados principais
        push!(results_main, Dict(
            :instancia => file, :a => a, :b => b, :p => p, :v => v, :N => nv, :M => na,
            :mean_chi => round(mean(chi), digits=2), 
            :se_chi => round(std(chi)/sqrt(N_REPETITIONS), digits=2),
            :mean_time => rmean(t_tot), 
            :se_time => se(t_tot)
        ))

        # 5. Resultados de Profiling (idêntico ao 2º arquivo)
        push!(results_prof, Dict(
            :instancia => file, :a => a, :b => b, :p => p, :v => v, :N => nv, :M => na,
            :mean_init_t => rmean(t_init), :se_init_t => se(t_init),
            :mean_upd_t  => rmean(t_upd),  :se_upd_t  => se(t_upd),
            :mean_fit_t  => rmean(t_fit),  :se_fit_t  => se(t_fit),
            :mean_sel_t  => rmean(t_sel),  :se_sel_t  => se(t_sel),
            :mean_cross_t => rmean(t_cross), :se_cross_t => se(t_cross),
            :mean_mut_t  => rmean(t_mut),  :se_mut_t  => se(t_mut)
        ))
    end

    # 6. Salvamento de dois arquivos CSV distintos
    df_m = DataFrame(results_main)
    select!(df_m, [:instancia, :a, :b, :p, :v, :N, :M, :mean_chi, :se_chi, :mean_time, :se_time])
    sort!(df_m, [:a, :b, :p, :v])
    CSV.write("results_GA_Bipartite_Final.csv", df_m)

    df_p = DataFrame(results_prof)
    select!(df_p, [:instancia, :a, :b, :p, :v, :N, :M, 
                   :mean_init_t, :se_init_t, :mean_upd_t, :se_upd_t, 
                   :mean_fit_t, :se_fit_t, :mean_sel_t, :se_sel_t,
                   :mean_cross_t, :se_cross_t, :mean_mut_t, :se_mut_t])
    sort!(df_p, [:a, :b, :p, :v])
    CSV.write("results_GA_Bipartite_Profiling.csv", df_p)
    
    println("\n--- Experimentos Concluídos! CSVs gerados. ---")
end

main()