using Metaheuristics
using LinearAlgebra
using Statistics
using Printf
using DataFrames 
using CSV        

include("../../colorGul.jl") 

global matriz_adj = Matrix{Int}(undef, 0, 0)
global num_vertices = 0

mutable struct CustomGAParams <: Metaheuristics.AbstractParameters
    N::Int
    p_crossover::Float64
    p_mutation::Float64
    stag_limit::Int
    last_best::Float64
    stag_iters::Int
end 

# Construtor com calibração exata 
CustomGAParams(; N = 1000, p_crossover = 0.5, p_mutation = 0.5,
                 stag_limit = 50, last_best = Inf, stag_iters = 0) =
    CustomGAParams(N, p_crossover, p_mutation, stag_limit, last_best, stag_iters)

# Crossover de filho único 
function crossover_media_simples(parents::Matrix{Float64})
    p1 = parents[1, :]
    p2 = parents[2, :]
    child = (p1 .+ p2) ./ 2.0
    return [child'] 
end

function graph_swap_mutation!(x::Vector{Float64}, m_adj)
    n = length(x)
    if n < 2 return x end
    v1 = rand(1:n)
    vizinhos = [i for i in 1:n if m_adj[v1, i] == 1]
    if isempty(vizinhos) return x end
    v2 = rand(vizinhos)
    x[v1], x[v2] = x[v2], x[v1]
    return x
end 

function Metaheuristics.initialize!(status, parameters::CustomGAParams, problem, information, options)
    initial_state = Metaheuristics.gen_initial_state(problem, parameters, information, options, status)
    parameters.last_best = minimum([s.f for s in initial_state.population])
    parameters.stag_iters = 0
    return initial_state
end 

function Metaheuristics.update_state!(state::Metaheuristics.State, parameters::CustomGAParams, problem, information, options)
    global matriz_adj
    
    function tournament_select(pop)
        candidates = rand(pop, 2)
        return candidates[argmin([s.f for s in candidates])].x
    end 

    p1 = tournament_select(state.population)
    p2 = tournament_select(state.population)
    parents = vcat(p1', p2')

    offsprings = crossover_media_simples(parents)

    for i in axes(offsprings, 1)
        graph_swap_mutation!(offsprings[i, :], matriz_adj)
        vec = offsprings[i, :]
        fit = float(problem.f(vec))
        push!(state.population, Metaheuristics.xf_solution(vec, fit))
    end 

    # Elitismo e Manutenção da População 
    sort!(state.population, by = s -> s.f)
    state.population = state.population[1:min(parameters.N, length(state.population))]

    melhor_atual = state.population[1].f
    if melhor_atual < parameters.last_best
        parameters.last_best = melhor_atual
        parameters.stag_iters = 0
    else
        parameters.stag_iters += 1
    end 

    return parameters.stag_iters < parameters.stag_limit
end 

function Metaheuristics.final_stage!(state, parameters::CustomGAParams, problem, information, options)
    return state
end

# Lógica do experimento

function run_ga_experiment(file_name::String, k_limit::Int, N_pop::Int)
    global matriz_adj
    num_v, num_a = leInfo!(file_name)
    matriz_adj = zeros(Int, num_v, num_v)
    leArestas!(file_name, matriz_adj) 

    fitness_func(ind) = float(maximum(NOVOcoloracaoHarmonicaGuloso!(matriz_adj, sortperm(ind, rev=true))))
    bounds = [zeros(num_v) ones(num_v)]'
    
    algo = Metaheuristics.Algorithm(CustomGAParams(N = N_pop, stag_limit = k_limit))
    
    time_taken = @elapsed begin
        result = Metaheuristics.optimize(fitness_func, bounds, algo)
    end
    
    return time_taken, Int(Metaheuristics.minimum(result)), num_v, num_a
end

function main()
    N_REPETITIONS = 5
    K_LIMIT = 50
    N_POP = 1000 

    # POSSIVELMENTE REMOVER A FILTRAGEM DE TAMANHO PARA EXPERIMENTOS COM TODOS OS ARQUIVOS
    all_files = filter(f -> startswith(f, "bi_") && endswith(f, ".col"), readdir())
    filtered = filter(all_files) do f
        m = match(r"bi_a(\d+)_b(\d+)", f)
        m !== nothing && parse(Int, m.captures[1]) <= 100 && parse(Int, m.captures[2]) <= 100
    end
    sort!(filtered)

    detailed_results = []

    for file in filtered
        println("Processando $file...")
        t_rep, c_rep = Float64[], Int[]
        
        # Extração de parâmetros 
        m = match(r"bi_a(\d+)_b(\d+)_p(\d+)%_v(\d+)\.col", file)
        a, b, p, v = parse(Int, m[1]), parse(Int, m[2]), parse(Int, m[3]), parse(Int, m[4])

        local nv, na
        for i in 1:N_REPETITIONS
            t, c, nv, na = run_ga_experiment(file, K_LIMIT, N_POP)
            push!(t_rep, t); push!(c_rep, c)
        end 

        # Armazenamento com arredondamento 
        push!(detailed_results, Dict(
            :instancia => file, :a => a, :b => b, :p => p, :v => v, :N => nv, :M => na,
            :mean_chi => round(mean(c_rep), digits=2),
            :se_chi => round(std(c_rep)/sqrt(N_REPETITIONS), digits=2),
            :mean_time => round(mean(t_rep), digits=5),
            :se_time => round(std(t_rep)/sqrt(N_REPETITIONS), digits=5)
        ))
    end

    # Geração do CSV Ordenado
    df = DataFrame(detailed_results)
    select!(df, [:instancia, :a, :b, :p, :v, :N, :M, :mean_chi, :se_chi, :mean_time, :se_time])
    CSV.write("results_GA_Bipartite_Final.csv", df)
    
    println("\n--- Experimentos Concluídos! CSV gerado. ---")
end

main()