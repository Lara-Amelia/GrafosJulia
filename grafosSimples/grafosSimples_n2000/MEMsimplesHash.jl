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

global ADJ = Vector{Vector{Int}}()
global V = 0

# estrutura para memoização de valores de fitness (via tabela hash)
const FITNESS_CACHE = Dict{UInt64, Float64}()

function hash_perm(p::Vector{Int})
    h = UInt64(1469598103934665603)
    for x in p
        h ⊻= UInt64(x)
        h *= 1099511628211
    end
    return h
end

function fitness_harmonious_coloring(individual::Vector{Float64})
    lista_prioridade = sortperm(individual, rev = true)
    h = hash_perm(lista_prioridade)

    if haskey(FITNESS_CACHE, h)
        return FITNESS_CACHE[h]
    end

    cores_vertices = coloracaoHarmonicaAdjVetAux!(ADJ, lista_prioridade)
    fitness = maximum(cores_vertices)
    FITNESS_CACHE[h] = fitness
    return fitness
end

function graph_swap_mutation!(x::Vector{Float64})
    n = length(x)
    while true
        v1 = rand(1:n)
        vizinhos = ADJ[v1]
        if !isempty(vizinhos)
            v2 = rand(vizinhos)
            x[v1], x[v2] = x[v2], x[v1]
            return x
        end
    end
end

mutable struct CustomGAParams <: Metaheuristics.AbstractParameters
    N::Int
    p_mutation::Float64
    stag_limit::Int
    last_best::Float64
    stag_iters::Int
end

CustomGAParams(; N = 200, p_mutation = 0.5, stag_limit = 50, 
                last_best = -1, stag_iters = 0) =
    CustomGAParams(N, p_mutation, stag_limit, last_best, stag_iters)

function Metaheuristics.initialize!(
    status,
    parameters::CustomGAParams,
    problem::Metaheuristics.Problem,
    information::Metaheuristics.Information,
    options::Metaheuristics.Options
)
    parameters.stag_iters = 0
    parameters.last_best = Inf
    state = Metaheuristics.gen_initial_state(problem, parameters, information, options, status)
    sort!(state.population, by = s -> s.f)
    parameters.last_best = state.population[1].f
    state.best_sol = deepcopy(state.population[1])
    return state
end

function Metaheuristics.update_state!(state, parameters::CustomGAParams, problem, information, options)
    pop = state.population
    N = length(pop)
    offspring = Metaheuristics.xf_solution[]

    function tournament_select(pop)
        candidates = rand(pop, 2)
        return candidates[argmin([c.f for c in candidates])].x
    end

    for _ in 1:div(N, 2)
        parent1 = tournament_select(state.population)
        parent2 = tournament_select(state.population)
        child_x = (parent1 .+ parent2) ./ 2.0
        if rand() < parameters.p_mutation
            graph_swap_mutation!(child_x)
        end
        fit = Float64(fitness_harmonious_coloring(child_x))
        push!(offspring, Metaheuristics.xf_solution(child_x, fit))
    end

    combined = vcat(pop, offspring)
    sort!(combined, by = s -> s.f)
    state.population = combined[1:N]
    gen_best = state.population[1]

    if state.best_sol === nothing || gen_best.f < state.best_sol.f
        state.best_sol = deepcopy(gen_best)
    end

    if gen_best.f <= parameters.last_best
        parameters.last_best = gen_best.f
        parameters.stag_iters = 0
    else
        parameters.stag_iters += 1
    end
    return parameters.stag_iters < parameters.stag_limit
end

function Metaheuristics.final_stage!(state, parameters::CustomGAParams, problem, information, options)
    return state
end

function run_ga_experiment(file_name::String, k_limit::Int, N_pop::Int)
    bounds = [zeros(V) ones(V)]'
    params = CustomGAParams(N=N_pop, p_mutation=0.5, stag_limit=k_limit, last_best=Inf)
    opt_settings = Metaheuristics.Options(f_calls_limit = typemax(Int), iterations = 1000, store_convergence = true)
    my_ga = Metaheuristics.Algorithm(params, options = opt_settings)
    result = Metaheuristics.optimize(fitness_harmonious_coloring, bounds, my_ga)
    return Int(Metaheuristics.minimum(result))
end

function main()
    N_REPETITIONS = 5
    K_STAG = 50
    N_POP = 200

    # para os testes de memória usaremos somente os grafos versão 1
    all_files = filter(f -> endswith(f, "v1.col"), readdir())
    sort!(all_files)

    results_main = []

    @showprogress 1 "Processando: " for (idx, file) in enumerate(all_files)
        println("\n[$idx/$(length(all_files))] Processing: $file")
        
        # extração de parâmetros do nome do arquivo
        m_n = match(r"n(\d+)", file); n_param = m_n !== nothing ? parse(Int, m_n.captures[1]) : 0
        m_p = match(r"p(\d+)", file); p_param = m_p !== nothing ? parse(Int, m_p.captures[1]) : 0
        m_v = match(r"v(\d+)", file); v_param = m_v !== nothing ? parse(Int, m_v.captures[1]) : 0

        t_tot, chi, mem_history = Float64[], Int[], Float64[]
        empty!(FITNESS_CACHE) 
        
        global ADJ, V
        num_v, num_a = leInfo!(file)
        ADJ = [Int[] for _ in 1:num_v]
        leArestasLista!(file, ADJ)
        V = num_v

        for i in 1:N_REPETITIONS
            GC.gc() 
            
            elapsed = @elapsed begin
                ch = run_ga_experiment(file, K_STAG, N_POP)
            end
            
            # footprint da run atual em megabytes
            current_mem = Base.gc_live_bytes() / 1024^2
            
            push!(t_tot, elapsed) 
            push!(chi, ch)
            push!(mem_history, current_mem)
            @printf("   Run %d Done (Chi: %d, Mem: %.2f MB)\n", i, ch, current_mem)
        end

        se(v) = round(std(v) / sqrt(length(v)), digits=5)
        meanr(v) = round(mean(v), digits=5)

        push!(results_main, Dict(
            :instancia => file,
            :n => n_param,
            :p => p_param,
            :v => v_param,
            :N => num_v,
            :M => num_a,
            :mean_chi => mean(chi),
            :se_chi => se(chi),
            :mean_time => meanr(t_tot),
            :se_time => se(t_tot),
            :peak_mem_mb => maximum(mem_history), # pico de uso de memória
            :mean_mem_mb => meanr(mem_history)    # média de uso em todas as runs
        ))
    end

    df_main = DataFrame(results_main)
    sort!(df_main, [:n, :p, :v])

    CSV.write("results_GA_MEM_n2000.csv", df_main)
    println("\n--- Experimentos concluídos. CSV gerado com métricas de memória. ---")
end

main()