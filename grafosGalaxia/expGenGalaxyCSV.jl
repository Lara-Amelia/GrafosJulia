using Metaheuristics
using LinearAlgebra
using Statistics
using Printf
using DataFrames
using CSV

import Metaheuristics: initialize!, update_state!, final_stage!
import Metaheuristics: AbstractParameters, gen_initial_state, Algorithm

# Inclui as funções de coloração do seu repositório 
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

# Crossover de filho único (média simples) 
function crossover_media_simples(parents::Matrix{Float64})
    p1 = parents[1, :]
    p2 = parents[2, :]
    child = (p1 .+ p2) ./ 2.0
    return [child'] 
end

function graph_swap_mutation!(x::Vector{Float64}, matriz_adj_local)
    n = length(x)
    if n < 2 return x end
    v1 = rand(1:n)
    vizinhos = [i for i in 1:n if matriz_adj_local[v1, i] == 1]
    if isempty(vizinhos) return x end
    v2 = rand(vizinhos)
    x[v1], x[v2] = x[v2], x[v1]
    return x
end

# Ciclo de vida do GA
function Metaheuristics.update_state!(
    state::Metaheuristics.State,
    parameters::CustomGAParams,
    problem::Metaheuristics.Problem,
    information::Metaheuristics.Information,
    options::Metaheuristics.Options
)
    function tournament_select(pop)
        candidates = rand(pop, 2)
        return candidates[argmin([c.f for c in candidates])].x
    end

    p1 = tournament_select(state.population)
    p2 = tournament_select(state.population)
    parents = vcat(p1', p2')  

    offsprings = crossover_media_simples(parents)

    for i in axes(offsprings, 1)  
        vec = offsprings[i, :]
        fit = Float64(problem.f(vec))
        push!(state.population, Metaheuristics.xf_solution(vec, fit))
    end

    sort!(state.population, by = s -> s.f)
    state.population = state.population[1:min(parameters.N, length(state.population))]

    melhor_atual = state.population[1].f
    if melhor_atual < parameters.last_best || parameters.last_best == -1
        parameters.last_best = Int(melhor_atual)
        parameters.stag_iters = 0
    else
        parameters.stag_iters += 1
    end

    return parameters.stag_iters < parameters.stag_limit
end

function Metaheuristics.initialize!(status, parameters::CustomGAParams, problem, information, options)
    state = Metaheuristics.gen_initial_state(problem, parameters, information, options, status)
    parameters.last_best = Int(minimum([s.f for s in state.population]))
    parameters.stag_iters = 0
    return state
end

function Metaheuristics.final_stage!(state, p, prob, info, opts) return state end

# lógica do experimento

function run_ga_experiment(file_name, k_limit, n_pop)
    nv, na = leInfo!(file_name)
    matriz_adj = zeros(Int, nv, nv)
    leArestas!(file_name, matriz_adj)

    fitness_func(ind) = maximum(NOVOcoloracaoHarmonicaGuloso!(matriz_adj, sortperm(ind, rev=true)))
    bounds = [zeros(nv) ones(nv)]'
    
    algo = Metaheuristics.Algorithm(CustomGAParams(N = n_pop, stag_limit = k_limit))
    
    # Medição robusta de tempo
    time_taken = @elapsed begin
        result = Metaheuristics.optimize(fitness_func, bounds, algo)
    end
    
    return time_taken, Int(Metaheuristics.minimum(result))
end

function main()
    all_files = filter(f -> startswith(f, "galaxy_") && endswith(f, ".col"), readdir())
    filtered = filter(f -> (m = match(r"a(\d+)", f); m !== nothing && parse(Int, m.captures[1]) <= 25), all_files)
    sort!(filtered)

    results_list = []

    for file in filtered
        println("Processando $file...")
        t_rep, c_rep = Float64[], Int[]
        nv, na = leInfo!(file)
        
        # Extração manual de parâmetros para o CSV
        m_a = match(r"a(\d+)", file); a = m_a !== nothing ? parse(Int, m_a.captures[1]) : 0
        m_b = match(r"b(\d+)", file); b = m_b !== nothing ? parse(Int, m_b.captures[1]) : 0
        m_c = match(r"c(\d+)", file); c = m_c !== nothing ? parse(Int, m_c.captures[1]) : 0
        m_v = match(r"v(\d+)", file); v = m_v !== nothing ? parse(Int, m_v.captures[1]) : 0

        for i in 1:5
            t, c = run_ga_experiment(file, 50, 1000)
            push!(t_rep, t); push!(c_rep, c)
        end

        
        push!(results_list, Dict(
            :instancia => file, :a => a, :b => b, :c => c, :v => v, :N => nv, :M => na,
            :mean_chi => round(mean(c_rep), digits=2),
            :se_chi => round(std(c_rep)/sqrt(5), digits=2),
            :mean_time => round(mean(t_rep), digits=5), # 5 casas decimais
            :se_time => round(std(t_rep)/sqrt(5), digits=4)
        ))
    end

    df = DataFrame(results_list)
    select!(df, [:instancia, :a, :b, :c, :v, :N, :M, :mean_chi, :se_chi, :mean_time, :se_time])
    CSV.write("results_GA_Final.csv", df)
    println("--- CSV Gerado com tempos arredondados ---")
end

main()

# NOTA 1: verificar a formatação dos campos de tempo ao importar para excel ou planilhas
# NOTA 2: remover filtros de a <= 25 para rodar testes