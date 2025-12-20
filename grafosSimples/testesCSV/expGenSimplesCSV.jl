# ===========================
#   HARMONIOUS COLORING GA - SIMPLE GRAPHS (ALL FILES - CSV)
# ===========================

using Metaheuristics
using LinearAlgebra
using Statistics
using Printf
using DataFrames
using CSV

# --- Carregar Módulo de Coloração (colorGul.jl) ---
include("../../colorGul.jl")

# --- Variáveis Globais ---
global matriz_adj = Matrix{Int}(undef, 0, 0)
global num_vertices = 0

# ===========================
#   Genetic Algorithm Setup
# ===========================

mutable struct CustomGAParams <: Metaheuristics.AbstractParameters
    N::Int
    p_crossover::Float64
    p_mutation::Float64
    stag_limit::Int
    last_best::Float64
    stag_iters::Int
end

# Construtor com calibração exata (N=1000 conforme gaUpdate.jl)
CustomGAParams(; N = 1000, p_crossover = 0.5, p_mutation = 0.5,
                 stag_limit = 50, last_best = Inf, stag_iters = 0) =
    CustomGAParams(N, p_crossover, p_mutation, stag_limit, last_best, stag_iters)

# Crossover de filho único (média simples) - SEPARADO conforme solicitado
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

# --- Metaheuristics Hooks ---

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
        return candidates[argmin([c.f for c in candidates])].x
    end

    p1 = tournament_select(state.population)
    p2 = tournament_select(state.population)
    parents = vcat(p1', p2')

    # Uso da função de filho único
    offsprings = crossover_media_simples(parents)

    for i in axes(offsprings, 1)
        graph_swap_mutation!(offsprings[i, :], matriz_adj)
        vec = offsprings[i, :]
        fit = float(problem.f(vec))
        push!(state.population, Metaheuristics.xf_solution(vec, fit))
    end

    # Elitismo: Mantém os N melhores
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

function Metaheuristics.final_stage!(state, parameters, problem, information, options)
    return state
end

# ===========================
#   Lógica do Experimento
# ===========================

function run_ga_experiment(file_name::String, k_limit::Int, N_pop::Int)
    global matriz_adj
    num_v, num_a = leInfo!(file_name)
    matriz_adj = zeros(Int, num_v, num_v)
    leArestas!(file_name, matriz_adj)

    fitness_func(ind) = float(maximum(NOVOcoloracaoHarmonicaGuloso!(matriz_adj, sortperm(ind, rev=true))))
    bounds = [zeros(num_v) ones(num_v)]'
    
    algo = Metaheuristics.Algorithm(CustomGAParams(N = N_pop, stag_limit = k_limit))
    
    # Medição robusta de tempo decorrido
    time_taken = @elapsed begin
        result = Metaheuristics.optimize(fitness_func, bounds, algo)
    end
    
    return time_taken, Int(Metaheuristics.minimum(result)), num_v, num_a
end

function extract_simple_params(file_name)
    m = match(r"simples_n(\d+)_p(\d+)%_v(\d+)\.col", file_name)
    if m !== nothing
        return (n=parse(Int, m[1]), p=parse(Int, m[2]), v=parse(Int, m[3]))
    else
        return (n=0, p=0, v=0)
    end
end

function main()
    N_REPETITIONS = 5
    K_LIMIT = 50
    N_POP = 1000

    # Coleta TODOS os arquivos sem filtros de N ou V
    all_files = filter(f -> startswith(f, "simples_") && endswith(f, ".col"), readdir())
    sort!(all_files)

    if isempty(all_files)
        println("Nenhum arquivo 'simples_*.col' encontrado.")
        return
    end

    num_files = length(all_files)
    println("--- Iniciando Experimentos GA: $num_files arquivos encontrados ---")

    csv_results = []

    for (i, file) in enumerate(all_files)
        println("($i/$num_files) Processando $file...")
        t_rep, c_rep = Float64[], Int[]
        
        params = extract_simple_params(file)
        local nv, na

        for rep in 1:N_REPETITIONS
            try
                t, c, nv, na = run_ga_experiment(file, K_LIMIT, N_POP)
                push!(t_rep, t); push!(c_rep, c)
            catch e
                @error "Erro no arquivo $file, repetição $rep: $e"
            end
        end

        if !isempty(c_rep)
            push!(csv_results, Dict(
                :instancia => file,
                :n => params.n, :p => params.p, :v => params.v,
                :N => nv, :M => na,
                :mean_chi => round(mean(c_rep), digits=2),
                :se_chi => round(std(c_rep)/sqrt(length(c_rep)), digits=2),
                :mean_time => round(mean(t_rep), digits=5),
                :se_time => round(std(t_rep)/sqrt(length(t_rep)), digits=5)
            ))
        end
    end

    # Geração e Ordenação Final do CSV
    df = DataFrame(csv_results)
    select!(df, [:instancia, :n, :p, :v, :N, :M, :mean_chi, :se_chi, :mean_time, :se_time])
    
    output_filename = "results_GA_Simple_All_Final.csv"
    CSV.write(output_filename, df)
    
    println("\n--- Concluído! Dados salvos em $output_filename ---")
end

main()