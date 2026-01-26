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

# hash simples para a permutação (associamos permutações a valores na tabela hash
# de forma a recuperar o valor de fitness, evitando recálcu-los que seriam caros)
# listas de prioridade que são essencialmente iguais geram o mesmo hash e armazenamos seu fitness
# usamos o offset base do algoritmo FNV-1a (para 64 bits) - álgebra A
function hash_perm(p::Vector{Int})
    h = UInt64(1469598103934665603)
    for x in p
        h ⊻= UInt64(x)     # XOR bit a bit
        h *= 1099511628211 # multiplicação pelo primo FNV
    end
    return h
end

# função de fitness utilizando tabela hash (recuperar valores para ordens já calculadas)
# e lista de adjacência (ao invés de matriz de adjacência)
function fitness_harmonious_coloring(individual::Vector{Float64})
    lista_prioridade = sortperm(individual, rev = true)
    h = hash_perm(lista_prioridade)

    # se a chave associada a ordem já está na tabela hash, apenas recuperamos o fitness
    if haskey(FITNESS_CACHE, h)
        return FITNESS_CACHE[h]
    end

    # fazemos a coloração para a ordem sse não encontramos associada a chave na tabela hash
    cores_vertices = coloracaoHarmonicaAdjVetAux!(ADJ, lista_prioridade)
    fitness = maximum(cores_vertices)
    FITNESS_CACHE[h] = fitness
    return fitness
end

# tentamos até encontrar um vértice com pelo menos um vizinho para a mutação
# caso contrário teríamos puramente repetição de indivíduos e, com isso, menos diversidade 
function graph_swap_mutation!(x::Vector{Float64})
    n = length(x)

    while true
        v1 = rand(1:n)                  # vértice escolhido aleatoriamente 
        vizinhos = ADJ[v1]              # recuperação dos vizinhos do primeiro vértice

        if !isempty(vizinhos)           # se há algum vizinho, selecionamos random para swap
            v2 = rand(vizinhos)
            x[v1], x[v2] = x[v2], x[v1] # swap de prioridades entre vértices vizinhos
            return x
        end
    end
end

#crossover retornando somente 1 filho, média simples dos pais
function crossover_media_simples(parents::Matrix{Float64})
    p1 = parents[1, :]
    p2 = parents[2, :]
    child = (p1 .+ p2) ./2.0
    return [child']
    #return child
end 

# NOTA: implementar os demais métodos personalizados adaptados se for necessário

# definição do GA personalizado - poderíamos incluir outros campos, como 
# probabilidade associada a crossover, a mutação, etc (por enquanto não será aplicado)
mutable struct CustomGAParams <: Metaheuristics.AbstractParameters
    N::Int
    p_mutation::Float64
    stag_limit::Int
    last_best::Float64 # provavelmente será Int p/ nro. de cores = fitness
    stag_iters::Int
end

# o valor k (limite de iterações estagnadas será lido da entrada)
CustomGAParams(; N = 200, p_mutation = 0.5, stag_limit = k, 
                 last_best = -1, stag_iters = 0) =
    CustomGAParams(N, p_mutation, stag_limit, last_best, stag_iters)

# método para gerar um indivíduo do genético a partir da ordem por grau definida
function create_greedy_individual(permutation::Vector{Int})
    n = length(permutation)
    x = zeros(Float64, n)
    
    for (rank, node_idx) in enumerate(permutation)
        # quanto maior o rank do vértice, maior é o valor em Float64 recebido 
        # (ou seja, mantemos a mesma ordem)
        x[node_idx] = (n - rank + 1) / n
    end
    return x
end

# inicialização do estado da população
# NOTA: alteração para incluir na pop. inicial as ordenações por grau máximo e grau mínimo
# (ainda não acrescentamos por grau de saturação devido à ordem dinâmica - pode ser complexo demais)
# deve haver diversificação na população para que o genetico consiga evlouir a partir das ordens gulosas
function Metaheuristics.initialize!(
    status,
    parameters::CustomGAParams,
    problem::Metaheuristics.Problem,
    information::Metaheuristics.Information,
    options::Metaheuristics.Options
)
    state = nothing
    parameters.stag_iters = 0
    parameters.last_best = Inf

    # pop. aleatória base 
    state = Metaheuristics.gen_initial_state(
            problem, parameters, information, options, status)

    # "injeção" de ordens gulosas
    #=greedy_orders = [
        obtemPrioridadePorGrauAdj(ADJ, true),
        obtemPrioridadePorGrauAdj(ADJ, false) ]

    for (i, order) in enumerate(greedy_orders)
        if i <= length(state.population)
            greedy_x = create_greedy_individual(order)
            greedy_f = fitness_harmonious_coloring(greedy_x)
            state.population[i] = Metaheuristics.xf_solution(greedy_x, Float64(greedy_f))
        end
    end=#

    sort!(state.population, by = s -> s.f)

    parameters.last_best = state.population[1].f
    state.best_sol = deepcopy(state.population[1])

    return state
end

function Metaheuristics.update_state!(
    state,
    parameters::CustomGAParams,
    problem,
    information,
    options
)
    pop = state.population
    N = length(pop)
    offspring = Metaheuristics.xf_solution[]

    # debug: print initial population fitness distribution
    #@info "Population fitnesses: " [s.f for s in pop]

    function tournament_select(pop)
        k = 2
        candidates = rand(pop, k)
        fitnesses = [c.f for c in candidates]
        return candidates[argmin(fitnesses)].x
    end

    # geração de filhos (geramos vários)
    for _ in 1:div(N, 2)
        parent1 = tournament_select(state.population)
        parent2 = tournament_select(state.population)

        #parents = vcat(parent1', parent2')  
        #child_x = crossover_media_simples(parents)
        child_x = (parent1 .+ parent2) ./2.0

        # mutação baseada na prob. de mutação definida
        if rand() < parameters.p_mutation
            graph_swap_mutation!(child_x)
        end

        fit = Float64(fitness_harmonious_coloring(child_x))
        #@info "Child fitness: $fit, child sample priorities: $(child_x[1:5])"
        push!(offspring, Metaheuristics.xf_solution(child_x, fit)) # objeto armazenando filhos
    end

    # elitist replacement
    combined = vcat(pop, offspring)    # une pop. atual e filhos produzidos
    sort!(combined, by = s -> s.f)     # ordena pelo valor do fitness
    state.population = combined[1:N]   # seleciona somente os de melhor fitness na pop. combinada

    
    # atualização de parâmetros e estagnação
    gen_best = state.population[1]

    #@info "Generation best fitness: $(gen_best.f)"

    # atualização da melhor solução
    if state.best_sol === nothing || gen_best.f < state.best_sol.f
        state.best_sol = deepcopy(gen_best)
        #@info "New GLOBAL best: $(state.best_sol.f)"
    end

    # estagnação
    if gen_best.f <= parameters.last_best
        parameters.last_best = gen_best.f
        parameters.stag_iters = 0
    else
        parameters.stag_iters += 1
    end

    if parameters.stag_iters >= parameters.stag_limit
        #@info "Parada por estagnação"
        return false
    end

    return true

end

# caso queiramos fazer algum tratamento extra sobre o estado final do GA (por agora não)
function Metaheuristics.final_stage!(
    state,
    parameters::CustomGAParams,
    problem,
    information,
    options
)
    return state
end

function run_ga_experiment(file_name::String, k_limit::Int, N_pop::Int)
    bounds = [zeros(V) ones(V)]'

    params = CustomGAParams(N=200, p_mutation=0.5, stag_limit=50, last_best=Inf)

    opt_settings = Metaheuristics.Options(f_calls_limit = typemax(Int), iterations = 1000, store_convergence = true)

    my_ga = Metaheuristics.Algorithm(params, options = opt_settings)
    problem = Metaheuristics.Problem(fitness_harmonious_coloring, bounds)

    result = Metaheuristics.optimize(fitness_harmonious_coloring, bounds, my_ga)

    return Int(Metaheuristics.minimum(result))
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
    
    TARGET_A = 1000
    TARGET_B = 5000

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
        empty!(FITNESS_CACHE)

        global ADJ, V
        num_v, num_a = leInfo!(file)
        ADJ = [Int[] for _ in 1:num_v]
        leArestasLista!(file, ADJ)
        V = num_v

        for i in 1:N_REPETITIONS
            elapsed = @elapsed begin
                ch = run_ga_experiment(file, K_LIMIT, N_POP)
            end
            push!(t_tot, elapsed) 
            push!(chi, ch)
            @printf("   Run %d Done (Chi: %d)\n", i, ch)
        end

        se(v) = std(v) / sqrt(length(v))
        meanr(v) = mean(v)

        # dicionário na ordem desejada para a saída do .csv final
        push!(results_main, Dict(
            :a => params[:a_param],
            :b => params[:b_param],
            :N => num_v,
            :p => params[:p_param],
            :M => num_a,
            :v => params[:v_param],
            :mean_time => meanr(t_tot),
            :se_time => se(t_tot),
            :mean_chi => meanr(chi),
            :se_chi => se(chi),
            :instancia => file
        ))

    end

    df_main = DataFrame(results_main)
    col_order = [:a, :b, :N, :p, :M, :v, :mean_time, :se_time, :mean_chi, :se_chi, :instancia]
    select!(df_main, col_order)

    CSV.write("results_GA_Final_a1000_b5000_final.csv", df_main)

    println("\n--- Experimentos concluídos ---")
end

main()