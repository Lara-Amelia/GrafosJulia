# versão do GA utilizando ao máximo os métodos do framework
using Metaheuristics
using LinearAlgebra
using Random

import Metaheuristics: initialize!, update_state!, final_stage!, gen_initial_state
import Metaheuristics: SBX_crossover, polynomial_mutation!, create_solution, is_better
import Metaheuristics: reset_to_violated_bounds!

#inclui o arquivo com as funções de coloração gulosa
include("colorGul.jl")

# leitura do arquivo de entrada e população da lista de adjacência 
println("Insira o nome do arquivo a ser lido: ")
nome_arquivo = readline()

num_vertices, num_arestas = leInfo!(nome_arquivo)
println("nro de vertices: $num_vertices")
println("nro de arestas: $num_arestas")

# cria lista de adjacência
lista_adj = [Int[] for _ in 1:num_vertices]
leArestasLista!(nome_arquivo, lista_adj)

# referência global constante para utilizar no GA
const ADJ = lista_adj
const V = num_vertices

# probabilidade associada a crossover, a mutação, etc (por enquanto não será aplicado)
mutable struct CustomGAParams <: Metaheuristics.AbstractParameters
    N::Int
    p_mutation::Float64
    stag_limit::Int
    last_best::Float64 # provavelmente será Int p/ nro. de cores = fitness
    stag_iters::Int
    selection::Metaheuristics.TournamentSelection
end

# o valor k (limite de iterações estagnadas será lido da entrada)
#=CustomGAParams(; N = 100, p_mutation = 0.5, stag_limit = k, 
                 last_best = -1, stag_iters = 0) =
    CustomGAParams(N, p_mutation, stag_limit, last_best, stag_iters)=#

function CustomGAParams(; N = 100, p_mutation = 0.5, stag_limit = 50, 
                         last_best = Inf, stag_iters = 0)
    selection_strategy = Metaheuristics.TournamentSelection(K=2, N=N)
    return CustomGAParams(N, p_mutation, stag_limit, last_best, stag_iters, selection_strategy)
end

# fitness da coloração sem utilizar hash (verificar depois se isso ajuda ou atrapalha)
function fitness_harmonious_coloring(individual::Vector{Float64})
    lista_prioridade = sortperm(individual, rev = true)

    cores_vertices = coloracaoHarmonicaAdjVetAux!(ADJ, lista_prioridade)
    fitness = maximum(cores_vertices)
    return fitness
end

# dessa vez tentaremos chamar os métodos personalizados pelo wrapper da própria biblioteca,
# possivelmente teremos de alterar alguns deles (para adequar às formas e tipos esperados)
# aparentemente, não será necessário sobrescrever o método initialize! (ou podemos apenas chamar gen_initial_state)
# population será um array simples de valores 

function crossover_simple_mean!(population)
    # seleciona os elementos/pais (os que estão lado a lado são pareados)
    Q = positions(population) # matriz para representar indivíduos e genes

    Q1 = Q[1:2:end-1, :]
    Q2 = Q[2:2:end,   :]

    children = (Q1 .+ Q2) ./ 2.0

    # children é uma matriz 50 x n_vertices
    return children 
end

# um jeito mais fácil de realmente integrar os novos métodos pode ser reescrever update_state!,
# mas apenas substituindo as chamadas a crossover e mutação por nossos próprios métodos
function graph_swap_mutation!(Q::AbstractMatrix{Float64})
    n_individuals, n_genes = size(Q)
    adj_list = ADJ
    p = 0.5 # probabilidade de ocorência de mutação nos filhos do crossover

    # "pré-seleção" dos filhos que serão mutados
    to_mutate = findall(rand(n_individuals) .< p)

    for i in to_mutate
       v1 = rand(1:n_genes)
        vizinhos = ADJ[v1]
        
        # até que um vizinho aleatório tenha outros vizinhos
        if isempty(vizinhos)
            continue 
        end

        v2 = rand(vizinhos)

        # swap otimizado
        @inbounds begin
            tmp = Q[i, v1]
            Q[i, v1] = Q[i, v2]
            Q[i, v2] = tmp
        end 
    end
    # como alteramos a matriz "in-place", não é necessário retornar nada
end

function replacement_elitism(population, offsprings, N)
    combined = append!(population, offsprings)

    sort!(combined, alg=PartialQuickSort(N), by = s -> s.f)

    deleteat!(population, (N+1):length(population))
end

# sobre-escrita de métodos do próprio GA (wrappers para as suas funções)
function Metaheuristics.initialize!(status,
    parameters::CustomGAParams,
    problem::Metaheuristics.Problem,
    information::Metaheuristics.Information,
    options::Metaheuristics.Options)

    state = nothing
    parameters.stag_iters = 0
    parameters.last_best = Inf

    state = Metaheuristics.gen_initial_state(
            problem, parameters, information, options, status)
end 

function Metaheuristics.update_state!(state,
    parameters::CustomGAParams,
    problem, # como avaliamos o fitness das soluções
    information,
    options)
    pop = state.population

    parent_mask = Metaheuristics.selection(pop, parameters.selection)

    Q = crossover_simple_mean!(pop[parent_mask])

    graph_swap_mutation!(Q)

    offsprings = Metaheuristics.create_solutions(Q, problem)

    replacement_elitism(pop, offsprings, params.N)
    current_best = Metaheuristics.get_best(pop)

    if Metaheuristics.is_better(current_best, state.best_sol)
        state.best_sol = deepcopy(current_best)
    end
end

function Metaheuristics.final_stage!(
    state,
    parameters::CustomGAParams,
    problem,
    information,
    options
)
    return state
end

# definição do genético personalizado utilizando os novos métodos
bounds = [zeros(V) ones(V)]'

params = CustomGAParams(N=100, p_mutation=0.5, stag_limit=50, last_best=Inf)

opt_settings = Metaheuristics.Options(f_calls_limit = typemax(Int), iterations = 1000, store_convergence = true)

my_ga = Metaheuristics.Algorithm(params, options = opt_settings)
problem = Metaheuristics.Problem(fitness_harmonious_coloring, bounds)

result = Metaheuristics.optimize(fitness_harmonious_coloring, bounds, my_ga)

# resultados da otimização
@show Metaheuristics.minimum(result)
@show result