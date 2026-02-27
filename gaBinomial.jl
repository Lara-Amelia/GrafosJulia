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
    crossover::Metaheuristics.BinomialCrossover
end

# o valor k (limite de iterações estagnadas será lido da entrada)
#=CustomGAParams(; N = 100, p_mutation = 0.5, stag_limit = k, 
                 last_best = -1, stag_iters = 0) =
    CustomGAParams(N, p_mutation, stag_limit, last_best, stag_iters)=#

function CustomGAParams(; N = 100, p_mutation = 0.5, stag_limit = 50, last_best = Inf, stag_iters = 0)
    selection_strategy = Metaheuristics.TournamentSelection(K=2, N=N)
    cross_op = Metaheuristics.BinomialCrossover(p = 0.5, n_offsprings = 2)
    
    return CustomGAParams(N, p_mutation, stag_limit, last_best, 
                          stag_iters, selection_strategy, cross_op)
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
    # adj_list = ADJ
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

# função para critério de parada desejado
function stop_on_stagnation(parameters)
    return parameters.stag_iters >= parameters.stag_limit
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

    # crossover binomial default com p = 0.5 e n_offsprings = 2
    Q = Metaheuristics.crossover(pop[parent_mask], parameters.crossover)

    graph_swap_mutation!(Q)

    offsprings = Metaheuristics.create_solutions(Q, problem)

    replacement_elitism(pop, offsprings, params.N)
    current_best = Metaheuristics.get_best(pop)

    # lógica para critério de parada por gerações em estagnação
    if parameters.last_best == Inf
        parameters.last_best = current_best.f
    end

    if current_best.f < parameters.last_best
        # se houve melhora, resetamos o contador de stag iters
        parameters.last_best = current_best.f
        parameters.stag_iters = 0
    else
        # se não houve melhora na geração, incrementamos o contador
        parameters.stag_iters += 1
    end

    if parameters.stag_iters < parameters.stag_limit
        println("nro. de iterações estagnadas: $(parameters.stag_iters)")
        println("AINDA NÃO ATINGIMOS O LIMITE DE ITERAÇÕES SEM MELHORA")
    end

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

# redefinição do método de critério de parada para garantir o uso de nosso critério
function Metaheuristics.stop_criteria!(
        status,
        parameters::CustomGAParams,
        problem::Metaheuristics.Problem,
        information::Metaheuristics.Information,
        options::Metaheuristics.Options,
    )

    status.stop = status.stop || Metaheuristics.call_limit_stop_check(status, information, options) ||
                  Metaheuristics.iteration_stop_check(status, information, options)  ||
                  Metaheuristics.time_stop_check(status, information, options)       ||
                  Metaheuristics.accuracy_stop_check(status, information, options)

    # inclusão do critério de parada personalizado - enviamos a mensagem de other_limit nesse caso
    # NOTA: calibrar o valor de stag_limit?
    if parameters.stag_iters >= parameters.stag_limit
        status.stop = true
        status.termination_status_code = Metaheuristics.OTHER_LIMIT
        return
    end  

    return
end

# definição do genético personalizado utilizando os novos métodos
bounds = [zeros(V) ones(V)]'

params = CustomGAParams(N=100, p_mutation=0.5, stag_limit=50, last_best=Inf)

#opt_settings = Metaheuristics.Options(f_calls_limit = typemax(Int), iterations = 1000, store_convergence = true)
# novos parâmetros definidos em Options para evitar convergências "automáticas" da biblioteca
opt_settings = Metaheuristics.Options(
    f_calls_limit = typemax(Int), 
    iterations = 10000, 
    store_convergence = true,
    f_tol = -1, 
    x_tol = -1
)

my_ga = Metaheuristics.Algorithm(params, options = opt_settings)
problem = Metaheuristics.Problem(fitness_harmonious_coloring, bounds)

checkall = checkany = [stop_on_stagnation]
result = Metaheuristics.optimize(fitness_harmonious_coloring, bounds, my_ga)

# resultados da otimização
@show Metaheuristics.minimum(result)
@show result
@show result.stop
@show result.termination_status_code