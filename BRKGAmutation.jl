using Metaheuristics
using LinearAlgebra
using Random 

include("colorGul.jl")

# entrada de dados
println("Insira o nome do arquivo a ser lido: ")
nome_arquivo = readline()

num_vertices, num_arestas = leInfo!(nome_arquivo)
println("nro de vertices: $num_vertices")
println("nro de arestas: $num_arestas")

lista_adj = [Int[] for _ in 1:num_vertices]
leArestasLista!(nome_arquivo, lista_adj)

const ADJ = lista_adj
const V = num_vertices

const P_BUFFER = zeros(Int, V)

function fitness_harmonious_coloring(individual::Vector{Float64})
    sortperm!(P_BUFFER, individual, rev = true)
    cores_vertices = coloracaoHarmonicaAdjVetAux!(ADJ, B_BUFFER)
    return Float64(maximum(cores_vertices))
end

# wrapper para adicionar novo critério de parada e fazer os "semi-overrides necessários"
mutable struct StagnationWrapper <: Metaheuristics.AbstractParameters
    internal_params::Metaheuristics.AbstractParameters
    stag_limit::Int
    stag_iters::Int
    last_best::Float64
end

# delegação de inicialização
function Metaheuristics.initialize!(status, parameters::StagnationWrapper, problem, information, options)
    # inicialização dos contadores do wrapper
    parameters.stag_iters = 0
    parameters.last_best = Inf
    
    return Metaheuristics.initialize!(status, parameters.internal_params, problem, information, options)
end 

# eveolução e contagem de iterações em estagnação
function Metaheuristics.update_state!(state, parameters::StagnationWrapper, problem, information, options)
    # executa um "passo" do BRKGA usando os parâmetros e operadores dados
    Metaheuristics.update_state!(state, parameters.internal_params, problem, information, options)

    # monitoramento após evolução
    current_best = state.best_sol.f
    
    if current_best < parameters.last_best
        parameters.last_best = current_best
        parameters.stag_iters = 0
    else
        parameters.stag_iters += 1
    end
    println("melhor fitness da iteração: $current_best")
end

# novo método para critério de parada
function Metaheuristics.stop_criteria!(status, parameters::StagnationWrapper, problem, information, options)
    # checks padrão
    status.stop = status.stop || 
                  Metaheuristics.iteration_stop_check(status, information, options) ||
                  Metaheuristics.time_stop_check(status, information, options)

    # check para limite de iterações estagnadas
    if parameters.stag_iters >= parameters.stag_limit
        status.stop = true
        status.termination_status_code = Metaheuristics.OTHER_LIMIT
        # println("\n > Parada por estagnação: $(parameters.stag_limit) gerações (Best: $(status.best_sol.f))")
    end
    return status.stop
end

# override de estágio final
function Metaheuristics.final_stage!(state, parameters::StagnationWrapper, problem, information, options)
    return Metaheuristics.final_stage!(state, parameters.internal_params, problem, information, options)
end

# operador de mutação customizado e wrapper correspondente
struct GraphSwapMutation
    p::Float64
end

# mutação personalizada com swaps no grafo
# para o BRKGA, a probabilidade será de 1 (os selecionados com certeza são mutados)
function Metaheuristics.mutation!(Q::AbstractMatrix{Float64}, mut::GraphSwapMutation)
    n_individuals, n_genes = size(Q)
    to_mutate = findall(rand(n_individuals) .< mut.p)

    for i in to_mutate
        v1 = rand(1:n_genes)
        vizinhos = ADJ[v1]
        if isempty(vizinhos) continue end

        v2 = rand(vizinhos)
        @inbounds begin
            tmp = Q[i, v1]
            Q[i, v1] = Q[i, v2]
            Q[i, v2] = tmp
        end 
    end
end

# construção do BRKGA com métodos personalizados (similar ao que é feito na biblioteca)
function MyStagnatedBRKGA(;
        num_elites = 20,
        num_mutants = 10,
        num_offsprings = 70,
        bias = 0.7,
        stag_limit = 50,
        kargs...
    )

    N = num_elites + num_mutants + num_offsprings
    
    # criação do GA base - BRKGA com mutação e critério de parada personalizados
    base_ga = Metaheuristics.GA(;
        initializer = Metaheuristics.RandomInBounds(; N),
        selection   = Metaheuristics.BiasedSelection(num_elites, num_offsprings),
        crossover   = Metaheuristics.BinomialCrossover(p = bias, n_offsprings = 1),
        mutation    = GraphSwapMutation(1.0),
        environmental_selection = Metaheuristics.ElitistReplacement(),
        kargs...
    )

    # parâmetros do novo GA no wrapper
    wrapped_params = StagnationWrapper(base_ga.parameters, stag_limit, 0, Inf)

    # retorna o objeto algorithm com os parâmetros definidos
    return Metaheuristics.Algorithm(wrapped_params, options = base_ga.options)
end

# instancia o algoritmo
genetico = MyStagnatedBRKGA(
    num_elites = 20, 
    num_mutants = 10, 
    num_offsprings = 70, 
    bias = 0.7,
    stag_limit = 50 # valor limite para quantidade de iterações estagnadas
)

# inicializa o objeto options
genetico.options = Metaheuristics.Options(
    f_calls_limit = typemax(Int), 
    iterations = 10000, 
    store_convergence = true,
    f_tol = -1, # desativa outros critérios de convergência
    x_tol = -1
)

bounds = [zeros(V) ones(V)]'

result = optimize(fitness_harmonious_coloring, bounds, genetico)

# Resultados finais
println("\n" * "="^20)
@show result
@show Metaheuristics.minimum(result)
@show result.iteration
@show result.termination_status_code