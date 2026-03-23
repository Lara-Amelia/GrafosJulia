using Metaheuristics
using LinearAlgebra
using Random 

include("colorGul.jl")

println("Insira o nome do arquivo a ser lido: ")
nome_arquivo = readline()

num_vertices, num_arestas = leInfo!(nome_arquivo)
println("nro de vertices: $num_vertices")
println("nro de arestas: $num_arestas")

lista_adj = [Int[] for _ in 1:num_vertices]
leArestasLista!(nome_arquivo, lista_adj)

const ADJ = lista_adj
const V = num_vertices

function fitness_harmonious_coloring(individual::Vector{Float64})
    lista_prioridade = sortperm(individual, rev = true)
    cores_vertices = coloracaoHarmonicaAdjVetAux!(ADJ, lista_prioridade)
    return Float64(maximum(cores_vertices))
end

# wrapper para adicionar o novo critério de parada no alg.
# optamos por não reescrever os métodos de fato, mas sim usar wrappers para chamá-los
mutable struct StagnationWrapper <: Metaheuristics.AbstractParameters
    internal_params::Metaheuristics.AbstractParameters
    stag_limit::Int
    stag_iters::Int
    last_best::Float64
end

# delegações para o framework 
function Metaheuristics.initialize!(status, parameters::StagnationWrapper, problem, information, options)
    parameters.stag_iters = 0
    parameters.last_best = Inf
    return Metaheuristics.initialize!(status, parameters.internal_params, problem, information, options)
end 

# faz uma geração do ga e as operações para observar a quantidade de iterações estagnadas
function Metaheuristics.update_state!(state, parameters::StagnationWrapper, problem, information, options)
    Metaheuristics.update_state!(state, parameters.internal_params, problem, information, options)
    
    current_best = state.best_sol.f
    if current_best < parameters.last_best
        parameters.last_best = current_best
        parameters.stag_iters = 0
    else
        parameters.stag_iters += 1
    end
    println("melhor fitness da iteração: $current_best")
end

# adiciona o novo critério de parada por quantidade de iterações estagnadas
function Metaheuristics.stop_criteria!(status, parameters::StagnationWrapper, problem, information, options)
    status.stop = status.stop || 
                  Metaheuristics.iteration_stop_check(status, information, options) ||
                  Metaheuristics.time_stop_check(status, information, options)

    if parameters.stag_iters >= parameters.stag_limit
        status.stop = true
        status.termination_status_code = Metaheuristics.OTHER_LIMIT
        # println("\n > Parada por estagnação: $(parameters.stag_limit) iterações (Fitness: $(status.best_sol.f))")
    end
    return status.stop
end

function Metaheuristics.final_stage!(state, parameters::StagnationWrapper, problem, information, options)
    return Metaheuristics.final_stage!(state, parameters.internal_params, problem, information, options)
end

# configuração do BRKGA com critério de parada personalizado
function StagnatedBRKGA(; num_elites=20, num_mutants=10, num_offsprings=70, bias=0.7, stag_limit=50)
    # instâmcia de um BRKGA regular da biblioteca
    brkga_base = BRKGA(
        num_elites = num_elites, 
        num_mutants = num_mutants, 
        num_offsprings = num_offsprings, 
        bias = bias
    )
    
    # parâmetros envolvidos no wrapper
    wrapped_params = StagnationWrapper(brkga_base.parameters, stag_limit, 0, Inf)
    
    return Metaheuristics.Algorithm(wrapped_params, options = brkga_base.options)
end

# instancia o BRKGA baseline com o limite de estagnação desejado
genetico = StagnatedBRKGA(
    num_elites = 20, 
    num_mutants = 10, 
    num_offsprings = 70, 
    bias = 0.7, 
    stag_limit = 50 # limite de gerações sem melhora
)

# ajuste de opções para desativar outros critérios de parada
genetico.options = Metaheuristics.Options(
    f_calls_limit = typemax(Int), 
    iterations = 10000, 
    store_convergence = true,
    f_tol = -1, # desativa outros critérios de convergência
    x_tol = -1
)

bounds = [zeros(V) ones(V)]'

result = optimize(fitness_harmonious_coloring, bounds, genetico)

# Resultados
println("\n" * "="^20)
@show result
@show Metaheuristics.minimum(result)
@show result.iteration
@show result.termination_status_code