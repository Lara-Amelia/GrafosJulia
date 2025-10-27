#CONSERTAR IMPRESSÃO DO ERRO
using Metaheuristics
using LinearAlgebra
# base methods
using Metaheuristics
import Metaheuristics: initialize!, update_state!, final_stage!
import Metaheuristics: AbstractParameters, gen_initial_state, Algorithm, get_position
# genetic operators
import Metaheuristics: SBX_crossover, polynomial_mutation!, create_solution, is_better
import Metaheuristics: reset_to_violated_bounds!

#inclui o arquivo com as funções de coloração gulosa
include("colorGul.jl")

#leitura do grafo
println("Insira o nome do arquivo a ser lido: ")
nome_arquivo = readline()

num_vertices, num_arestas = leInfo!(nome_arquivo)
println("nro de vertices: ", num_vertices)
println("nro de arestas: ", num_arestas)

matriz_adj = zeros(Int, num_vertices, num_vertices)
leArestas!(nome_arquivo, matriz_adj)

print("Insira o número máximo (inteiro) de iterações sem melhora do fitness (k): ")
input_string = readline()
k = parse(Int, input_string)
println("O parametro escolhido foi: ", k)

#coloração harmônica utilizando somente o algoritmo guloso
cores_aresta_distinguivel = fill(-1, num_vertices)
coloracaoHarmonica!(matriz_adj, cores_aresta_distinguivel, num_vertices)
nro_cores_aresta_dist = maximum(cores_aresta_distinguivel)
println("\n--- Coloração Distinguível por Arestas ---")
#=for i in 1:num_vertices
    println("cor usada pelo vertice $i: $(cores_aresta_distinguivel[i])")
end=#
println("O número total de cores usado somente com o guloso foi: $nro_cores_aresta_dist") 

#funções de operadores personalizados
#função para crossover (gerando 2 filhos)
function crossover_harmonious_coloring(parents::Matrix{Float64})
    p1 = parents[1, :]
    p2 = parents[2, :]
    child1 = (p1 .+ p2) ./ 2.0
    #randomização para não gerar 2 filhos exatamente iguais
    child2 = (p1 .+ p2) ./ 2.0 .+ 0.05 .* randn(length(p1))
    return [child1'; child2']
end

#crossover retornando somente 1 filho, média simples dos pais
function crossover_media_simples(parents::Matrix{Float64})
    p1 = parents[1, :]
    p2 = parents[2, :]
    child = (p1 .+ p2) ./2.0
    return [child']
end 

#ou talvez tenha 1 peso para cada (checar, porque o cálculo dos valores sera diferente)
#possivelmente checagem se permanece nos limites do espaço de busca (0 e 1)
function crossover_media_ponderada(parents::Matrix{Float64}, alpha::Float64)
    p1 = parents[1, :]
    p2 = parents[2, :]

    fit1 = fitness_harmonious_coloring(p1)
    fit2 = fitness_harmonious_coloring(p2)

    if alpha >= 0.5
        maior = alpha
        menor = 1.0 - alpha
    else
        maior = 1.0 - alpha
        menor = alpha
    end

    if fit1 >= fit2
        melhor = p2
        pior = p1
    else
        melhor = p1
        pior = p2
    end

    child = melhor .* maior .+ pior .* menor

    return[child']
end 

#swap p - p será passado como parametro e entrará em nosso GAParams
function graph_swap_p(x::Vector{Float64}, p::Int)
    n = length(x)
    if n < 2 || p < 1
        return x
    end

    for k in 1:p
        # escolhe um vértice/posição aleatoriamente
        v1 = rand(1:n)

        # obtem os vizinhos do vértice escolhido
        vizinhos = [i for i in 1:n if matriz_adj[v1, i] == 1]

        # se não há vizinhos, pula para a próxima iteração
        if isempty(vizinhos)
            continue
        end

        # seleciona aleatoriamente algum dos vizinhos e troca
        v2 = rand(vizinhos)
        x[v1], x[v2] = x[v2], x[v1]
    end
end

#função para mutação
function graph_swap_mutation!(x::Vector{Float64})
    n = length(x)
    if n < 2 return x end
    #escolhe um vértice/posição aleatoriamente
    v1 = rand(1:n)
    #obtem os vizinhos do 1° obtido aleatoriamente
    vizinhos = [i for i in 1:n if matriz_adj[v1, i] == 1]
    #se não há vizinhos, não é feito nada (talvez mudar para obter outro vértice)
    if isempty(vizinhos)
        return x
    end
    #seleciona aleatoriamente algum dos vizinhos
    v2 = rand(vizinhos)
    x[v1], x[v2] = x[v2], x[v1]
    return x
end

#definição do Algoritmo Genético personalizado ---
mutable struct CustomGAParams <: Metaheuristics.AbstractParameters
    N::Int
    p_crossover::Float64
    p_mutation::Float64
    stag_limit::Int        
    last_best::Int          
    stag_iters::Int 
end

#adicionamos o limite de iterações estagnadas (sem melhora do fitness)
#falta adicionar os parametros para crossover e para swap com p elementos
#após isso, poderemos escolher qual será usado e integrá-los à lógica do GA
CustomGAParams(; N = 100, p_crossover = 0.5, p_mutation = 0.5, stag_limit = k, 
                 last_best = -1, stag_iters = 0) =
    CustomGAParams(N, p_crossover, p_mutation, stag_limit, last_best, stag_iters)

#inicialização do estado
function Metaheuristics.initialize!(
    status,
    parameters::CustomGAParams,
    problem::Metaheuristics.Problem,
    information::Metaheuristics.Information,
    options::Metaheuristics.Options
)
    return Metaheuristics.gen_initial_state(problem, parameters, information, options, status)
    parameters.last_best = Metaheuristics.best_alternative(population)
    parameters.stagnant_iters = 0
end

function Metaheuristics.update_state!(
    state::Metaheuristics.State,
    parameters::CustomGAParams,
    problem::Metaheuristics.Problem,
    information::Metaheuristics.Information,
    options::Metaheuristics.Options
)
    #implementa torunament selection
    function tournament_select(pop)
        k = 2  # tournament size
        candidates = rand(pop, k)                      #k soluções aleatórias 
        candidate_vectors = [c.x for c in candidates]  #extrai vetores dos objetos solução
        fitnesses = [fitness_harmonious_coloring(v) for v in candidate_vectors]
        return candidates[argmin(fitnesses)].x         #retorna vetor com os melhores
    end

    #seleciona 2 pais na população
    p1 = tournament_select(state.population)
    p2 = tournament_select(state.population)

    parents = vcat(p1', p2')  

    #crossover - pode ser o original ou os de 1 filho, com média simples ou ponderada
    offsprings = crossover_harmonious_coloring(parents)

    #mutação - graph_swap ou graph_swap_p (falta incorporar e testar)
    for i in axes(offsprings, 1)  
        graph_swap_mutation!(offsprings[i, :])
    end

    #criar objetos solução com a offspring produzida
    offspring_solutions = Metaheuristics.xf_solution[]
    for i in axes(offsprings, 1)
        vec = offsprings[i, :]
        fit = Float64(fitness_harmonious_coloring(vec))
        push!(offspring_solutions, Metaheuristics.xf_solution(vec, fit))
    end

    #replace! "não existe"
    #=Metaheuristics.replace!(
        Metaheuristics.ElitistReplacement(),
        state,
        offspring_solutions,  # <-- vector of xf_solution
        options
    )=#

    #mantem o melhor indivíduo da população atual
    best_ind = findmin((s.f for s in state.population))[2]  #índice do melhor
    best_ind = state.population[best_ind]                   #melhor solução

    #troca os demais com a offspring produzida
    state.population = [best_ind; offspring_solutions[2:end]]


    #checa "estagnação"
    melhor_atual = minimum([s.f for s in state.population]) 
    if melhor_atual < parameters.last_best
        parameters.last_best = melhor_atual
        parameters.stag_iters = 0
    else
        parameters.stag_iters += 1
    end

    #condição de Parada
    if parameters.stag_iters >= parameters.stag_limit
        @info "Parada: $(parameters.stag_limit) iterações sem melhora"
        return false  #para otimização
    end

    return true  #continua otimização
end

#define a função final_stage!(ficará vazia por agora)
#"filtragem" dos resultados finais obtidos
function Metaheuristics.final_stage!(
    state,
    parameters::CustomGAParams,
    problem,
    information,
    options
)
    return state
end

#função de fitness
function fitness_harmonious_coloring(individual::Vector{Float64})
    lista_prioridade = sortperm(individual, rev = true)
    cores_vertices = coloracaoHarmonicaGuloso!(matriz_adj, lista_prioridade)
    return maximum(cores_vertices) #talvez seja melhor contar o nro de cores distintas usadas
end

#limites de busca e execução
bounds = [zeros(num_vertices) ones(num_vertices)]'

#cria uma instância do GA
opts = Metaheuristics.Options(f_calls_limit = typemax(Int), iterations = typemax(Int))
params = CustomGAParams(N = 1000, p_crossover = 0.5, p_mutation = 0.5, stag_limit = k, last_best = -1, stag_iters = 0)
my_custom_ga = Metaheuristics.Algorithm(params, options = opts) 

problem = Metaheuristics.Problem(fitness_harmonious_coloring, bounds)

#executa a otimização
result = Metaheuristics.optimize(fitness_harmonious_coloring, bounds, my_custom_ga)

#mostra os resultados da otimização
@show result

melhor_individuo = Metaheuristics.minimizer(result)
lista_prioridade = sortperm(melhor_individuo, rev = true)
cores_vertices = coloracaoHarmonicaGuloso!(matriz_adj, lista_prioridade)

println("\n--- Coloração harmônica final ---")
for i in 1:num_vertices
    println("Vértice $i: cor $(cores_vertices[i])")
end
println("Número total de cores usadas: ", maximum(cores_vertices))