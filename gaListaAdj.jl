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
    @info "Population fitnesses: " [s.f for s in pop]

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

# parâmetros para geração dos indivíduos - são permutações de valores reais entre 0 e 1
# representando listas de prioridade na coloração (tentamos encontrar prioridades melhores)
bounds = [zeros(V) ones(V)]'

params = CustomGAParams(N=100, p_mutation=0.5, stag_limit=50, last_best=Inf)

opt_settings = Metaheuristics.Options(f_calls_limit = typemax(Int), iterations = 1000, store_convergence = true)

my_ga = Metaheuristics.Algorithm(params, options = opt_settings)
problem = Metaheuristics.Problem(fitness_harmonious_coloring, bounds)

result = Metaheuristics.optimize(fitness_harmonious_coloring, bounds, my_ga)

# resultados da otimização
@show Metaheuristics.minimum(result)
@show result

#=melhor_x = Metaheuristics.minimizer(result)
println("Melhor número de cores: ",
        fitness_harmonious_coloring(melhor_x))=#

# NOTA: tentativas de melhoria para o desempenho do GA em termos de qualidade da coloração obtida - 
# inserir ordens gulosas na população inicial, alterar operadores (especialmente o crossover, 
# que parece ser fraco/gerar pouca diversidade na população - nas últimas gerações do GA, temos indivíduos 
# com basicamente a exata mesma ordem, ou seja, podemos extar explorando pouco o espaço de busca), 
# implementação de uma busca local (GRASP) que poderia ser utilizada como operador de mutação 
# (apesar de que a mutação atual parece ser boa para a diversificação da população)

# NOTA: verificar questão com grafos simples densos que sempre caem no pior caso (possivelmente pq
# a população inicial em si sempre começa com todos os elementos no pior caso, e o algoritmo não consegue
# realmente evoluir a partir daí - nesses casos os 3 gulosos também resultam sempre no pior caso)