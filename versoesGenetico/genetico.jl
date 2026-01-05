using Metaheuristics
using LinearAlgebra

# Inclui o arquivo com as funções de coloração gulosa
include("../colorGul.jl")

# --- Leitura do grafo ---
println("Insira o nome do arquivo a ser lido: ")
nome_arquivo = readline()

num_vertices, num_arestas = leInfo!(nome_arquivo)
println("nro de vertices: $num_vertices")
println("nro de arestas: $num_arestas")

matriz_adj = zeros(Int, num_vertices, num_vertices)
leArestas!(nome_arquivo, matriz_adj)

#Coloração harmônica utilizando somente o algoritmo guloso
cores_aresta_distinguivel = fill(-1, num_vertices)
coloracaoHarmonica!(matriz_adj, cores_aresta_distinguivel, num_vertices)
nro_cores_aresta_dist = maximum(cores_aresta_distinguivel)
println("\n--- Coloração Distinguível por Arestas ---")
#=for i in 1:num_vertices
    println("cor usada pelo vertice $i: $(cores_aresta_distinguivel[i])")
end=#
println("O número total de cores usado somente com o guloso foi: $nro_cores_aresta_dist") 

# --- Função de fitness ---
function fitness_harmonious_coloring(individual::Vector{Float64})
    lista_prioridade = sortperm(individual, rev=true)
    cores_vertices = coloracaoHarmonicaPrioridade(matriz_adj, lista_prioridade)
    return maximum(cores_vertices)
end

# --- Operador de crossover ---
function crossover_harmonious_coloring(parents::Matrix{Float64})
    p1 = parents[1, :]
    p2 = parents[2, :]
    child1 = (p1 .+ p2) ./ 2.0
    child2 = (p1 .+ p2) ./ 2.0 .+ 0.05 .* randn(length(p1))
    return [child1'; child2']
end

# --- Operador de mutação ---
function graph_swap_mutation(x::Vector{Float64})
    n = length(x)
    if n < 2 return x end
    v1 = rand(1:n)
    vizinhos = [i for i in 1:n if matriz_adj[v1, i] == 1]
    if isempty(vizinhos)
        return x
    end
    v2 = rand(vizinhos)
    x[v1], x[v2] = x[v2], x[v1]
    return x
end

# --- Limites de busca ---
bounds = [zeros(num_vertices) ones(num_vertices)]'

# --- GA ---
#=ga = GA(; N=100, p_mutation=0.5, p_crossover = 0.5, initializer = RandomInBounds(),
          selection = TournamentSelection(), crossover = crossover_harmonious_coloring(), 
          mutation = graph_swap_mutation(), environmental_selection = ElitistReplacement()
       )=#

# --- Execução ---
result = optimize(fitness_harmonious_coloring, bounds, ga)

@show result

melhor_individuo = Metaheuristics.minimizer(result)
lista_prioridade = sortperm(melhor_individuo, rev=true)
cores_vertices = coloracaoHarmonicaPrioridade!(matriz_adj, lista_prioridade)

println("\n--- Coloração harmônica final ---")
for i in 1:num_vertices
    println("Vértice $i: cor $(cores_vertices[i])")
end
total_cores = maximum(cores_vertices)
println("Número total de cores usadas: $total_cores")

#=using Metaheuristics
using LinearAlgebra

include("colorGul.jl")

# --- Leitura do grafo ---
println("Insira o nome do arquivo a ser lido: ")
nome_arquivo = readline()

num_vertices, num_arestas = leInfo!(nome_arquivo)
println("nro de vertices: ", num_vertices)
println("nro de arestas: ", num_arestas)

matriz_adj = zeros(Int, num_vertices, num_vertices)
leArestas!(nome_arquivo, matriz_adj)

# --- Função de fitness ---
function fitness_harmonious_coloring(individual::Vector{Float64})
    lista_prioridade = sortperm(individual, rev=true)
    cores_vertices = coloracaoHarmonicaPrioridade!(matriz_adj, lista_prioridade)
    return maximum(cores_vertices)
end

# --- Operador de crossover ---
struct CrossoverHarmonious <: AbstractRecombination end

function Metaheuristics.recombine!(::CrossoverHarmonious, parents::Matrix{Float64})
    p1 = parents[1, :]
    p2 = parents[2, :]
    child1 = (p1 .+ p2) ./ 2.0
    child2 = (p1 .+ p2) ./ 2.0 .+ 0.05 .* randn(length(p1))
    return [child1'; child2']
end

# --- Operador de mutação ---
struct GraphSwapMutation <: AbstractMutation end

function Metaheuristics.mutation!(::GraphSwapMutation, x::Vector{Float64})
    n = length(x)
    if n < 2 return x end
    v1 = rand(1:n)
    vizinhos = [i for i in 1:n if matriz_adj[v1, i] == 1]
    if isempty(vizinhos)
        return x
    end
    v2 = rand(vizinhos)
    x[v1], x[v2] = x[v2], x[v1]
    return x
end

# --- Limites de busca ---
bounds = [zeros(num_vertices) ones(num_vertices)]'

# --- GA ---
ga = GA(
    populationSize = 100,
    recombination = CrossoverHarmonious(),
    mutation = GraphSwapMutation()
)

# --- Execução ---
result = optimize(fitness_harmonious_coloring, bounds, ga)

@show result

melhor_individuo = result.best_solution
lista_prioridade = sortperm(melhor_individuo, rev=true)
cores_vertices = coloracaoHarmonicaPrioridade!(matriz_adj, lista_prioridade)

println("\nColoração harmônica final:")
for i in 1:num_vertices
    println("Vértice $i: cor $(cores_vertices[i])")
end
println("Número total de cores usadas: ", maximum(cores_vertices))=#

#=using Metaheuristics
using LinearAlgebra

# Inclui o arquivo com as funções de coloração gulosa
include("colorGul.jl")

# O processamento/leitura dos dados do grafo para a matriz de adjacência é feito aqui
println("Insira o nome do arquivo a ser lido: ")
nome_arquivo = readline()

num_vertices, num_arestas = leInfo!(nome_arquivo)
println("nro de vertices: ", num_vertices)
println("nro de arestas: ", num_arestas)

# Declara a matriz de adjacência como uma constante global
matriz_adj = zeros(Int, num_vertices, num_arestas)
leArestas!(nome_arquivo, matriz_adj)

# Função de fitness (avalia a qualidade da solução)
function fitness_harmonious_coloring(individual::Vector{Float64})
    lista_prioridade = sortperm(individual, rev=true)
    cores_vertices = coloracaoHarmonicaPrioridade(matriz_adj, lista_prioridade)
    num_cores = maximum(cores_vertices)
    return num_cores
end

# Define os limites de busca para os valores de prioridade (entre 0 e 1)
lower_bounds = [0.0 for _ in 1:num_vertices]
upper_bounds = [1.0 for _ in 1:num_vertices]
bounds = [lower_bounds upper_bounds]'

# Cria o otimizador GA com os tipos dos seus operadores personalizados
ga_optimizer = GA(
    populationSize=100,
    recombination=UniformCrossover(),
    mutation=InsertRandomMutation()
)

# Executa a otimização
result = optimize(fitness_harmonious_coloring, bounds, ga_optimizer)

# Mostra o resultado final
@show result=#


#=using Metaheuristics
using LinearAlgebra

# Inclui o arquivo com as funções de coloração gulosa
include("colorGul.jl")

# --- Leitura do grafo ---
println("Insira o nome do arquivo a ser lido: ")
nome_arquivo = readline()

num_vertices, num_arestas = leInfo!(nome_arquivo)
println("nro de vertices: ", num_vertices)
println("nro de arestas: ", num_arestas)

# Cria a matriz de adjacência
matriz_adj = zeros(Int, num_vertices, num_vertices)
leArestas!(nome_arquivo, matriz_adj)

# --- Função de fitness ---
function fitness_harmonious_coloring(individual::Vector{Float64})
    # Converte o vetor de prioridade do GA em uma ordem de vértices
    lista_prioridade = sortperm(individual, rev=true)
    
    # Chama o algoritmo guloso adaptado
    cores_vertices = coloracaoHarmonicaPrioridade!(matriz_adj, lista_prioridade)
    
    # Fitness = número máximo de cores usadas
    return maximum(cores_vertices)
end

# --- Função de cruzamento ---
function crossover_harmonious(parents::Matrix{Float64})
    parent1 = parents[1, :]
    parent2 = parents[2, :]
    
    # Dois filhos diferentes: média + média + ruído pequeno
    child1 = (parent1 .+ parent2) ./ 2
    child2 = (parent1 .+ parent2) ./ 2 .+ 0.01 .* rand(length(parent1))
    
    return [child1'; child2']
end

# --- Função de mutação ---
function mutate_harmonious(x::Vector{Float64})
    num_vertices = length(x)
    if num_vertices < 2
        return x
    end
    
    # Seleciona aleatoriamente um vértice
    v1_idx = rand(1:num_vertices)
    
    # Encontra vizinhos
    neighbor_indices = [i for i in 1:num_vertices if matriz_adj[v1_idx, i] == 1]
    
    if isempty(neighbor_indices)
        return x
    end
    
    # Seleciona aleatoriamente um vizinho
    v2_idx = rand(neighbor_indices)
    
    # Troca prioridades
    x[v1_idx], x[v2_idx] = x[v2_idx], x[v1_idx]
    
    return x
end

# --- Espaço de busca ---
lower_bounds = zeros(num_vertices)
upper_bounds = ones(num_vertices)
bounds = [lower_bounds upper_bounds]'

# --- Cria o GeneticAlgorithm com operadores customizados ---
ga_optimizer = GeneticAlgorithm(
    selection = TournamentSelection(),
    recombination = crossover_harmonious,
    mutation = mutate_harmonious,
    replacement = ElitistReplacement(),
    populationSize = 100
)

# --- Executa a otimização ---
result = optimize(fitness_harmonious_coloring, bounds, ga_optimizer)

# --- Mostra o resultado ---
println("Melhor fitness encontrado (menor número de cores): ", result.minimum)
println("Melhor vetor de prioridade: ", result.minimizer)=#


#=using Metaheuristics
using LinearAlgebra

# Inclui o arquivo com as funções de coloração gulosa
include("colorGul.jl")

# --- Leitura do grafo ---
println("Insira o nome do arquivo a ser lido: ")
nome_arquivo = readline()

num_vertices, num_arestas = leInfo!(nome_arquivo)
println("nro de vertices: ", num_vertices)
println("nro de arestas: ", num_arestas)

# Cria a matriz de adjacência
matriz_adj = zeros(Int, num_vertices, num_vertices)
leArestas!(nome_arquivo, matriz_adj)

# --- Função de fitness ---
function fitness_harmonious_coloring(individual::Vector{Float64})
    # Converte o vetor de prioridade do GA em uma ordem de vértices
    lista_prioridade = sortperm(individual, rev=true)
    
    # Chama o algoritmo guloso adaptado
    cores_vertices = coloracaoHarmonicaPrioridade!(matriz_adj, lista_prioridade)
    
    # Fitness = número máximo de cores usadas
    return maximum(cores_vertices)
end

# --- Definição do operador de cruzamento customizado ---
struct MyCrossover <: AbstractRecombination end

function Metaheuristics.recombine!(::MyCrossover, parents::Matrix{Float64})
    parent1 = parents[1, :]
    parent2 = parents[2, :]
    
    # Dois filhos diferentes: média + média + ruído pequeno
    child1 = (parent1 .+ parent2) ./ 2
    child2 = (parent1 .+ parent2) ./ 2 .+ 0.01 .* rand(length(parent1))
    
    return [child1'; child2']
end

# --- Definição do operador de mutação customizado ---
struct MyMutation <: AbstractMutation end

function Metaheuristics.mutate!(::MyMutation, x::Vector{Float64})
    num_vertices = length(x)
    if num_vertices < 2
        return x
    end
    
    # Seleciona aleatoriamente um vértice
    v1_idx = rand(1:num_vertices)
    
    # Encontra vizinhos
    neighbor_indices = [i for i in 1:num_vertices if matriz_adj[v1_idx, i] == 1]
    
    if isempty(neighbor_indices)
        return x
    end
    
    # Seleciona aleatoriamente um vizinho
    v2_idx = rand(neighbor_indices)
    
    # Troca prioridades
    x[v1_idx], x[v2_idx] = x[v2_idx], x[v1_idx]
    
    return x
end

# --- Espaço de busca ---
lower_bounds = zeros(num_vertices)
upper_bounds = ones(num_vertices)
bounds = [lower_bounds upper_bounds]'

# --- Cria o GeneticAlgorithm com operadores customizados ---
ga_optimizer = GeneticAlgorithm(
    selection = TournamentSelection(),
    recombination = MyCrossover(),
    mutation = MyMutation(),
    replacement = ElitistReplacement(),
    populationSize = 100
)

# --- Executa a otimização ---
result = optimize(fitness_harmonious_coloring, bounds, ga_optimizer)

# --- Mostra o resultado ---
println("Melhor fitness encontrado (menor número de cores): ", result.minimum)
println("Melhor vetor de prioridade: ", result.minimizer)=#


#=using Metaheuristics
using LinearAlgebra

# Inclui o arquivo com as funções de coloração gulosa
include("colorGul.jl")

# --- Leitura do grafo ---
println("Insira o nome do arquivo a ser lido: ")
nome_arquivo = readline()

num_vertices, num_arestas = leInfo!(nome_arquivo)
println("nro de vertices: ", num_vertices)
println("nro de arestas: ", num_arestas)

# Matriz de adjacência
matriz_adj = zeros(Int, num_vertices, num_vertices)
leArestas!(nome_arquivo, matriz_adj)

# --- Função de fitness ---
function fitness_harmonious_coloring(individual::Vector{Float64})
    lista_prioridade = sortperm(individual, rev=true)
    cores_vertices = coloracaoHarmonicaPrioridade!(matriz_adj, lista_prioridade)
    return maximum(cores_vertices)
end

# --- Operador de crossover ---
struct CrossoverHarmonious end

function Metaheuristics.recombine!(::CrossoverHarmonious, parents::Matrix{Float64})
    parent1 = parents[1, :]
    parent2 = parents[2, :]

    # cria dois filhos diferentes aplicando pequenas perturbações
    child1 = (parent1 .+ parent2) ./ 2.0
    child2 = (parent1 .+ parent2) ./ 2.0 .+ 0.01 .* randn(length(parent1))
    
    return [child1'; child2']
end

# --- Operador de mutação ---
struct MutationHarmonious end

function Metaheuristics.mutation!(::MutationHarmonious, x::Vector{Float64})
    n = length(x)
    if n < 2
        return x
    end

    # seleciona um vértice aleatório
    v1 = rand(1:n)
    neighbors = [i for i in 1:n if matriz_adj[v1, i] == 1]

    if isempty(neighbors)
        return x
    end

    v2 = rand(neighbors)
    x[v1], x[v2] = x[v2], x[v1]
    return x
end

# --- Limites de busca ---
lower_bounds = zeros(num_vertices)
upper_bounds = ones(num_vertices)
bounds = [lower_bounds upper_bounds]'

# --- Otimizador GA ---
ga_optimizer = GeneticAlgorithm(
    populationSize=100,
    recombination=CrossoverHarmonious(),
    mutation=MutationHarmonious()
)

# --- Execução ---
result = optimize(fitness_harmonious_coloring, bounds, ga_optimizer)

@show result=#


#=using Metaheuristics
using LinearAlgebra # Necessário para o operador '

# Inclui o arquivo com as funções de coloração gulosa
include("colorGul.jl")

# O processamento/leitura dos dados do grafo para a matriz de adjacência é feito aqui
println("Insira o nome do arquivo a ser lido: ")
nome_arquivo = readline()

num_vertices, num_arestas = leInfo!(nome_arquivo)
println("nro de vertices: ", num_vertices)
println("nro de arestas: ", num_arestas)

# Declara a matriz de adjacência como uma constante global
#const matriz_adj = zeros(Int, num_vertices, num_vertices)
matriz_adj = zeros(Int, num_vertices, num_arestas)
leArestas!(nome_arquivo, matriz_adj)

# --- Definições do Algoritmo Genético ---

# Função de fitness (avalia a qualidade da solução)
function fitness_harmonious_coloring(individual::Vector{Float64})
    # Converte o vetor de prioridade do GA em uma ordem de vértices
    lista_prioridade = sortperm(individual, rev=true)
    
    # Chama o algoritmo guloso adaptado com a lista de prioridade de inteiros
    cores_vertices = coloracaoHarmonicaPrioridade(matriz_adj, lista_prioridade)
    
    # O fitness é o número máximo de cores usadas, que o GA vai tentar minimizar
    num_cores = maximum(cores_vertices)
    
    return num_cores
end

# Define um tipo para o seu operador de cruzamento
struct CrossoverHarmonious end

# Diz à biblioteca Metaheuristics como recombinar indivíduos com o seu tipo
function Metaheuristics.recombine!(::CrossoverHarmonious, parents::Matrix{Float64})
    parent1 = parents[1, :]
    parent2 = parents[2, :]
    
    child1 = (parent1 .+ parent2) ./ 2.0
    child2 = (parent1 .+ parent2) ./ 2.0
    
    return [child1'; child2']
end

# Define um tipo para o seu operador de mutação
struct MutationHarmonious end

# Diz à biblioteca Metaheuristics como mutar indivíduos com o seu tipo
function Metaheuristics.mutation!(::MutationHarmonious, x::Vector{Float64})
    num_vertices = length(x)
    if num_vertices < 2
        return x
    end
    
    # Seleciona aleatoriamente um vértice (posição)
    v1_idx = rand(1:num_vertices)
    
    # Encontra os vizinhos do vértice na matriz de adjacência global
    neighbor_indices = Int[]
    for i in 1:num_vertices
        if matriz_adj[v1_idx, i] == 1
            push!(neighbor_indices, i)
        end
    end
    
    if isempty(neighbor_indices)
        return x
    end
    
    v2_idx = rand(neighbor_indices)
    
    # Faz a troca dos valores de prioridade
    x[v1_idx], x[v2_idx] = x[v2_idx], x[v1_idx]
    
    return x
end

# --- Execução do Algoritmo Genético ---

# Define os limites de busca para os valores de prioridade (entre 0 e 1)
lower_bounds = [0.0 for _ in 1:num_vertices]
upper_bounds = [1.0 for _ in 1:num_vertices]
bounds = [lower_bounds upper_bounds]'

# Cria o otimizador GA com os tipos dos seus operadores personalizados
ga_optimizer = GA(
    populationSize=100,
    recombination=CrossoverHarmonious(),
    mutation=MutationHarmonious()
)

# Executa a otimização
result = optimize(fitness_harmonious_coloring, bounds, ga_optimizer)

# Mostra o resultado final
@show result=#

#=using Metaheuristics
include("colorGul.jl")

#o processamento/leitura dos dados do grafo para a matriz de adjcaencia pode ser feito aqui,
#uma vez que também teremos de usar a matriz para aplicar o guloso ao longo do alg genético
println("Insira o nome do arquivo a ser lido: ")
nome_arquivo = readline()

num_vertices, num_arestas = leInfo!(nome_arquivo)
println("nro de vertices: $num_vertices")
println("nro de arestas: $num_arestas")

matriz_adj = zeros(Int, num_vertices, num_vertices)
leArestas!(nome_arquivo, matriz_adj)

#função que obtem o fitness de uma dada solução/indivíduo
#a "função" de fitness a ser minimizadas consistirá no número de cores
#utilizado na coloração harmônica gerada por lagoritmo guloso
function fitness_harmonious_coloring(individual::Vector{Float64})
    #ordena por vértices de acordo com a ordem de prioridade definida na lista
    lista_prioridade = sortperm(individual, rev=true)
    num_cores = 0
    #cores_vertices = fill(-1, num_vertices)
    #método para coloração harmônica deve ser adaptado para produzir a coloração
    #segundo a ordem de prioridade estabelecida pelo indivíduo
    cores_vertices = coloracaoHarmonicaPrioridade!(matriz_adj, lista_prioridade)
    num_cores = length(unique(cores_vertices))
    #o número de cores utilizado pelo guloso será nosso fitness
    return num_cores
end

function crossover_harmonious_coloring(parents::Matrix{Float64})
    if size(parents, 1) != 2
        error("São necessários 2 parents")
    end

    parent1 = parents[1, :]
    parent2 = parents[2, :]
    
    #gera 2 filhos idênticos - possivelmente alterar
    child1 = (parent1 .+ parent2) ./ 2.0
    child2 = (parent1 .+ parent2) ./ 2.0
    
    return [child1'; child2']
end

function graph_swap_mutation(x::Vector{Float64})
    num_vertices = length(x)
    if num_vertices < 2
        return x
    end
    
    #seleciona aleatoriamente algum vértice
    v1_idx = rand(1:num_vertices)
    
    #encontra os vizinhos do vértice selcionado
    neighbor_indices = Int[]
    for i in 1:num_vertices
        #encontra registra vizinhos se for necessário
        if graph_matrix[v1_idx, i] == 1
            push!(neighbor_indices, i)
        end
    end
    
    #checa se há algum vizinho at all
    if isempty(neighbor_indices)
        #se não há vizinhos, ou não fazemos o swap ou retornamos ao primeiro passo
        #e selecinamos outro vértice para que ocorra a troca
        return x
    end
    
    #seleciona aleatoriamente algum dos vizinhos para fazer o swap
    v2_idx = rand(neighbor_indices)
    
    #faz o swap dos valores de prioridade
    x[v1_idx], x[v2_idx] = x[v2_idx], x[v1_idx]
    #retorna o vetor com as alterações realizadas
    return x
end

#define os limites de busca para os valores de prioridade (entre 0 e 1)
lower_bounds = [0.0 for _ in 1:num_vertices]
upper_bounds = [1.0 for _ in 1:num_vertices]
bounds = [lower_bounds upper_bounds]'

ga_optimizer = GA(
    populationSize=100,
    recombination=crossover_harmonious_coloring,
    mutation=graph_swap_mutation
)

result = optimize(fitness_harmonious_coloring, bounds, ga_optimizer)

@show result =#