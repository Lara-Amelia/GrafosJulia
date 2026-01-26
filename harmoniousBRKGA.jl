using Metaheuristics
using LinearAlgebra
using Random 

include("colorGul.jl")

# NOTA: seguirá com o mesmo problema dos grafos da base de stanford, que têm arestas repetidas

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

# BRKGA usando parâmetros propostos por Mauricio et. aleatoriamente
# o BRKGA poderá ser usado como baseline ou poderemos analisar o seu desempenho
# para tentar novos operadores (em especial de crossover e mutação) para o ga personalizado
genetico = BRKGA(num_elites = 20, num_mutants = 10, num_offsprings = 70, bias = 0.7)

result = optimize(fitness_harmonious_coloring, [zeros(V) ones(V)], genetico)

# resultados da otimização
@show Metaheuristics.minimum(result)
@show result