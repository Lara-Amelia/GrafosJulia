using Metaheuristics
using LinearAlgebra
using Random 

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

function fitness_harmonious_coloring(individual::Vector{Float64})
    lista_prioridade = sortperm(individual, rev = true)

    cores_vertices = coloracaoHarmonicaAdjVetAux!(ADJ, lista_prioridade)
    fitness = Float64(maximum(cores_vertices))
    return fitness
end

# struct para que possamos aplicar procedimentos da biblioteca diretamente
struct GraphSwapMutation
    p::Float64 # probabilidade de mutação
end

# implementação da lógica no wrapper da mutação
function Metaheuristics.mutation!(Q::AbstractMatrix{Float64}, mut::GraphSwapMutation)
    n_individuals, n_genes = size(Q)
    
    # "pré-seleção" dos filhos que serão mutados usando a probabilidade do struct
    to_mutate = findall(rand(n_individuals) .< mut.p)

    for i in to_mutate
        v1 = rand(1:n_genes)
        vizinhos = ADJ[v1]
        
        if isempty(vizinhos)
            continue 
        end

        v2 = rand(vizinhos)

        @inbounds begin
            tmp = Q[i, v1]
            Q[i, v1] = Q[i, v2]
            Q[i, v2] = tmp
        end 
    end
end

# setup do BRKGA com operador de mutação personalizado
function MyCustomBRKGA(;
        num_elites = 20,
        num_mutants = 10,
        num_offsprings = 70,
        N = num_elites + num_mutants + num_offsprings,
        bias = 0.7,
        kargs...
    )

    # configuração manual dos componentes - o único que muda é a mutação em si
    initializer = Metaheuristics.RandomInBounds(; N)
    selection   = Metaheuristics.BiasedSelection(num_elites, num_offsprings)
    crossover   = Metaheuristics.BinomialCrossover(p = bias, n_offsprings = 1)
    
    # uso da mutação customizada
    mutation    = GraphSwapMutation(0.5) 
    
    environmental_selection = Metaheuristics.ElitistReplacement()

    # retorn o objeto GA modificado (exatamente como no src da biblioteca)
    return Metaheuristics.GA(;
        initializer,
        selection,
        crossover,
        mutation,
        environmental_selection,
        kargs...
      )
end

# instância do BRKGA customizado
genetico = MyCustomBRKGA(num_elites = 20, num_mutants = 10, num_offsprings = 70, bias = 0.7)

bounds = [zeros(V) ones(V)]'

result = optimize(fitness_harmonious_coloring, bounds, genetico)

# resultados
println("\n" * "="^20)
@show Metaheuristics.minimum(result)
@show result