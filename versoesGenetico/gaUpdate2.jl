using Metaheuristics
using LinearAlgebra

#inclui o arquivo com as funções de coloração gulosa
include("../colorGul.jl")

#leitura do grafo
println("Insira o nome do arquivo a ser lido: ")
nome_arquivo = readline()

num_vertices, num_arestas = leInfo!(nome_arquivo)
println("nro de vertices: ", num_vertices)
println("nro de arestas: ", num_arestas)

matriz_adj = zeros(Int, num_vertices, num_vertices)
leArestas!(nome_arquivo, matriz_adj)

#coloração harmônica utilizando somente o algoritmo guloso
cores_aresta_distinguivel = fill(-1, num_vertices)
coloracaoHarmonica!(matriz_adj, cores_aresta_distinguivel, num_vertices)
nro_cores_aresta_dist = maximum(cores_aresta_distinguivel)
println("\n--- Coloração Distinguível por Arestas ---")
println("O número total de cores usado somente com o guloso foi: $nro_cores_aresta_dist") 

#funções de operadores personalizados
#função para crossover
function crossover_harmonious_coloring(parents::Matrix{Float64})
    p1 = parents[1, :]
    p2 = parents[2, :]
    child1 = (p1 .+ p2) ./ 2.0
    child2 = (p1 .+ p2) ./ 2.0 .+ 0.05 .* randn(length(p1))
    return [child1'; child2']
end

#função para mutação
function graph_swap_mutation!(x::Vector{Float64})
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

#definição do Algoritmo Genético personalizado ---
mutable struct CustomGAParams <: Metaheuristics.AbstractParameters
    N::Int
    p_crossover::Float64
    p_mutation::Float64
    max_iter_sem_melhoria::Int
    
    # NOVO CAMPO PARA OPÇÕES
    options::Metaheuristics.Options
end

# NOVO CONSTRUTOR COM O CAMPO DE OPÇÕES
CustomGAParams(; N = 100, p_crossover = 0.5, p_mutation = 0.5, max_iter_sem_melhoria = 20,
                 options = Metaheuristics.Options()) =
    CustomGAParams(N, p_crossover, p_mutation, max_iter_sem_melhoria, options)

#inicialização do estado
function Metaheuristics.initialize!(
    status,
    parameters::CustomGAParams,
    problem::Metaheuristics.Problem,
    information::Metaheuristics.Information,
    options::Metaheuristics.Options
)
    return Metaheuristics.gen_initial_state(problem, parameters, information, options, status)
end

#lógica do algoritmo por iteração (como atualizamos o estado)
function Metaheuristics.update_state!(
    state::Metaheuristics.State,
    parameters::CustomGAParams,
    problem::Metaheuristics.Problem,
    information::Metaheuristics.Information,
    options::Metaheuristics.Options
)
    #seleção de pais
    parents = Metaheuristics.select_parents(
        state.population,
        Metaheuristics.TournamentSelection(),
        options
    )

    #crossover
    offsprings = crossover_harmonious_coloring(parents)

    #mutação
    for i in axes(offsprings, 1)
        graph_swap_mutation!(offsprings[i, :])
    end

    #avaliação
    Metaheuristics.evaluate_population!(offsprings, problem, options)

    #substituição
    Metaheuristics.replace!(
        Metaheuristics.GreedyReplacement(),
        state,
        offsprings,
        options
    )

    # LÓGICA DE PARADA
    if information.iterations_without_improvement >= parameters.max_iter_sem_melhoria
        return false
    end

    return true
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
    return maximum(cores_vertices)
end

#limites de busca e execução
bounds = [zeros(num_vertices) ones(num_vertices)]'

#cria uma instância do GA com o novo parâmetro
# NOVO CÓDIGO A SER USADO
options = Metaheuristics.Options(f_calls_limit=Inf)
params = CustomGAParams(N=100000, p_crossover=0.5, p_mutation=0.5, max_iter_sem_melhoria=20, options=options)

my_custom_ga = Metaheuristics.Algorithm(params) 

problem = Metaheuristics.Problem(fitness_harmonious_coloring, bounds)

# Não precisa mais passar 'options' separadamente, pois está dentro de 'params'
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

#=CONSERTANDO METHOD ERROR EM OPTIMIZE
using Metaheuristics
using LinearAlgebra

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

#coloração harmônica utilizando somente o algoritmo guloso
cores_aresta_distinguivel = fill(-1, num_vertices)
coloracaoHarmonica!(matriz_adj, cores_aresta_distinguivel, num_vertices)
nro_cores_aresta_dist = maximum(cores_aresta_distinguivel)
println("\n--- Coloração Distinguível por Arestas ---")
println("O número total de cores usado somente com o guloso foi: $nro_cores_aresta_dist") 

#funções de operadores personalizados
#função para crossover
function crossover_harmonious_coloring(parents::Matrix{Float64})
    p1 = parents[1, :]
    p2 = parents[2, :]
    child1 = (p1 .+ p2) ./ 2.0
    child2 = (p1 .+ p2) ./ 2.0 .+ 0.05 .* randn(length(p1))
    return [child1'; child2']
end

#função para mutação
function graph_swap_mutation!(x::Vector{Float64})
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

#definição do Algoritmo Genético personalizado ---
mutable struct CustomGAParams <: Metaheuristics.AbstractParameters
    N::Int
    p_crossover::Float64
    p_mutation::Float64
    max_iter_sem_melhoria::Int
end

CustomGAParams(; N = 100, p_crossover = 0.5, p_mutation = 0.5, max_iter_sem_melhoria = 20) =
    CustomGAParams(N, p_crossover, p_mutation, max_iter_sem_melhoria)

#inicialização do estado
function Metaheuristics.initialize!(
    status,
    parameters::CustomGAParams,
    problem::Metaheuristics.Problem,
    information::Metaheuristics.Information,
    options::Metaheuristics.Options
)
    return Metaheuristics.gen_initial_state(problem, parameters, information, options, status)
end

#lógica do algoritmo por iteração (como atualizamos o estado)
function Metaheuristics.update_state!(
    state::Metaheuristics.State,
    parameters::CustomGAParams,
    problem::Metaheuristics.Problem,
    information::Metaheuristics.Information,
    options::Metaheuristics.Options
)
    #seleção de pais
    parents = Metaheuristics.select_parents(
        state.population,
        Metaheuristics.TournamentSelection(),
        options
    )

    #crossover
    offsprings = crossover_harmonious_coloring(parents)

    #mutação
    for i in axes(offsprings, 1)
        graph_swap_mutation!(offsprings[i, :])
    end

    #avaliação
    Metaheuristics.evaluate_population!(offsprings, problem, options)

    #substituição
    Metaheuristics.replace!(
        Metaheuristics.GreedyReplacement(),
        state,
        offsprings,
        options
    )

    # NOVO CÓDIGO DE PARADA
    if information.iterations_without_improvement >= parameters.max_iter_sem_melhoria
        return false # Sinaliza que o algoritmo deve parar
    end

    return true
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
    return maximum(cores_vertices)
end

#limites de busca e execução
bounds = [zeros(num_vertices) ones(num_vertices)]'

#cria uma instância do GA
params = CustomGAParams(N = 100, p_crossover = 0.5, p_mutation = 0.5)
my_custom_ga = Metaheuristics.Algorithm(params) 

problem = Metaheuristics.Problem(fitness_harmonious_coloring, bounds)

# NOVO CÓDIGO
# Define as opções de execução, desabilitando o limite de iterações padrão
options = Metaheuristics.Options(f_calls_limit=Inf)

# Agora, execute a otimização, passando o objeto 'options'
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
println("Número total de cores usadas: ", maximum(cores_vertices))=#

#=using Metaheuristics
using LinearAlgebra

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
#função para crossover
function crossover_harmonious_coloring(parents::Matrix{Float64})
    p1 = parents[1, :]
    p2 = parents[2, :]
    child1 = (p1 .+ p2) ./ 2.0
    #randomização para não gerar 2 filhos exatamente iguais
    child2 = (p1 .+ p2) ./ 2.0 .+ 0.05 .* randn(length(p1))
    return [child1'; child2']
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
end

CustomGAParams(; N = 100, p_crossover = 0.5, p_mutation = 0.5) =
    CustomGAParams(N, p_crossover, p_mutation)

#inicialização do estado
function Metaheuristics.initialize!(
    status,
    parameters::CustomGAParams,
    problem::Metaheuristics.Problem,
    information::Metaheuristics.Information,
    options::Metaheuristics.Options
)
    return Metaheuristics.gen_initial_state(problem, parameters, information, options, status)
    #Metaheuristics.evaluate_population!(population, problem, options)
    #best_solution = Metaheuristics.best(population)
    #return Metaheuristics.State(best_solution, population)
end

#lógica do algoritmo por iteração (como atualizamos o estado)
function Metaheuristics.update_state!(
    state::Metaheuristics.State,
    parameters::CustomGAParams,
    problem::Metaheuristics.Problem,
    information::Metaheuristics.Information,
    options::Metaheuristics.Options
)
    #seleção de pais
    parents = Metaheuristics.select_parents(
        state.population,
        Metaheuristics.TournamentSelection(),
        options
    )

    #crossover
    offsprings = crossover_harmonious_coloring(parents)

    #mutação
    for i in axes(offsprings, 1)
        graph_swap_mutation!(offsprings[i, :])
    end

    #avaliação
    Metaheuristics.evaluate_population!(offsprings, problem, options)

    #substituição
    Metaheuristics.replace!(
        Metaheuristics.GreedyReplacement(),
        state,
        offsprings,
        options
    )

    return true
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
params = CustomGAParams(N = 100, p_crossover = 0.5, p_mutation = 0.5)
my_custom_ga = Metaheuristics.Algorithm(params) 

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
println("Número total de cores usadas: ", maximum(cores_vertices))=#