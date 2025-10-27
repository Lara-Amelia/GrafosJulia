
#= TENTAR DA FORMA MAIS SIMPLES AGORA
using Metaheuristics
using LinearAlgebra

# Inclui o arquivo com as funções de coloração gulosa
include("colorGul.jl")

# --- Leitura do grafo ---
println("Insira o nome do arquivo a ser lido: ")
nome_arquivo = readline()

num_vertices, num_arestas = leInfo!(nome_arquivo)
println("nro de vertices: ", num_vertices)
println("nro de arestas: ", num_arestas)

matriz_adj = zeros(Int, num_vertices, num_vertices)
leArestas!(nome_arquivo, matriz_adj)

# --- Funções de operadores personalizados ---
function crossover_harmonious_coloring(parents::Matrix{Float64})
    p1 = parents[1, :]
    p2 = parents[2, :]
    child1 = (p1 .+ p2) ./ 2.0
    child2 = (p1 .+ p2) ./ 2.0 .+ 0.05 .* randn(length(p1))
    return [child1'; child2']
end

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

# --- Algoritmo genético personalizado ---
mutable struct CustomGAParams <: Metaheuristics.AbstractParameters
    N::Int                 # population size
    p_crossover::Float64
    p_mutation::Float64
end

# Default constructor
CustomGAParams() = CustomGAParams(100, 0.5, 0.5)

# Inicialização do estado
function Metaheuristics.initialize!(
    status,
    parameters::CustomGAParams,
    problem::Metaheuristics.Problem,
    information::Metaheuristics.Information,
    options::Metaheuristics.Options
)
    population = Metaheuristics.gen_initial_state(problem, parameters.N, information, options, status)
    # Framework now evaluates automatically, no evaluate_population! needed
    best_solution = population[1]  # placeholder; best is tracked internally
    return Metaheuristics.State(best_solution, population)
end

# Lógica do algoritmo por iteração
function Metaheuristics.update_state!(
    state::Metaheuristics.State,
    parameters::CustomGAParams,
    problem::Metaheuristics.Problem,
    information::Metaheuristics.Information,
    options::Metaheuristics.Options
)
    # 1. Seleção de pais
    parents = Metaheuristics.select_parents(
        state.population,
        Metaheuristics.TournamentSelection(),
        options
    )

    # 2. Crossover
    offsprings = crossover_harmonious_coloring(parents)

    # 3. Mutação
    for child in eachrow(offsprings)
        graph_swap_mutation!(child)
    end

    # 4. Substituir no estado
    Metaheuristics.replace!(
        Metaheuristics.GreedyReplacement(),
        state,
        offsprings,
        options
    )

    return true
end

# --- Função de Fitness ---
function fitness_harmonious_coloring(individual::Vector{Float64})
    lista_prioridade = sortperm(individual, rev = true)
    cores_vertices = coloracaoHarmonicaPrioridade!(matriz_adj, lista_prioridade)
    return maximum(cores_vertices)
end

# --- Limites de busca ---
bounds = [zeros(num_vertices) ones(num_vertices)]'

# --- Instancia parâmetros do GA ---
params = CustomGAParams(100, 0.5, 0.5)

# --- Cria o algoritmo genético ---
my_custom_ga = Metaheuristics.Algorithm(params)

# --- Cria o problema ---
problem = Metaheuristics.Problem(fitness_harmonious_coloring, bounds)

# --- Executa otimização ---
result = Metaheuristics.optimize(fitness_harmonious_coloring, bounds, my_custom_ga)

@show result

# --- Extrai melhor indivíduo e cores ---
melhor_individuo = Metaheuristics.minimizer(result)
lista_prioridade = sortperm(melhor_individuo, rev = true)
cores_vertices = coloracaoHarmonicaPrioridade!(matriz_adj, lista_prioridade)

println("\n--- Coloração harmônica final ---")
for i in 1:num_vertices
    println("Vértice $i: cor $(cores_vertices[i])")
end
println("Número total de cores usadas: ", maximum(cores_vertices))=#


#=ERRO COM evaluate_population
using Metaheuristics
using LinearAlgebra

# Inclui o arquivo com as funções de coloração gulosa
include("colorGul.jl")

# --- Leitura do grafo ---
println("Insira o nome do arquivo a ser lido: ")
nome_arquivo = readline()

num_vertices, num_arestas = leInfo!(nome_arquivo)
println("nro de vertices: ", num_vertices)
println("nro de arestas: ", num_arestas)

matriz_adj = zeros(Int, num_vertices, num_vertices)
leArestas!(nome_arquivo, matriz_adj)

# --- Funções de operadores personalizados ---
function crossover_harmonious_coloring(parents::Matrix{Float64})
    p1 = parents[1, :]
    p2 = parents[2, :]
    child1 = (p1 .+ p2) ./ 2.0
    child2 = (p1 .+ p2) ./ 2.0 .+ 0.05 .* randn(length(p1))
    return [child1'; child2']
end

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

# --- Algoritmo genético personalizado ---
mutable struct CustomGAParams <: Metaheuristics.AbstractParameters
    N::Int                 # population size
    p_crossover::Float64
    p_mutation::Float64
end

# Default constructor for Metaheuristics
CustomGAParams() = CustomGAParams(100, 0.5, 0.5)

# Inicialização do estado
function Metaheuristics.initialize!(
    status,
    parameters::CustomGAParams,
    problem::Metaheuristics.Problem,
    information::Metaheuristics.Information,
    options::Metaheuristics.Options
)
    population = Metaheuristics.gen_initial_state(problem, parameters, information, options, status)
    Metaheuristics.evaluate_population(population, problem, options)
    best_solution = Metaheuristics.best(population)
    return Metaheuristics.State(best_solution, population)
end

# Lógica do algoritmo por iteração
function Metaheuristics.update_state!(
    state::Metaheuristics.State,
    parameters::CustomGAParams,
    problem::Metaheuristics.Problem,
    information::Metaheuristics.Information,
    options::Metaheuristics.Options
)
    # 1. Seleção de pais
    parents = Metaheuristics.select_parents(
        state.population,
        Metaheuristics.TournamentSelection(),
        options
    )

    # 2. Crossover
    offsprings = crossover_harmonious_coloring(parents)

    # 3. Mutação
    for child in eachrow(offsprings)
        graph_swap_mutation!(child)
    end

    # 4. Avaliar
    Metaheuristics.evaluate_population(offsprings, problem, options)

    # 5. Substituir
    Metaheuristics.replace!(
        Metaheuristics.GreedyReplacement(),
        state,
        offsprings,
        options
    )

    return true
end

# --- Função de Fitness ---
function fitness_harmonious_coloring(individual::Vector{Float64})
    lista_prioridade = sortperm(individual, rev = true)
    cores_vertices = coloracaoHarmonicaPrioridade!(matriz_adj, lista_prioridade)
    return maximum(cores_vertices)
end

# --- Limites de busca ---
bounds = [zeros(num_vertices) ones(num_vertices)]' 

# --- Cria o algoritmo genético ---
#=my_custom_ga = Metaheuristics.Algorithm(CustomGAParams)  # pass the type, not an instance

# --- Cria o problema ---
problem = Metaheuristics.Problem(fitness_harmonious_coloring, bounds)

# --- Executa otimização ---
result = Metaheuristics.optimize(fitness_harmonious_coloring, bounds, my_custom_ga)=#
# --- Instancia parâmetros do GA ---
params = CustomGAParams(100, 0.5, 0.5)

# --- Cria o algoritmo genético ---
my_custom_ga = Metaheuristics.Algorithm(params)

# --- Executa otimização ---
result = Metaheuristics.optimize(fitness_harmonious_coloring, bounds, my_custom_ga)


@show result

# --- Extrai melhor indivíduo e cores ---
melhor_individuo = Metaheuristics.minimizer(result)
lista_prioridade = sortperm(melhor_individuo, rev = true)
cores_vertices = coloracaoHarmonicaPrioridade!(matriz_adj, lista_prioridade)

println("\n--- Coloração harmônica final ---")
for i in 1:num_vertices
    println("Vértice $i: cor $(cores_vertices[i])")
end
println("Número total de cores usadas: ", maximum(cores_vertices))=#


#=ERRO COM CHECK_COMPACT
using Metaheuristics
using LinearAlgebra

# Inclui o arquivo com as funções de coloração gulosa
include("colorGul.jl")

# --- Leitura do grafo ---
println("Insira o nome do arquivo a ser lido: ")
nome_arquivo = readline()

num_vertices, num_arestas = leInfo!(nome_arquivo)
println("nro de vertices: ", num_vertices)
println("nro de arestas: ", num_arestas)

matriz_adj = zeros(Int, num_vertices, num_vertices)
leArestas!(nome_arquivo, matriz_adj)

# --- Funções de operadores personalizados ---
function crossover_harmonious_coloring(parents::Matrix{Float64})
    p1 = parents[1, :]
    p2 = parents[2, :]
    child1 = (p1 .+ p2) ./ 2.0
    child2 = (p1 .+ p2) ./ 2.0 .+ 0.05 .* randn(length(p1))
    return [child1'; child2']
end

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

# --- Algoritmo genético personalizado ---
mutable struct CustomGAParams <: Metaheuristics.AbstractParameters
    N::Int                 # population size
    p_crossover::Float64
    p_mutation::Float64
end

CustomGAParams(; N = 100, p_crossover = 0.5, p_mutation = 0.5) =
    CustomGAParams(N, p_crossover, p_mutation)

# Inicialização do estado
function Metaheuristics.initialize!(
    status,
    parameters::CustomGAParams,
    problem::Metaheuristics.Problem,
    information::Metaheuristics.Information,
    options::Metaheuristics.Options
)
    population = Metaheuristics.gen_initial_state(problem, parameters, information, options, status)
    Metaheuristics.evaluate_population!(population, problem, options)
    best_solution = Metaheuristics.best(population)
    return Metaheuristics.State(best_solution, population)
end

# Lógica do algoritmo por iteração
function Metaheuristics.update_state!(
    state::Metaheuristics.State,
    parameters::CustomGAParams,
    problem::Metaheuristics.Problem,
    information::Metaheuristics.Information,
    options::Metaheuristics.Options
)
    # 1. Seleção de pais
    parents = Metaheuristics.select_parents(
        state.population,
        Metaheuristics.TournamentSelection(),
        options
    )

    # 2. Crossover
    offsprings = crossover_harmonious_coloring(parents)

    # 3. Mutação
    for child in eachrow(offsprings)
        graph_swap_mutation!(child)
    end

    # 4. Avaliar
    Metaheuristics.evaluate_population!(offsprings, problem, options)

    # 5. Substituir
    Metaheuristics.replace!(
        Metaheuristics.GreedyReplacement(),
        state,
        offsprings,
        options
    )

    return true
end

# --- Função de Fitness ---
function fitness_harmonious_coloring(individual::Vector{Float64})
    lista_prioridade = sortperm(individual, rev = true)
    cores_vertices = coloracaoHarmonicaPrioridade(matriz_adj, lista_prioridade)
    return maximum(cores_vertices)
end

# --- Limites de busca ---
bounds = [zeros(num_vertices) ones(num_vertices)]'

# --- Instancia parâmetros do GA ---
#=params = CustomGAParams(N = 100, p_crossover = 0.5, p_mutation = 0.5)

# --- Executa otimização ---
result = Metaheuristics.optimize(fitness_harmonious_coloring, bounds, params)=#
# --- Instancia parâmetros do GA e wrap em Algorithm ---
my_custom_ga = Metaheuristics.Algorithm(CustomGAParams)

# --- Executa otimização ---
result = Metaheuristics.optimize(fitness_harmonious_coloring, bounds, my_custom_ga)

@show result

# --- Extrai melhor indivíduo e cores ---
melhor_individuo = Metaheuristics.minimizer(result)
lista_prioridade = sortperm(melhor_individuo, rev = true)
cores_vertices = coloracaoHarmonicaPrioridade(matriz_adj, lista_prioridade)

println("\n--- Coloração harmônica final ---")
for i in 1:num_vertices
    println("Vértice $i: cor $(cores_vertices[i])")
end
println("Número total de cores usadas: ",  maximum(cores_vertices))=#


#=PROBLEMA COM PARAMETRO N PARA USAR gen_initial_state
using Metaheuristics
using LinearAlgebra

# Inclui o arquivo com as funções de coloração gulosa
include("colorGul.jl")

# --- Leitura do grafo ---
println("Insira o nome do arquivo a ser lido: ")
nome_arquivo = readline()

num_vertices, num_arestas = leInfo!(nome_arquivo)
println("nro de vertices: ", num_vertices)
println("nro de arestas: ", num_arestas)

matriz_adj = zeros(Int, num_vertices, num_vertices)
leArestas!(nome_arquivo, matriz_adj)

# --- Funções de operadores personalizados ---
function crossover_harmonious_coloring(parents::Matrix{Float64})
    p1 = parents[1, :]
    p2 = parents[2, :]
    child1 = (p1 .+ p2) ./ 2.0
    child2 = (p1 .+ p2) ./ 2.0 .+ 0.05 .* randn(length(p1))
    return [child1'; child2']
end

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

# --- Definição do Algoritmo Personalizado ---
mutable struct CustomGAParams <: Metaheuristics.AbstractParameters
    N::Int
    p_crossover::Float64
    p_mutation::Float64
end

CustomGAParams(; N = 100, p_crossover = 0.5, p_mutation = 0.5) =
    CustomGAParams(N, p_crossover, p_mutation)

# Inicialização do estado
function Metaheuristics.initialize!(
    status,
    parameters::CustomGAParams,
    problem::Metaheuristics.Problem,
    information::Metaheuristics.Information,
    options::Metaheuristics.Options
)
    population = Metaheuristics.gen_initial_state(problem, parameters.populationSize, information, options, status)
    Metaheuristics.evaluate_population!(population, problem, options)
    best_solution = Metaheuristics.best(population)
    return Metaheuristics.State(best_solution, population)
end

# Lógica do algoritmo por iteração
function Metaheuristics.update_state!(
    state::Metaheuristics.State,
    parameters::CustomGAParams,
    problem::Metaheuristics.Problem,
    information::Metaheuristics.Information,
    options::Metaheuristics.Options
)
    # 1. Seleção de pais
    parents = Metaheuristics.select_parents(
        state.population,
        Metaheuristics.TournamentSelection(),
        options
    )

    # 2. Crossover
    offsprings = crossover_harmonious_coloring(parents)

    # 3. Mutação
    for i in axes(offsprings, 1)
        graph_swap_mutation!(offsprings[i, :])
    end

    # 4. Avaliar
    Metaheuristics.evaluate_population!(offsprings, problem, options)

    # 5. Substituir
    Metaheuristics.replace!(
        Metaheuristics.GreedyReplacement(),
        state,
        offsprings,
        options
    )

    return true
end

# --- Função de Fitness ---
function fitness_harmonious_coloring(individual::Vector{Float64})
    lista_prioridade = sortperm(individual, rev = true)
    cores_vertices = coloracaoHarmonicaPrioridade(matriz_adj, lista_prioridade)
    return maximum(cores_vertices)
end

# --- Limites de busca e Execução ---
bounds = [zeros(num_vertices) ones(num_vertices)]'

# Instância do seu GA
params = CustomGAParams(N = 100, p_crossover = 0.5, p_mutation = 0.5)
my_custom_ga = Metaheuristics.Algorithm(params)  # wrap parameters

problem = Metaheuristics.Problem(fitness_harmonious_coloring, bounds)

# Executa otimização
result = Metaheuristics.optimize(fitness_harmonious_coloring, bounds, my_custom_ga)

@show result

melhor_individuo = Metaheuristics.minimizer(result)
lista_prioridade = sortperm(melhor_individuo, rev = true)
cores_vertices = coloracaoHarmonicaPrioridade(matriz_adj, lista_prioridade)

println("\n--- Coloração harmônica final ---")
for i in 1:num_vertices
    println("Vértice $i: cor $(cores_vertices[i])")
end
println("Número total de cores usadas: ", maximum(cores_vertices))=#


#=CHAT COM ERRO EM INTIALIZE!
using Metaheuristics
using LinearAlgebra

# Inclui o arquivo com as funções de coloração gulosa
include("colorGul.jl")

# --- Leitura do grafo ---
println("Insira o nome do arquivo a ser lido: ")
nome_arquivo = readline()

num_vertices, num_arestas = leInfo!(nome_arquivo)
println("nro de vertices: ", num_vertices)
println("nro de arestas: ", num_arestas)

matriz_adj = zeros(Int, num_vertices, num_vertices)
leArestas!(nome_arquivo, matriz_adj)

# --- Funções de operadores personalizados ---
function crossover_harmonious_coloring(parents::Matrix{Float64})
    p1 = parents[1, :]
    p2 = parents[2, :]
    child1 = (p1 .+ p2) ./ 2.0
    child2 = (p1 .+ p2) ./ 2.0 .+ 0.05 .* randn(length(p1))
    return [child1'; child2']
end

function graph_swap_mutation!(x::Vector{Float64})
    n = length(x)
    if n < 2
        return x
    end
    v1 = rand(1:n)
    vizinhos = [i for i in 1:n if matriz_adj[v1, i] == 1]
    if isempty(vizinhos)
        return x
    end
    v2 = rand(vizinhos)
    x[v1], x[v2] = x[v2], x[v1]
    return x
end

# --- Definição do Algoritmo Personalizado ---
mutable struct CustomGAParams <: Metaheuristics.AbstractParameters
    populationSize::Int
    p_crossover::Float64
    p_mutation::Float64
end

CustomGAParams(; N = 100, p_crossover = 0.5, p_mutation = 0.5) =
    CustomGAParams(N, p_crossover, p_mutation)

# Inicialização do estado
function Metaheuristics.initial_state!(
    parameters::CustomGAParams,
    problem::Metaheuristics.Problem,
    information::Metaheuristics.Information,
    options::Metaheuristics.Options
)
    population = Metaheuristics.gen_initial_state(problem, parameters.populationSize, information, options)
    Metaheuristics.evaluate_population!(population, problem, options)
    best_solution = Metaheuristics.best(population)
    return Metaheuristics.State(best_solution, population)
end

# Lógica do algoritmo por iteração
function Metaheuristics.update_state!(
    state::Metaheuristics.State,
    parameters::CustomGAParams,
    problem::Metaheuristics.Problem,
    information::Metaheuristics.Information,
    options::Metaheuristics.Options
)
    # 1. Seleção de pais
    parents = Metaheuristics.select_parents(
        state.population,
        Metaheuristics.TournamentSelection(),
        options
    )

    # 2. Crossover
    offsprings = crossover_harmonious_coloring(parents)

    # 3. Mutação
    for eachindex in offsprings
        graph_swap_mutation!(offsprings[i, :])
    end

    # 4. Avaliar
    Metaheuristics.evaluate_population!(offsprings, problem, options)

    # 5. Substituir
    Metaheuristics.replace!(
        Metaheuristics.GreedyReplacement(),
        state,
        offsprings,
        options
    )

    return true
end

# --- Função de Fitness ---
function fitness_harmonious_coloring(individual::Vector{Float64})
    lista_prioridade = sortperm(individual, rev = true)
    cores_vertices = coloracaoHarmonicaPrioridade(matriz_adj, lista_prioridade)
    return maximum(cores_vertices)
end

# --- Limites de busca e Execução ---
bounds = [zeros(num_vertices) ones(num_vertices)]'

# Instância do seu GA
params = CustomGAParams(N = 100, p_crossover = 0.5, p_mutation = 0.5)
my_custom_ga = Metaheuristics.Algorithm(params)  # <-- wrap aqui

problem = Metaheuristics.Problem(fitness_harmonious_coloring, bounds)

# Executa otimização
result = Metaheuristics.optimize(problem, my_custom_ga)

@show result

melhor_individuo = Metaheuristics.minimizer(result)
lista_prioridade = sortperm(melhor_individuo, rev = true)
cores_vertices = coloracaoHarmonicaPrioridade(matriz_adj, lista_prioridade)

println("\n--- Coloração harmônica final ---")
for i in 1:num_vertices
    println("Vértice $i: cor $(cores_vertices[i])")
end
println("Número total de cores usadas: ", maximum(cores_vertices))=#


#O QUE ESTÁ DANDO PROBLEMA COM PARAMETERS
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

matriz_adj = zeros(Int, num_vertices, num_vertices)
leArestas!(nome_arquivo, matriz_adj)

# --- Funções de operadores personalizados (as "ferramentas") ---
# Nota: A matriz_adj é acessada globalmente nessas funções
function crossover_harmonious_coloring(parents::Matrix{Float64})
    p1 = parents[1, :]
    p2 = parents[2, :]
    child1 = (p1 .+ p2) ./ 2.0
    child2 = (p1 .+ p2) ./ 2.0 .+ 0.05 .* randn(length(p1))
    return [child1'; child2']
end

function graph_swap_mutation(x::Vector{Float64})
    n = length(x)
    if n < 2 return x end
    #seleciona o primeiro vértice a ser alterado
    v1 = rand(1:n)
    vizinhos = [i for i in 1:n if matriz_adj[v1, i] == 1]
    #se não há vizinhos, não faz a troca 
    #talvez mudar para selecionar algum outro vértice
    if isempty(vizinhos)
        return x
    end
    v2 = rand(vizinhos)
    x[v1], x[v2] = x[v2], x[v1]
    return x
end

# --- Definição do Algoritmo Personalizado ---
# Este é o container para as regras do seu GA

mutable struct CustomGA <: Metaheuristics.AbstractAlgorithm
    populationSize::Int
    p_crossover::Float64
    p_mutation::Float64
    information::Metaheuristics.Information
    options::Metaheuristics.Options
end

function CustomGA(; N = 100, p_crossover = 0.5, p_mutation = 0.5, 
                    information = Metaheuristics.Information(), options = Metaheuristics.Options())
    return CustomGA(N, p_crossover, p_mutation, information, options)
    #Algorithm(parameters, information = information, options = options,)
end
#CustomGA(; populationSize = 100) = CustomGA(populationSize)

# Define o estado inicial do seu algoritmo
function Metaheuristics.initialize!(status, parameters::CustomGA, problem, 
                                    information, options)
    population = Metaheuristics.gen_initial_state(problem, parameters, information, options, status)
    #=    parameters.populationSize,
        problem,
        Metaheuristics.RandomInBounds()
    )=#
    return Metaheuristics.State(best_solution, population)
end

# Define a lógica do seu algoritmo a cada iteração
function Metaheuristics.update_state!(
    state::Metaheuristics.State,
    parameters::CustomGA,
    problem::Metaheuristics.Problem,
    information::Metaheuristics.Information,
    options::Metaheuristics.Options
)
    # 1. Seleção de pais
    parents = Metaheuristics.select_parents(
        state.population,
        Metaheuristics.TournamentSelection(),
        options
    )
    
    # 2. Geração de filhos (crossover) usando o SEU método
    # Nota: A sua função é chamada diretamente aqui.
    offsprings = crossover_harmonious_coloring(parents)
    
    # 3. Mutação nos filhos usando o SEU método
    # Nota: A sua função é chamada diretamente aqui.
    Metaheuristics.mutation!(
        graph_swap_mutation,
        offsprings[1]
    )
    
    # 4. Avaliar o fitness dos novos indivíduos
    Metaheuristics.evaluate_population!(offsprings, problem, options)
    
    # 5. Substituir os indivíduos da população antiga
    Metaheuristics.replace!(
        Metaheuristics.GreedyReplacement(),
        state,
        offsprings,
        options
    )
    
    # Retorna 'true' para continuar a otimização
    return true
end

# --- Função de Fitness ---
function fitness_harmonious_coloring(individual::Vector{Float64})
    lista_prioridade = sortperm(individual, rev=true)
    cores_vertices = coloracaoHarmonicaPrioridade(matriz_adj, lista_prioridade)
    return maximum(cores_vertices)
end

# --- Limites de busca e Execução ---
bounds = [zeros(num_vertices) ones(num_vertices)]'

my_custom_ga = CustomGA()
problem = Metaheuristics.Problem(fitness_harmonious_coloring, bounds)

#result = optimize(problem, my_custom_ga)
result = optimize(fitness_harmonious_coloring, bounds, my_custom_ga)

@show result

melhor_individuo = Metaheuristics.minimizer(result)
lista_prioridade = sortperm(melhor_individuo, rev=true)
cores_vertices = coloracaoHarmonicaPrioridade(matriz_adj, lista_prioridade)

println("\n--- Coloração harmônica final ---")
for i in 1:num_vertices
    println("Vértice $i: cor $(cores_vertices[i])")
end
println("Número total de cores usadas: ", maximum(cores_vertices))=#

#=using Metaheuristics
import Metaheuristics: initialize!, update_state!, final_stage!
import Metaheuristics: AbstractParameters, gen_initial_state, Algorithm, get_position

import Metaheuristics: SBX_crossover, polynomial_mutation!, create_solution, is_better
import Metaheuristics: reset_to_violated_bounds!

include("colorGul.jl")

mutable struct MeuGA <: AbstractParameters
    N::Int # tamanho da população
    p_crossover::Float64 # probabilidade de crossover
    p_mutation::Float64 # probabilidade de mutação
end

#construtor para a struct criada
function MeuGA(;N = 100,
                p_crossover = 0.9,
                p_mutation = 0.1,
                information = Information(),
                options = Options()
              )
    parameters = MeuGA(N, p_crossover, p_mutation)

    Algorithm(
        parameters,
        information = information,
        options = options,)
end

function initialize!(
        status,
        parameters::MeuGA,
        problem,
        information,
        options,
        args...;
        kargs...
    )

    if options.iterations == 0
        options.iterations = 500
    end

    if options.f_calls_limit == 0
        options.f_calls_limit = options.iterations * parameters.N + 1
    end

    # gen_initial_state require that `parameters.N` is defined.
    return gen_initial_state(problem,parameters,information,options,status)

end

function update_state!(
        status,
        parameters::MeuGA,
        problem,
        information,
        options,
        args...;
        kargs...
    )

    population = status.population
    N = parameters.N

    for i in 1:N
        # selection
        parent_1 = get_position(rand(population))
        parent_2 = get_position(rand(population))

        # generate offspring  via SBX crossover
        c,_ = crossover_harmonious_coloring(parent_1, parent_2, problem.bounds, 20, parameters.p_crossover)

        # Mutate solution
        polynomial_mutation!(c, problem.bounds, 15, parameters.p_mutation)
        # Fix solution if necessary
        reset_to_violated_bounds!(c, problem.bounds)

        # crate the solution and evaluate fitness (x, f(x))
        offspring = create_solution(c, problem)

        push!(population, offspring)
    end

    # environmental selection
    sort!(population, lt = is_better, alg=PartialQuickSort(N))
    deleteat!(population, N+1:length(population))

end=#

#=ERRO NA VERSÃO SIMPLES - NÃO USA OS PARAMETROS DADOS
using Metaheuristics
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
    lista_prioridade = sortperm(individual, rev = true)
    cores_vertices = coloracaoHarmonicaPrioridade!(matriz_adj, lista_prioridade)
    return maximum(cores_vertices)
end

# --- Limites de busca ---
bounds = [zeros(num_vertices) ones(num_vertices)]  # [lower, upper] for each variable

# --- Operadores personalizados ---
function crossover_harmonious(parents::Matrix{Float64})
    p1, p2 = parents[1, :], parents[2, :]
    child1 = (p1 .+ p2) ./ 2
    child2 = (p1 .+ p2) ./ 2 .+ 0.05 .* randn(length(p1))
    return [child1'; child2']
end

function mutation_swap!(x::Vector{Float64})
    n = length(x)
    v1 = rand(1:n)
    neighbors = [i for i in 1:n if matriz_adj[v1, i] == 1]
    if !isempty(neighbors)
        v2 = rand(neighbors)
        x[v1], x[v2] = x[v2], x[v1]
    end
    return x
end

# --- Parâmetros do GA ---
mutable struct CustomGAParams <: Metaheuristics.AbstractParameters
    N::Int
    crossover_rate::Float64
    mutation_rate::Float64
end

CustomGAParams(; N=100, crossover_rate=0.5, mutation_rate=0.5) =
    CustomGAParams(N, crossover_rate, mutation_rate)

params = CustomGAParams()

# --- Cria o algoritmo genético usando seus operadores ---
ga = Metaheuristics.Algorithm(params;
    crossover=crossover_harmonious,
    mutation=mutation_swap!
)

# --- Cria o problema ---
problem = Metaheuristics.Problem(fitness_harmonious_coloring, bounds)

# --- Executa otimização ---
result = Metaheuristics.optimize(fitness_harmonious_coloring, bounds, ga)

# --- Mostra resultado ---
melhor_individuo = Metaheuristics.minimizer(result)
lista_prioridade = sortperm(melhor_individuo, rev = true)
cores_vertices = coloracaoHarmonicaPrioridade!(matriz_adj, lista_prioridade)

println("\n--- Coloração harmônica final ---")
for i in 1:num_vertices
    println("Vértice $i: cor $(cores_vertices[i])")
end
println("Número total de cores usadas: ", maximum(cores_vertices))=#


#=TENTANDO ENTENDER O QUE FAZER NO MÉTODO "SIMPLES"
using Metaheuristics
using LinearAlgebra

# Inclui o arquivo com as funções de coloração gulosa
include("colorGul.jl")

# --- Leitura do grafo ---
println("Insira o nome do arquivo a ser lido: ")
nome_arquivo = readline()

num_vertices, num_arestas = leInfo!(nome_arquivo)
println("nro de vertices: ", num_vertices)
println("nro de arestas: ", num_arestas)

matriz_adj = zeros(Int, num_vertices, num_vertices)
leArestas!(nome_arquivo, matriz_adj)

# --- Funções de operadores personalizados ---
function crossover_harmonious_coloring(parents::Matrix{Float64})
    p1 = parents[1, :]
    p2 = parents[2, :]
    child1 = (p1 .+ p2) ./ 2.0
    child2 = (p1 .+ p2) ./ 2.0 .+ 0.05 .* randn(length(p1))
    return [child1'; child2']
end

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

meuGA = GA(N=100, p_mutation = 0.5, p_crossover = 0.5, initializer = RandomInBounds(), 
           selection=TournamentSelection(), crossover = crossover_harmonious_coloring(), 
           mutation = graph_swap_mutation!(), environmental_selection = ElitistReplacement())=#

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

matriz_adj = zeros(Int, num_vertices, num_vertices)
leArestas!(nome_arquivo, matriz_adj)

# --- Funções de operadores personalizados ---
function crossover_harmonious_coloring(parents::Matrix{Float64})
    p1 = parents[1, :]
    p2 = parents[2, :]
    child1 = (p1 .+ p2) ./ 2.0
    child2 = (p1 .+ p2) ./ 2.0 .+ 0.05 .* randn(length(p1))
    return [child1'; child2']
end

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

# --- Função de Fitness ---
function fitness_harmonious_coloring(individual::Vector{Float64})
    lista_prioridade = sortperm(individual, rev=true)
    cores_vertices = coloracaoHarmonicaPrioridade!(matriz_adj, lista_prioridade)
    return maximum(cores_vertices)
end

# --- Limites de busca e Execução ---
bounds = [zeros(num_vertices) ones(num_vertices)]'

# A forma mais simples e robusta: o GA é criado passando as funções diretamente
ga_optimizer = GA(
    populationSize=100,
    recombination=crossover_harmonious_coloring,
    mutation=graph_swap_mutation
)

result = optimize(fitness_harmonious_coloring, bounds, ga_optimizer)

@show result

melhor_individuo = Metaheuristics.minimizer(result)
lista_prioridade = sortperm(melhor_individuo, rev=true)
cores_vertices = coloracaoHarmonicaPrioridade!(matriz_adj, lista_prioridade)

println("\n--- Coloração harmônica final ---")
for i in 1:num_vertices
    println("Vértice $i: cor $(cores_vertices[i])")
end
println("Número total de cores usadas: ", maximum(cores_vertices))=#

#=versão simples que não funciona sla pq
using Metaheuristics
using LinearAlgebra

# Inclui o arquivo com as funções de coloração gulosa
include("colorGul.jl")

# --- Leitura do grafo ---
println("Insira o nome do arquivo a ser lido: ")
nome_arquivo = readline()

num_vertices, num_arestas = leInfo!(nome_arquivo)
println("nro de vertices: ", num_vertices)
println("nro de arestas: ", num_arestas)

matriz_adj = zeros(Int, num_vertices, num_vertices)
leArestas!(nome_arquivo, matriz_adj)

# --- Funções de operadores personalizados ---
function crossover_harmonious_coloring(parents::Matrix{Float64})
    p1 = parents[1, :]
    p2 = parents[2, :]
    child1 = (p1 .+ p2) ./ 2.0
    child2 = (p1 .+ p2) ./ 2.0 .+ 0.05 .* randn(length(p1))
    return [child1'; child2']
end

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

# --- Função de Fitness ---
function fitness_harmonious_coloring(individual::Vector{Float64})
    lista_prioridade = sortperm(individual, rev = true)
    cores_vertices = coloracaoHarmonicaPrioridade!(matriz_adj, lista_prioridade)
    return maximum(cores_vertices)
end

# --- Limites de busca ---
bounds = [zeros(num_vertices) ones(num_vertices)]

# --- Cria o algoritmo genético usando seus operadores ---
ga = GA(
    population_size = 100,
    crossover_rate  = 0.5,
    mutation_rate   = 0.5,
    crossover       = crossover_harmonious_coloring,
    mutation        = graph_swap_mutation!
)

# --- Executa otimização ---
result = optimize(fitness_harmonious_coloring, bounds, ga)

@show result

# --- Extrai melhor indivíduo e cores ---
melhor_individuo = result.minimizer
lista_prioridade = sortperm(melhor_individuo, rev = true)
cores_vertices = coloracaoHarmonicaPrioridade!(matriz_adj, lista_prioridade)

println("\n--- Coloração harmônica final ---")
for i in 1:num_vertices
    println("Vértice $i: cor $(cores_vertices[i])")
end
println("Número total de cores usadas: ", maximum(cores_vertices))=#