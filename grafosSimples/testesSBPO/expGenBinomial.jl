# experimentos SBPO GA crossover Binomial com Operadores Customizados e Checkpoints
using Metaheuristics
using LinearAlgebra
using Statistics
using Printf
using DataFrames
using CSV
using ProgressMeter
using Random

import Metaheuristics: initialize!, update_state!, final_stage!, gen_initial_state

# inclui o arquivo com as funções de coloração gulosa
include("../../colorGul.jl")

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

function CustomGAParams(; N = 100, p_mutation = 0.5, stag_limit = 50, last_best = Inf, stag_iters = 0)
    selection_strategy = Metaheuristics.TournamentSelection(K=2, N=N)
    cross_op = Metaheuristics.BinomialCrossover(p = 0.5, n_offsprings = 2)
    
    return CustomGAParams(N, p_mutation, stag_limit, last_best, 
                          stag_iters, selection_strategy, cross_op)
end

# fitness da coloração sem utilizar hash (verificar depois se isso ajuda ou atrapalha)
function fitness_harmonious_coloring(individual::Vector{Float64})
    sortperm!(P_BUFFER, individual, rev = true)

    cores_vertices = coloracaoHarmonicaAdjVetAux!(ADJ, P_BUFFER)
    fitness = maximum(cores_vertices)
    return fitness
end

# dessa vez tentaremos chamar os métodos personalizados pelo wrapper da própria biblioteca,
# possivelmente teremos de alterar alguns deles (para adequar às formas e tipos esperados)
# aparentemente, não será necessário sobrescrever o método initialize! (ou podemos apenas chamar gen_initial_state)
# population será um array simples de valores 

# um jeito mais fácil de realmente integrar os novos métodos pode ser reescrever update_state!,
# mas apenas substituindo as chamadas a crossover e mutação por nossos próprios métodos
function graph_swap_mutation!(Q::AbstractMatrix{Float64})
    n_individuals, n_genes = size(Q)
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
    #=combined = append!(population, offsprings)
    sort!(combined, alg=PartialQuickSort(N), by = s -> s.f)
    deleteat!(population, (N+1):length(population))=#

    N = length(population) 
    append!(population, offsprings) 
    #sort!(population, by = s -> s.f) 
    sort!(population, lt=Metaheuristics.is_better)
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

# sobre-escrita do método de atualização de estado com os operadores e critério de parada próprios
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

    replacement_elitism(pop, offsprings, parameters.N)
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

    #=if parameters.stag_iters < parameters.stag_limit
        println("nro. de iterações estagnadas: $(parameters.stag_iters)")
        println("AINDA NÃO ATINGIMOS O LIMITE DE ITERAÇÕES SEM MELHORA")
    end=#

    if Metaheuristics.is_better(current_best, state.best_sol)
        state.best_sol = current_best 
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

# execução do experimento
function run_ga_experiment(k_limit::Int, N_pop::Int)
    bounds = [zeros(V) ones(V)]'
    params = CustomGAParams(N=N_pop, p_mutation=0.5, stag_limit=k_limit)
    
    # f_tol = -1 garante que o algoritmo não pare por precisão numérica
    opt_settings = Metaheuristics.Options(
    f_calls_limit = typemax(Int), 
    iterations = 10000, 
    store_convergence = true,
    f_tol = -1, 
    x_tol = -1)
    my_ga = Metaheuristics.Algorithm(params, options = opt_settings)
    
    result = Metaheuristics.optimize(fitness_harmonious_coloring, bounds, my_ga)
    @show result
    return Int(Metaheuristics.minimum(result))
end

function main()
    # configuração de filtros para instâncias 
    FILTER_N = [100, 500, 1000, 1500, 2000] 
    # limitamos a probabilidade porque sabemos que, a partir de um certo ponto, sempre cai no pior caso
    FILTER_P = [1, 3, 5, 10, 20, 30, 40] 
    FILTER_V = [1, 2]

    # configuração dos testes
    N_REPETITIONS = 5
    K_STAG = 50
    N_POP = 100
    SAVE_EVERY = 10  # frequência de limpeza da memória e salvamento no disco
    
    csv_path = "results_GABinomial_progresso.csv"

    # localização e filtragem dos arquivos a serem lidos e processados
    raiz_grafos = dirname(@__DIR__) 
    all_files_raw = String[]
    for (root, dirs, files) in walkdir(raiz_grafos)
        if occursin("testesSBPO", root) || occursin("testesCSV", root) || occursin("LATEX", root)
            continue
        end
        for f in files
            if endswith(f, ".col")
                push!(all_files_raw, joinpath(root, f))
            end
        end
    end

    # aplica os filtros definidos 
    all_files = filter(file_path -> begin
        fname = basename(file_path)
        m_n = match(r"n(\d+)", fname); n = m_n !== nothing ? parse(Int, m_n.captures[1]) : -1
        m_p = match(r"p(\d+)", fname); p = m_p !== nothing ? parse(Int, m_p.captures[1]) : -1
        m_v = match(r"v(\d+)", fname); v = m_v !== nothing ? parse(Int, m_v.captures[1]) : -1

        cond_n = isnothing(FILTER_N) || n in FILTER_N
        cond_p = isnothing(FILTER_P) || p in FILTER_P
        cond_v = isnothing(FILTER_V) || v in FILTER_V

        return cond_n && cond_p && cond_v
    end, all_files_raw)

    if isempty(all_files)
        println("AVISO: Nenhum arquivo .col corresponde aos filtros definidos em $raiz_grafos.")
        return
    end

    sort!(all_files)
    println("Filtro aplicado! Total de instâncias para processar: ", length(all_files))

    # preparação do CSV
    header_df = DataFrame(
        instancia=String[], n=Int[], p=Int[], v=Int[], 
        N=Int[], M=Int[], mean_chi=Float64[], se_chi=Float64[], 
        mean_time=Float64[], se_time=Float64[]
    )
    CSV.write(csv_path, header_df)

    results_buffer = Dict{Symbol, Any}[]

    # loop de experimentos
    @showprogress 1 "Processando: " for (idx, file_path) in enumerate(all_files)
        file_name = basename(file_path)
        
        m_n = match(r"n(\d+)", file_name); n_param = m_n !== nothing ? parse(Int, m_n.captures[1]) : 0
        m_p = match(r"p(\d+)", file_name); p_param = m_p !== nothing ? parse(Int, m_p.captures[1]) : 0
        m_v = match(r"v(\d+)", file_name); v_param = m_v !== nothing ? parse(Int, m_v.captures[1]) : 0

        # carregamento de novo grafo e atualização de variáveis globais
        num_v, num_a = leInfo!(file_path)
        global V = num_v
        global ADJ = [Int[] for _ in 1:V]
        leArestasLista!(file_path, ADJ)
        global P_BUFFER = zeros(Int, V) 

        println("\n[$idx/$(length(all_files))] Instância: $file_name (V=$V, A=$num_a)")

        t_tot, chi = Float64[], Int[]

        for i in 1:N_REPETITIONS
            elapsed = @elapsed begin
                ch = run_ga_experiment(K_STAG, N_POP)
            end
            push!(t_tot, elapsed) 
            push!(chi, ch)
            @printf("   Run %d -> Chi: %d | Tempo: %.4fs\n", i, ch, elapsed)
        end

        # armazenamento no buffer temporário
        push!(results_buffer, Dict(
            :instancia => file_name, :n => n_param, :p => p_param, :v => v_param,
            :N => num_v, :M => num_a, :mean_chi => mean(chi),
            :se_chi => std(chi)/sqrt(N_REPETITIONS),
            :mean_time => mean(t_tot), :se_time => std(t_tot)/sqrt(N_REPETITIONS)
        ))

        # lógica de flush (escrita e liberação de memória)
        if idx % SAVE_EVERY == 0 || idx == length(all_files)
            df_to_append = DataFrame(results_buffer)
            CSV.write(csv_path, df_to_append, append=true)
            
            empty!(results_buffer)
            GC.gc()
            
            println("\n[SISTEMA] Checkpoint salvo e RAM liberada (Instância $idx).")
        end
    end

    println("\n" * "="^40)
    println("EXPERIMENTOS CONCLUÍDOS COM SUCESSO!")
    println("Resultados acumulados em: $csv_path")
    println("="^40)
end

main()