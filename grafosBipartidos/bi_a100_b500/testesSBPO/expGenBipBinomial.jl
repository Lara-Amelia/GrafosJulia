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
include("../../../colorGul.jl")

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

function CustomGAParams(; N = 100, p_mutation = 0.5, stag_limit = 20, last_best = Inf, stag_iters = 0)
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

# experimentação
function run_ga_experiment(k_limit::Int, N_pop::Int)
    bounds = [zeros(V) ones(V)]'

    params = CustomGAParams(N=N_pop, p_mutation=0.5, stag_limit=k_limit)

    opt_settings = Metaheuristics.Options(
    f_calls_limit = typemax(Int), 
    iterations = 10000, 
    store_convergence = true,
    f_tol = -1, 
    x_tol = -1)

    my_ga = Metaheuristics.Algorithm(params, options = opt_settings)
    
    result = Metaheuristics.optimize(fitness_harmonious_coloring, bounds, my_ga)
    #@show result
    return Int(Metaheuristics.minimum(result))
end

function extract_bipartite_params(file_name)
    # Regex ajustada para o padrão bi_a100_b100_p...
    match_result = match(r"bi_a(\d+)_b(\d+)_p(\d+)%_v(\d+)\.col", file_name)
    if match_result !== nothing
        return (
            a = parse(Int, match_result.captures[1]),
            b = parse(Int, match_result.captures[2]),
            p = parse(Int, match_result.captures[3]),
            v = parse(Int, match_result.captures[4])
        )
    end
    return nothing
end

function main()
    # filtros de instâncias
    TARGET_A = [100, 500, 1000]        
    TARGET_B = [100, 500]        
    TARGET_P = nothing      # processa todas as probabilidades
    TARGET_V = [1, 2] # processa versões de 1 a 5

    # configs. do experimento
    N_REPETITIONS = 5
    K_STAG = 20
    N_POP = 100
    SAVE_EVERY = 10
    
    csv_path = "results_GABinomial_Bi_pmutation05_stag20_progresso1.csv"

    # busca e filtragem de arquivos de entrada
    raiz_busca = dirname(@__DIR__) 
    
    println("Buscando grafos em: ", raiz_busca)

    # Lista apenas os arquivos .col na pasta pai
    all_files_raw = filter(f -> endswith(f, ".col"), readdir(raiz_busca))
    
    filtered_files = filter(f -> begin
        p_params = extract_bipartite_params(f)
        if isnothing(p_params) return false end
        
        cond_a = isnothing(TARGET_A) || p_params.a in TARGET_A
        cond_b = isnothing(TARGET_B) || p_params.b in TARGET_B
        cond_p = isnothing(TARGET_P) || p_params.p in TARGET_P
        cond_v = isnothing(TARGET_V) || p_params.v in TARGET_V
        
        return cond_a && cond_b && cond_p && cond_v
    end, all_files_raw)

    if isempty(filtered_files)
        println("AVISO: Nenhum arquivo bi_*.col corresponde aos filtros em: ", raiz_busca)
        return
    end

    sort!(filtered_files)
    println("Filtro aplicado! Total de instâncias: ", length(filtered_files))

    # preparação do csv
    header_df = DataFrame(
        a=Int[], b=Int[], N=Int[], p=Int[], M=Int[], v=Int[], 
        mean_time=Float64[], se_time=Float64[], 
        mean_chi=Float64[], se_chi=Float64[], instancia=String[]
    )
    CSV.write(csv_path, header_df)

    results_buffer = Dict{Symbol, Any}[]

    # loop de processamento de instâncias
    @showprogress 1 "Bipartidos: " for (idx, file) in enumerate(filtered_files)
        p_info = extract_bipartite_params(file)
        
        full_path = joinpath(raiz_busca, file)
        
        # atualização de globais para novo grafo/instância
        num_v, num_a = leInfo!(full_path)
        global V = num_v
        global ADJ = [Int[] for _ in 1:V]
        leArestasLista!(full_path, ADJ)
        global P_BUFFER = zeros(Int, V) 

        println("\n[$idx/$(length(filtered_files))] Processando: $file")

        t_tot, chi = Float64[], Int[]

        for i in 1:N_REPETITIONS
            elapsed = @elapsed begin
                ch = run_ga_experiment(K_STAG, N_POP)
            end
            push!(t_tot, elapsed) 
            push!(chi, ch)
            @printf("   Run %d -> Chi: %d | %.4fs\n", i, ch, elapsed)
        end

        # Adiciona ao buffer com os dados mapeados corretamente
        push!(results_buffer, Dict(
            :a => p_info.a, 
            :b => p_info.b, 
            :N => num_v, 
            :p => p_info.p, 
            :M => num_a, 
            :v => p_info.v,
            :mean_time => mean(t_tot),
            :se_time   => std(t_tot)/sqrt(N_REPETITIONS),
            :mean_chi  => mean(chi),
            :se_chi    => std(chi)/sqrt(N_REPETITIONS),
            :instancia => file
        ))

        # flush e limpeza de RAM a cada 10 arquivos
        if idx % SAVE_EVERY == 0 || idx == length(filtered_files)
            df_to_save = DataFrame(results_buffer)
            # Reordena colunas para bater com o cabeçalho CSV
            col_order = [:a, :b, :N, :p, :M, :v, :mean_time, :se_time, :mean_chi, :se_chi, :instancia]
            select!(df_to_save, col_order)
            
            CSV.write(csv_path, df_to_save, append=true)
            
            empty!(results_buffer)
            GC.gc()
            println("\n[SISTEMA] Checkpoint salvo e RAM liberada.")
        end
    end
    println("\n--- Experimentos bipartidos concluídos! ---")
end

main()