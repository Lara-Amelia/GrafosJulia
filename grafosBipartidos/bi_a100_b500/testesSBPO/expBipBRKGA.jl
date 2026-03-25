# experimentos SBPO BRKGA Puro
using Metaheuristics
using LinearAlgebra
using Statistics
using Printf
using DataFrames
using CSV
using ProgressMeter
using Random

import Metaheuristics: initialize!, update_state!, final_stage!, gen_initial_state

# Inclui o arquivo com as funções de coloração gulosa
include("../../../colorGul.jl")

# vars. globais
global V = 0
global ADJ = Vector{Int}[]
global P_BUFFER = Int[]

# wrapper
mutable struct StagnationWrapper <: Metaheuristics.AbstractParameters
    internal_params::Metaheuristics.AbstractParameters
    stag_limit::Int
    stag_iters::Int
    last_best::Float64
end

function Metaheuristics.initialize!(status, parameters::StagnationWrapper, problem, information, options)
    parameters.stag_iters = 0
    parameters.last_best = Inf
    return Metaheuristics.initialize!(status, parameters.internal_params, problem, information, options)
end 

function Metaheuristics.update_state!(state, parameters::StagnationWrapper, problem, information, options)
    old_best = state.best_sol.f

    Metaheuristics.update_state!(state, parameters.internal_params, problem, information, options)

    if state.best_sol.f < old_best
        parameters.last_best = state.best_sol.f
        parameters.stag_iters = 0
    else
        parameters.stag_iters += 1
    end
end

function Metaheuristics.stop_criteria!(status, parameters::StagnationWrapper, problem, information, options)
    status.stop = status.stop || 
                  Metaheuristics.iteration_stop_check(status, information, options) ||
                  Metaheuristics.time_stop_check(status, information, options)

    if parameters.stag_iters >= parameters.stag_limit
        status.stop = true
        status.termination_status_code = Metaheuristics.OTHER_LIMIT
    end
    return status.stop
end

function Metaheuristics.final_stage!(state, parameters::StagnationWrapper, problem, information, options)
    return Metaheuristics.final_stage!(state, parameters.internal_params, problem, information, options)
end

# configuração do BRKGA com critério de parada personalizado
function StagnatedBRKGA(; num_elites=20, num_mutants=10, num_offsprings=70, bias=0.7, stag_limit=50)
    brkga_base = BRKGA(
        num_elites = num_elites, 
        num_mutants = num_mutants, 
        num_offsprings = num_offsprings, 
        bias = bias
    )
    wrapped_params = StagnationWrapper(brkga_base.parameters, stag_limit, 0, Inf)
    return Metaheuristics.Algorithm(wrapped_params, options = brkga_base.options)
end

# fitness (problem-dependent part)
function fitness_harmonious_coloring(individual::Vector{Float64})
    sortperm!(P_BUFFER, individual, rev = true)
    cores_vertices = coloracaoHarmonicaAdjVetAux!(ADJ, P_BUFFER)
    return Float64(maximum(cores_vertices))
end

#instanciação e execução experimento com BRKGA
function run_brkga_experiment(k_limit::Int, N_pop::Int)
    bounds = [zeros(V) ones(V)]'
    
    # Define as proporções baseadas no N_pop (ex: 20% elite, 10% mutante)
    n_elite = Int(floor(0.20 * N_pop))
    n_mutant = Int(floor(0.10 * N_pop))
    n_offspring = N_pop - n_elite - n_mutant

    # Instancia o BRKGA com o Wrapper de estagnação
    alg = StagnatedBRKGA(
        num_elites = n_elite, 
        num_mutants = n_mutant, 
        num_offsprings = n_offspring, 
        bias = 0.7, 
        stag_limit = k_limit
    )

    # Configurações de precisão e limites
    alg.options = Metaheuristics.Options(
        f_calls_limit = typemax(Int), 
        iterations = 10000, 
        f_tol = -1.0, 
        x_tol = -1.0
    )

    result = Metaheuristics.optimize(fitness_harmonious_coloring, bounds, alg)
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
    K_STAG = 50
    N_POP = 100
    SAVE_EVERY = 10
    
    csv_path = "results_BRKGA_Bi_stag50_progresso1.csv"

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
                ch = run_brkga_experiment(K_STAG, N_POP)
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