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
include("../../colorGul.jl")

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

# instanciação e execução experimento com BRKGA
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

function main()
    # configuração de filtros para instâncias 
    FILTER_N = [100, 500, 1000, 1500, 2000] 
    FILTER_P = [1, 3, 5, 10, 20, 30, 40] 
    FILTER_V = [1, 2]

    # configuração dos testes
    N_REPETITIONS = 5
    K_STAG = 50
    N_POP = 100
    SAVE_EVERY = 10  # frequência de limpeza da memória e salvamento no disco
    
    csv_path = "results_BRKGA_progresso.csv"

    # localização e filtragem
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
        println("AVISO: Nenhum arquivo .col corresponde aos filtros.")
        return
    end

    sort!(all_files)
    println("Filtro aplicado! Total de instâncias: ", length(all_files))

    # cabeçalho corrigido (dados se alinham corretamente com as infos do cabeçalho)
    header_df = DataFrame(
        instancia=String[], n=Int[], p=Int[], v=Int[], 
        N=Int[], M=Int[], mean_chi=Float64[], se_chi=Float64[], 
        mean_time=Float64[], se_time=Float64[]
    )
    CSV.write(csv_path, header_df)

    results_buffer = Dict{Symbol, Any}[]

    @showprogress 1 "Processando: " for (idx, file_path) in enumerate(all_files)
        file_name = basename(file_path)
        
        m_n = match(r"n(\d+)", file_name); n_param = m_n !== nothing ? parse(Int, m_n.captures[1]) : 0
        m_p = match(r"p(\d+)", file_name); p_param = m_p !== nothing ? parse(Int, m_p.captures[1]) : 0
        m_v = match(r"v(\d+)", file_name); v_param = m_v !== nothing ? parse(Int, m_v.captures[1]) : 0

        num_v, num_a = leInfo!(file_path)
        global V = num_v
        global ADJ = [Int[] for _ in 1:V]
        leArestasLista!(file_path, ADJ)
        global P_BUFFER = zeros(Int, V) 

        println("\n[$idx/$(length(all_files))] Instância: $file_name")

        t_tot, chi = Float64[], Int[]

        for i in 1:N_REPETITIONS
            elapsed = @elapsed begin
                ch = run_brkga_experiment(K_STAG, N_POP)
            end
            push!(t_tot, elapsed) 
            push!(chi, ch)
            @printf("   Run %d -> Chi: %d | Tempo: %.4fs\n", i, ch, elapsed)
        end

        push!(results_buffer, Dict(
            :instancia => file_name,
            :n         => n_param,
            :p         => p_param,
            :v         => v_param,
            :N         => num_v,
            :M         => num_a,
            :mean_chi  => mean(chi),
            :se_chi    => std(chi)/sqrt(N_REPETITIONS),
            :mean_time => mean(t_tot),
            :se_time   => std(t_tot)/sqrt(N_REPETITIONS)
        ))

        if idx % SAVE_EVERY == 0 || idx == length(all_files)
            df_to_append = DataFrame(results_buffer)
           
            select!(df_to_append, [:instancia, :n, :p, :v, :N, :M, :mean_chi, :se_chi, :mean_time, :se_time])
            
            CSV.write(csv_path, df_to_append, append=true)
            
            empty!(results_buffer)
            GC.gc()
            println("\n[SISTEMA] Checkpoint salvo e RAM liberada.")
        end
    end
end

main()