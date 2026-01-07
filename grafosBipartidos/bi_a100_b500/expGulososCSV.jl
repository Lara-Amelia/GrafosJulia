using Statistics
using Printf
using DelimitedFiles
using LinearAlgebra
using DataFrames  
using CSV          
using ProgressMeter

# Inclui as funções de coloração
include("../../colorGul.jl")

#=global matriz_adj = Matrix{Int}(undef, 0, 0)
global num_vertices = 0
global num_arestas = 0=#

global ADJ = Vector{Vector{Int}}()
global V = 0

function extract_bipartite_params(file_name)
    match_result = match(r"bi_a(\d+)_b(\d+)_p(\d+)%_v(\d+)\.col", file_name)
    if match_result !== nothing
        return Dict(
            :a_param => parse(Int, match_result.captures[1]),
            :b_param => parse(Int, match_result.captures[2]),
            :p_param => parse(Int, match_result.captures[3]),
            :v_param => parse(Int, match_result.captures[4])
        )
    end
    return nothing
end

function run_greedy_experiment(file_name::String, heuristic::Function)
    #global matriz_adj
    global ADJ, V
    try
        #=num_v, num_a = leInfo!(file_name)
        matriz_adj = zeros(Int, num_v, num_v)
        leArestas!(file_name, matriz_adj)=#

        global ADJ, V
        num_v, num_a = leInfo!(file_name)
        ADJ = [Int[] for _ in 1:num_v]
        leArestasLista!(file_name, ADJ)
        V = num_v
        
        local cores_resultado
        elapsed_time = @elapsed begin
            #cores_resultado = heuristic(matriz_adj)
            cores_resultado = heuristic(ADJ)
        end
        return maximum(cores_resultado), elapsed_time, num_v, num_a
    catch e
        @error "Erro no arquivo $file_name: $e"
        return 0, 0.0, 0, 0
    end
end

function main()
    # 1. Definição das Heurísticas
    heuristics = Dict(
        #=:max_deg => coloracaoHarmonicaGrauMax!,
        :min_deg => coloracaoHarmonicaGrauMin!,
        :sat_deg => NOVOcoloracaoHarmonicaSaturacao!=#
        :max_deg => coloracaoHarmonicaGrauMaxAdj!, 
        :min_deg => coloracaoHarmonicaGrauMinAdj!,
        :sat_deg => coloracaoHarmonicaSaturacaoAdj!
    )

    # 2. Filtragem de Arquivos
    all_files = filter(f -> startswith(f, "bi_") && endswith(f, ".col"), readdir())
    filtered_file_names = String[]

    limite_a = 1000
    limite_b = 1000

    for file_name in all_files
        m = match(r"bi_a(\d+)_b(\d+)_p(\d+)%_v(\d+)\.col", file_name)
        if m === nothing continue end

        a_p = parse(Int, m.captures[1])
        b_p = parse(Int, m.captures[2])

        if a_p <= limite_a && b_p <= limite_b
            push!(filtered_file_names, file_name)
        end
    end

    if isempty(filtered_file_names)
        println("AVISO: Nenhum arquivo corresponde aos critérios.")
        return
    end

    sort!(filtered_file_names)
    num_files = length(filtered_file_names)
    println("--- Starting Greedy Experiments: $num_files files ---")

    detailed_results = []
    summary_results = []

    # 3. Loop de Execução com ProgressMeter
    @showprogress 1 "Processando: " for (idx, file_name) in enumerate(filtered_file_names)
        # Extração de parâmetros para o CSV
        params = extract_bipartite_params(file_name)

        # Execução das 3 heurísticas
        max_c, max_t, nv, na = run_greedy_experiment(file_name, heuristics[:max_deg])
        min_c, min_t, _, _   = run_greedy_experiment(file_name, heuristics[:min_deg])
        sat_c, sat_t, _, _   = run_greedy_experiment(file_name, heuristics[:sat_deg])

        # Lógica para encontrar a melhor entre as gulosas
        res_list = [
            (name = "MaxGrau", colors = max_c, time = max_t),
            (name = "MinGrau", colors = min_c, time = min_t),
            (name = "Satur",   colors = sat_c, time = sat_t)
        ]
        
        valid = filter(r -> r.colors > 0, res_list)
        best_c, best_t, best_name = 0, 0.0, "Falha"

        if !isempty(valid)
            sort!(valid, by = r -> (r.colors, r.time))
            best_c, best_t, best_name = valid[1].colors, valid[1].time, valid[1].name
        end

        @printf("\n[%d/%d] %s | Best Greedy: %d (%s em %.4fs)\n", 
                idx, num_files, file_name, best_c, best_name, best_t)

        # 4. Armazenamento para CSV
        push!(detailed_results, Dict(
            :instancia => file_name,
            :a => params[:a_param], :b => params[:b_param], :p => params[:p_param], :v => params[:v_param],
            :N => nv, :M => na,
            :chi_sat => sat_c, :tempo_sat => round(sat_t, digits=4),
            :chi_min => min_c, :tempo_min => round(min_t, digits=4),
            :chi_max => max_c, :tempo_max => round(max_t, digits=4)
        ))

        push!(summary_results, Dict(
            :instancia => file_name,
            :a => params[:a_param], :b => params[:b_param], :p => params[:p_param], :v => params[:v_param],
            :N => nv, :M => na,
            :melhor_chi => best_c,
            :tempo_melhor => round(best_t, digits=4),
            :heuristica_vencedora => best_name
        ))
    end

    # 5. Geração dos CSVs
    df_detailed = DataFrame(detailed_results)
    sort!(df_detailed, [:a, :b, :p, :v])
    CSV.write("resultsGulADJ.csv", df_detailed)

    df_summary = DataFrame(summary_results)
    sort!(df_summary, [:a, :b, :p, :v])
    CSV.write("resultsGulADJ_sumario.csv", df_summary)

    println("\n--- Experimentos Concluídos! ---")
end

main()