using Statistics
using Printf
using DelimitedFiles
using LinearAlgebra
using DataFrames  # Essencial para organizar os dados em tabelas
using CSV         # Essencial para exportar arquivos .csv

# Inclui as funções de coloração (Busca no diretório pai)
include("../../colorGul.jl")

# Variáveis globais (usadas pelas funções de leitura no colorGul.jl)
global matriz_adj = Matrix{Int}(undef, 0, 0)
global num_vertices = 0
global num_arestas = 0

"""
Extrai parâmetros (a, b, p, v) do nome do arquivo Bipartido.
"""
function extract_bipartite_params(file_name)
    match_result = match(r"bi_a(\d+)_b(\d+)_p(\d+)%_v(\d+)\.col", file_name)
    
    if match_result !== nothing
        return Dict(
            :a_param => parse(Int, match_result.captures[1]),
            :b_param => parse(Int, match_result.captures[2]),
            :p_param => parse(Int, match_result.captures[3]),
            :v_param => parse(Int, match_result.captures[4])
        )
    else
        return nothing
    end
end

"""
Executa uma heurística gulosa e mede o tempo.
Retorna o número de cores usadas e o tempo de execução.
"""
function run_greedy_experiment(file_name::String, heuristic::Function)
    global matriz_adj
    global num_vertices
    global num_arestas

    try
        num_vertices, num_arestas = leInfo!(file_name)
        matriz_adj = zeros(Int, num_vertices, num_vertices)
        leArestas!(file_name, matriz_adj)
    catch e
        @error "Falha ao ler grafo $file_name: $e"
        return 0, 0.0
    end

    local cores_resultado
    local elapsed_time

    try
        elapsed_time = @elapsed begin
            cores_resultado = heuristic(matriz_adj)
        end
        return maximum(cores_resultado), elapsed_time
    catch e
        @error "Erro na heurística no arquivo $file_name: $e"
        return 0, 0.0
    end
end

function main()
    # 1. Definição das Heurísticas
    heuristics = Dict(
        :max_deg => coloracaoHarmonicaGrauMax!,
        :min_deg => coloracaoHarmonicaGrauMin!,
        :sat_deg => NOVOcoloracaoHarmonicaSaturacao! 
    )

    # 2. Coleta e filtragem (a <= 100, b <= 100)
    all_files = filter(f -> startswith(f, "bi_") && endswith(f, ".col"), readdir())
    filtered_names = String[]

    for file_name in all_files
        m = extract_bipartite_params(file_name)
        if m !== nothing && m[:a_param] <= 100 && m[:b_param] <= 100
            push!(filtered_names, file_name)
        end
    end

    if isempty(filtered_names)
        println("AVISO: Nenhum arquivo bipartido encontrado para o filtro. Saindo.")
        return
    end

    sort!(filtered_names)
    num_files = length(filtered_names)
    println("--- Starting Experiments: $num_files Bipartite files ---")

    detailed_results = []
    summary_results = []

    # 3. Execução dos experimentos
    for (i, file_name) in enumerate(filtered_names)
        println("($i/$num_files) Processing $file_name...")

        # Coleta dimensões
        nv, na = leInfo!(file_name)
        params = extract_bipartite_params(file_name)

        # Roda as 3 heurísticas
        max_c, max_t = run_greedy_experiment(file_name, heuristics[:max_deg])
        min_c, min_t = run_greedy_experiment(file_name, heuristics[:min_deg])
        sat_c, sat_t = run_greedy_experiment(file_name, heuristics[:sat_deg])

        # Lógica do melhor resultado
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

        # 4. Organização dos dados
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

    # 5. Geração e Ordenação dos CSVs
    df_detailed = DataFrame(detailed_results)
    select!(df_detailed, [:instancia, :a, :b, :p, :v, :N, :M, :chi_sat, :tempo_sat, :chi_min, :tempo_min, :chi_max, :tempo_max])

    df_summary = DataFrame(summary_results)
    select!(df_summary, [:instancia, :a, :b, :p, :v, :N, :M, :melhor_chi, :tempo_melhor, :heuristica_vencedora])

    CSV.write("results_Greedy_Bipartite_Detailed.csv", df_detailed)
    CSV.write("results_Greedy_Bipartite_Summary.csv", df_summary)

    println("\n--- Geração Finalizada: Arquivos CSV salvos com sucesso. ---")
end

main()