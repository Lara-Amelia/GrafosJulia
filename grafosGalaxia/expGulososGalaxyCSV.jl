using Metaheuristics
using LinearAlgebra
using Statistics
using Printf
using DataFrames  
using CSV         

# Inclui as funções de coloração do arquivo pai (colorGul.jl) 
include("../colorGul.jl")

# Variáveis globais para processamento 
global matriz_adj = Matrix{Int}(undef, 0, 0)
global num_vertices = 0
global num_arestas = 0

"""
Extrai parâmetros (a, b, c, v) do nome do arquivo Galaxy[cite: 22].
"""
function extract_galaxy_params(file_name)
    match_result = match(r"galaxy_a(\d+)_b(\d+)_c(\d+)_v(\d+)\.col", file_name)
    
    if match_result !== nothing
        return Dict(
            :a_param => parse(Int, match_result.captures[1]),
            :b_param => parse(Int, match_result.captures[2]),
            :c_param => parse(Int, match_result.captures[3]),
            :v_param => parse(Int, match_result.captures[4])
        )
    else
        return nothing
    end
end

"""
Executa uma heurística gulosa e mede o tempo
"""
function run_greedy_experiment(file_name, heuristic::Function)
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
    
    local elapsed_time
    local cores_resultado
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
    heuristics = Dict(
        :max_deg => coloracaoHarmonicaGrauMax!,
        :min_deg => coloracaoHarmonicaGrauMin!,
        :sat_deg => NOVOcoloracaoHarmonicaSaturacao!
    )

    all_files = filter(f -> startswith(f, "galaxy_") && endswith(f, ".col"), readdir())
    
    # Filtragem a <= 25
    filtered_files = filter(all_files) do f
        p = extract_galaxy_params(f)
        p !== nothing && p[:a_param] <= 25
    end
    
    sort!(filtered_files)
    num_files = length(filtered_files)
    println("--- Iniciando Experimento em $num_files arquivos ---")

    detailed_results = []
    summary_results = []

    for (index, file_name) in enumerate(filtered_files)
        println("($index/$num_files) Processando $file_name...")

        # Execução das 3 heurísticas 
        max_c, max_t = run_greedy_experiment(file_name, heuristics[:max_deg])
        min_c, min_t = run_greedy_experiment(file_name, heuristics[:min_deg])
        sat_c, sat_t = run_greedy_experiment(file_name, heuristics[:sat_deg])

        params = extract_galaxy_params(file_name)
        nv, na = leInfo!(file_name)

        # Lógica para encontrar o melhor resultado 
        results_list = [
            (name = "MaxGrau", colors = max_c, time = max_t),
            (name = "MinGrau", colors = min_c, time = min_t),
            (name = "Satur",   colors = sat_c, time = sat_t)
        ]
        
        valid = filter(r -> r.colors > 0, results_list)
        best_c, best_t, best_name = 0, 0.0, "N/A"
        
        if !isempty(valid)
            sort!(valid, by = r -> (r.colors, r.time))
            best_c, best_t, best_name = valid[1].colors, valid[1].time, valid[1].name
        end

        # Adiciona ao conjunto detalhado 
        push!(detailed_results, Dict(
            :instancia => file_name,
            :a => params[:a_param], :b => params[:b_param], :c => params[:c_param], :v => params[:v_param],
            :N => nv, :M => na,
            :chi_sat => sat_c, :tempo_sat => sat_t,
            :chi_min => min_c, :tempo_min => min_t,
            :chi_max => max_c, :tempo_max => max_t
        ))

        # Adiciona ao conjunto de resumo 
        push!(summary_results, Dict(
            :instancia => file_name,
            :a => params[:a_param], :b => params[:b_param], :c => params[:c_param], :v => params[:v_param],
            :N => nv, :M => na,
            :melhor_chi => best_c,
            :tempo_melhor => best_t,
            :heuristica_vencedora => best_name
        ))
    end

    # Criação dos DataFrames e Ordenação de Colunas
    df_detailed = DataFrame(detailed_results)
    select!(df_detailed, [:instancia, :a, :b, :c, :v, :N, :M, :chi_sat, :tempo_sat, :chi_min, :tempo_min, :chi_max, :tempo_max])

    df_summary = DataFrame(summary_results)
    select!(df_summary, [:instancia, :a, :b, :c, :v, :N, :M, :melhor_chi, :tempo_melhor, :heuristica_vencedora])

    # Geração dos arquivos CSV
    CSV.write("results_Greedy_todos.csv", df_detailed)
    CSV.write("results_Greedy_melhor.csv", df_summary)

    println("\n--- Concluído! CSVs gerados e ordenados com sucesso. ---")
end

main()