using Metaheuristics
using LinearAlgebra
using Statistics
using Printf
using DataFrames  
using CSV         

# Inclui as funções de coloração do arquivo pai (colorGul.jl) 
include("../../colorGul.jl")

# Variáveis globais para processamento 
#=global matriz_adj = Matrix{Int}(undef, 0, 0)
global num_vertices = 0
global num_arestas = 0=#

global ADJ = Vector{Vector{Int}}()
global V = 0

"""
Extrai parâmetros (a, b, c, v) do nome do arquivo Galaxy.
"""
function extract_simples_params(file_name)
    # Regex ajustada para capturar n (vértices), p (probabilidade %) e v (versão/seed)
    match_result = match(r"simples_n(\d+)_p(\d+)%_v(\d+)\.col", file_name)
    
    if match_result !== nothing
        return Dict(
            :n_param => parse(Int, match_result.captures[1]),
            :p_param => parse(Int, match_result.captures[2]),
            :v_param => parse(Int, match_result.captures[3])
        )
    else
        return nothing
    end
end

"""
Executa uma heurística gulosa e mede o tempo
"""
function run_greedy_experiment(file_name, heuristic::Function)
    #=global matriz_adj
    global num_vertices
    global num_arestas

    # A leitura do grafo e alocação da matriz deve ser feita uma única vez por arquivo no loop principal,
    # mas para manter a compatibilidade com sua estrutura de 'run_greedy_experiment',
    # garantimos que os dados globais estejam carregados.
    
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
    end=#

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
    heuristics = Dict(
        :max_deg => coloracaoHarmonicaGrauMaxAdj!, 
        :min_deg => coloracaoHarmonicaGrauMinAdj!,
        :sat_deg => coloracaoHarmonicaSaturacaoAdj!
    )

    # Coleta todos os arquivos com o novo padrão "simples_"
    all_files = filter(f -> startswith(f, "simples_") && endswith(f, ".col"), readdir())
    
    # Ordenação natural dos arquivos
    sort!(all_files)
    num_files = length(all_files)
    
    println("--- Iniciando Experimento em $num_files arquivos (Padrão Simples) ---")

    detailed_results = []
    summary_results = []

    for (index, file_name) in enumerate(all_files)
        println("($index/$num_files) Processando $file_name...")

        # Execução das 3 heurísticas
        max_c, max_t, nv, na = run_greedy_experiment(file_name, heuristics[:max_deg])
        min_c, min_t, _, _    = run_greedy_experiment(file_name, heuristics[:min_deg])
        sat_c, sat_t, _, _    = run_greedy_experiment(file_name, heuristics[:sat_deg])

        params = extract_simples_params(file_name)
        
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

        # Metadados extraídos do nome do arquivo
        n_file = params !== nothing ? params[:n_param] : 0
        p_file = params !== nothing ? params[:p_param] : 0
        v_file = params !== nothing ? params[:v_param] : 0

        # Adiciona ao conjunto detalhado 
        push!(detailed_results, Dict(
            :instancia => file_name,
            :n_config => n_file, :p_config => p_file, :v_config => v_file,
            :N_real => nv, :M_real => na,
            :chi_sat => sat_c, :tempo_sat => round(sat_t, digits=6),
            :chi_min => min_c, :tempo_min => round(min_t, digits=6),
            :chi_max => max_c, :tempo_max => round(max_t, digits=6)
        ))

        # Adiciona ao conjunto de resumo 
        push!(summary_results, Dict(
            :instancia => file_name,
            :n_config => n_file, :p_config => p_file, :v_config => v_file,
            :N_real => nv, :M_real => na,
            :melhor_chi => best_c,
            :tempo_melhor => round(best_t, digits=6),
            :heuristica_vencedora => best_name
        ))
    end

    # Criação dos DataFrames e Seleção de Colunas
    if !isempty(detailed_results)
        df_detailed = DataFrame(detailed_results)
        # Ordenação por N e P para facilitar leitura do CSV
        sort!(df_detailed, [:n_config, :p_config, :v_config])
        
        df_summary = DataFrame(summary_results)
        sort!(df_summary, [:n_config, :p_config, :v_config])

        # Geração dos arquivos CSV
        CSV.write("resultsGulSimplesADJ.csv", df_detailed)
        CSV.write("resultsGulSimplesADJ_sumario.csv", df_summary)

        println("\n--- Concluído! CSVs gerados para $num_files arquivos. ---")
    else
        println("Nenhum resultado processado.")
    end
end

main()