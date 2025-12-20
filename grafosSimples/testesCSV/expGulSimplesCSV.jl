using Statistics
using Printf
using DelimitedFiles
using LinearAlgebra 
using DataFrames  #
using CSV         #

# Inclui as funções de coloração (Busca no diretório pai)
include("../../colorGul.jl") 

# Variáveis globais para processamento
global matriz_adj = Matrix{Int}(undef, 0, 0)
global num_vertices = 0
global num_arestas = 0

"""
Extrai parâmetros (n, p, v) do nome do arquivo Simples.
"""
function extract_simple_params(file_name)
    # Padrão: simples_n(N)_p(P)%_v(V).col
    match_result = match(r"simples_n(\d+)_p(\d+)%_v(\d+)\.col", file_name)
    
    if match_result !== nothing
        return Dict(
            :n_param => parse(Int, match_result.captures[1]),
            :p_param => parse(Int, match_result.captures[2]),
            :v_param => parse(Int, match_result.captures[3])
        )
    else
        # Caso o arquivo não siga o padrão, retorna valores padrão para evitar erro
        return Dict(:n_param => 0, :p_param => 0, :v_param => 0)
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
        @error "Falha ao ler info do grafo $file_name: $e"
        return 0, 0.0
    end

    local elapsed_time
    try
        # Medição robusta utilizando @elapsed
        elapsed_time = @elapsed begin
            cores_resultado = heuristic(matriz_adj)
        end
        return maximum(cores_resultado), elapsed_time
    catch e
        @error "Erro na execução da heurística no arquivo $file_name: $e"
        return 0, 0.0
    end
end

function main()
    # Definição das Heurísticas
    heuristics = Dict(
        :max_deg => coloracaoHarmonicaGrauMax!,
        :min_deg => coloracaoHarmonicaGrauMin!,
        :sat_deg => NOVOcoloracaoHarmonicaSaturacao! 
    )

    # Coleta de TODOS os arquivos que seguem o padrão de nome simples
    all_files = filter(f -> startswith(f, "simples_") && endswith(f, ".col"), readdir())

    if isempty(all_files)
        println("AVISO: Nenhum arquivo Simples (.col) encontrado no diretório. Saindo.")
        return
    end

    sort!(all_files)
    num_files = length(all_files)
    println("--- Iniciando Experimentos Gulosos (Sem Filtros) ---")
    println("Processando $num_files arquivos Simples encontrados...")

    detailed_results = []
    summary_results = []

    # Execução dos Experimentos
    for (i, file_name) in enumerate(all_files)
        # Extração de parâmetros para o CSV e para o log no terminal
        params = extract_simple_params(file_name)
        
        println("($i/$num_files) [N=$(params[:n_param]), p=$(params[:p_param]), v=$(params[:v_param])] Processando $file_name...")

        nv, na = leInfo!(file_name) 

        # Executa as 3 heurísticas
        max_c, max_t = run_greedy_experiment(file_name, heuristics[:max_deg])
        min_c, min_t = run_greedy_experiment(file_name, heuristics[:min_deg])
        sat_c, sat_t = run_greedy_experiment(file_name, heuristics[:sat_deg])

        # Lógica do melhor resultado (Cores primeiro, depois Tempo como desempate)
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

        # Armazenamento detalhado
        push!(detailed_results, Dict(
            :instancia => file_name,
            :n => params[:n_param], :p => params[:p_param], :v => params[:v_param],
            :N => nv, :M => na,
            :chi_sat => sat_c, :tempo_sat => round(sat_t, digits=4),
            :chi_min => min_c, :tempo_min => round(min_t, digits=4),
            :chi_max => max_c, :tempo_max => round(max_t, digits=4)
        ))

        # Armazenamento resumo
        push!(summary_results, Dict(
            :instancia => file_name,
            :n => params[:n_param], :p => params[:p_param], :v => params[:v_param],
            :N => nv, :M => na,
            :melhor_chi => best_c, 
            :tempo_melhor => round(best_t, digits=4),
            :heuristica_vencedora => best_name
        ))
    end

    # Geração dos DataFrames e Ordenação de Colunas
    df_detailed = DataFrame(detailed_results)
    select!(df_detailed, [:instancia, :n, :p, :v, :N, :M, :chi_sat, :tempo_sat, :chi_min, :tempo_min, :chi_max, :tempo_max])

    df_summary = DataFrame(summary_results)
    select!(df_summary, [:instancia, :n, :p, :v, :N, :M, :melhor_chi, :tempo_melhor, :heuristica_vencedora])

    # Exportação Final
    CSV.write("results_Greedy_Simple_Detailed_All.csv", df_detailed)
    CSV.write("results_Greedy_Simple_Summary_All.csv", df_summary)

    println("\n--- Concluído! Todos os $(num_files) arquivos processados. ---")
    println("Arquivos CSV gerados com sucesso.")
end

main()