using Metaheuristics
using LinearAlgebra
using Statistics
using Printf
using DataFrames  
using CSV         

# Inclui as funções de coloração do arquivo pai (colorGul.jl) 
include("../colorGul.jl")

# Variáveis globais para processamento 
#=global matriz_adj = Matrix{Int}(undef, 0, 0)
global num_vertices = 0
global num_arestas = 0=#

global ADJ = Vector{Vector{Int}}()
global V = 0

"""
Extrai parâmetros (a, b, c, v) do nome do arquivo Galaxy.
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
        #=:max_deg => coloracaoHarmonicaGrauMax!,
        :min_deg => coloracaoHarmonicaGrauMin!,
        :sat_deg => NOVOcoloracaoHarmonicaSaturacao!=#
        :max_deg => coloracaoHarmonicaGrauMaxAdj!, 
        :min_deg => coloracaoHarmonicaGrauMinAdj!,
        :sat_deg => coloracaoHarmonicaSaturacaoAdj!
    )

    # Coleta todos os arquivos sem filtros de parâmetros
    all_files = filter(f -> startswith(f, "galaxy_") && endswith(f, ".col"), readdir())
    
    # Ordenação natural dos arquivos
    sort!(all_files)
    num_files = length(all_files)
    
    println("--- Iniciando Experimento em $num_files arquivos ---")

    detailed_results = []
    summary_results = []

    for (index, file_name) in enumerate(all_files)
        println("($index/$num_files) Processando $file_name...")

        # Carrega o grafo uma única vez para este arquivo para economizar I/O
        #=global matriz_adj, num_vertices, num_arestas
        try
            num_vertices, num_arestas = leInfo!(file_name)
            matriz_adj = zeros(Int, num_vertices, num_vertices)
            leArestas!(file_name, matriz_adj)
        catch e
            @warn "Falha ao carregar $file_name, pulando... Erro: $e"
            continue
        end=#

        # Execução das 3 heurísticas 
        max_c, max_t, num_vertices, num_arestas = run_greedy_experiment(file_name, heuristics[:max_deg])
        min_c, min_t, num_vertices, num_arestas = run_greedy_experiment(file_name, heuristics[:min_deg])
        sat_c, sat_t, num_vertices, num_arestas = run_greedy_experiment(file_name, heuristics[:sat_deg])

        params = extract_galaxy_params(file_name)
        
        # Lógica para encontrar o melhor resultado 
        results_list = [
            (name = "MaxGrau", colors = max_c, time = max_t),
            (name = "MinGrau", colors = min_c, time = min_t),
            (name = "Satur",   colors = sat_c, time = sat_t)
        ]
        
        valid = filter(r -> r.colors > 0, results_list)
        best_c, best_t, best_name = 0, 0.0, "N/A"
        
        if !isempty(valid)
            # Ordena por cor (menor melhor) e tempo (menor melhor)
            sort!(valid, by = r -> (r.colors, r.time))
            best_c, best_t, best_name = valid[1].colors, valid[1].time, valid[1].name
        end

        # Metadados básicos para o dicionário
        # Se params for nothing (arquivo fora do padrão galaxy), usamos valores padrão
        a = params !== nothing ? params[:a_param] : 0
        b = params !== nothing ? params[:b_param] : 0
        c = params !== nothing ? params[:c_param] : 0
        v = params !== nothing ? params[:v_param] : 0

        # Adiciona ao conjunto detalhado 
        push!(detailed_results, Dict(
            :instancia => file_name,
            :a => a, :b => b, :c => c, :v => v,
            :N => num_vertices, :M => num_arestas,
            :chi_sat => sat_c, :tempo_sat => sat_t,
            :chi_min => min_c, :tempo_min => min_t,
            :chi_max => max_c, :tempo_max => max_t
        ))

        # Adiciona ao conjunto de resumo 
        push!(summary_results, Dict(
            :instancia => file_name,
            :a => a, :b => b, :c => c, :v => v,
            :N => num_vertices, :M => num_arestas,
            :melhor_chi => best_c,
            :tempo_melhor => best_t,
            :heuristica_vencedora => best_name
        ))
    end

    # Criação dos DataFrames e Seleção de Colunas
    if !isempty(detailed_results)
        df_detailed = DataFrame(detailed_results)
        select!(df_detailed, [:instancia, :a, :b, :c, :v, :N, :M, :chi_sat, :tempo_sat, :chi_min, :tempo_min, :chi_max, :tempo_max])

        df_summary = DataFrame(summary_results)
        select!(df_summary, [:instancia, :a, :b, :c, :v, :N, :M, :melhor_chi, :tempo_melhor, :heuristica_vencedora])

        # Geração dos arquivos CSV
        CSV.write("resultsGulGalaxyADJ.csv", df_detailed)
        CSV.write("resultsGulGalaxyADJ_sumario.csv", df_summary)

        println("\n--- Concluído! CSVs gerados para $num_files arquivos. ---")
    else
        println("Nenhum resultado processado.")
    end
end

main()