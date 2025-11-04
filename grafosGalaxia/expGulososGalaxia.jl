# Carrega os pacotes essenciais
using Statistics
using Printf
using DelimitedFiles

# Inclui as funções de coloração do arquivo pai (colorGul.jl)
# Certifique-se de que o colorGul.jl está um nível acima deste diretório.
include("../colorGul.jl")

# Variáveis globais para armazenar temporariamente o grafo durante o processamento
global matriz_adj = Matrix{Int}(undef, 0, 0)
global num_vertices = 0
global num_arestas = 0

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
            :v_param => parse(Int, match_result.captures[4]),
            :graph_type => "Galaxy"
        )
    else
        return nothing
    end
end

"""
Função para executar uma heurística gulosa (Max, Min, Sat) e medir o tempo.
(Revertido para 1 execução, conforme solicitado).
Retorna o número de cores usadas e o tempo de execução.
"""
function run_greedy_experiment(file_name, heuristic::Function)
    global matriz_adj
    global num_vertices
    global num_arestas

    try
        # leInfo! e leArestas! estão em colorGul.jl
        num_vertices, num_arestas = leInfo!(file_name)
    catch e
        @error "Falha ao ler info do grafo $file_name: $e"
        return 0, 0.0
    end
    
    matriz_adj = zeros(Int, num_vertices, num_vertices)
    try
        leArestas!(file_name, matriz_adj)
    catch e
        @error "Falha ao ler arestas do grafo $file_name: $e"
        return 0, 0.0
    end
    
    local cores_resultado
    local elapsed_time

    try
        elapsed_time = @elapsed begin
            cores_resultado = heuristic(matriz_adj)
        end
        
        max_colors = maximum(cores_resultado)

        return max_colors, elapsed_time
    catch e
        @error "Erro durante a execução da heurística no arquivo $file_name: $e"
        return 0, 0.0
    end
end

# --- FUNÇÃO DE GERAÇÃO DE TABELA DE RESUMO (MODIFICADA) ---

"""
Gera a tabela LaTeX de Resumo (Melhor de 3), indicando a heurística.
"""
function generate_best_summary_table(results_data, filename, graph_type)
    
    title = "Melhor Resultado Global das Heurísticas Gulosas - Grafos $graph_type (a \\le 25, v=1)"
    
    # Preâmbulo LaTeX
    latex_output = """
\\documentclass[11pt, a4paper]{article}
\\usepackage[a4paper, top=2.5cm, bottom=2.5cm, left=2cm, right=2cm]{geometry}
\\usepackage{booktabs}
\\usepackage{amsmath}
\\usepackage{amssymb}
\\usepackage{fontspec}
\\usepackage[english, bidi=basic, provide=*]{babel}

\\babelfont{rm}{Noto Sans}

\\title{$title}
\\author{Experimentos de Coloração Harmônica}
\\date{\\today}

\\begin{document}
\\maketitle

\\begin{table}[htbp]
  \\centering
  \\caption{Resumo: Melhor Coloração Harmônica (\\min \\chi_h) Encontrada entre as 3 Heurísticas}
  \\label{tab:greedy_summary_$graph_type}
  \\begin{tabular}{@{} l ccc l @{}}
    \\toprule
"""
    
    # --- CABEÇALHO MODIFICADO ---
    # Adicionada a coluna "Heurística"
    header_row = "Instância & N & M & " * "\\chi_h" * " (Melhor) & Tempo (s) & Heurística \\\\\n"
    
    latex_output *= header_row * "\\midrule\n"
    
    # Preenche as linhas da tabela
    for row in results_data
        # --- SPRINTF MODIFICADO ---
        latex_output *= @sprintf(
            "%s & %d & %d & %d & %.4f & %s \\\\\n",
            row[:instance_id],
            row[:N_vertices],
            row[:N_edges],
            row[:best_colors],
            row[:best_time],
            row[:best_heuristic] # Nova coluna
        )
    end

    latex_output *= """
    \\bottomrule
  \\end{tabular}
\\end{table}

\\end{document}
"""
    
    open(filename, "w") do f
        write(f, latex_output)
    end
end

# --- FUNÇÃO DE GERAÇÃO DE TABELA DETALHADA (Revertida) ---

"""
Gera a tabela LaTeX Detalhada (3 Heurísticas).
(Revertida para não mencionar "Tempo Médio")
"""
function generate_detailed_table(results_data, filename, graph_type)
    
    title = "Resultados Detalhados das Heurísticas Gulosas - Grafos $graph_type (a \\le 25, v=1)"

    # Preâmbulo LaTeX
    latex_output = """
\\documentclass[11pt, a4paper]{article}
\\usepackage[a4paper, top=2.5cm, bottom=2.5cm, left=2cm, right=2cm]{geometry}
\\usepackage{booktabs}
\\usepackage{amsmath}
\\usepackage{amssymb}
\\usepackage{fontspec}
\\usepackage[english, bidi=basic, provide=*]{babel}

\\babelfont{rm}{Noto Sans}

\\title{$title}
\\author{Experimentos de Coloração Harmônica}
\\date{\\today}

\\begin{document}
\\maketitle

\\begin{table}[htbp]
  \\centering
  \\caption{Comparativo Detalhado das Heurísticas Gulosas de Coloração Harmônica por Ordem de Vértices}
  \\label{tab:greedy_detailed_$graph_type}
  \\resizebox{\\textwidth}{!}{
  \\begin{tabular}{@{} l cc ccc ccc @{}}
    \\toprule
    \\multicolumn{3}{@{}l}{Parâmetros do Grafo} & \\multicolumn{2}{c}{Grau Máximo} & \\multicolumn{2}{c}{Grau Mínimo} & \\multicolumn{2}{c @{}}{Saturação} \\\\
    \\cmidrule(lr){1-3} \\cmidrule(lr){4-5} \\cmidrule(lr){6-7} \\cmidrule(lr){8-9}
"""
    
    # Linha de cabeçalho (Revertida para "Tempo (s)")
    header_row = "Instância (a, b, c, v) & N & M & " * "\\chi_h" * " & Tempo (s) & " * "\\chi_h" * " & Tempo (s) & " * "\\chi_h" * " & Tempo (s) \\\\\n"
    
    latex_output *= header_row * "\\midrule\n"
    
    # Preenche as linhas da tabela
    for row in results_data
        latex_output *= @sprintf(
            "%s & %d & %d & %d & %.4f & %d & %.4f & %d & %.4f \\\\\n",
            row[:instance_id],
            row[:N_vertices],
            row[:N_edges],
            row[:max_colors],
            row[:max_time],
            row[:min_colors],
            row[:min_time],
            row[:sat_colors],
            row[:sat_time]
        )
    end

    # Finaliza o documento LaTeX
    latex_output *= """
    \\bottomrule
  \\end{tabular}
  } % end resizebox
\\end{table}

\\end{document}
"""
    
    open(filename, "w") do f
        write(f, latex_output)
    end
end

# --- FUNÇÃO PRINCIPAL (LÓGICA DO "MELHOR" ATUALIZADA) ---

function main()
    # 1. Definição das Heurísticas
    heuristics = Dict(
        :max_deg => coloracaoHarmonicaGrauMax!,
        :min_deg => coloracaoHarmonicaGrauMin!,
        :sat_deg => NOVOcoloracaoHarmonicaSaturacao!
    )

    # 2. Filtragem de Arquivos
    all_files = readdir()
    galaxy_files = filter(f -> startswith(f, "galaxy_") && endswith(f, ".col"), all_files)
    
    filtered_file_names = String[]
    
    println("--- Filtrando Arquivos Galaxy (a <= 25 e v=1) ---")
    for file_name in galaxy_files
        data = extract_galaxy_params(file_name)
        if data !== nothing
            a_param = data[:a_param]
            v_param = data[:v_param]

            if a_param <= 25 && v_param == 1
                push!(filtered_file_names, file_name)
            else
                println("  Pulando $file_name (a=$a_param ou v=$v_param)")
            end
        else
            println("  Ignorando $file_name (formato inválido)")
        end
    end
    
    sort!(filtered_file_names)
    num_files = length(filtered_file_names)
    println("\n--- Processando $num_files Arquivos Galaxy Filtrados ---")

    # 3. Execução dos Experimentos
    detailed_results = []
    summary_results = []

    for (index, file_name) in enumerate(filtered_file_names)
        println("($index/$num_files) Processando $file_name...")

        # Executa as 3 heurísticas (1x cada)
        max_colors, max_time = run_greedy_experiment(file_name, heuristics[:max_deg])
        min_colors, min_time = run_greedy_experiment(file_name, heuristics[:min_deg])
        sat_colors, sat_time = run_greedy_experiment(file_name, heuristics[:sat_deg])

        data = extract_galaxy_params(file_name)
        
        # --- LÓGICA ATUALIZADA PARA ENCONTRAR O MELHOR ---
        
        # 1. Armazena os resultados em uma lista de tuplas para facilitar a ordenação
        results_list = [
            (name = "MaxGrau", colors = max_colors, time = max_time),
            (name = "MinGrau", colors = min_colors, time = min_time),
            (name = "Satur",   colors = sat_colors, time = sat_time)
        ]

        # 2. Filtra resultados inválidos (onde a heurística falhou e retornou 0 cores)
        valid_results = filter(r -> r.colors > 0, results_list)
        
        local best_colors
        local best_time
        local best_heuristic_name

        if isempty(valid_results)
            best_colors = 0
            best_time = 0.0
            best_heuristic_name = "Nenhum (Falha)"
            println("  AVISO: Todas as heurísticas falharam para $file_name.")
        else
            # 3. Ordena os resultados:
            #    Critério 1: Menor número de cores (r.colors)
            #    Critério 2 (Desempate): Menor tempo de execução (r.time)
            sort!(valid_results, by = r -> (r.colors, r.time))
            
            # O melhor resultado é o primeiro da lista ordenada
            best_result = valid_results[1]
            
            best_colors = best_result.colors
            best_time = best_result.time
            best_heuristic_name = best_result.name
        end

        # 4. Armazena para as Duas Tabelas
        
        # A. Detalhada (Exatamente como antes)
        push!(detailed_results, Dict(
            :instance_id => @sprintf("a%d\\_b%d\\_c%d\\_v%d", data[:a_param], data[:b_param], data[:c_param], data[:v_param]),
            :N_vertices => num_vertices,
            :N_edges => num_arestas,
            :max_colors => max_colors,
            :max_time => max_time,
            :min_colors => min_colors,
            :min_time => min_time,
            :sat_colors => sat_colors,
            :sat_time => sat_time
        ))

        # B. Resumo (Com a nova informação :best_heuristic)
        push!(summary_results, Dict(
            :instance_id => @sprintf("a%d\\_b%d\\_c%d\\_v%d", data[:a_param], data[:b_param], data[:c_param], data[:v_param]),
            :N_vertices => num_vertices,
            :N_edges => num_arestas,
            :best_colors => best_colors,
            :best_time => best_time,
            :best_heuristic => best_heuristic_name
        ))

        println("  Cores Max/Min/Sat: $max_colors / $min_colors / $sat_colors")
        println("  Tempo Max/Min/Sat: $(@sprintf("%.4f", max_time))s / $(@sprintf("%.4f", min_time))s / $(@sprintf("%.4f", sat_time))s")
        println("  Melhor Resultado (\\chi_h): $best_colors ($best_heuristic_name) em $(@sprintf("%.4f", best_time))s\n")
    end

    # 5. Geração dos Arquivos LaTeX
    
    # Tabela 1: Detalhada
    output_detailed = "results_Greedy_Galaxy_a01-a25_v1_Detailed.tex"
    generate_detailed_table(detailed_results, output_detailed, "Galaxy")
    println("--- Tabela Detalhada Gerada ---")
    println("Resultados salvos em: $output_detailed")

    # Tabela 2: Resumo
    output_summary = "results_Greedy_Galaxy_a01-a25_v1_Summary.tex"
    generate_best_summary_table(summary_results, output_summary, "Galaxy")
    println("--- Tabela de Resumo Gerada ---")
    println("Resultados salvos em: $output_summary")
end

# Executa a função principal
main()