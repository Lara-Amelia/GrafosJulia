# Importa pacotes necessários
using Statistics
using Printf
using DelimitedFiles
using LinearAlgebra # Adicionado para garantir a matriz_adj zeroed

# Inclui as funções de coloração (Busca no diretório pai)
include("../../colorGul.jl")

# Variáveis globais (usadas pelas funções de leitura no colorGul.jl)
global matriz_adj = Matrix{Int}(undef, 0, 0)
global num_vertices = 0
global num_arestas = 0

# --- FUNÇÕES AUXILIARES DE LEITURA (Assumidas no arquivo pai) ---

# Implementação de 'leInfo!' e 'leArestas!' assumidas no ../colorGul.jl
# e a função 'coloracaoHarmonicaSaturacao!' e suas variantes também.

function extract_bipartite_params(file_name)
    match_result = match(r"bi_a(\d+)_b(\d+)_p(\d+)%_v(\d+)\.col", file_name)
    
    if match_result !== nothing
        return Dict(
            :a_param => parse(Int, match_result.captures[1]),
            :b_param => parse(Int, match_result.captures[2]),
            :p_param => parse(Int, match_result.captures[3]),
            :v_param => parse(Int, match_result.captures[4]),
            :graph_type => "Bipartite"
        )
    else
        return nothing
    end
end

"""
Função para executar UMA HEURÍSTICA GULOSA (Max, Min, Sat) e medir o tempo.
Retorna o número de cores usadas e o tempo de execução.
"""
function run_greedy_experiment(file_name::String, heuristic::Function)
    global matriz_adj
    global num_vertices
    global num_arestas

    # 1. Leitura do Grafo (Feita uma vez para o arquivo)
    try
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

    # 2. Execução Única e Medição
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

# --- FUNÇÃO DE GERAÇÃO DE TABELA DE RESUMO (ÚNICO RESULTADO) ---

"""
Gera a tabela LaTeX de Resumo (Melhor entre as 3 heurísticas).
"""
function generate_best_summary_table(results_data, filename, graph_type)
    
    title = "Melhor Resultado Global das Heurísticas Gulosas - Grafos $graph_type (a \\le 100, b \\le 100)"
    
    latex_output = """
\\documentclass[11pt, a4paper]{article}
\\usepackage[a4paper, top=2.5cm, bottom=2.5cm, left=2cm, right=2cm]{geometry}
\\usepackage{booktabs}
\\usepackage{amsmath}
\\usepackage{amssymb}
\\usepackage{siunitx}

\\sisetup{round-mode = places, round-precision = 4, table-align-text-pre = false}

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
    Instância & N & M & \$\\chi_h\$ (Melhor) & Tempo (s) & Heurística \\\\
    \\midrule
"""
    
    # Preenche as linhas da tabela
    for row in results_data
        
        instance_id = @sprintf("a%d\\_b%d\\_p%d\\%%\\_v%d", row[:a_param], row[:b_param], row[:p_param], row[:v_param])
        
        latex_output *= @sprintf(
            "%s & %d & %d & %d & \\num{%.4f} & %s \\\\\n",
            instance_id,
            row[:N_vertices],
            row[:N_edges],
            row[:best_colors],
            row[:best_time],
            row[:best_heuristic] 
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

# --- FUNÇÃO DE GERAÇÃO DE TABELA DETALHADA (ÚNICO RESULTADO) ---

"""
Gera a tabela LaTeX Detalhada (3 Heurísticas - 1 execução).
"""
function generate_detailed_table(results_data, filename, graph_type)
    
    title = "Resultados Detalhados das Heurísticas Gulosas - Grafos $graph_type (a \\le 100, b \\le 100)"

    latex_output = """
\\documentclass[11pt, a4paper]{article}
\\usepackage[a4paper, top=2.5cm, bottom=2.5cm, left=1.5cm, right=1.5cm]{geometry}
\\usepackage{booktabs}
\\usepackage{amsmath}
\\usepackage{amssymb}
\\usepackage{siunitx}
\\usepackage{longtable}

\\sisetup{round-mode = places, round-precision = 4, table-align-text-pre = false}

\\title{$title}
\\author{Experimentos de Coloração Harmônica}
\\date{\\today}

\\begin{document}
\\maketitle

\\begin{longtable}{@{} l cc ccc ccc @{}}
    \\toprule
    \\multicolumn{3}{@{}l}{Parâmetros do Grafo} & \\multicolumn{2}{c}{Grau Máximo} & \\multicolumn{2}{c}{Grau Mínimo} & \\multicolumn{2}{c @{}}{Saturação} \\\\
    \\cmidrule(lr){1-3} \\cmidrule(lr){4-5} \\cmidrule(lr){6-7} \\cmidrule(lr){8-9}
    Instância (a,b,p,v) & N & M & \$\\chi_h\$ & Tempo (s) & \$\\chi_h\$ & Tempo (s) & \$\\chi_h\$ & Tempo (s) \\\\
    \\midrule
    \\endfirsthead

    \\toprule
    Instância (a,b,p,v) & N & M & \$\\chi_h\$ & Tempo (s) & \$\\chi_h\$ & Tempo (s) & \$\\chi_h\$ & Tempo (s) \\\\
    \\midrule
    \\endhead

    \\midrule
    \\multicolumn{9}{r}{\\textit{Continua na próxima página}} \\\\
    \\bottomrule
    \\endfoot

    \\bottomrule
    \\endlastfoot
"""
    
    # Preenche as linhas da tabela
    for row in results_data
        instance_id = @sprintf("a%d\\_b%d\\_p%d\\%%\\_v%d", row[:a_param], row[:b_param], row[:p_param], row[:v_param])
        
        # Usando \num{} para formatar os tempos com a precisão do siunitx (4 casas)
        latex_output *= @sprintf(
            "%s & %d & %d & %d & \\num{%.4f} & %d & \\num{%.4f} & %d & \\num{%.4f} \\\\\n",
            instance_id,
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

    # Finaliza o ambiente longtable e o documento
    latex_output *= """
\\end{longtable}

\\end{document}
"""
    
    open(filename, "w") do f
        write(f, latex_output)
    end
end


# --- FUNÇÃO PRINCIPAL (EXECUÇÃO) ---

function main()
    # 1. Definição das Heurísticas (Assume que estão no arquivo pai)
    heuristics = Dict(
        :max_deg => coloracaoHarmonicaGrauMax!,
        :min_deg => coloracaoHarmonicaGrauMin!,
        :sat_deg => coloracaoHarmonicaSaturacao! 
    )

    # 2. COLETA DE ARQUIVOS E FILTRAGEM (a <= 100, b <= 100)
    all_bipartite_files = filter(f -> startswith(f, "bi_") && endswith(f, ".col"), readdir())

    filtered_file_names = String[]
    for file_name in all_bipartite_files
        m = match(r"bi_a(\d+)_b(\d+)_p(\d+)%_v(\d+)\.col", file_name)
        
        if m === nothing continue end

        a_param = parse(Int, m.captures[1])
        b_param = parse(Int, m.captures[2])

        # CRITÉRIO DE FILTRAGEM: a <= 100 e b <= 100 (incluindo todos os p e v)
        if a_param <= 100 && b_param <= 100
            push!(filtered_file_names, file_name)
        end
    end

    if isempty(filtered_file_names)
        println("AVISO: Nenhum arquivo bipartido encontrado que satisfaça a<=100 e b<=100. Saindo.")
        return
    end

    sort!(filtered_file_names)
    num_files = length(filtered_file_names)

    println("--- Starting Greedy Algorithm Experiments (1 execution per instance) ---")
    println("Processing $num_files Bipartite files...")

    # 3. EXECUÇÃO DOS EXPERIMENTOS E COLETA DE RESULTADOS
    detailed_results = []
    summary_results = []

    for (i, file_name) in enumerate(filtered_file_names)
        
        @printf "\nProcessing file (%d/%d): %s\n" i num_files file_name

        # Dicionário para coletar os resultados desta instância
        instance_metrics = Dict{Symbol, Any}()
        
        # Armazena N e M para a linha da tabela
        N_v, N_e = leInfo!(file_name)
        
        data_params = extract_bipartite_params(file_name)
        
        best_colors = typemax(Int)
        best_time = typemax(Float64)
        best_heuristic_name = ""

        # Roda CADA heurística apenas UMA VEZ
        for (key, heuristic_func) in heuristics
            
            colors, time = run_greedy_experiment(file_name, heuristic_func)
            
            # Armazenamento detalhado
            instance_metrics[Symbol(string(key) * "_colors")] = colors
            instance_metrics[Symbol(string(key) * "_time")] = time
            
            # Atualização do Melhor Resultado
            if colors < best_colors
                best_colors = colors
                best_time = time
                best_heuristic_name = string(key)
            elseif colors == best_colors && time < best_time
                # Desempate pelo tempo mais rápido
                best_time = time
                best_heuristic_name = string(key)
            end
            
            @printf "  %-10s -> Chi: %d | Time: %.4fs\n" string(key) colors time
        end
        
        # 4. Armazena para as Duas Tabelas
        
        # A. Tabela Detalhada
        push!(detailed_results, Dict(
            :a_param => data_params[:a_param], :b_param => data_params[:b_param], 
            :p_param => data_params[:p_param], :v_param => data_params[:v_param],
            :N_vertices => N_v, :N_edges => N_e,

            :max_colors => instance_metrics[:max_deg_colors], :max_time => instance_metrics[:max_deg_time],
            :min_colors => instance_metrics[:min_deg_colors], :min_time => instance_metrics[:min_deg_time],
            :sat_colors => instance_metrics[:sat_deg_colors], :sat_time => instance_metrics[:sat_deg_time],
        ))

        # B. Tabela Resumo
        push!(summary_results, Dict(
            :a_param => data_params[:a_param], :b_param => data_params[:b_param], 
            :p_param => data_params[:p_param], :v_param => data_params[:v_param],
            :N_vertices => N_v, :N_edges => N_e,
            :best_colors => best_colors,
            :best_time => best_time,
            :best_heuristic => best_heuristic_name
        ))
    end

    # 5. Geração dos Arquivos LaTeX
    output_detailed = "results_Greedy_Bipartite_a100_b100_Detailed.tex"
    generate_detailed_table(detailed_results, output_detailed, "Bipartite")

    output_summary = "results_Greedy_Bipartite_a100_b100_Summary.tex"
    generate_best_summary_table(summary_results, output_summary, "Bipartite")

    println("\n--- Geração Finalizada ---")
    println("Resultados detalhados salvos em: $output_detailed")
    println("Resultados resumo salvos em: $output_summary")
end

# Executa a função principal quando o script é chamado
main()