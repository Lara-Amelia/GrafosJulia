#=
# Importa pacotes necessários
using Statistics
using Printf
using DelimitedFiles
using LinearAlgebra 

# Inclui as funções de coloração (Busca no diretório pai)
# **NOTE:** Certifique-se de que o caminho abaixo está correto.
include("../../colorGul.jl") 

# Variáveis globais (usadas pelas funções de leitura no colorGul.jl)
global matriz_adj = Matrix{Int}(undef, 0, 0)
global num_vertices = 0
global num_arestas = 0

# --- FUNÇÕES AUXILIARES DE LEITURA E PARSE ---

function extract_simple_params(file_name)
    # Padrão: simples_n(N)_p(P)%_v(V).col
    match_result = match(r"simples_n(\d+)_p(\d+)%_v(\d+)\.col", file_name)
    
    if match_result !== nothing
        return Dict(
            :n_param => parse(Int, match_result.captures[1]),
            :p_param => parse(Int, match_result.captures[2]),
            :v_param => parse(Int, match_result.captures[3]),
            :graph_type => "Simples"
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
            # Heurística recebe a matriz de adjacência
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
function generate_best_summary_table(results_data, filename, N_LIMIT_MIN, N_LIMIT_MAX)
    
    title = "Melhor Resultado Global das Heurísticas Gulosas - Grafos Simples (N: $N_LIMIT_MIN - $N_LIMIT_MAX, v1)"
    
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
  \\label{tab:greedy_summary_simple}
  \\begin{tabular}{@{} l ccc l @{}}
    \\toprule
    Instância & N & M & \$\\chi_h\$ (Melhor) & Tempo (s) & Heurística \\\\
    \\midrule
"""
    
    # Preenche as linhas da tabela
    for row in results_data
        
        # Formata o ID da instância
        instance_id = @sprintf("n%d\\_p%d\\%%\\_v%d", row[:n_param], row[:p_param], row[:v_param])
        
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
function generate_detailed_table(results_data, filename, N_LIMIT_MIN, N_LIMIT_MAX)
    
    title = "Resultados Detalhados das Heurísticas Gulosas - Grafos Simples (N: $N_LIMIT_MIN - $N_LIMIT_MAX, v1)"

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
    \\caption{Resultados Detalhados da Coloração Harmônica por Heurística (1 Execução)}\\label{tab:greedy_detailed_simple}\\\\
    \\toprule
    \\multicolumn{3}{@{}l}{Parâmetros do Grafo} & \\multicolumn{2}{c}{Grau Máximo} & \\multicolumn{2}{c}{Grau Mínimo} & \\multicolumn{2}{c @{}}{Saturação} \\\\
    \\cmidrule(lr){1-3} \\cmidrule(lr){4-5} \\cmidrule(lr){6-7} \\cmidrule(lr){8-9}
    Instância (N, p, v) & N & M & \$\\chi_h\$ & Tempo (s) & \$\\chi_h\$ & Tempo (s) & \$\\chi_h\$ & Tempo (s) \\\\
    \\midrule
    \\endfirsthead

    \\toprule
    Instância (N, p, v) & N & M & \$\\chi_h\$ & Tempo (s) & \$\\chi_h\$ & Tempo (s) & \$\\chi_h\$ & Tempo (s) \\\\
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
        instance_id = @sprintf("n%d\\_p%d\\%%\\_v%d", row[:n_param], row[:p_param], row[:v_param])
        
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
    # Limites do filtro de parâmetros
    N_LIMIT_MIN = 100
    N_LIMIT_MAX = 1000
    V_FILTER = 1

    # 1. Definição das Heurísticas (Assume que estão no arquivo pai)
    heuristics = Dict(
        :max_deg => coloracaoHarmonicaGrauMax!,
        :min_deg => coloracaoHarmonicaGrauMin!,
        :sat_deg => NOVOcoloracaoHarmonicaSaturacao! 
    )

    # 2. COLETA DE ARQUIVOS E FILTRAGEM
    all_simple_files = filter(f -> startswith(f, "simples_") && endswith(f, ".col"), readdir())

    filtered_file_names = String[]
    for file_name in all_simple_files
        m = match(r"simples_n(\d+)_p(\d+)%_v(\d+)\.col", file_name)
        
        if m === nothing continue end

        n_param = parse(Int, m.captures[1])
        v_param = parse(Int, m.captures[3])

        # CRITÉRIO DE FILTRAGEM: N entre 100 e 1000 E v = 1
        if (N_LIMIT_MIN <= n_param <= N_LIMIT_MAX) && (v_param == V_FILTER)
            push!(filtered_file_names, file_name)
        end
    end

    if isempty(filtered_file_names)
        @printf "AVISO: Nenhum arquivo Simples encontrado que satisfaça N=[%d,%d] e v=%d. Saindo.\n" N_LIMIT_MIN N_LIMIT_MAX V_FILTER
        return
    end

    sort!(filtered_file_names)
    num_files = length(filtered_file_names)

    @printf "--- Iniciando Experimentos de Algoritmos Gulosos (1 execução por instância) ---\n"
    @printf "Processando %d arquivos Simples (N: %d-%d, v%d)...\n" num_files N_LIMIT_MIN N_LIMIT_MAX V_FILTER

    # 3. EXECUÇÃO DOS EXPERIMENTOS E COLETA DE RESULTADOS
    detailed_results = []
    summary_results = []

    for (i, file_name) in enumerate(filtered_file_names)
        
        @printf "\nProcessando arquivo (%d/%d): %s\n" i num_files file_name

        instance_metrics = Dict{Symbol, Any}()
        
        # Para fins de extração de N e M para a tabela (leInfo! atualiza globais N e M)
        N_v, N_e = leInfo!(file_name) 
        
        data_params = extract_simple_params(file_name)
        
        best_colors = typemax(Int)
        best_time = typemax(Float64)
        best_heuristic_name = ""

        # Roda CADA heurística apenas UMA VEZ
        for (key, heuristic_func) in heuristics
            
            # Garante que o grafo seja lido novamente para cada execução, 
            # já que as funções gulosas modificam a matriz de adjacência (ou estruturas internas).
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
            :n_param => data_params[:n_param], :p_param => data_params[:p_param], 
            :v_param => data_params[:v_param],
            :N_vertices => N_v, :N_edges => N_e,

            :max_colors => instance_metrics[:max_deg_colors], :max_time => instance_metrics[:max_deg_time],
            :min_colors => instance_metrics[:min_deg_colors], :min_time => instance_metrics[:min_deg_time],
            :sat_colors => instance_metrics[:sat_deg_colors], :sat_time => instance_metrics[:sat_deg_time],
        ))

        # B. Tabela Resumo
        push!(summary_results, Dict(
            :n_param => data_params[:n_param], :p_param => data_params[:p_param], 
            :v_param => data_params[:v_param],
            :N_vertices => N_v, :N_edges => N_e,
            :best_colors => best_colors,
            :best_time => best_time,
            :best_heuristic => best_heuristic_name
        ))
    end

    # 5. Geração dos Arquivos LaTeX
    output_detailed = "results_Greedy_Simple_N$(N_LIMIT_MIN)-$(N_LIMIT_MAX)_v$(V_FILTER)_Detailed.tex"
    generate_detailed_table(detailed_results, output_detailed, N_LIMIT_MIN, N_LIMIT_MAX)

    output_summary = "results_Greedy_Simple_N$(N_LIMIT_MIN)-$(N_LIMIT_MAX)_v$(V_FILTER)_Summary.tex"
    generate_best_summary_table(summary_results, output_summary, N_LIMIT_MIN, N_LIMIT_MAX)

    println("\n--- Geração Finalizada ---")
    println("Resultados detalhados salvos em: $output_detailed")
    println("Resultados resumo salvos em: $output_summary")
end

# Executa a função principal quando o script é chamado
main()=#

# Importa pacotes necessários
using Statistics
using Printf
using DelimitedFiles
using LinearAlgebra 

# Inclui as funções de coloração (Busca no diretório pai)
# **NOTE:** Certifique-se de que o caminho abaixo está correto.
include("../../colorGul.jl") 

# Variáveis globais (usadas pelas funções de leitura no colorGul.jl)
global matriz_adj = Matrix{Int}(undef, 0, 0)
global num_vertices = 0
global num_arestas = 0

# --- FUNÇÕES AUXILIARES DE LEITURA E PARSE ---

function extract_simple_params(file_name)
    # Padrão: simples_n(N)_p(P)%_v(V).col
    match_result = match(r"simples_n(\d+)_p(\d+)%_v(\d+)\.col", file_name)
    
    if match_result !== nothing
        return Dict(
            :n_param => parse(Int, match_result.captures[1]),
            :p_param => parse(Int, match_result.captures[2]),
            :v_param => parse(Int, match_result.captures[3]),
            :graph_type => "Simples"
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
            # Heurística recebe a matriz de adjacência
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
function generate_best_summary_table(results_data, filename, N_LIMITS)
    
    N_STR = join(N_LIMITS, ", ") # "100, 500, 1000"
    title = "Melhor Resultado Global das Heurísticas Gulosas - Grafos Simples (N=$N_STR, v1)"
    
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
  \\label{tab:greedy_summary_simple}
  % Usando formato c para N, M, Chi_h para garantir que sejam inteiros limpos.
  \\begin{tabular}{@{} l ccc l @{}}
    \\toprule
    Instância & N & M & \$\\chi_h\$ (Melhor) & Tempo (s) & Heurística \\\\
    \\midrule
"""
    
    # Preenche as linhas da tabela
    for row in results_data
        
        # Formata o ID da instância
        instance_id = @sprintf("n%d\\_p%d\\%%\\_v%d", row[:n_param], row[:p_param], row[:v_param])
        
        # N, M e Chi_h formatados como inteiros (%d)
        # O Tempo é formatado com \num{%.4f}
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
function generate_detailed_table(results_data, filename, N_LIMITS)
    
    N_STR = join(N_LIMITS, ", ") # "100, 500, 1000"
    title = "Resultados Detalhados das Heurísticas Gulosas - Grafos Simples (N=$N_STR, v1)"

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
    \\caption{Resultados Detalhados da Coloração Harmônica por Heurística (1 Execução)}\\label{tab:greedy_detailed_simple}\\\\
    \\toprule
    \\multicolumn{3}{@{}l}{Parâmetros do Grafo} & \\multicolumn{2}{c}{Grau Máximo} & \\multicolumn{2}{c}{Grau Mínimo} & \\multicolumn{2}{c @{}}{Saturação} \\\\
    \\cmidrule(lr){1-3} \\cmidrule(lr){4-5} \\cmidrule(lr){6-7} \\cmidrule(lr){8-9}
    Instância (N, p, v) & N & M & \$\\chi_h\$ & Tempo (s) & \$\\chi_h\$ & Tempo (s) & \$\\chi_h\$ & Tempo (s) \\\\
    \\midrule
    \\endfirsthead

    \\toprule
    Instância (N, p, v) & N & M & \$\\chi_h\$ & Tempo (s) & \$\\chi_h\$ & Tempo (s) & \$\\chi_h\$ & Tempo (s) \\\\
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
        instance_id = @sprintf("n%d\\_p%d\\%%\\_v%d", row[:n_param], row[:p_param], row[:v_param])
        
        # N, M e todas as cores formatadas como inteiros (%d)
        # Os tempos formatados com \num{%.4f}
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
    # Limites do filtro de parâmetros
    # Filtra explicitamente por N=100, 500, 1000
    N_LIMITS = [100, 500, 1000] 
    V_FILTER = 1

    # 1. Definição das Heurísticas (Assume que estão no arquivo pai)
    heuristics = Dict(
        :max_deg => coloracaoHarmonicaGrauMax!,
        :min_deg => coloracaoHarmonicaGrauMin!,
        :sat_deg => NOVOcoloracaoHarmonicaSaturacao! 
    )

    # 2. COLETA DE ARQUIVOS E FILTRAGEM
    # Filtra arquivos que começam com "simples_" e terminam com ".col"
    all_simple_files = filter(f -> startswith(f, "simples_") && endswith(f, ".col"), readdir())

    filtered_file_names = String[]
    for file_name in all_simple_files
        # Padrão: simples_n(N)_p(P)%_v(V).col
        m = match(r"simples_n(\d+)_p(\d+)%_v(\d+)\.col", file_name)
        
        if m === nothing 
            # O nome do arquivo não corresponde ao padrão esperado. 
            continue 
        end

        n_param = parse(Int, m.captures[1])
        v_param = parse(Int, m.captures[3])

        # CRITÉRIO DE FILTRAGEM: N pertence à lista N_LIMITS E v = 1
        if (n_param in N_LIMITS) && (v_param == V_FILTER)
            push!(filtered_file_names, file_name)
        end
    end

    if isempty(filtered_file_names)
        N_STR = join(N_LIMITS, ", ")
        @printf "AVISO: Nenhum arquivo Simples encontrado que satisfaça N={%s} e v=%d. Saindo.\n" N_STR V_FILTER
        return
    end

    sort!(filtered_file_names)
    num_files = length(filtered_file_names)

    N_STR = join(N_LIMITS, ", ")
    @printf "--- Iniciando Experimentos de Algoritmos Gulosos (1 execução por instância) ---\n"
    @printf "Processando %d arquivos Simples (N={%s}, v%d)...\n" num_files N_STR V_FILTER

    # 3. EXECUÇÃO DOS EXPERIMENTOS E COLETA DE RESULTADOS
    detailed_results = []
    summary_results = []

    for (i, file_name) in enumerate(filtered_file_names)
        
        @printf "\nProcessando arquivo (%d/%d): %s\n" i num_files file_name

        instance_metrics = Dict{Symbol, Any}()
        
        N_v, N_e = leInfo!(file_name) 
        
        data_params = extract_simple_params(file_name)
        
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
                best_time = time
                best_heuristic_name = string(key)
            end
            
            @printf "  %-10s -> Chi: %d | Time: %.4fs\n" string(key) colors time
        end
        
        # 4. Armazena para as Duas Tabelas
        
        # A. Tabela Detalhada
        push!(detailed_results, Dict(
            :n_param => data_params[:n_param], :p_param => data_params[:p_param], 
            :v_param => data_params[:v_param],
            :N_vertices => N_v, :N_edges => N_e,

            :max_colors => instance_metrics[:max_deg_colors], :max_time => instance_metrics[:max_deg_time],
            :min_colors => instance_metrics[:min_deg_colors], :min_time => instance_metrics[:min_deg_time],
            :sat_colors => instance_metrics[:sat_deg_colors], :sat_time => instance_metrics[:sat_deg_time],
        ))

        # B. Tabela Resumo
        push!(summary_results, Dict(
            :n_param => data_params[:n_param], :p_param => data_params[:p_param], 
            :v_param => data_params[:v_param],
            :N_vertices => N_v, :N_edges => N_e,
            :best_colors => best_colors,
            :best_time => best_time,
            :best_heuristic => best_heuristic_name
        ))
    end

    # 5. Geração dos Arquivos LaTeX
    N_FILE_STR = join(N_LIMITS, "-") # Para o nome do arquivo: "100-500-1000"
    output_detailed = "results_Greedy_Simple_N$(N_FILE_STR)_v$(V_FILTER)_Detailed.tex"
    generate_detailed_table(detailed_results, output_detailed, N_LIMITS)

    output_summary = "results_Greedy_Simple_N$(N_FILE_STR)_v$(V_FILTER)_Summary.tex"
    generate_best_summary_table(summary_results, output_summary, N_LIMITS)

    println("\n--- Geração Finalizada ---")
    println("Resultados detalhados salvos em: $output_detailed")
    println("Resultados resumo salvos em: $output_summary")
end

# Executa a função principal quando o script é chamado
main()