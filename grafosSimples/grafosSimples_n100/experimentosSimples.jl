# Importa pacotes necessários
using Metaheuristics
using LinearAlgebra
using Statistics
using Printf

# Inclui as funções de coloração do arquivo pai (colorGul.jl)
# Certifique-se de que o colorGul.jl está um nível acima deste diretório.
include("../../colorGul.jl")

# Variáveis globais para armazenar a matriz de adjacência e o número de vértices
# Estas são necessárias porque as funções Metaheuristics chamam fitness_harmonious_coloring
# no escopo global e precisam acessar a matriz e o tamanho do grafo.
global matriz_adj = Matrix{Int}(undef, 0, 0)
global num_vertices = 0

# --- ESTRUTURAS DO ALGORITMO GENÉTICO (GA) ---

# Implementação do CustomGAParams (parâmetros do GA)
mutable struct CustomGAParams <: Metaheuristics.AbstractParameters
    N::Int
    p_crossover::Float64
    p_mutation::Float64
    stag_limit::Int        
    last_best::Float64     # Alterado para Float64
    stag_iters::Int 
end

# Construtor com valores padrão, usando 'k' (stag_limit) da função principal
CustomGAParams(; N = 200, p_crossover = 0.5, p_mutation = 0.5, stag_limit = 50, 
                 last_best = Inf, stag_iters = 0) =
    CustomGAParams(N, p_crossover, p_mutation, stag_limit, last_best, stag_iters)

# Função de fitness: usa a ordem do indivíduo para colorir harmonicamente o grafo global
function fitness_harmonious_coloring(individual::Vector{Float64})
    # O GA usa o valor do indivíduo como prioridade. Sortperm decodifica essa prioridade.
    lista_prioridade = sortperm(individual, rev = true)
    
    # É CRUCIAL que coloracaoHarmonicaGuloso! utilize a matriz_adj global
    global matriz_adj
    cores_vertices = coloracaoHarmonicaGuloso!(matriz_adj, lista_prioridade)
    
    # Retorna o número máximo de cores (o objetivo é minimizar isso)
    return float(maximum(cores_vertices)) 
end

# Função de mutação (graph_swap_mutation! adaptada do seu gaUpdate.jl)
function graph_swap_mutation!(x::Vector{Float64})
    global matriz_adj
    n = length(x)
    if n < 2 return x end
    
    v1 = rand(1:n)
    # Obtém vizinhos (usando a matriz_adj global)
    vizinhos = [i for i in 1:n if matriz_adj[v1, i] == 1]
    
    if isempty(vizinhos)
        return x
    end
    
    v2 = rand(vizinhos)
    x[v1], x[v2] = x[v2], x[v1]
    return x
end

# Hook de inicialização do Metaheuristics
function Metaheuristics.initialize!(
    status,
    parameters::CustomGAParams,
    problem::Metaheuristics.Problem,
    information::Metaheuristics.Information,
    options::Metaheuristics.Options
)
    return Metaheuristics.gen_initial_state(problem, parameters, information, options, status)
end

# Hook de atualização de estado do Metaheuristics (lógica do GA)
function Metaheuristics.update_state!(
    state::Metaheuristics.State,
    parameters::CustomGAParams,
    problem::Metaheuristics.Problem,
    information::Metaheuristics.Information,
    options::Metaheuristics.Options
)
    # Implementa tournament selection (seleção de pais)
    function tournament_select(pop)
        k = 2
        candidates = rand(pop, k)
        candidate_vectors = [c.x for c in candidates]
        fitnesses = [fitness_harmonious_coloring(v) for v in candidate_vectors]
        return candidates[argmin(fitnesses)].x
    end

    p1 = tournament_select(state.population)
    p2 = tournament_select(state.population)
    parents = vcat(p1', p2')  

    # Crossover (média simples e randomização)
    p1 = parents[1, :]
    p2 = parents[2, :]
    child1 = (p1 .+ p2) ./ 2.0
    child2 = (p1 .+ p2) ./ 2.0 .+ 0.05 .* randn(length(p1))
    offsprings = [child1'; child2']

    # Mutação (swap nos vizinhos)
    for i in axes(offsprings, 1)  
        graph_swap_mutation!(offsprings[i, :])
    end

    # Cria objetos solução para os filhos
    offspring_solutions = Metaheuristics.xf_solution[]
    for i in axes(offsprings, 1)
        vec = offsprings[i, :]
        fit = fitness_harmonious_coloring(vec)
        push!(offspring_solutions, Metaheuristics.xf_solution(vec, fit))
    end

    # Elitismo: Mantém o melhor indivíduo da população atual
    best_ind = findmin((s.f for s in state.population))[2]
    best_ind = state.population[best_ind]

    # Substitui o restante da população pelos novos filhos
    state.population = [best_ind; offspring_solutions[2:end]]

    # Checa estagnação (parada)
    melhor_atual = minimum([s.f for s in state.population]) 
    if melhor_atual < parameters.last_best
        parameters.last_best = melhor_atual
        parameters.stag_iters = 0
    else
        parameters.stag_iters += 1
    end

    # Condição de Parada
    if parameters.stag_iters >= parameters.stag_limit
        @info "Parada: $(parameters.stag_limit) iterações sem melhora"
        return false
    end

    return true
end

# Hook de estágio final
function Metaheuristics.final_stage!(state, parameters::CustomGAParams, problem, information, options)
    return state
end


# --- FUNÇÃO CENTRAL DE EXECUÇÃO DO EXPERIMENTO ---

"""
run_ga_experiment(file_name, k_limit, N_pop)
Executa o GA 1 vez para um dado arquivo e retorna o melhor fitness e tempo.
"""
function run_ga_experiment(file_name::String, k_limit::Int, N_pop::Int)
    global matriz_adj
    global num_vertices

    # 1. LEITURA E PREPARAÇÃO DO GRAFO (Atualiza globais)
    num_vertices, num_arestas = leInfo!(file_name)
    matriz_adj = zeros(Int, num_vertices, num_vertices)
    leArestas!(file_name, matriz_adj)
    
    # 2. CONFIGURAÇÃO DO GA
    bounds = [zeros(num_vertices) ones(num_vertices)]'

    opts = Metaheuristics.Options(f_calls_limit = typemax(Int), iterations = typemax(Int))
    
    params = CustomGAParams(N = N_pop, p_crossover = 0.5, p_mutation = 0.5, 
                            stag_limit = k_limit, last_best = Inf, stag_iters = 0)
    
    my_custom_ga = Metaheuristics.Algorithm(params, options = opts) 
    
    problem = Metaheuristics.Problem(fitness_harmonious_coloring, bounds)

    # 3. EXECUÇÃO DA OTIMIZAÇÃO E MEDIÇÃO DO TEMPO
    start_time = time()
    result = Metaheuristics.optimize(fitness_harmonious_coloring, bounds, my_custom_ga)
    end_time = time()
    
    time_taken = end_time - start_time
    best_colors = Metaheuristics.minimum(result)

    return time_taken, best_colors
end


# --- FUNÇÃO PRINCIPAL ---

function main()
    # --- PARÂMETROS DO EXPERIMENTO ---
    N_REPETITIONS = 5
    K_LIMIT = 50      # Limite de estagnação (k)
    N_POP = 200       # Tamanho da população (N)

    # --- COLETA DE ARQUIVOS E FILTRAGEM ---
    
    # 1. Busca por todos os arquivos simples na pasta atual
    all_simple_files = filter(f -> startswith(f, "simples_") && endswith(f, ".col"), readdir())

    if isempty(all_simple_files)
        println("ERRO: Nenhum arquivo 'simples_*.col' encontrado no diretório atual.")
        return
    end

    # 2. Filtragem: N (vertices) entre 100 e 1000 E versão v1
    filtered_file_names = String[]
    for file_name in all_simple_files
        # Padrão: simples_nXXX_pYY%_vZ.col
        m = match(r"simples_n(\d+)_p(\d+)%_v(\d+)\.col", file_name)
        
        if m === nothing
            continue
        end

        num_vertices_in_name = parse(Int, m.captures[1])
        v_param = parse(Int, m.captures[3])

        if (100 <= num_vertices_in_name <= 1000) && (v_param == 1)
            push!(filtered_file_names, file_name)
        end
    end

    if isempty(filtered_file_names)
        println("AVISO: Nenhum arquivo Simples encontrado que satisfaça 100<=N<=1000 e v=1.")
        return
    end

    # 3. Estrutura para armazenar todos os resultados
    simple_results = []
    num_files = length(filtered_file_names)

    println("--- Starting Genetic Algorithm Experiments (k=$(K_LIMIT), N=$(N_POP), $(N_REPETITIONS) repetitions) ---")
    println("Processing $(num_files) Simple files (N: 100-1000, v1)...")

    # --- EXECUÇÃO DOS EXPERIMENTOS ---
    for (i, file_name) in enumerate(filtered_file_names)
        
        # Array para armazenar 5 resultados para este arquivo
        colors = Float64[]
        times = Float64[]
        
        println("\nProcessing file ($(i)/$(num_files)): $(file_name)")

        # Executa 5 repetições
        for rep in 1:N_REPETITIONS
            try
                time_taken, best_colors = run_ga_experiment(file_name, K_LIMIT, N_POP)
                push!(colors, best_colors)
                push!(times, time_taken)
                @printf "  Run %d/%d... Done. Colors: %.2f, Time: %.3fs\n" rep N_REPETITIONS best_colors time_taken
            catch e
                @printf "  Run %d/%d... ERROR: %s\n" rep N_REPETITIONS sprint(showerror, e)
            end
        end

        # --- CÁLCULO DAS ESTATÍSTICAS ---
        if isempty(colors)
            println("AVISO: Nenhuma execução bem-sucedida para $(file_name). Pulando cálculo de estatísticas.")
            continue
        end

        mean_colors = mean(colors)
        std_err_colors = std(colors) / sqrt(length(colors)) # Erro padrão (SE)
        mean_time = mean(times)
        std_time = std(times) # Desvio Padrão (sigma)
        
        # Extrai parâmetros para a tabela
        m = match(r"simples_n(\d+)_p(\d+)%_v(\d+)\.col", file_name)
        num_vertices_curr, num_arestas_curr = leInfo!(file_name) 
        
        data = (
            file_name = file_name,
            n_param = parse(Int, m.captures[1]),
            p_param = parse(Int, m.captures[2]),
            v_param = parse(Int, m.captures[3]),
            num_vertices = num_vertices_curr,
            num_arestas = num_arestas_curr,
            mean_colors = mean_colors,
            std_err_colors = std_err_colors,
            mean_time = mean_time,
            std_time = std_time,
        )

        # Imprime o resumo no terminal
        @printf "\nRESULTS for %s:\n" file_name
        @printf "  Mean Colors (\\bar{\\chi}_h): %.2f \\pm %.2f\n" mean_colors std_err_colors
        @printf "  Mean Time (s): %.3f \\pm %.3f\n" mean_time std_time

        # Adiciona à lista de resultados
        push!(simple_results, data)
    end
    
    # --- GERAÇÃO DA TABELA LATEX ---
    
    # Ordena os resultados para agrupar instâncias semelhantes
    sort!(simple_results, by = r -> (r.n_param, r.p_param, r.v_param))

    output_file = "results_Simple_N100-N1000_v1.tex"
    generate_latex_table(simple_results, output_file)

    println("\n--- Generating LaTeX Table for Simple Graphs ---")
    println("LaTeX table saved to $(output_file)")
end

"""
generate_latex_table(results, filepath)
Cria o arquivo LaTeX (.tex) com os resultados formatados.
"""
function generate_latex_table(results, filepath)
    
    # 1. Cabeçalho LaTeX (CORRIGIDO PARA O FILTRO N100-N1000)
    latex_output = """
\\documentclass{article}
\\usepackage[a4paper, margin=1in]{geometry}
\\usepackage{booktabs}
\\usepackage{amsmath}
\\usepackage{siunitx}
\\sisetup{round-mode = places, round-precision = 2}
\\title{Resultados do Algoritmo Genético em Grafos Simples (100 \\\\le N \\\\le 1000, v1)}
\\author{}
\\date{}

\\begin{document}
\\maketitle

\\begin{table}[htbp]
    \\centering
    \\caption{Desempenho da Coloração Harmônica (Média de 5 Execuções)}
    \\label{tab:simple_results}
    \\begin{tabular}{@{} l ccc ccc @{}}
        \\toprule
        Instância (N, p, v) & \$\\text{N}_{total}\$ & \$\\text{M}\$ & \$\\bar{\\chi}_h \\pm \\text{SE}\$ & \$\\bar{\\text{Tempo (s)}} \\pm \\sigma\$ \\\\
        \\midrule
    """
    
    # 2. Linhas de Dados
    for r in results
        
        # Formata a string do nome da instância (nXXX_pYY%_vZ)
        instance_id = "n$(r.n_param)\\_p$(r.p_param)\\%\\_v$(r.v_param)"
        
        # Formatação dos resultados
        colors_fmt = @sprintf "\\num{%.2f} \\pm \\num{%.2f}" r.mean_colors r.std_err_colors
        time_fmt = @sprintf "\\num{%.3f} \\pm \\num{%.3f}" r.mean_time r.std_time
        
        # Note que N_total é r.num_vertices, que deve ser igual a r.n_param
        latex_output *= @sprintf "        %s & %d & %d & %s & %s \\\\\n" \
            instance_id r.num_vertices r.num_arestas colors_fmt time_fmt
    end
    
    # 3. Rodapé LaTeX
    latex_output *= """
        \\bottomrule
    \\end{tabular}
\\end{table}

\\end{document}
"""

    open(filepath, "w") do f
        write(f, latex_output)
    end
end
# Executa a função principal quando o script é chamado
main()