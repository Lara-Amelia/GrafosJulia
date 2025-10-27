using Metaheuristics
using LinearAlgebra
using Statistics
using Printf
using DelimitedFiles

# Inclui as funções de coloração (busca no diretório pai, como solicitado)
include("../colorGul.jl")

# --- ESTRUTURAS E MÉTODOS GA PERSONALIZADOS ---

# Definição dos parâmetros do Algoritmo Genético personalizado
mutable struct CustomGAParams <: Metaheuristics.AbstractParameters
    N::Int                  # Tamanho da População
    p_crossover::Float64
    p_mutation::Float64
    stag_limit::Int         # Limite de Estagnação (k)
    last_best::Int
    stag_iters::Int
end

CustomGAParams(; N = 100, p_crossover = 0.5, p_mutation = 0.5, stag_limit = 50,
                 last_best = -1, stag_iters = 0) =
    CustomGAParams(N, p_crossover, p_mutation, stag_limit, last_best, stag_iters)

# Função para crossover (gerando 2 filhos)
function crossover_harmonious_coloring(parents::Matrix{Float64})
    p1 = parents[1, :]
    p2 = parents[2, :]
    child1 = (p1 .+ p2) ./ 2.0
    # Randomização para não gerar 2 filhos exatamente iguais
    child2 = (p1 .+ p2) ./ 2.0 .+ 0.05 .* randn(length(p1))
    return [child1'; child2']
end

# Função para mutação (usando global matriz_adj, como no seu gaUpdate.jl)
function graph_swap_mutation!(x::Vector{Float64})
    global matriz_adj # Acessa a matriz global do grafo atual
    n = length(x)
    if n < 2 return x end
    
    # Escolhe um vértice/posição aleatoriamente
    v1 = rand(1:n)
    
    # Obtem os vizinhos do vértice escolhido (depende de matriz_adj)
    vizinhos = [i for i in 1:n if matriz_adj[v1, i] == 1]
    
    if isempty(vizinhos)
        return x
    end
    
    # Seleciona aleatoriamente algum dos vizinhos e troca
    v2 = rand(vizinhos)
    x[v1], x[v2] = x[v2], x[v1]
    return x
end

# Função de atualização de estado do Metaheuristics
function Metaheuristics.update_state!(
    state::Metaheuristics.State,
    parameters::CustomGAParams,
    problem::Metaheuristics.Problem,
    information::Metaheuristics.Information,
    options::Metaheuristics.Options
)
    global matriz_adj 
    
    function tournament_select(pop, fitness_func)
        k = 2  
        candidates = rand(pop, k)                      
        candidate_vectors = [c.x for c in candidates]  
        fitnesses = [fitness_func(v) for v in candidate_vectors]
        return candidates[argmin(fitnesses)].x         
    end

    p1 = tournament_select(state.population, problem.f)
    p2 = tournament_select(state.population, problem.f)
    parents = vcat(p1', p2')  
    offsprings = crossover_harmonious_coloring(parents)

    for i in axes(offsprings, 1)  
        # Chama a mutação sem argumentos, confiando no escopo global (como no seu código)
        graph_swap_mutation!(offsprings[i, :]) 
    end

    offspring_solutions = Metaheuristics.xf_solution[]
    for i in axes(offsprings, 1)
        vec = offsprings[i, :]
        fit = Float64(problem.f(vec))
        push!(offspring_solutions, Metaheuristics.xf_solution(vec, fit))
    end

    best_ind_idx = findmin((s.f for s in state.population))[2]
    best_ind = state.population[best_ind_idx]

    if isempty(offspring_solutions)
        state.population = [best_ind]
    else
        current_pop = vcat(state.population, offspring_solutions)
        sort!(current_pop, by = s -> s.f)
        N = parameters.N
        state.population = current_pop[1:min(N, length(current_pop))]
    end

    melhor_atual = minimum([s.f for s in state.population]) 
    if melhor_atual < parameters.last_best
        parameters.last_best = melhor_atual
        parameters.stag_iters = 0
    else
        parameters.stag_iters += 1
    end

    # Condição de Parada por Estagnação
    if parameters.stag_iters >= parameters.stag_limit
        # FIX: Removida a referência a information.f_calls que causava o FieldError
        @info "Stop condition: $(parameters.stag_limit) iterations without improvement." 
        return false
    end

    return true
end

function Metaheuristics.initialize!(
    status,
    parameters::CustomGAParams,
    problem::Metaheuristics.Problem,
    information::Metaheuristics.Information,
    options::Metaheuristics.Options
)
    initial_state = Metaheuristics.gen_initial_state(problem, parameters, information, options, status)
    
    population = initial_state.population
    parameters.last_best = minimum([s.f for s in population])
    parameters.stag_iters = 0

    return initial_state
end

function Metaheuristics.final_stage!(state, parameters::CustomGAParams, problem, information, options)
    return state
end


# --- LÓGICA DO EXPERIMENTO ---

# Variáveis globais necessárias para que as funções do GA (fitness e mutação) possam acessar o grafo.
global matriz_adj = Matrix{Int}(undef, 0, 0)
global num_vertices = 0


"""
    run_ga_experiment(file_name, k_limit, N_pop)

Executa o algoritmo genético uma vez para um dado grafo e parâmetros.
Retorna (tempo_decorrido, melhor_cores_encontradas).
"""
function run_ga_experiment(file_name, k_limit, N_pop)
    global matriz_adj
    global num_vertices

    try
        num_vertices, num_arestas = leInfo!(file_name) 
    catch e
        @error "Não foi possível ler as informações de $file_name: $e"
        return NaN, NaN
    end

    matriz_adj = zeros(Int, num_vertices, num_vertices)
    leArestas!(file_name, matriz_adj)
    
    # Função de fitness que usa matriz_adj e num_vertices globalmente
    function fitness_harmonious_coloring(individual::Vector{Float64})
        lista_prioridade = sortperm(individual, rev = true) 
        cores_vertices = coloracaoHarmonicaGuloso!(matriz_adj, lista_prioridade)
        return maximum(cores_vertices)
    end
    
    bounds = [zeros(num_vertices) ones(num_vertices)]'
    opts = Metaheuristics.Options(f_calls_limit = typemax(Int), iterations = typemax(Int))
    
    # Cria os parâmetros do GA
    params = CustomGAParams(N = N_pop, p_crossover = 0.5, p_mutation = 0.5, 
                            stag_limit = k_limit, last_best = -1, stag_iters = 0)
    my_custom_ga = Metaheuristics.Algorithm(params, options = opts) 
    
    problem = Metaheuristics.Problem(fitness_harmonious_coloring, bounds)

    start_time = time()
    result = Metaheuristics.optimize(fitness_harmonious_coloring, bounds, my_custom_ga)
    end_time = time()

    time_taken = end_time - start_time
    best_colors = Metaheuristics.minimum(result)

    return time_taken, Int(round(best_colors))
end


"""
    parse_filename_params(file_name, N_v, N_e)

Extrai os parâmetros do grafo (a, b, c) a partir do nome do arquivo.
"""
function parse_filename_params(file_name, N_v, N_e)
    
    match_a = match(r"a(\d+)", file_name)
    match_b = match(r"b(\d+)", file_name)
    match_c = match(r"c(\d+)", file_name)
    
    a_param = match_a !== nothing ? parse(Int, match_a.captures[1]) : "N/A"
    b_param = match_b !== nothing ? parse(Int, match_b.captures[1]) : "N/A"
    c_param = match_c !== nothing ? parse(Int, match_c.captures[1]) : "N/A"
    
    size_str = "$N_v \\times $N_e" 
    
    return a_param, b_param, c_param, size_str
end


"""
    generate_latex_table(results_subset)

Gera a string da tabela LaTeX.
"""
function generate_latex_table(results_subset)
    
    sort!(results_subset, by = x -> (x[:a_param], x[:b_param], x[:c_param]))
    
    header = """
\\documentclass[a4paper]{article}
\\usepackage[utf8]{inputenc}
\\usepackage[T1]{fontenc}
\\usepackage{amsmath}
\\usepackage{booktabs}
\\usepackage{siunitx}
\\usepackage[margin=1in]{geometry}
\\sisetup{round-mode=figures, round-precision=2}

\\begin{document}

\\section*{Resultados do Algoritmo Genético para Coloração Harmônica - Grafos Tipo \\textbf{Galaxy}}

\\begin{table}[h!]
    \\centering
    \\caption{ Resultados da Execução do GA (N_{pop}=100, k=50) em Grafos Tipo \\textbf{Galaxy} (\$a \\le 25, v=1\$) (5 Repetições) }
    \\label{tab:ga_results_galaxy_v1}
    \\begin{tabular}{@{} l c S[table-format=2.2] S[table-format=1.2] S[table-format=1.2] @{}}
        \\toprule
        Instância (a, b, c) & Parâmetros (N, M) & {Média de Cores (\\si{\\bar{\\chi}_h})} & {Erro Padrão (SE)} & {Tempo Médio (s)} \\\\
        \\midrule
"""
    
    body = ""
    for data in results_subset
        instance_id = "a$(data[:a_param])\\_b$(data[:b_param])\\_c$(data[:c_param])\\_v$(data[:v_param])"
        size_str = data[:size_str]
        mean_colors = data[:mean_colors]
        se_colors = data[:se_colors]
        mean_time = data[:mean_time]
        
        row = Printf.format(
            Printf.Format("        %s & %s & %.2f & %.2f & %.2f \\\\ \n"),
            instance_id,
            size_str,
            mean_colors,
            se_colors,
            mean_time
        )
        
        body *= row
    end
    
    footer = """
        \\bottomrule
    \\end{tabular}
\\end{table}

\\end{document}
"""
    
    return header * body * footer
end


# --- EXECUÇÃO PRINCIPAL ---
function main()
    # 1. ENCONTRA E FILTRA ARQUIVOS
    
    all_galaxy_files = filter(f -> startswith(f, "galaxy_") && endswith(f, ".col"), readdir())
    
    filtered_file_names = String[]
    
    println("--- Filtering Galaxy Files (Stopping at a=25, only v1) ---")

    for file_name in all_galaxy_files
        match_a = match(r"a(\d+)", file_name)
        match_v = match(r"v(\d+)", file_name) 
        
        if match_a !== nothing && match_v !== nothing
            a_param = parse(Int, match_a.captures[1])
            v_param = parse(Int, match_v.captures[1])
            
            if a_param <= 25 && v_param == 1
                push!(filtered_file_names, file_name)
            else
                if a_param > 25
                    println("Skipping $file_name (a=$a_param > 25)")
                elseif v_param != 1
                    println("Skipping $file_name (v=$v_param, only v1 is requested)")
                else
                    println("Skipping $file_name (Unknown reason/combination)")
                end
            end
        else
            println("Warning: Could not extract 'a' or 'v' parameter from $file_name. Skipping.")
        end
    end
    
    if isempty(filtered_file_names)
        println("\nNo Galaxy files found or none matched the a <= 25 and v1 filter. Exiting.")
        return
    end

    sort!(filtered_file_names)

    # 2. PARÂMETROS DO EXPERIMENTO
    N_REPETITIONS = 5
    K_LIMIT = 50      # Limite de Estagnação (k)
    N_POP = 100       # Tamanho da População (N)

    galaxy_results = []
    
    println("\n--- Starting Genetic Algorithm Experiments (k=$K_LIMIT, N=$N_POP, $N_REPETITIONS repetitions) ---")
    println("Processing $(length(filtered_file_names)) Galaxy files...")

    for file_name in filtered_file_names
        println("\nProcessing file: $file_name")
        
        color_results = Int[]
        time_results = Float64[]

        num_vertices, num_arestas = leInfo!(file_name)
        a, b, c, size_str = parse_filename_params(file_name, num_vertices, num_arestas)
        match_v = match(r"v(\d+)", file_name)
        v_param = match_v !== nothing ? parse(Int, match_v.captures[1]) : "N/A"
        
        for i in 1:N_REPETITIONS
            print("  Run $i/$N_REPETITIONS...")
            
            # Executa o experimento
            time_taken, best_colors = run_ga_experiment(file_name, K_LIMIT, N_POP)
            
            if isnan(time_taken)
                println(" ERROR. Skipping file.")
                break
            end
            
            push!(time_results, time_taken)
            push!(color_results, best_colors)
            println(" Done. Colors: $best_colors, Time: $(@sprintf("%.2f", time_taken))s")
        end
        
        if length(color_results) == N_REPETITIONS
            # Calcula métricas
            mean_colors = mean(color_results)
            mean_time = mean(time_results)
            std_dev_colors = std(color_results)
            se_colors = std_dev_colors / sqrt(N_REPETITIONS)
            
            # Armazena resultados
            file_result_data = Dict(
                :a_param => a,
                :b_param => b,
                :c_param => c,
                :v_param => v_param,
                :size_str => size_str,
                :mean_colors => mean_colors,
                :se_colors => se_colors,
                :mean_time => mean_time
            )
            push!(galaxy_results, file_result_data)
            
            println("\nRESULTS for $file_name:")
            println("  Mean Colors (\\bar{\\chi}_h): $(@sprintf("%.2f", mean_colors)) ± $(@sprintf("%.2f", se_colors))")
            println("  Mean Time (s): $(@sprintf("%.2f", mean_time))")
        end
    end

    # 3. GERA E SALVA A TABELA
    println("\n--- Generating LaTeX Table for Galaxy Graphs ---")
    latex_output = generate_latex_table(galaxy_results)
    
    output_filename = "results_Galaxy_a01-a25_v1.tex"
    open(output_filename, "w") do f
        write(f, latex_output)
    end
    println("LaTeX table saved to $output_filename")
end

# Chamada da função principal
main()