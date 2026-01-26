using Metaheuristics
using Statistics
using Printf
using DataFrames
using CSV
using ProgressMeter
using Random

include("../../colorGul.jl")

# referência global constante para utilizar no GA
global ADJ = Vector{Vector{Int}}()
global V = 0

# estrutura para memoização de valores de fitness (via tabela hash)
const FITNESS_CACHE = Dict{UInt64, Float64}()

# hash simples para a permutação (associamos permutações a valores na tabela hash
# de forma a recuperar o valor de fitness, evitando recálcu-los que seriam caros)
# listas de prioridade que são essencialmente iguais geram o mesmo hash e armazenamos seu fitness
# usamos o offset base do algoritmo FNV-1a (para 64 bits) - álgebra A
function hash_perm(p::Vector{Int})
    h = UInt64(1469598103934665603)
    for x in p
        h ⊻= UInt64(x)     # XOR bit a bit
        h *= 1099511628211 # multiplicação pelo primo FNV
    end
    return h
end

# função de fitness utilizando tabela hash (recuperar valores para ordens já calculadas)
# e lista de adjacência (ao invés de matriz de adjacência)
function fitness_harmonious_coloring(individual::Vector{Float64})
    lista_prioridade = sortperm(individual, rev = true)
    h = hash_perm(lista_prioridade)

    # se a chave associada a ordem já está na tabela hash, apenas recuperamos o fitness
    if haskey(FITNESS_CACHE, h)
        return FITNESS_CACHE[h]
    end

    # fazemos a coloração para a ordem sse não encontramos associada a chave na tabela hash
    cores_vertices = coloracaoHarmonicaAdjVetAux!(ADJ, lista_prioridade)
    fitness = maximum(cores_vertices)
    FITNESS_CACHE[h] = fitness
    return fitness
end

# BRKGA usando parâmetros propostos por Mauricio et. aleatoriamente
# o BRKGA poderá ser usado como baseline ou poderemos analisar o seu desempenho
# para tentar novos operadores (em especial de crossover e mutação) para o ga personalizado
#genetico = BRKGA(num_elites = 20, num_mutants = 10, num_offsprings = 70, bias = 0.7)

# NOTA: possíveis alterações são incluir critério de parada por estagnação, alterar tam. da população
function run_ga_experiment(file_name::String)
    genetico = BRKGA(num_elites = 20, num_mutants = 10, num_offsprings = 70, bias = 0.7)

    result = optimize(fitness_harmonious_coloring, [zeros(V) ones(V)], genetico)

    return Int(Metaheuristics.minimum(result))
end

# método para extração de parâmetros do nome do arquivo de entrada
function extract_bipartite_params(file_name)
    match_result = match(r"bi_a(\d+)_b(\d+)_p(\d+)%_v(\d+)\.col", file_name)
    if match_result !== nothing
        return Dict(
            :a_param => parse(Int, match_result.captures[1]),
            :b_param => parse(Int, match_result.captures[2]),
            :p_param => parse(Int, match_result.captures[3]),
            :v_param => parse(Int, match_result.captures[4])
        )
    end
    return nothing
end

function main()
    N_REPETITIONS = 5
    
    TARGET_A = 100
    TARGET_B = 100

    # filtragem de arquivos de entrada de acordo com os parâmetros desejados
    all_files = filter(f -> endswith(f, ".col"), readdir())
    filtered_files = String[]

    for file in all_files
        p = extract_bipartite_params(file)
        if p !== nothing && p[:a_param] == TARGET_A && p[:b_param] == TARGET_B
            push!(filtered_files, file)
        end
    end

    if isempty(filtered_files)
        println("AVISO: Nenhum arquivo corresponde aos critérios (a=$TARGET_A, b=$TARGET_B).")
        return
    end

    sort!(filtered_files)
    results_main = []

    @showprogress 1 "Processando: " for (idx, file) in enumerate(filtered_files)
        println("\n[$idx/$(length(filtered_files))] Processing: $file")

        params = extract_bipartite_params(file)
        t_tot, chi = Float64[], Int[]
        empty!(FITNESS_CACHE)

        global ADJ, V
        num_v, num_a = leInfo!(file)
        ADJ = [Int[] for _ in 1:num_v]
        leArestasLista!(file, ADJ)
        V = num_v

        for i in 1:N_REPETITIONS
            elapsed = @elapsed begin
                ch = run_ga_experiment(file)
            end
            push!(t_tot, elapsed) 
            push!(chi, ch)
            @printf("   Run %d Done (Chi: %d)\n", i, ch)
        end

        se(v) = std(v) / sqrt(length(v))
        meanr(v) = mean(v)

        # dicionário na ordem desejada para a saída do .csv final
        push!(results_main, Dict(
            :a => params[:a_param],
            :b => params[:b_param],
            :N => num_v,
            :p => params[:p_param],
            :M => num_a,
            :v => params[:v_param],
            :mean_time => meanr(t_tot),
            :se_time => se(t_tot),
            :mean_chi => meanr(chi),
            :se_chi => se(chi),
            :instancia => file
        ))

    end

    df_main = DataFrame(results_main)
    col_order = [:a, :b, :N, :p, :M, :v, :mean_time, :se_time, :mean_chi, :se_chi, :instancia]
    select!(df_main, col_order)

    CSV.write("results_BRKGA_a100_b100.csv", df_main)

    println("\n--- Experimentos concluídos ---")
end

main()