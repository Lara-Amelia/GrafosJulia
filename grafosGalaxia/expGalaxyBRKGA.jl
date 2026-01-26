using Metaheuristics
using Statistics
using Printf
using DataFrames
using CSV
using ProgressMeter
using Random

include("../colorGul.jl")

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

function main()
    all_files = filter(f -> startswith(f, "galaxy_") && endswith(f, ".col"), readdir())
    sort!(all_files)

    results_main= []
    n_rep = 5
    
    @showprogress 1 "Processando: " for (idx, file) in enumerate(all_files)
        println("\n[$idx/$(length(all_files))] Processing: $file")
        
        t_tot, chi = Float64[], Int[]
        empty!(FITNESS_CACHE)

        global ADJ, V
        num_v, num_a = leInfo!(file)
        ADJ = [Int[] for _ in 1:num_v]
        leArestasLista!(file, ADJ)
        V = num_v

        # Regex para extrair parâmetros do padrão galaxy_aX_bX_cX_vX
        m_a = match(r"a(\d+)", file); a = m_a !== nothing ? parse(Int, m_a.captures[1]) : 0
        m_b = match(r"b(\d+)", file); b = m_b !== nothing ? parse(Int, m_b.captures[1]) : 0
        m_c = match(r"c(\d+)", file); c = m_c !== nothing ? parse(Int, m_c.captures[1]) : 0
        m_v = match(r"v(\d+)", file); v = m_v !== nothing ? parse(Int, m_v.captures[1]) : 0

        for i in 1:n_rep
            elapsed = @elapsed begin
                ch = run_ga_experiment(file)
            end
                push!(t_tot, elapsed)
                push!(chi, ch)
                @printf("   Run %d Done (Chi: %d)\n", i, ch)
        end

        se(v) = round(std(v)/sqrt(n_rep), digits=4)
        mn(v) = round(mean(v), digits=5)

        push!(results_main, Dict(
            :instancia => file, :a => a, :b => b, :c => c, :v => v, :N => num_v, :M => num_a,
            :mean_chi => mean(chi), :se_chi => se(chi),
            :mean_time => mn(t_tot), :se_time => se(t_tot)
        ))
    end

    df = DataFrame(results_main)

    df = df[:, [
        :instancia,
        :a,
        :b,
        :c,
        :v,
        :N,
        :M,
        :mean_chi,
        :se_chi,
        :mean_time,
        :se_time
    ]]

    CSV.write("results_BRKGA.csv", df)

    println("\nCSVs Finalizados.")
end

main()