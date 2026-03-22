# versão do GA utilizando ao máximo os métodos do framework
using Metaheuristics
using LinearAlgebra
using Random
using Statistics
using Printf
using DataFrames
using CSV
using ProgressMeter

import Metaheuristics: initialize!, update_state!, final_stage!, gen_initial_state
import Metaheuristics: SBX_crossover, polynomial_mutation!, create_solution, is_better
import Metaheuristics: reset_to_violated_bounds!

include("../../colorGul.jl")

function fitness_harmonious_coloring(individual::Vector{Float64})
    lista_prioridade = sortperm(individual, rev = true)

    cores_vertices = coloracaoHarmonicaAdjVetAux!(ADJ, lista_prioridade)
    fitness = Float64(maximum(cores_vertices))
    return fitness
end

# struct para que possamos aplicar procedimentos da biblioteca diretamente
struct GraphSwapMutation
    p::Float64 # probabilidade de mutação
end

# implementação da lógica no wrapper da mutação
function Metaheuristics.mutation!(Q::AbstractMatrix{Float64}, mut::GraphSwapMutation)
    n_individuals, n_genes = size(Q)
    
    # "pré-seleção" dos filhos que serão mutados usando a probabilidade do struct
    to_mutate = findall(rand(n_individuals) .< mut.p)

    for i in to_mutate
        v1 = rand(1:n_genes)
        vizinhos = ADJ[v1]
        
        if isempty(vizinhos)
            continue 
        end

        v2 = rand(vizinhos)

        @inbounds begin
            tmp = Q[i, v1]
            Q[i, v1] = Q[i, v2]
            Q[i, v2] = tmp
        end 
    end
end

# setup do BRKGA com operador de mutação personalizado
function MyCustomBRKGA(;
        num_elites = 20,
        num_mutants = 10,
        num_offsprings = 70,
        N = num_elites + num_mutants + num_offsprings,
        bias = 0.7,
        kargs...
    )

    # configuração manual dos componentes - o único que muda é a mutação em si
    initializer = Metaheuristics.RandomInBounds(; N)
    selection   = Metaheuristics.BiasedSelection(num_elites, num_offsprings)
    crossover   = Metaheuristics.BinomialCrossover(p = bias, n_offsprings = 1)
    
    # uso da mutação customizada
    mutation    = GraphSwapMutation(0.5) 
    
    environmental_selection = Metaheuristics.ElitistReplacement()

    # retorn o objeto GA modificado (exatamente como no src da biblioteca)
    return Metaheuristics.GA(;
        initializer,
        selection,
        crossover,
        mutation,
        environmental_selection,
        kargs...
      )
end

# parâmetros para geração dos indivíduos - são permutações de valores reais entre 0 e 1
# representando listas de prioridade na coloração (tentamos encontrar prioridades melhores)
function run_ga_experiment(file_name::String, k_limit::Int, N_pop::Int)
    # instância do BRKGA customizado
    genetico = MyCustomBRKGA(num_elites = 20, num_mutants = 10, num_offsprings = 70, bias = 0.7)

    bounds = [zeros(V) ones(V)]'

    result = optimize(fitness_harmonious_coloring, bounds, genetico)

    # resultados
    #=
    println("\n" * "="^20)
    @show Metaheuristics.minimum(result)
    @show result
    =#

    return Int(Metaheuristics.minimum(result))
end

# talvez passar a leitura do grafo para esse loop, para evitar reler o mesmo grafo várias vezes
function main()
    all_files = filter(f -> startswith(f, "galaxy_") && endswith(f, ".col"), readdir(".."))
    sort!(all_files)

    results_main= []
    n_rep = 5
    K_STAG = 50
    N_POP = 100

    @showprogress 1 "Processando: " for (idx, file) in enumerate(all_files)
        println("\n[$idx/$(length(all_files))] Processing: $file")
        
        filepath = joinpath("..", file)
        t_tot, chi = Float64[], Int[]

        global ADJ, V
        num_v, num_a = leInfo!(filepath)
        ADJ = [Int[] for _ in 1:num_v]
        leArestasLista!(filepath, ADJ)
        V = num_v

        # Regex para extrair parâmetros do padrão galaxy_aX_bX_cX_vX
        m_a = match(r"a(\d+)", file); a = m_a !== nothing ? parse(Int, m_a.captures[1]) : 0
        m_b = match(r"b(\d+)", file); b = m_b !== nothing ? parse(Int, m_b.captures[1]) : 0
        m_c = match(r"c(\d+)", file); c = m_c !== nothing ? parse(Int, m_c.captures[1]) : 0
        m_v = match(r"v(\d+)", file); v = m_v !== nothing ? parse(Int, m_v.captures[1]) : 0

        for i in 1:n_rep
            elapsed = @elapsed begin
                ch = run_ga_experiment(filepath, K_STAG, N_POP)
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

    CSV.write("results_BRKGA_graphSwapMutation.csv", df)

    println("\nCSVs Finalizados.")
end

main()