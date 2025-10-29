# funções que alteram os argumentos recebidos tem ! no nome, por convenção
# algoritmo de coloração guloso - apenas não colocaremos a mesma cor para vértices adjacentes
# são usadas no máximo num_vertices cores


#println("rodando!!")
#struct para auxiliar na coloração por graus
mutable struct infoVertice
    id::Int
    grau::Int
    cor::Int
end 

#le as informações inciais do grafo: número de vértices e de arestas
function leInfo!(nome_arquivo)
    num_vertices = 0
    num_arestas = 0
    open(nome_arquivo, "r") do arquivo
        for linha in eachline(arquivo)
            linha = strip(linha)
            if startswith(linha, "p")
                partes = split(linha)
                num_vertices = parse(Int, partes[3])
                num_arestas = parse(Int, partes[4])
                break
            end
        end
    end
    return num_vertices, num_arestas
end

# le as linhas com infos sobre arestas e popula a matriz de adjacência
function leArestas!(nome_arquivo, matriz_adj)
    open(nome_arquivo, "r") do arquivo
        for linha in eachline(arquivo)
            linha = strip(linha)
            if startswith(linha, "e")
                partes = split(linha)
                u = parse(Int, partes[2])
                v = parse(Int, partes[3])
                matriz_adj[u, v] = 1 
                matriz_adj[v, u] = 1 
            end
        end
    end
end

#coloração gulosa com prioridade sequencial (1, 2, 3, ...)
function coloracao!(matriz_adj, cores_vertices, num_vertices)

    #o array 'cores_disponiveis_para_vertice' deve ser resetado para TRUE para CADA vértice que estamos colorindo
    #ele deve ter tamanho suficiente para as cores (de 1 a num_vertices, já que o máximo de cores é o nro de vértices em si)
    cores_disponiveis_para_vertice = Vector{Bool}(undef, num_vertices + 1)

    #itera sobre cada vértice para atribuir uma cor
    for i in 1:num_vertices #i é o vértice atual que estamos tentando colorir
        fill!(cores_disponiveis_para_vertice, true)

        #checamos quais cores já foram usadas por vértices adjacentes ao i 
        #e, portanto, não estão mais disponíveis
        for j in 1:num_vertices 
            if matriz_adj[i, j] == 1
                if cores_vertices[j] != -1 && cores_vertices[j] <= num_vertices
                    cor_do_adjacente = cores_vertices[j]
                    cores_disponiveis_para_vertice[cor_do_adjacente] = false
                end
            end
        end

        #atribui a menor cor disponível ao vértice i analisado
        for k in 1:num_vertices 
            if cores_disponiveis_para_vertice[k] == true 
                cores_vertices[i] = k 
                break
            end
        end
    end
end

#também faremos uma coloração seguindo o grau máximo e o grau mínimo dos vértices
#e uma coloração harmônica de fato

#a função recebe um vetor 'prioridade' que define a ordem dos vértices a serem coloridos
function coloracaoPrioridade!(matriz_adj, cores_vertices, num_vertices, prioridade)
    #cores_disponíveis rastreia as cores que podem ser usadas para o vértice atual
    cores_disponiveis = Vector{Bool}(undef, num_vertices + 1)

    #itera sobre o vetor 'prioridade' para colorir os vértices na ordem correta
    for i in 1:num_vertices
        #obtém o vértice a ser colorido com base na lista de prioridade
        vertice_prioritario = prioridade[i]
        
        #reseta o vetor de cores disponíveis para o novo vértice
        fill!(cores_disponiveis, true)

        #checa as cores dos vizinhos e marca as cores como indisponíveis
        for j in 1:num_vertices
            if matriz_adj[vertice_prioritario, j] == 1
                cor_do_adjacente = cores_vertices[j]
                if cor_do_adjacente != -1 && cor_do_adjacente <= num_vertices
                    cores_disponiveis[cor_do_adjacente] = false
                end
            end
        end

        #atribui a menor cor disponível ao vértice prioritário
        for k in 1:num_vertices
            if cores_disponiveis[k] == true 
                cores_vertices[vertice_prioritario] = k 
                break
            end
        end
    end 
end

#vertices será um vetor de structs infoVertice
#obtem os graus de cada vertice e os armazena em um vetor graus
#também popula os objetos infoVertice com a flag de sem cor atribuída
function obtemGrauVertice!(matriz_adj, vertices, num_vertices)
    for i in 1:num_vertices
        grauVert = 0 

        for j in 1:num_vertices
            vertices[i].id = i
            vertices[i].cor = -1 #flag para indicar que ainda não houve atribuição de cor
            if matriz_adj[i, j] == 1
                grauVert += 1 
            end
        end 
        vertices[i].grau = grauVert
    end 
end

#vertices é o vetor de structs infoVertice, que é populado anteriormente chamando
#o método obtemGrauVertice
function coloracaoGrauMax!(matriz_adj, vertices, num_vertices)
    #ordena os infoVertices pelo grau em ordem decrescente
    #o restante da lógica é similar à utilizada na coloração com grau mínimo
    sort!(vertices, by = v -> v.grau, rev = true)    
    cores_disponiveis = Vector{Bool}(undef, num_vertices + 1)

    cores_vertices_resultado = Vector{Int}(undef, num_vertices)
    fill!(cores_vertices_resultado, -1)
    
    for v_info in vertices
        vertex_id = v_info.id
        fill!(cores_disponiveis, true)

        for j in 1:num_vertices
            if matriz_adj[vertex_id, j] == 1 && cores_vertices_resultado[j] != -1
                cor_vizinho = cores_vertices_resultado[j]
                if cor_vizinho >= 1 && cor_vizinho <= num_vertices + 1
                    cores_disponiveis[cor_vizinho] = false
                end
            end
        end

        for k in 1:num_vertices
                if cores_disponiveis[k] == true 
                v_info.cor = k 
                cores_vertices_resultado[vertex_id] = k
                break
            end
        end
    end
    return cores_vertices_resultado
end

#realiza a coloração por ordem crescente dos graus dos vértices
#antes de chamarmos esse método, devemos obter as infos de graus para o vetor vértices
function coloracaoGrauMin!(matriz_adj, vertices, num_vertices)
    #ordena vértices em ordem crescente de grau
    sort!(vertices, by = v -> v.grau)
    cores_disponiveis = Vector{Bool}(undef, num_vertices + 1)

    #como ordenamos o vetor vértices, as cores podem não estar corretamente relacionados a seus vértices
    cores_vertices_resultado = Vector{Int}(undef, num_vertices)
    fill!(cores_vertices_resultado, -1)
    
    for v_info in vertices
        vertex_id = v_info.id
        fill!(cores_disponiveis, true)

        for j in 1:num_vertices
            if matriz_adj[vertex_id, j] == 1 && cores_vertices_resultado[j] != -1
                cor_vizinho = cores_vertices_resultado[j]
                if cor_vizinho >= 1 && cor_vizinho <= num_vertices + 1
                    cores_disponiveis[cor_vizinho] = false
                end
            end
        end

        for k in 1:num_vertices
                if cores_disponiveis[k] == true 
                v_info.cor = k 
                cores_vertices_resultado[vertex_id] = k
                break
            end
        end
    end
    return cores_vertices_resultado
end

#provavelmente precisa de correções
function coloracaoHarmonica!(matriz_adj, cores_vertices, num_vertices)
    cores_arestas_usadas = Set{Tuple{Int, Int}}()

    # Itera sobre cada vértice na ordem sequencial
    for i in 1:num_vertices
        
        # Encontra a menor cor válida para o vértice i
        for k in 1:num_vertices 
            eh_valida_para_aresta = true
            
            # Checa se a cor 'k' é válida para todas as arestas do vértice i
            for j in 1:num_vertices
                if matriz_adj[i, j] == 1
                    cor_vizinho = cores_vertices[j]
                    
                    if cor_vizinho != -1 # Aresta já colorida
                        novo_par_cores = (min(k, cor_vizinho), max(k, cor_vizinho))
                        if novo_par_cores in cores_arestas_usadas
                            eh_valida_para_aresta = false
                            break
                        end
                    end
                end
            end
            
            # Se a cor 'k' é válida, atribui e atualiza o conjunto de cores usadas
            if eh_valida_para_aresta
                cores_vertices[i] = k
                for j in 1:num_vertices
                    if matriz_adj[i, j] == 1 && cores_vertices[j] != -1
                        cor_vizinho = cores_vertices[j]
                        novo_par_cores = (min(k, cor_vizinho), max(k, cor_vizinho))
                        push!(cores_arestas_usadas, novo_par_cores)
                    end
                end
                break
            end
        end
    end
end

function coloracaoHarmonicaPrioridade!(matriz_adj, lista_prioridade::Vector{Int})
    num_vertices = size(matriz_adj, 1)
    cores_vertices = zeros(Int, num_vertices)
    
    #conjunto para armazenar pares de cores de arestas já usados
    cores_arestas_usadas = Set{Tuple{Int, Int}}()
    
    #itera sobre os vértices na ordem de prioridade definida pelo GA
    for v_id in lista_prioridade
        cor_valida = 1
        
        while true
            eh_valida = true
            
            #checa cada vizinho do vértice atual (v_id)
            for neighbor_id in 1:num_vertices
                if matriz_adj[v_id, neighbor_id] == 1
                    cor_vizinho = cores_vertices[neighbor_id]
                    
                    #se o vizinho já foi colorido
                    if cor_vizinho != 0
                        #cria o par de cores da aresta (sempre em ordem crescente)
                        novo_par_cores = (min(cor_valida, cor_vizinho), max(cor_valida, cor_vizinho))
                        
                        #se o par já foi usado, a cor não é válida
                        if novo_par_cores in cores_arestas_usadas
                            eh_valida = false
                            break
                        end
                    end
                end
            end
            
            #se a cor 'cor_valida' não formou nenhum par repetido, ela pode ser usada
            if eh_valida
                cores_vertices[v_id] = cor_valida
                
                #adiciona os novos pares de cores de arestas ao conjunto de usados
                for neighbor_id in 1:num_vertices
                    if matriz_adj[v_id, neighbor_id] == 1
                        cor_vizinho = cores_vertices[neighbor_id]
                        if cor_vizinho != 0
                            novo_par_cores = (min(cor_valida, cor_vizinho), max(cor_valida, cor_vizinho))
                            push!(cores_arestas_usadas, novo_par_cores)
                        end
                    end
                end
                break
            else
                #tenta a próxima cor
                cor_valida += 1
            end
        end
    end
    return cores_vertices
end

function coloracaoHarmonicaGuloso!(matriz_adj, lista_prioridade::Vector{Int})
    num_vertices = size(matriz_adj, 1)
    cores_vertices = zeros(Int, num_vertices)
    cores_arestas_usadas = Set{Tuple{Int, Int}}()
    
    # Itera sobre os vértices na ordem de prioridade
    for v_id in lista_prioridade
        cor_candidata = 1
        
        while true
            eh_valida = true
            
            # Checa se a cor candidata é válida para todos os vizinhos
            for neighbor_id in 1:num_vertices
                if matriz_adj[v_id, neighbor_id] == 1
                    cor_vizinho = cores_vertices[neighbor_id]
                    
                    if cor_vizinho != 0
                        novo_par_cores = (min(cor_candidata, cor_vizinho), max(cor_candidata, cor_vizinho))
                        
                        # Se o par de cores já foi usado, a cor candidata não é válida
                        if novo_par_cores in cores_arestas_usadas
                            eh_valida = false
                            break # Sai do loop de vizinhos
                        end
                    end
                end
            end
            
            # Se a cor candidata não gerou nenhum conflito, ela é a escolhida
            if eh_valida
                cores_vertices[v_id] = cor_candidata
                
                # Adiciona os novos pares de arestas ao conjunto de usados
                for neighbor_id in 1:num_vertices
                    if matriz_adj[v_id, neighbor_id] == 1
                        cor_vizinho = cores_vertices[neighbor_id]
                        if cor_vizinho != 0
                            novo_par_cores = (min(cor_candidata, cor_vizinho), max(cor_candidata, cor_vizinho))
                            push!(cores_arestas_usadas, novo_par_cores)
                        end
                    end
                end
                break # Sai do loop while e passa para o próximo vértice
            else
                # Se não for válida, tenta a próxima cor
                cor_candidata += 1
            end
        end
    end
    return cores_vertices
end

function obtemPrioridadePorGrau(matriz_adj, num_vertices, rev::Bool=false)
    vertices = [infoVertice(i, 0, -1) for i in 1:num_vertices]
    
    # popula os graus (usando a função existente)
    obtemGrauVertice!(matriz_adj, vertices, num_vertices)
    
    # ordena o vetor de structs
    sort!(vertices, by = v -> v.grau, rev = rev)
    
    # extrai apenas os IDs (a lista de prioridade)
    lista_prioridade = [v.id for v in vertices]
    return lista_prioridade
end

function coloracaoHarmonicaGrauMin!(matriz_adj)
    num_vertices = size(matriz_adj, 1)
    lista_prioridade = obtemPrioridadePorGrau(matriz_adj, num_vertices, false)
    return coloracaoHarmonicaGuloso!(matriz_adj, lista_prioridade)
end

function coloracaoHarmonicaGrauMax!(matriz_adj)
    num_vertices = size(matriz_adj, 1)
    lista_prioridade = obtemPrioridadePorGrau(matriz_adj, num_vertices, true)
    return coloracaoHarmonicaGuloso!(matriz_adj, lista_prioridade)
end

function coloracaoHarmonicaSaturacao!(matriz_adj)
    num_vertices = size(matriz_adj, 1)
    cores_vertices = zeros(Int, num_vertices)
    cores_arestas_usadas = Set{Tuple{Int, Int}}()
    
    saturation_degree = [sum(matriz_adj[v, :]) for v in 1:num_vertices]
    degree_orig = copy(saturation_degree) 
    
    for step in 1:num_vertices
        
        best_v = -1
        max_sat = -1
        max_deg = -1
        
        for v in 1:num_vertices
            if cores_vertices[v] == 0 
                current_sat = saturation_degree[v]
                current_deg = degree_orig[v]

                if current_sat > max_sat
                    max_sat = current_sat
                    max_deg = current_deg
                    best_v = v
                elseif current_sat == max_sat && current_deg > max_deg
                    max_deg = current_deg
                    best_v = v
                end
            end
        end
        
        v_id = best_v
        
        if v_id == -1
            break 
        end
        
        cor_candidata = 1
        while true
            eh_valida = true

            for neighbor_id in 1:num_vertices
                if matriz_adj[v_id, neighbor_id] == 1
                    cor_vizinho = cores_vertices[neighbor_id]
                    
                    if cor_candidata == cor_vizinho
                        eh_valida = false
                        break
                    end
                    
                    if cor_vizinho != 0
                        novo_par_cores = (min(cor_candidata, cor_vizinho), max(cor_candidata, cor_vizinho))
                        
                        if novo_par_cores in cores_arestas_usadas
                            eh_valida = false
                            break
                        end
                    end
                end
            end 
            
            if eh_valida
                cores_vertices[v_id] = cor_candidata
                
                # Adiciona os novos pares de arestas ao conjunto de usados
                for neighbor_id in 1:num_vertices
                    if matriz_adj[v_id, neighbor_id] == 1
                        cor_vizinho = cores_vertices[neighbor_id]
                        if cor_vizinho != 0
                            novo_par_cores = (min(cor_candidata, cor_vizinho), max(cor_candidata, cor_vizinho))
                            push!(cores_arestas_usadas, novo_par_cores)
                        end
                    end
                end
                break 
            else
                cor_candidata += 1
            end
        end
        
        for neighbor_id in 1:num_vertices
            if matriz_adj[v_id, neighbor_id] == 1 && cores_vertices[neighbor_id] == 0
                saturation_degree[neighbor_id] -= 1
            end
        end
    end # Fim loop step
    return cores_vertices
end
#println("Insira o nome do arquivo a ser lido: ")
#nome_arquivo = readline()

#num_vertices, num_arestas = leInfo!(nome_arquivo)
#println("nro de vertices: $num_vertices")
#println("nro de arestas: $num_arestas")

#=matriz_adj = zeros(Int, num_vertices, num_vertices)
leArestas!(nome_arquivo, matriz_adj)

#1. Coloração Gulosa Sequencial (1, 2, 3...)
cores_sequencial = fill(-1, num_vertices)
coloracao!(matriz_adj, cores_sequencial, num_vertices)
for i in 1:num_vertices
    println("cor do vertice $i: $(cores_sequencial[i])")
end
nro_cores_seq = maximum(cores_sequencial)
println("\n--- Coloração Sequencial ---")
println("O número total de cores usado foi: $nro_cores_seq")

#provavelmente mudaremos a forma de passar e receber os vetores com resultados das 
#colorações, para evitar uso de memória para fins redundantes
# 2. Coloração por Grau Máximo
#inicialização do vetor de structs infoVertice
vertices_grau_max = [infoVertice(i, 0, -1) for i in 1:num_vertices]
#preenche o vetor com os graus
obtemGrauVertice!(matriz_adj, vertices_grau_max, num_vertices)
#realiza a coloração em si
cores_grau_max = coloracaoGrauMax!(matriz_adj, vertices_grau_max, num_vertices)
nro_cores_max = maximum(cores_grau_max)
println("\n--- Coloração por Grau Máximo ---")
println("O número total de cores usado foi: $nro_cores_max")

# 3. Coloração por Grau Mínimo
#mesmo processo do realizado para obter a coloração por grau máximo
#(criar o vetor, popular, etc)
vertices_grau_min = [infoVertice(i, 0, -1) for i in 1:num_vertices]
obtemGrauVertice!(matriz_adj, vertices_grau_min, num_vertices)
cores_grau_min = coloracaoGrauMin!(matriz_adj, vertices_grau_min, num_vertices)
nro_cores_min = maximum(cores_grau_min)
println("\n--- Coloração por Grau Mínimo ---")
println("O número total de cores usado foi: $nro_cores_min")

#4. Coloração harmônica
cores_aresta_distinguivel = fill(-1, num_vertices)
coloracaoHarmonica!(matriz_adj, cores_aresta_distinguivel, num_vertices)
nro_cores_aresta_dist = maximum(cores_aresta_distinguivel)
println("\n--- Coloração Distinguível por Arestas ---")
for i in 1:num_vertices
    println("cor usada pelo vertice $i: $(cores_aresta_distinguivel[i])")
end
println("O número total de cores usado foi: $nro_cores_aresta_dist") =#