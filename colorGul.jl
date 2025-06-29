# funções que alteram os argumentos recebidos tem ! no nome, por convenção
# algoritmo de coloração guloso - apenas não colocaremos a mesma cor para vértices adjacentes
# são usadas no máximo num_vertices cores

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

#coloração é feita segundo uma "lista de prioridade" definida em um vetor (vetor aleatório)
function coloracaoPrioridade!(matriz_adj, cores_vertices, num_vertices, prioridade)
    cores_disponíveis = Vector{Bool}(undef, num_vertices + 1)

    for i in 1:num_vertices
         #originalmente, todas as cores estão disponíveis e se tornarão indisponíveis á medida 
         #em que "descobrimos" que já foram utilizadas e violam alguam restrição
         fill!(cores_disponiveis, true)
         for i in 1:num_vertices
            #obtem o elemento/vértice de prioridade i do vetor de prioridades
            prioritario = prioridade[i]
            for j in 1:num_vertices
                #se encontramos uma adjacência
                if matriz_adj[prioritario, j] == 1
                    if cores_vertices[j] != -1 && cores_vertices[j] <= num_vertices
                        corAdj = cores_vertices[j]
                        #a cor não está mais disponível
                        cores_disponíveis[corAdj] = false
                    end
                end
            end
            for k in 1:num_vertices
                if cores_disponiveis[k] == true 
                cores_vertices[i] = k 
                break
            end
        end
        #CHECAR
        cores_vertices[prioritario] = k
    end
end 

#vertices será um vetor de structs infoVertice
#obtem os graus de cada vertice e os armazena em um vetor graus
#também popula os objetos infoVertice com a flag de sem cor atribuída
function obtemGrauVertice!(matriz_adj, vertices, num_vertices)
    for i in 1:num_vertices
        grauVert = 0 

        for j in 1:num_vertices
            vertices.id = i
            vertices.cor = -1 #flag para indicar que ainda não houve atribuição de cor
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
end

function coloracaoHarmonica(matriz_adj, cores_vertices, num_vertices)

end

println("Insira o nome do arquivo a ser lido: ")
nome_arquivo = readline()

num_vertices, num_arestas = leInfo!(nome_arquivo)
println("nro de vertices: $num_vertices")
println("nro de arestas: $num_arestas")

#inicializa a matriz de adjacência de dimensões num_vertices com zeros
#e depois a inicializa com as infos lida do arquivo sobre arestas
matriz_adj = zeros(Int, num_vertices, num_vertices)
leArestas!(nome_arquivo, matriz_adj)
#println("Matriz de adjacência inicializada com arestas:")
#println(matriz_adj) 

#cria o array com entradas -1 e tam = num_vertices
#-1 indica que ainda não foi atribuída cor para o vértice
cores_vertices = fill(-1, num_vertices)

coloracao!(matriz_adj, cores_vertices, num_vertices)

println("\nResultados da Coloração:")
for i in 1:num_vertices
    println("Cor do vértice $i: $(cores_vertices[i])") 
end

nro_cores = maximum(cores_vertices)
println("O número total de cores usado foi: $nro_cores")
end