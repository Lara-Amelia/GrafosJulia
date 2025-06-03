#funções que alteram os argumentos recebidos tem ! no nome, por convenção
#algoritmo de coloração guloso - apenas não colocaremos a mesma cor para vértices adjacentes

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

#le as linhas com infos sobre arestas e popula a matriz de adjacência
function leArestas!(nome_arquivo, matriz_adj)
    u = v = 0
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

# "main"
println("Insira o nome do arquivo a ser lido: ")
nome_arquivo = readline()

num_vertices, num_arestas = leInfo!(nome_arquivo)
println("nro de vertices: $num_vertices")
println("nro de arestas: $num_arestas")

#inicializa a matriz de adjacência de dimensões num_vertices com zeros 
#e depois a inicializa com as infos lida do arquivo sobre arestas 
matriz_adj = zeros(Int, num_vertices, num_vertices)
leArestas!(nome_arquivo, matriz_adj)
println("Matriz de adjacência inicializada com arestas:")
println(matriz_adj)
