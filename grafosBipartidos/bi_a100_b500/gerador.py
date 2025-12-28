import networkx as nx
import os

# --- Parâmetros de Execução ---
SIZES = [100, 500, 1000, 5000]
P_VALUES = [0.01, 0.03, 0.05, 0.10, 0.20, 0.30, 0.40, 0.50, 0.60, 0.70, 0.80, 0.90]
NUM_VERSIONS = 5 

print(f"Iniciando a verificação e geração de grafos bipartidos...")

for i in range(len(SIZES)):
    a = SIZES[i]
    for j in range(i, len(SIZES)):
        b = SIZES[j]
        for p_value in P_VALUES:
            for v in range(1, NUM_VERSIONS + 1):
                
                # 1. Formata o nome do arquivo exatamente como o padrão existente
                p_percent = int(p_value * 100) 
                file_name = f"bi_a{a}_b{b}_p{p_percent:02d}%_v{v}.col"

                # 2. VERIFICAÇÃO DE EXISTÊNCIA
                if os.path.exists(file_name):
                    # Opcional: print(f"Pulando: {file_name} já existe.")
                    continue

                # 3. Geração apenas se não existir
                try:
                    G = nx.bipartite.random_graph(a, b, p_value)
                    
                    with open(file_name, 'w') as f:
                        num_nodes = G.number_of_nodes()
                        num_edges = G.number_of_edges()
                        f.write(f"p edge {num_nodes} {num_edges}\n")
                        for u, v in G.edges():
                            f.write(f"e {u+1} {v+1}\n")

                    print(f"Gerado: {file_name} (Nós: {num_nodes}, Arestas: {num_edges})")
                
                except MemoryError:
                    print(f"ERRO DE MEMÓRIA: Grafo bi_a{a}_b{b} com p={p_value} é muito grande para a RAM.")
                except Exception as e:
                    print(f"Erro inesperado em {file_name}: {e}")

print("Processo concluído.")

"""
# import networkx as nx

# --- Parâmetros de Execução para Grafos Bipartidos ---

# Tamanhos das partições (a e b) a serem testados
SIZES = [100, 500, 1000, 5000]
# Lista de probabilidades de aresta a serem testadas
P_VALUES = [0.01, 0.03, 0.05, 0.10, 0.20, 0.30, 0.40, 0.50, 0.60, 0.70, 0.80, 0.90]
NUM_VERSIONS = 5 # Gerar 5 versões para cada combinação (a, b, p)

print(f"Iniciando a geração de grafos bipartidos...")

# LAÇO 1: Itera sobre todos os possíveis tamanhos para a partição 'a'
for i in range(len(SIZES)):
    a = SIZES[i]
    
    # LAÇO 2: Itera sobre todos os tamanhos b >= a (evita repetição de (b, a))
    for j in range(i, len(SIZES)):
        b = SIZES[j]
        
        # LAÇO 3: Itera sobre cada probabilidade (p) (AQUI ESTÁ A MUDANÇA)
        for p_value in P_VALUES:
            
            # LAÇO 4: Gera as 5 versões
            for v in range(1, NUM_VERSIONS + 1):
                
                # 1. Gera o grafo bipartido aleatório (G_a_b)
                G = nx.bipartite.random_graph(a, b, p_value)

                # 2. Formata o nome do arquivo
                # Ex: bi_a100_b100_p01%_v1.col
                p_percent = int(p_value * 100) 
                file_name = f"bi_a{a}_b{b}_p{p_percent:02d}%_v{v}.col"

                # 3. Escreve o grafo no arquivo no formato DIMACS
                with open(file_name, 'w') as f:
                    num_nodes = G.number_of_nodes()
                    num_edges = G.number_of_edges()
                    
                    # Cabeçalho: p edge [nós] [arestas]
                    f.write(f"p edge {num_nodes} {num_edges}\n")

                    # Corpo: Lista de arestas 'e u v'
                    for u, v in G.edges():
                        # Adiciona 1 para que os nós comecem em 1 (padrão DIMACS)
                        f.write(f"e {u+1} {v+1}\n")

                print(f"Gerado: {file_name} (Nós: {num_nodes}, Arestas: {num_edges})")

print("Geração de todos os arquivos de grafos bipartidos concluída.")
"""