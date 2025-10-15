import networkx as nx

# --- Parâmetros de Execução ---
N = 2000 
P_VALUES = [0.01, 0.03, 0.05, 0.10, 0.20, 0.30, 0.40, 0.50, 0.60, 0.70, 0.80, 0.90]
NUM_VERSIONS = 5 # Gerar 5 versões para cada probabilidade

print(f"Iniciando a geração de grafos para N={N}...")

# Loop principal para iterar sobre cada probabilidade (p)
for p_value in P_VALUES:
    
    # Loop interno para gerar as 5 versões de cada probabilidade
    for v in range(1, NUM_VERSIONS + 1):
        
        # 1. Gera o grafo aleatório
        G = nx.erdos_renyi_graph(N, p_value)

        # 2. Formata o nome do arquivo
        # Ex: simples_n2000_p01%_v1.col (p_value * 100 para % formatada)
        p_percent = int(p_value * 100) 
        file_name = f"simples_n{N}_p{p_percent:02d}%_v{v}.col"

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

        print(f"Gerado: {file_name} (Arestas: {num_edges})")

print("Geração de todos os arquivos concluída.")