import networkx as nx
import random

# --- DEFINIÇÃO DA FUNÇÃO GERADORA DE GRAFOS GALÁXIA ---

def galaxy(a, b_min, c_max):
    """
    Cria um grafo "galáxia" composto por 'a' grafos estrela, que permanecem desconectados.

    Args:
        a (int): Número de centros (grafos estrela).
        b_min (int): Número mínimo de satélites por estrela.
        c_max (int): Número máximo de satélites por estrela.

    Returns:
        nx.Graph: O grafo galáxia gerado (com componentes desconectados).
    """
    if a < 1:
        return nx.Graph()
        
    G = nx.Graph()
    current_node_offset = 0
    
    # 1. GERAÇÃO DOS GRAFOS ESTRELA
    for _ in range(a):
        # Sorteia uniformemente o número de vértices satélites para a estrela
        num_satellites = random.randint(b_min, c_max)
        
        # Cria a estrela com o número de satélites sorteado
        star = nx.star_graph(num_satellites)
        
        # Renomeia os nós para garantir IDs exclusivos (começando em current_node_offset)
        mapping = {node: node + current_node_offset for node in star.nodes()}
        star_relabel = nx.relabel_nodes(star, mapping)
        
        # Adiciona o grafo estrela ao grafo principal G.
        G = nx.union(G, star_relabel)
        
        # Atualiza o offset para a próxima estrela (num_satellites + 1 centro)
        current_node_offset += num_satellites + 1
    
    # O grafo G agora é composto por 'a' componentes estrela desconectados.
    return G

# --- PARÂMETROS DE EXECUÇÃO PARA GRAFOS GALÁXIA (NOVOS PARÂMETROS) ---

# Geração de 'a' variando de 1 a 100 de 5 em 5, corrigida.
# Inclui o 1, e depois gera de 5 até 100 de 5 em 5.
NUM_CENTERS = [1] + list(range(5, 101, 5))

# Geração de pares (b_min, c_max) onde b_min = 10 e c_max varia de 20 a 100 de 10 em 10.
MIN_B = 10
MAX_C_VALUES = list(range(20, 101, 10))
MIN_MAX_PAIRS = [(MIN_B, c_max) for c_max in MAX_C_VALUES]

NUM_VERSIONS = 5 # Gerar 5 versões para cada combinação

print(f"Iniciando a geração de grafos galáxia (desconectados) com novos parâmetros de variação...")

# LAÇO 1: Itera sobre o número de centros (a)
for a in NUM_CENTERS:
    
    # LAÇO 2: Itera sobre os pares de mínimo e máximo de satélites (b e c)
    for b_min, c_max in MIN_MAX_PAIRS:
        
        # LAÇO 3: Gera as 5 versões
        for v in range(1, NUM_VERSIONS + 1):
            
            # 1. Gera o grafo GALÁXIA
            G = galaxy(a, b_min, c_max)

            # 2. Formata o nome do arquivo
            # Ex: galaxy_a01_b10_c20_v1.col
            file_name = f"galaxy_a{a:02d}_b{b_min:02d}_c{c_max:02d}_v{v}.col"

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

            print(f"Gerado: {file_name} (Centros: {a}, Satélites: [{b_min}-{c_max}], Nós: {num_nodes}, Arestas: {num_edges})")

print("Geração de todos os arquivos de grafos galáxia concluída.")