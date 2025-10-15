import networkx as nx

n = 100  # número de vértices
p = 0.90  # probabilidade de existência de aresta

G = nx.erdos_renyi_graph(n, p)

file_name = "simples_n100_p90%_v5.col"

# escreve as infos do grafo em arquivo no formato DIMACS
with open(file_name, 'w') as f:
    num_nodes = G.number_of_nodes()
    num_edges = G.number_of_edges()
    f.write(f"p edge {num_nodes} {num_edges}\n")

    for u, v in G.edges():
        f.write(f"e {u} {v}\n")