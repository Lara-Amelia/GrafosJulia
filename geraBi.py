import networkx as nx

a = 100
b = 100
p = 0.01

G = nx.bipartite.random_graph(a, b, p)

file_name = "bi_a100_b100_p1%_v1.col"

# escreve as infos do grafo em arquivo no formato DIMACS
with open(file_name, 'w') as f:
    num_nodes = G.number_of_nodes()
    num_edges = G.number_of_edges()
    f.write(f"p edge {num_nodes} {num_edges}\n")

    for u, v in G.edges():
        f.write(f"e {u} {v}\n")