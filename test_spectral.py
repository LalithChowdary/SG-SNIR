import networkx as nx
from src.data_loader import load_dataset
from src.sg_snir import get_kscc

G = load_dataset("p2p-Gnutella")
S_star = get_kscc(G)
subgraph = G.subgraph(S_star)

intra_degrees = {}
for n in S_star:
    intra_degrees[n] = (subgraph.in_degree(n), subgraph.out_degree(n))

vol_S_star = sum(d_out for _, d_out in intra_degrees.values())
deg_product_sum = sum(d_in * d_out for d_in, d_out in intra_degrees.values())
rho_S_star = deg_product_sum / vol_S_star

drops = []
for u, v in subgraph.edges():
    d_in_u = intra_degrees[u][0]
    d_out_v = intra_degrees[v][1]
    new_numerator = deg_product_sum - d_in_u - d_out_v
    new_rho = new_numerator / (vol_S_star - 1)
    drop = rho_S_star - new_rho
    drops.append(drop)

print(f"rho={rho_S_star}")
print(f"min drop: {min(drops)}")
print(f"max drop: {max(drops)}")
