import sys
import os
import json
import time
import networkx as nx
import numpy as np
import random

# Add paper-2 to path to import DINO
paper2_path = os.path.abspath(os.path.join(os.path.dirname(__file__), '../../paper-2'))
if paper2_path not in sys.path:
    sys.path.append(paper2_path)
    
from measures.dino import Dino

from src.data_loader import load_dataset, assign_fixed_weights
from src.snir_model import SNIRParams, compute_influence_range
from src.sg_snir import sg_snir_blocking

def run_dino_edge_equivalent():
    dataset_name = 'p2p-Gnutella'
    print(f"Loading {dataset_name}...")
    G = load_dataset(dataset_name)
    G = assign_fixed_weights(G)
    
    # SNIR params
    params = SNIRParams(
        alpha=0.033, beta=0.022, delta=0.020,
        eta=0.140, gamma=0.315, xi=0.300
    )
    T = 10
    n_seeds = 5
    
    random.seed(42)
    nodes_list = list(G.nodes())
    initial_I = set(random.sample(nodes_list, n_seeds))
    initial_N = set()
    initial_R = set()
    initial_S = set(nodes_list) - initial_I
    
    print("Preparing adjacency matrix for DINO...")
    nodes_list_sorted = sorted(nodes_list)
    A = nx.to_numpy_array(G, nodelist=nodes_list_sorted)
    
    k_nodes = 5
    print(f"Running DINO to select {k_nodes} nodes...")
    start_dino = time.time()
    dino = Dino(None)
    # Get immunized nodes (returns integer indices of A)
    dino_indices = dino.get_immunized_nodes(A, k_nodes)
    dino_nodes = [nodes_list_sorted[i] for i in dino_indices]
    dino_time = time.time() - start_dino
    
    print(f"DINO selected nodes: {dino_nodes} in {dino_time:.1f}s")
    
    # Count how many edges this removes
    dino_edges_to_remove = set()
    for u in dino_nodes:
        # All outgoing edges
        for v in G.successors(u):
            dino_edges_to_remove.add((u, v))
        # All incoming edges
        for pre in G.predecessors(u):
            dino_edges_to_remove.add((pre, u))
            
    num_edges_budget = len(dino_edges_to_remove)
    print(f"These {k_nodes} nodes correspond to {num_edges_budget} total edges.")
    
    # Evaluate DINO's graph
    G_dino = G.copy()
    G_dino.remove_edges_from(dino_edges_to_remove)
    print("Evaluating DINO intervention via SNIR simulation...")
    H_dino, _ = compute_influence_range(G_dino, initial_S, initial_N, initial_I, initial_R, params, T)
    print(f"DINO H_final: {H_dino:.4f}")
    
    # Run SG-SNIR for exactly num_edges_budget
    print(f"\nRunning SG-SNIR for {num_edges_budget} edges...")
    start_sg = time.time()
    _, sg_H_history, sg_evals = sg_snir_blocking(G, initial_S, initial_N, initial_I, initial_R, params, num_edges_budget, T)
    sg_time = time.time() - start_sg
    H_sg = sg_H_history[-1]
    
    print(f"SG-SNIR H_final: {H_sg:.4f} in {sg_time:.1f}s")
    print(f"Quality Gap (SG-SNIR improvement): {H_dino - H_sg:.4f}")
    
    # Save results
    results = {
        'dataset': dataset_name,
        'k_nodes': k_nodes,
        'dino': {
            'selected_nodes': [int(n) if isinstance(n, np.integer) else n for n in dino_nodes],
            'edges_removed': num_edges_budget,
            'H_final': float(H_dino),
            'time_s': dino_time
        },
        'sg_snir': {
            'edges_removed': num_edges_budget,
            'H_final': float(H_sg),
            'time_s': sg_time,
            'evals': sg_evals
        }
    }
    os.makedirs('results/experiment9', exist_ok=True)
    with open(f'results/experiment9/{dataset_name}_baseline_comparison.json', 'w') as f:
        json.dump(results, f, indent=4)
        
    print(f"\nSaved to results/experiment9/{dataset_name}_baseline_comparison.json")

if __name__ == '__main__':
    run_dino_edge_equivalent()
