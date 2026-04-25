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
    
    k_nodes_max = 5
    print(f"Running DINO to select {k_nodes_max} nodes...")
    start_dino = time.time()
    dino = Dino(None)
    # Get immunized nodes (returns integer indices of A)
    dino_indices = dino.get_immunized_nodes(A, k_nodes_max)
    dino_nodes = [nodes_list_sorted[i] for i in dino_indices]
    dino_time = time.time() - start_dino
    
    print(f"DINO selected nodes: {dino_nodes} in {dino_time:.1f}s")
    
    # Pre-calculate the cumulative edges removed by DINO at each k step
    dino_cumulative_edges = []
    edges_so_far = set()
    for u in dino_nodes:
        # All outgoing edges
        for v in G.successors(u):
            edges_so_far.add((u, v))
        # All incoming edges
        for pre in G.predecessors(u):
            edges_so_far.add((pre, u))
        dino_cumulative_edges.append(len(edges_so_far))
        
    max_budget = dino_cumulative_edges[-1]
    
    # Run SG-SNIR once up to max_budget
    print(f"\nRunning SG-SNIR for up to {max_budget} edges...")
    start_sg = time.time()
    _, sg_H_history, sg_evals = sg_snir_blocking(G, initial_S, initial_N, initial_I, initial_R, params, max_budget, T)
    sg_time = time.time() - start_sg
    
    # Evaluate DINO at each k step
    results = []
    edges_so_far_eval = set()
    
    for i, u in enumerate(dino_nodes):
        k = i + 1
        for v in G.successors(u):
            edges_so_far_eval.add((u, v))
        for pre in G.predecessors(u):
            edges_so_far_eval.add((pre, u))
            
        budget = len(edges_so_far_eval)
        
        # Evaluate DINO's graph
        G_dino = G.copy()
        G_dino.remove_edges_from(edges_so_far_eval)
        H_dino, _ = compute_influence_range(G_dino, initial_S, initial_N, initial_I, initial_R, params, T)
        
        # Get SG-SNIR at exact budget (or its last value if budget exceeds total possible transmitting edges)
        sg_idx = min(budget, len(sg_H_history) - 1)
        H_sg = sg_H_history[sg_idx]
        
        print(f"k={k} nodes | Budget={budget} edges | DINO H={H_dino:.4f} | SG-SNIR H={H_sg:.4f}")
        
        results.append({
            'k_nodes': k,
            'budget_edges': budget,
            'dino_H': float(H_dino),
            'sg_H': float(H_sg)
        })

    # Save results
    os.makedirs('results/experiment9', exist_ok=True)
    with open(f'results/experiment9/{dataset_name}_baseline_comparison.json', 'w') as f:
        json.dump({'dataset': dataset_name, 'data': results}, f, indent=4)
        
    print(f"\nSaved to results/experiment9/{dataset_name}_baseline_comparison.json")

if __name__ == '__main__':
    run_dino_edge_equivalent()
