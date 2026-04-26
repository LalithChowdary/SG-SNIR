import sys
import os
import json
import time
import networkx as nx
import numpy as np
import random
import logging

# Add paper-2 to path to import DINO
paper2_path = os.path.abspath(os.path.join(os.path.dirname(__file__), '../../paper-2'))
if paper2_path not in sys.path:
    sys.path.append(paper2_path)
    
from measures.dino import Dino

from src.data_loader import load_dataset, assign_fixed_weights
from src.snir_model import SNIRParams, compute_influence_range
from src.sg_snir import sg_snir_blocking

def setup_logger():
    log_dir = 'results/experiment9/logs'
    os.makedirs(log_dir, exist_ok=True)
    
    logger = logging.getLogger('exp9')
    logger.setLevel(logging.INFO)
    
    # File handler
    fh = logging.FileHandler(os.path.join(log_dir, 'exp9_run.log'), mode='w')
    fh.setLevel(logging.INFO)
    
    # Console handler
    ch = logging.StreamHandler()
    ch.setLevel(logging.INFO)
    
    formatter = logging.Formatter('%(asctime)s - %(message)s', datefmt='%Y-%m-%d %H:%M:%S')
    fh.setFormatter(formatter)
    ch.setFormatter(formatter)
    
    logger.addHandler(fh)
    logger.addHandler(ch)
    return logger

def run_dino_edge_equivalent():
    logger = setup_logger()
    dataset_name = 'p2p-Gnutella'
    logger.info(f"Loading {dataset_name}...")
    G = load_dataset(dataset_name)
    G = assign_fixed_weights(G)
    
    # SNIR params
    params = SNIRParams(
        alpha=0.033, beta=0.022, delta=0.020,
        eta=0.140, gamma=0.315, xi=0.300
    )
    T = 10
    
    nodes_list = list(G.nodes())
    nodes_list_sorted = sorted(nodes_list)
    
    logger.info("Preparing adjacency matrix for DINO...")
    A = nx.to_numpy_array(G, nodelist=nodes_list_sorted)
    
    configs = [
        (300, [5, 10, 15]),
        (600, [10, 15, 20]),
        (900, [15, 20, 25])
    ]
    
    all_results = {}
    
    for seeds, dino_k_list in configs:
        logger.info(f"\n{'='*50}\nStarting Configuration: {seeds} Seeds\n{'='*50}")
        
        random.seed(42)
        initial_I = set(random.sample(nodes_list, seeds))
        initial_N = set()
        initial_R = set()
        initial_S = set(nodes_list) - initial_I
        
        k_nodes_max = max(dino_k_list)
        logger.info(f"Running DINO to select {k_nodes_max} nodes...")
        start_dino = time.time()
        dino = Dino(None)
        # Get immunized nodes
        dino_indices = dino.get_immunized_nodes(A, k_nodes_max)
        dino_nodes = [nodes_list_sorted[i] for i in dino_indices]
        dino_time = time.time() - start_dino
        logger.info(f"DINO finished in {dino_time:.1f}s")
        
        # Pre-calculate the cumulative edges removed by DINO at each k step
        dino_cumulative_edges = []
        edges_so_far = set()
        for u in dino_nodes:
            for v in G.successors(u):
                edges_so_far.add((u, v))
            for pre in G.predecessors(u):
                edges_so_far.add((pre, u))
            dino_cumulative_edges.append(len(edges_so_far))
            
        max_budget = dino_cumulative_edges[-1]
        
        logger.info(f"Running SG-SNIR for up to {max_budget} edges...")
        start_sg = time.time()
        # Pass verbose=False or capture it
        _, sg_H_history, sg_evals = sg_snir_blocking(G, initial_S, initial_N, initial_I, initial_R, params, max_budget, T)
        sg_time = time.time() - start_sg
        logger.info(f"SG-SNIR finished in {sg_time:.1f}s. History length: {len(sg_H_history)}")
        
        results = []
        
        for k in dino_k_list:
            idx = k - 1
            budget = dino_cumulative_edges[idx]
            
            # Nodes for this k
            k_nodes = dino_nodes[:k]
            edges_to_remove = set()
            for u in k_nodes:
                for v in G.successors(u):
                    edges_to_remove.add((u, v))
                for pre in G.predecessors(u):
                    edges_to_remove.add((pre, u))
                    
            # Evaluate DINO's graph
            G_dino = G.copy()
            G_dino.remove_edges_from(edges_to_remove)
            H_dino, _ = compute_influence_range(G_dino, initial_S, initial_N, initial_I, initial_R, params, T)
            
            # Get SG-SNIR at exact budget (or its last value)
            sg_idx = min(budget, len(sg_H_history) - 1)
            H_sg = sg_H_history[sg_idx]
            
            logger.info(f"  --> k={k} nodes | Budget={budget} edges | DINO H={H_dino:.4f} | SG-SNIR H={H_sg:.4f}")
            
            results.append({
                'k_nodes': k,
                'budget_edges': budget,
                'dino_H': float(H_dino),
                'sg_H': float(H_sg)
            })
            
        all_results[f"{seeds}_seeds"] = results
        
        # Save incremental
        with open(f'results/experiment9/logs/results_{seeds}_seeds.json', 'w') as f:
            json.dump({'dataset': dataset_name, 'seeds': seeds, 'data': results}, f, indent=4)
            
    logger.info("\nAll configurations completed successfully!")

if __name__ == '__main__':
    run_dino_edge_equivalent()
