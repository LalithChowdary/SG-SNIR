"""
src/data_loader.py
==================
Dataset loaders for SG-SNIR experiments.

HIV Transmission Network
------------------------
Source: ICPSR #22140 (DS0002 — dyadic contact data)
The dyadic file encodes directed contact pairs between individuals.
Columns used: ID1 (ego/source), ID2 (alter/target).
Each row represents a directed edge ID1 → ID2.

Directed Erdős–Rényi Synthetic Networks
----------------------------------------
Three fixed configurations matching HIV network density (avg degree ≈ 4):
  ER-Small  (n=1000,  p=0.004)
  ER-Medium (n=3000,  p=0.0013)
  ER-Large  (n=5000,  p=0.0008)
p = 4/n in each case, maintaining constant expected average degree of 4.
"""

import random
from pathlib import Path
from typing import Dict, Set, Tuple

import networkx as nx
import numpy as np


# ---------------------------------------------------------------------------
# HIV Transmission Network
# ---------------------------------------------------------------------------

# Path to the dyadic data file relative to the project root
HIV_TSV_PATH = Path(__file__).parent.parent / "data" / "ICPSR_22140" / "DS0002" / "22140-0002-Data.tsv"


def load_hiv_network(
    tsv_path: Path = HIV_TSV_PATH,
    default_weight: float = 0.05,
) -> nx.DiGraph:
    """
    Parse ICPSR #22140 DS0002 and build a directed contact network.

    Returns a NetworkX DiGraph where each edge (u→v) represents a
    directed contact from individual ID1=u to individual ID2=v with
    a uniform transmission probability of `default_weight`.

    Self-loops are removed (no epidemiological meaning).
    """
    G = nx.DiGraph()

    with open(tsv_path, "r") as f:
        header = f.readline().strip().split("\t")
        try:
            id1_col = header.index("ID1")
            id2_col = header.index("ID2")
        except ValueError:
            raise ValueError("Expected columns 'ID1' and 'ID2' in the TSV header.")

        for line in f:
            parts = line.strip().split("\t")
            if len(parts) <= max(id1_col, id2_col):
                continue
            try:
                u = int(parts[id1_col])
                v = int(parts[id2_col])
            except ValueError:
                continue

            if u != v:  # drop self-loops
                G.add_edge(u, v, weight=default_weight)

    print(f"[HIV] Loaded network: {G.number_of_nodes()} nodes, {G.number_of_edges()} edges")
    return G


# ---------------------------------------------------------------------------
# Directed Erdős–Rényi Synthetic Networks
# ---------------------------------------------------------------------------

ER_CONFIGS = {
    "ER-Small":  {"n": 1000, "p": 0.004},
    "ER-Medium": {"n": 3000, "p": 0.0013},
    "ER-Large":  {"n": 5000, "p": 0.0008},
}
# p = 4/n in all three cases → constant expected average degree ≈ 4


def generate_er_network(
    n: int,
    p: float,
    seed: int = 42,
    default_weight: float = 0.05,
) -> nx.DiGraph:
    """
    Generate a directed Erdős–Rényi graph G(n, p) with uniform edge weights.
    """
    G = nx.erdos_renyi_graph(n, p, seed=seed, directed=True)
    nx.set_edge_attributes(G, default_weight, "weight")
    print(f"[ER] n={n}, p={p}: {G.number_of_nodes()} nodes, {G.number_of_edges()} edges "
          f"(avg degree ≈ {G.number_of_edges() / G.number_of_nodes():.2f})")
    return G


# ---------------------------------------------------------------------------
# Initial state seeding
# ---------------------------------------------------------------------------

def seed_initial_states(
    G: nx.DiGraph,
    frac_I: float = 0.02,
    frac_N: float = 0.02,
    seed: int = 42,
) -> Tuple[Set, Set, Set, Set]:
    """
    Randomly assign initial epidemic states to graph nodes.

    frac_I : Fraction of nodes placed in Infectious (I) state.
    frac_N : Fraction of nodes placed in Non-symptomatic (N) state.
    Remaining nodes are Susceptible (S). R = ∅ at t=0.

    Returns (initial_S, initial_N, initial_I, initial_R).
    """
    rng = random.Random(seed)
    nodes = list(G.nodes())
    rng.shuffle(nodes)

    n = len(nodes)
    n_I = max(1, int(frac_I * n))
    n_N = max(1, int(frac_N * n))

    initial_I = set(nodes[:n_I])
    initial_N = set(nodes[n_I:n_I + n_N])
    initial_R: Set = set()
    initial_S = set(nodes[n_I + n_N:])

    print(f"[Seed] S={len(initial_S)}, N={len(initial_N)}, I={len(initial_I)}, R=0")
    return initial_S, initial_N, initial_I, initial_R
