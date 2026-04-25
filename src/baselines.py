"""
src/baselines.py
================
Baseline edge blocking algorithms for SG-SNIR comparison.

Three baselines:
  1. MaxExpectedH  — Greedy exhaustive simulation over all Γ(W) edges.
  2. DegreeProduct — Fast structural heuristic: rank by d_out(u)·d_in(v).
                     No theoretical grounding for edge-level containment;
                     included purely as a computationally cheap contrast.
  3. Random        — Uniformly random edge selection from Γ(W).
"""

import random
from typing import Dict, List, Set, Tuple

import networkx as nx

from src.snir_model import SNIRParams, compute_influence_range


# ---------------------------------------------------------------------------
# Baseline 1: MaxExpectedH (exhaustive greedy)
# ---------------------------------------------------------------------------

def max_expected_h_blocking(
    G: nx.DiGraph,
    initial_S: Set,
    initial_N: Set,
    initial_I: Set,
    initial_R: Set,
    params: SNIRParams,
    k: int,
    T: int = 10,
    default_weight: float = 0.05,
    verbose: bool = False,
) -> Tuple[List[Tuple], List[float], List[int]]:
    """
    Greedy edge blocking: at each iteration evaluate ALL edges in Γ(W)
    with a full SNIR simulation and remove the one minimising H.

    Γ(W) = {(u,v) | u ∈ W, v ∈ S, (u,v) ∈ E}
    (same definition as SG-SNIR Step 2)

    This is the accuracy gold standard; SG-SNIR should match it closely
    while requiring far fewer SNIR evaluations per iteration.

    Returns
    -------
    selected_edges : Edges removed in order.
    H_history      : H after each removal.
    eval_counts    : SNIR evaluations per iteration (= |Γ(W)| each time).
    """
    G = G.copy()
    W: Set = set(initial_N) | set(initial_I)

    selected_edges: List[Tuple] = []
    H_history:      List[float] = []
    eval_counts:    List[int]   = []

    for iteration in range(k):
        gamma_W = [
            (u, v) for u in W if u in G
            for v in G.successors(u) if v in initial_S
        ]

        if not gamma_W:
            if verbose:
                print(f"[MaxExpH iter {iteration}] Γ(W) empty.")
            break

        best_edge = None
        best_H    = float("inf")

        for (u, v) in gamma_W:
            # Save original edge attributes before removal so they are
            # restored exactly — crucial for heterogeneous edge weights.
            original_attrs = G[u][v].copy()
            G.remove_edge(u, v)
            H_e, _ = compute_influence_range(  # _ = per-timestep state history, not needed here
                G, initial_S, initial_N, initial_I, initial_R,
                params, T, default_weight
            )
            G.add_edge(u, v, **original_attrs)
            if H_e < best_H:
                best_H    = H_e
                best_edge = (u, v)

        if best_edge is None:
            break

        G.remove_edge(*best_edge)
        selected_edges.append(best_edge)
        H_history.append(best_H)
        eval_counts.append(len(gamma_W))

        if verbose:
            print(f"[MaxExpH iter {iteration}] Removed {best_edge} → H={best_H:.4f} (evals={len(gamma_W)})")

    return selected_edges, H_history, eval_counts


# ---------------------------------------------------------------------------
# Baseline 2: DegreeProduct heuristic
# ---------------------------------------------------------------------------

def degree_product_blocking(
    G: nx.DiGraph,
    initial_S: Set,
    initial_N: Set,
    initial_I: Set,
    initial_R: Set,
    k: int,
    verbose: bool = False,
) -> Tuple[List[Tuple], List[None], List[int]]:
    """
    Fast structural baseline: rank Γ(W) edges by d_out(u) · d_in(v)
    (full-graph degrees) and greedily remove the top-k.

    This is a computationally cheap O(|E|) heuristic with no theoretical
    guarantee for epidemic containment. It is included purely to show the
    cost of ignoring epidemic dynamics.

    No SNIR simulations are run; eval_counts = [0, ..., 0].
    H values are not tracked (None placeholders).
    """
    G = G.copy()
    W: Set = set(initial_N) | set(initial_I)

    selected_edges: List[Tuple] = []

    for iteration in range(k):
        gamma_W = [
            (u, v) for u in W if u in G
            for v in G.successors(u) if v in initial_S
        ]

        if not gamma_W:
            break

        # Score = d_out(u) * d_in(v) — full-graph degrees
        best_edge = max(
            gamma_W,
            key=lambda e: G.out_degree(e[0]) * G.in_degree(e[1]),
        )

        G.remove_edge(*best_edge)
        selected_edges.append(best_edge)

        if verbose:
            u, v = best_edge
            score = G.out_degree(u) * G.in_degree(v)
            print(f"[DegProd iter {iteration}] Removed {best_edge} (score={score})")

    return selected_edges, [None] * len(selected_edges), [0] * len(selected_edges)


# ---------------------------------------------------------------------------
# Baseline 3: Random
# ---------------------------------------------------------------------------

def random_blocking(
    G: nx.DiGraph,
    initial_S: Set,
    initial_N: Set,
    initial_I: Set,
    initial_R: Set,
    params: SNIRParams,
    k: int,
    T: int = 10,
    default_weight: float = 0.05,
    seed: int = 42,
    verbose: bool = False,
) -> Tuple[List[Tuple], List[float], List[int]]:
    """
    Random baseline: uniformly sample one edge from Γ(W) per iteration.
    H is computed once after each removal for fair comparison.
    """
    rng = random.Random(seed)
    G   = G.copy()
    W: Set = set(initial_N) | set(initial_I)

    selected_edges: List[Tuple] = []
    H_history:      List[float] = []

    for iteration in range(k):
        gamma_W = [
            (u, v) for u in W if u in G
            for v in G.successors(u) if v in initial_S
        ]

        if not gamma_W:
            break

        chosen = rng.choice(gamma_W)
        G.remove_edge(*chosen)
        H_e, _ = compute_influence_range(
            G, initial_S, initial_N, initial_I, initial_R,
            params, T, default_weight
        )
        selected_edges.append(chosen)
        H_history.append(H_e)

        if verbose:
            print(f"[Random iter {iteration}] Removed {chosen} → H={H_e:.4f}")

    return selected_edges, H_history, [1] * len(selected_edges)
