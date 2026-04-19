"""
src/sg_snir.py
==============
SG-SNIR: Spectral-Guided SNIR Edge Blocking Algorithm.

Implements the 6-step algorithm from the paper:
  Step 1 — Initialize SNIR (set W = N ∪ I)
  Step 2 — Compute candidate edges Γ(W) restricted to Susceptible targets
  Step 3 — KSCC identification using intra-SCC degrees + bridge detection
  Step 4 — Adaptive spectral filtering (bimodal: bridges bypass threshold)
  Step 5 — SNIR simulation for each candidate edge
  Step 6 — Select and remove minimum-H edge (static W throughout)
"""

import networkx as nx
from typing import Dict, FrozenSet, List, Optional, Set, Tuple

from src.snir_model import SNIRParams, compute_influence_range


# ---------------------------------------------------------------------------
# Step 2 — Candidate edge generation
# ---------------------------------------------------------------------------

def get_candidates(
    G: nx.DiGraph,
    W: Set,
    initial_S: Set,
) -> List[Tuple]:
    """
    Γ(W) = {(u,v) | u ∈ W, v ∈ S, (u,v) ∈ E}

    This refines the original MaxExpectedH definition (v ∉ W) by
    further restricting to v ∈ S, since edges to Recovered nodes
    carry zero transmission probability under the SNIR model.

    Returns a list of (u, v) edge tuples.
    """
    candidates = []
    for u in W:
        if u not in G:
            continue
        for v in G.successors(u):
            if v in initial_S:
                candidates.append((u, v))
    return candidates


# ---------------------------------------------------------------------------
# Step 3 — KSCC identification
# ---------------------------------------------------------------------------

def get_kscc(
    G: nx.DiGraph,
) -> Tuple[Optional[FrozenSet], float, Dict, float]:
    """
    Identify the Key Strongly Connected Component S* using intra-SCC degrees.

    ρ(S_i) ≈ Σ_{u∈S_i} d_in^{S_i}(u) · d_out^{S_i}(u) / vol(S_i)

    Only non-trivial SCCs (|S_i| > 1) contribute; single-node SCCs have
    spectral radius 0 (no self-loops assumed).

    Returns
    -------
    S_star        : frozenset of nodes in the KSCC, or None if no non-trivial SCC.
    rho_S_star    : Approximated spectral radius of KSCC.
    intra_degrees : {node: (d_in_intra, d_out_intra)} for all nodes in S_star.
    vol_S_star    : vol(S*) = Σ d_in^{S*}(u) for u ∈ S*.
    """
    best_scc: Optional[FrozenSet] = None
    best_rho: float = -1.0
    best_intra: Dict = {}
    best_vol: float = 0.0

    sccs = list(nx.strongly_connected_components(G))

    for scc in sccs:
        if len(scc) < 2:
            continue  # trivial SCC, ρ = 0

        sub = G.subgraph(scc)
        # Intra-SCC in/out degrees (only edges within the SCC count)
        d_in  = dict(sub.in_degree())
        d_out = dict(sub.out_degree())

        vol = sum(d_in[u] for u in scc)
        if vol == 0:
            continue

        deg_product_sum = sum(d_in[u] * d_out[u] for u in scc)
        rho = deg_product_sum / vol

        if rho > best_rho:
            best_rho  = rho
            best_scc  = frozenset(scc)
            best_intra = {u: (d_in[u], d_out[u]) for u in scc}
            best_vol   = vol

    return best_scc, best_rho, best_intra, best_vol


# ---------------------------------------------------------------------------
# Step 4 — Adaptive spectral filtering
# ---------------------------------------------------------------------------

def spectral_filter(
    gamma_prime: List[Tuple],
    S_star: Optional[FrozenSet],
    rho_S_star: float,
    intra_degrees: Dict,
    vol_S_star: float,
    deg_product_sum: float,
    eps: float = 0.01,
) -> List[Tuple]:
    """
    Build candidate set C via bimodal filtering.

    For INTERNAL edges (u,v) ∈ E_{S*}:
        SpectralDrop(u→v) ≈ ρ(S*) − (deg_product_sum − d_in(u) − d_out(v)) / (vol(S*)−1)
        Edge enters C iff SpectralDrop > ε · ρ(S*)

        NOTE: Equation 17 is defined ONLY for edges where u,v ∈ S*.
              Bridge edges MUST NOT be passed to this formula.

    For BRIDGE edges (u∉S*, v∈S*):
        SpectralDrop = 0 by construction (u's degrees do not appear in the
        KSCC degree-product sum), so they are admitted directly without
        spectral threshold check.

    Fallback: If C = ∅ after filtering, relax to C = Γ'(W).

    Parameters
    ----------
    gamma_prime    : Filtered candidate list Γ'(W) = Γ(W) ∩ (E_{S*} ∪ B_{S*}).
    S_star         : Frozenset of KSCC nodes.
    rho_S_star     : Spectral radius of KSCC.
    intra_degrees  : {node: (d_in_intra, d_out_intra)} for nodes in S*.
    vol_S_star     : vol(S*).
    deg_product_sum: Σ_{i∈S*} d_in(i) · d_out(i), pre-computed.
    eps            : Spectral threshold coefficient ε.

    Returns
    -------
    C : List of candidate edges after bimodal spectral filtering.
    """
    if S_star is None or vol_S_star <= 1:
        return list(gamma_prime)

    C = []
    threshold = eps * rho_S_star

    for (u, v) in gamma_prime:
        u_internal = u in S_star
        v_internal = v in S_star

        if u_internal and v_internal:
            # Internal KSCC edge — apply SpectralDrop filter (Eq. 17)
            d_in_u  = intra_degrees[u][0]
            d_out_v = intra_degrees[v][1]
            new_numerator = deg_product_sum - d_in_u - d_out_v
            new_rho = new_numerator / (vol_S_star - 1)
            spectral_drop = rho_S_star - new_rho
            if spectral_drop > threshold:
                C.append((u, v))
        else:
            # Bridge edge: u ∉ S*, v ∈ S* — bypass threshold (Eq. 18)
            # SpectralDrop = 0 for bridge edges; importance is epidemiological.
            C.append((u, v))

    # Fallback: if all internal edges are below threshold, relax to Γ'(W)
    if len(C) == 0:
        C = list(gamma_prime)

    return C


# ---------------------------------------------------------------------------
# Master Algorithm: sg_snir_blocking
# ---------------------------------------------------------------------------

def sg_snir_blocking(
    G: nx.DiGraph,
    initial_S: Set,
    initial_N: Set,
    initial_I: Set,
    initial_R: Set,
    params: SNIRParams,
    k: int,
    T: int = 10,
    eps: float = 0.01,
    default_weight: float = 0.05,
    verbose: bool = False,
) -> Tuple[List[Tuple], List[float], List[int]]:
    """
    SG-SNIR: Spectral-Guided SNIR Edge Blocking (Algorithm from paper).

    Static budget allocation: W = N ∪ I is fixed at t=0 and held constant
    throughout all k iterations, modeling a 'Day 1' preemptive intervention
    where all k edges must be identified and blocked upfront.

    Parameters
    ----------
    G            : Directed graph (modified in-place during selection,
                   then edges restored after each eval).
    initial_S/N/I/R : Initial state node sets.
    params       : SNIRParams.
    k            : Edge blocking budget.
    T            : SNIR simulation horizon.
    eps          : Spectral threshold ε for SpectralDrop filter.
    default_weight: Default edge transmission probability.
    verbose      : Print per-iteration diagnostics if True.

    Returns
    -------
    selected_edges : List of edges removed (in order of selection).
    H_history      : H value after each removal (length k).
    eval_counts    : Number of SNIR evaluations per iteration (for efficiency
                     comparison against MaxExpectedH baseline).
    """
    G = G.copy()

    # Step 1 — Static W: fixed at t=0, never updated
    W: Set = set(initial_N) | set(initial_I)

    selected_edges: List[Tuple] = []
    H_history:     List[float] = []
    eval_counts:   List[int]   = []

    for iteration in range(k):
        # Step 2 — Compute Γ(W) restricted to Susceptible targets
        gamma_W = get_candidates(G, W, initial_S)

        if not gamma_W:
            if verbose:
                print(f"[iter {iteration}] Γ(W) is empty — budget exhausted early.")
            break

        # Step 3 — KSCC identification using intra-SCC degrees
        S_star, rho_S_star, intra_degrees, vol_S_star = get_kscc(G)

        # Build Γ'(W) = Γ(W) ∩ (E_{S*} ∪ B_{S*})
        if S_star is not None:
            gamma_prime = [
                (u, v) for (u, v) in gamma_W
                if v in S_star  # v in S* covers both internal and bridge edges
            ]
        else:
            gamma_prime = []

        # Fallback: if epidemic disconnected from KSCC, use full Γ(W)
        if len(gamma_prime) == 0 and len(gamma_W) > 0:
            if verbose:
                print(f"[iter {iteration}] Γ'(W) empty — epidemic disconnected from KSCC. Fallback to Γ(W).")
            gamma_prime = gamma_W

        # Step 4 — Adaptive spectral filtering
        if S_star is not None and vol_S_star > 1:
            deg_product_sum = sum(
                intra_degrees[u][0] * intra_degrees[u][1] for u in S_star
            )
        else:
            deg_product_sum = 0.0

        C = spectral_filter(
            gamma_prime, S_star, rho_S_star,
            intra_degrees, vol_S_star, deg_product_sum, eps
        )

        if verbose:
            print(f"[iter {iteration}] |Γ(W)|={len(gamma_W)}, |Γ'(W)|={len(gamma_prime)}, |C|={len(C)}, ρ(S*)={rho_S_star:.4f}")

        # Step 5 — SNIR simulation for each edge in C
        best_edge = None
        best_H    = float("inf")
        evals     = 0

        for (u, v) in C:
            G.remove_edge(u, v)
            H_e, _ = compute_influence_range(
                G, initial_S, initial_N, initial_I, initial_R,
                params, T, default_weight
            )
            G.add_edge(u, v, **_edge_attrs(G, u, v, default_weight))
            evals += 1

            if H_e < best_H:
                best_H   = H_e
                best_edge = (u, v)

        if best_edge is None:
            break

        # Step 6 — Permanently remove best edge
        G.remove_edge(*best_edge)
        selected_edges.append(best_edge)
        H_history.append(best_H)
        eval_counts.append(evals)

        if verbose:
            print(f"[iter {iteration}] Removed {best_edge} → H={best_H:.4f} (evals={evals})")

    return selected_edges, H_history, eval_counts


def _edge_attrs(G: nx.DiGraph, u, v, default_weight: float) -> Dict:
    """Return edge attributes to restore after temporary removal."""
    # Edge doesn't exist anymore after removal; return default attrs
    return {"weight": default_weight}
