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

    Complexity: O(|W| × avg_degree) with O(1) set lookup for initial_S.

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
) -> Tuple[Optional[FrozenSet], float, Dict, float, float]:
    """
    Identify the Key Strongly Connected Component S* using intra-SCC degrees.

    ρ(S_i) ≈ Σ_{u∈S_i} d_in^{S_i}(u) · d_out^{S_i}(u) / vol(S_i)

    Only non-trivial SCCs (|S_i| > 1) contribute; single-node SCCs have
    spectral radius 0 (no self-loops assumed).

    Returns
    -------
    S_star            : frozenset of nodes in the KSCC, or None if no non-trivial SCC.
    rho_S_star        : Approximated spectral radius of KSCC.
    intra_degrees     : {node: (d_in_intra, d_out_intra)} for all nodes in S_star.
    vol_S_star        : vol(S*) = Σ d_in^{S*}(u) for u ∈ S*.
    deg_product_sum   : Σ_{u∈S*} d_in(u)·d_out(u), returned directly to avoid
                        redundant recomputation in the caller.
    """
    best_scc: Optional[FrozenSet] = None
    best_rho: float = -1.0
    best_intra: Dict = {}
    best_vol: float = 0.0
    best_deg_product_sum: float = 0.0

    sccs = list(nx.strongly_connected_components(G))

    for scc in sccs:
        if len(scc) < 2:
            # Paper Step 3: "Trivial SCCs (single nodes with no self-loop) are discarded."
            # A single node WITH a self-loop is non-trivial and must not be skipped.
            node = next(iter(scc))
            if not G.has_edge(node, node):
                continue  # truly trivial — no self-loop, spectral radius = 0
            # else: single node with self-loop — fall through to compute ρ

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
            best_rho            = rho
            best_scc            = frozenset(scc)
            best_intra          = {u: (d_in[u], d_out[u]) for u in scc}
            best_vol            = vol
            best_deg_product_sum = deg_product_sum

    return best_scc, best_rho, best_intra, best_vol, best_deg_product_sum


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

    Fallback: Bridge edges are unconditionally added above. Therefore, C is
    only empty if there are NO bridge edges AND all internal edges fail the
    spectral threshold. In this strict case, we relax the threshold and
    return C = Γ'(W).

    Parameters
    ----------
    gamma_prime    : Filtered candidate list Γ'(W) = Γ(W) ∩ (E_{S*} ∪ B_{S*}).
    S_star         : Frozenset of KSCC nodes.
    rho_S_star     : Spectral radius of KSCC.
    intra_degrees  : {node: (d_in_intra, d_out_intra)} for nodes in S*.
    vol_S_star     : vol(S*).
    deg_product_sum: Σ_{i∈S*} d_in(i) · d_out(i), pre-computed by get_kscc.
    eps            : Spectral threshold coefficient ε.

    Returns
    -------
    C : List of candidate edges after bimodal spectral filtering.
    """
    if S_star is None or vol_S_star <= 1:
        return list(gamma_prime)

    # Defensive assertion: all edges in gamma_prime must have v ∈ S*.
    # This is guaranteed by the Γ'(W) construction above, but assert
    # defensively to catch any future construction errors immediately.
    # Use `if __debug__:` so this is stripped in production (`python -O`).
    if __debug__ and S_star is not None:
        assert all(v in S_star for (_, v) in gamma_prime), (
            "gamma_prime contains edge with v ∉ S* — Γ'(W) construction error"
        )

    C = []
    threshold = eps * rho_S_star

    for (u, v) in gamma_prime:
        # All edges in gamma_prime are guaranteed to have v ∈ S* by construction.
        # Only u_internal varies — True for internal edges, False for bridge edges.
        u_internal = u in S_star

        if u_internal:
            # Internal KSCC edge — apply SpectralDrop filter (Eq. 17)
            #
            # Derivation of Δ_sum when removing edge u→v:
            #   u loses one outgoing edge: d_out(u) → d_out(u)-1
            #     Δ_u = d_in(u)·d_out(u) − d_in(u)·(d_out(u)−1) = d_in(u)
            #   v loses one incoming edge: d_in(v) → d_in(v)-1
            #     Δ_v = d_in(v)·d_out(v) − (d_in(v)−1)·d_out(v) = d_out(v)
            # Paper Eq. 17: subtract d_in(u) + d_out(v) directly.
            #

            # Defensive guard MUST come before any intra_degrees[v] access.
            # v ∈ S* is guaranteed by gamma_prime construction, but if
            # get_kscc ever has a bug this prevents a KeyError crash.
            if v not in intra_degrees:
                C.append((u, v))  # conservative: treat as bridge
                continue

            d_in_u  = intra_degrees[u][0]   # d_in^{S*}(u)
            d_out_v = intra_degrees[v][1]   # d_out^{S*}(v)

            # Δ_sum = d_in(u) + d_out(v) — direct paper formula (Eq. 17)
            new_numerator = deg_product_sum - d_in_u - d_out_v
            new_rho = new_numerator / (vol_S_star - 1)
            spectral_drop = rho_S_star - new_rho
            if spectral_drop > threshold:
                C.append((u, v))
        else:
            # Bridge edge: u ∉ S*, v ∈ S* — bypass threshold (Eq. 18)
            # SpectralDrop = 0 for bridge edges; importance is epidemiological.
            C.append((u, v))

    # Fallback: bridge edges are unconditionally added above, so C is only
    # empty when there are zero bridge edges AND all internal edges are below
    # the threshold. In this case relax the filter to Γ'(W).
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
    G            : Directed graph. G.copy() is made internally so the
                   caller's graph is never modified.
    initial_S    : Susceptible nodes at t=0. Held fixed throughout all
                   iterations — intentional design choice per static budget
                   allocation (paper Section IV, Step 6).
    initial_N/I/R: Initial state node sets.
    params       : SNIRParams.
    k            : Edge blocking budget.
    T            : SNIR simulation horizon.
    eps          : Spectral threshold ε for SpectralDrop filter.
    default_weight: Fallback transmission probability when an edge carries
                   no explicit 'weight' attribute. A warning is logged if
                   this fallback is triggered, since silent use of 0.05 on
                   a dataset with heterogeneous probabilities would silently
                   corrupt simulation results.
    verbose      : Print per-iteration diagnostics if True.

    Degradation note
    ----------------
    If no non-trivial SCC exists (S_star = None), the KSCC-based filtering
    is skipped and the algorithm degrades to vanilla MaxExpectedH over Γ(W).

    Returns
    -------
    selected_edges : List of edges removed (in order of selection).
    H_history      : H value after each removal (length k).
    eval_counts    : Number of SNIR evaluations per iteration (for efficiency
                     comparison against MaxExpectedH baseline).
    """
    G = G.copy()  # caller's graph is never modified

    # Warn once if any edge is missing a 'weight' attribute — using the
    # default_weight fallback silently on a heterogeneous dataset would
    # corrupt all H(e) computations.
    edges_without_weight = [(u, v) for u, v, d in G.edges(data=True) if "weight" not in d]
    if edges_without_weight:
        import warnings
        warnings.warn(
            f"{len(edges_without_weight)} edges have no 'weight' attribute; "
            f"defaulting to {default_weight}. Set explicit weights for correct results.",
            stacklevel=2,
        )

    # Step 1 — Static W: fixed at t=0, never updated.
    # W is intentionally static — see paper Section IV, Step 6.
    W: Set = set(initial_N) | set(initial_I)

    selected_edges: List[Tuple] = []
    H_history:     List[float] = []
    eval_counts:   List[int]   = []

    for iteration in range(k):
        # Step 2 — Compute Γ(W) restricted to Susceptible targets.
        # initial_S is static (same as W), consistent with static budget model.
        gamma_W = get_candidates(G, W, initial_S)

        if not gamma_W:
            if verbose:
                print(f"[iter {iteration}] Γ(W) is empty — budget exhausted early.")
            break

        # Step 3 — KSCC identification using intra-SCC degrees.
        # get_kscc returns deg_product_sum directly to avoid redundant recomputation.
        S_star, rho_S_star, intra_degrees, vol_S_star, deg_product_sum = get_kscc(G)

        # Build Γ'(W) = Γ(W) ∩ (E_{S*} ∪ B_{S*})
        # Both internal edges (u∈S*, v∈S*) and bridge edges (u∉S*, v∈S*)
        # share the same predicate: v ∈ S*. Outgoing KSCC edges (u∈S*, v∉S*)
        # are correctly excluded. The two-category split (internal vs bridge)
        # is handled inside spectral_filter where it matters for SpectralDrop.
        # Degradation: if S_star is None (no non-trivial SCC), gamma_prime is
        # left empty and the fallback below degrades to MaxExpectedH over Γ(W).
        if S_star is not None:
            gamma_prime = [
                (u, v) for (u, v) in gamma_W
                if v in S_star
                # Covers both: internal edges (u ∈ S*, v ∈ S*)
                # and bridge edges (u ∉ S*, v ∈ S*)
                # Edges where v ∉ S* are excluded — they have no spectral impact on S*
            ]
        else:
            gamma_prime = []

        # Fallback: if W is entirely disconnected from KSCC (or S_star is None),
        # degrade gracefully to MaxExpectedH behaviour over the full Γ(W).
        if len(gamma_prime) == 0 and len(gamma_W) > 0:
            if verbose:
                print(f"[iter {iteration}] Γ'(W) empty — no KSCC reachable from W. Degrading to Γ(W).")
            gamma_prime = gamma_W

        # Step 4 — Adaptive spectral filtering
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
            # Save original edge attributes before removal so they are
            # restored exactly — crucial for heterogeneous transmission
            # probabilities that differ from default_weight.
            original_attrs = G[u][v].copy()
            G.remove_edge(u, v)
            # Pass initial state sets — static allocation means we always evaluate
            # from the t=0 configuration, consistent with W being fixed at t=0
            # (paper Step 6: "The infected set W is held fixed throughout all k iterations").
            H_e, _ = compute_influence_range(  # _ = per-timestep history, not needed here
                G, initial_S, initial_N, initial_I, initial_R,
                params, T, default_weight
            )
            G.add_edge(u, v, **original_attrs)
            evals += 1

            if H_e < best_H:
                best_H    = H_e
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
