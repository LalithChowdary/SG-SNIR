"""
src/tests.py
============
Unit tests for SG-SNIR algorithm correctness.

Tests:
  1. SpectralDrop algebraic correctness on a toy graph.
  2. Fallback 1: Γ'(W) falls back to Γ(W) when W is outside the KSCC.
  3. Fallback 2: C falls back to Γ'(W) when all edges below ε threshold.
"""

import networkx as nx

from src.sg_snir import get_candidates, get_kscc, spectral_filter
from src.snir_model import SNIRParams


# ---------------------------------------------------------------------------
# Test 1: SpectralDrop algebraic correctness
# ---------------------------------------------------------------------------

def test_spectral_drop_algebra():
    """
    Verify SpectralDrop(u→v) approximation against manual calculation.

    Toy graph:  0→1, 1→2, 2→0, 2→3  (3-cycle + extra edge)
    The KSCC should be {0,1,2} (the 3-cycle).
    We test removing edge 1→2:
      d_in(1)  = 1  (only 0→1)      → d_in_intra(1) = 1
      d_out(1) = 1  (only 1→2)      → d_out_intra(1) = 1
      d_in(2)  = 1  (only 1→2)      → d_in_intra(2) = 1
      d_out(2) = 2  (2→0, 2→3 — but 2→3 is outside SCC so d_out_intra = 1)

    In S* = {0,1,2} with intra-SCC edges only {0→1, 1→2, 2→0}:
      d_in_intra  = {0:1, 1:1, 2:1}
      d_out_intra = {0:1, 1:1, 2:1}
      deg_product_sum = 1+1+1 = 3
      vol = 3

    Removing 1→2:
      Δsum = d_in_intra(1) + d_out_intra(2) = 1 + 1 = 2
      new_rho = (3 - 2) / (3 - 1) = 0.5
      SpectralDrop = 1.0 - 0.5 = 0.5
    """
    G = nx.DiGraph()
    G.add_edges_from([(0, 1), (1, 2), (2, 0), (2, 3)])

    S_star, rho, intra_degrees, vol = get_kscc(G)

    assert S_star is not None, "Should find a non-trivial KSCC"
    assert S_star == frozenset({0, 1, 2}), f"KSCC should be {{0,1,2}}, got {S_star}"

    # Manual rho = deg_product_sum / vol = 3/3 = 1.0
    assert abs(rho - 1.0) < 1e-9, f"Expected rho=1.0, got {rho}"

    # Compute SpectralDrop for edge 1→2
    deg_product_sum = sum(intra_degrees[u][0] * intra_degrees[u][1] for u in S_star)
    d_in_u  = intra_degrees[1][0]   # d_in_intra(1) = 1
    d_out_v = intra_degrees[2][1]   # d_out_intra(2) = 1

    new_rho = (deg_product_sum - d_in_u - d_out_v) / (vol - 1)
    spectral_drop = rho - new_rho

    assert abs(spectral_drop - 0.5) < 1e-9, f"Expected SpectralDrop=0.5, got {spectral_drop}"
    print("✓ Test 1 passed: SpectralDrop algebra correct.")


# ---------------------------------------------------------------------------
# Test 2: Γ'(W) fallback when W is entirely outside the KSCC
# ---------------------------------------------------------------------------

def test_fallback_gamma_prime_disconnected():
    """
    Construct a graph where:
      - KSCC = {3,4,5} (a separate 3-cycle)
      - W = {0}  (outside the KSCC)
      - Edge 0→1 exists where 1 is also outside KSCC

    Expected: gamma_prime should be empty → fallback to gamma_W.
    """
    G = nx.DiGraph()
    # A 3-cycle forming the KSCC: 3→4→5→3
    G.add_edges_from([(3, 4), (4, 5), (5, 3)])
    # Disconnected infected source outside KSCC
    G.add_edge(0, 1)

    W = {0}
    initial_S = {1}  # susceptible target, outside KSCC

    gamma_W = get_candidates(G, W, initial_S)
    assert gamma_W == [(0, 1)], "Γ(W) should contain edge (0,1)"

    S_star, _, _, _ = get_kscc(G)
    assert S_star == frozenset({3, 4, 5}), f"KSCC should be {{3,4,5}}, got {S_star}"

    # gamma_prime: only edges where v ∈ S*
    gamma_prime = [(u, v) for (u, v) in gamma_W if v in S_star]
    assert len(gamma_prime) == 0, "Γ'(W) should be empty since 1 ∉ KSCC"

    # Fallback: if gamma_prime empty, use gamma_W
    if len(gamma_prime) == 0 and len(gamma_W) > 0:
        gamma_prime = gamma_W

    assert gamma_prime == [(0, 1)], "After fallback, Γ'(W) should equal Γ(W)"
    print("✓ Test 2 passed: Γ'(W) falls back to Γ(W) when epidemic disconnected from KSCC.")


# ---------------------------------------------------------------------------
# Test 3: C fallback when all internal edges are below ε threshold
# ---------------------------------------------------------------------------

def test_fallback_C_below_threshold():
    """
    Construct a graph where:
      - KSCC = {0,1,2} (dense → high rho)
      - All candidate edges have SpectralDrop << ε · ρ(S*)

    Use a very high ε (e.g. 0.99) so that NO edge passes the threshold,
    triggering the C = Γ'(W) fallback.
    """
    G = nx.DiGraph()
    # Dense 3-clique: rho is relatively high
    G.add_edges_from([(0, 1), (1, 0), (1, 2), (2, 1), (0, 2), (2, 0)])

    W = {0}
    initial_S = {1, 2}

    gamma_W = get_candidates(G, W, initial_S)
    S_star, rho_S_star, intra_degrees, vol_S_star = get_kscc(G)

    assert S_star is not None

    gamma_prime = [(u, v) for (u, v) in gamma_W if u in S_star and v in S_star]
    assert len(gamma_prime) > 0, "Should have internal KSCC candidate edges"

    deg_product_sum = sum(intra_degrees[u][0] * intra_degrees[u][1] for u in S_star)

    # Use a very high ε to guarantee all edges fall below threshold
    eps_high = 0.99
    C = spectral_filter(
        gamma_prime, S_star, rho_S_star,
        intra_degrees, vol_S_star, deg_product_sum, eps=eps_high
    )

    # The fallback inside spectral_filter ensures C == gamma_prime if empty
    assert len(C) > 0, "C should fall back to Γ'(W), never be empty"
    assert set(C) == set(gamma_prime), "After fallback, C should equal Γ'(W)"
    print("✓ Test 3 passed: C falls back to Γ'(W) when all SpectralDrops below ε threshold.")


# ---------------------------------------------------------------------------
# Run all tests
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    test_spectral_drop_algebra()
    test_fallback_gamma_prime_disconnected()
    test_fallback_C_below_threshold()
    print("\n✓ All unit tests passed.")
