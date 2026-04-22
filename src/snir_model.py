"""
src/snir_model.py
=================
Core SNIR propagation model for directed networks.

State Definitions:
    S  — Susceptible: not infected, not spreading
    N  — Non-symptomatic: infected but no symptoms, can spread
    I  — Infectious: infected with symptoms, can spread
    R  — Recovered: has antibodies, may return to S over time

Transition Probabilities (baseline):
    S → N  with prob α  (asymptomatic exposure)
    S → I  with prob β  (symptomatic exposure)
    N → I  with prob δ
    N → R  with prob η
    I → R  with prob γ
    R → S  with prob ξ

Directed extension:
    r(u, t+1) = 1 - ∏_{v:(v→u)∈E} [ 1 - p_{v,u} · a(v, t) ]
    Only predecessors (in-neighbours) can transmit to u.
"""

import networkx as nx
import numpy as np
from typing import Dict, Tuple, Set


class SNIRParams:
    """Container for all SNIR transition probabilities."""

    def __init__(
        self,
        alpha: float = 0.033,   # S → N (asymptomatic exposure)
        beta: float  = 0.022,   # S → I (symptomatic exposure)
        delta: float = 0.020,   # N → I (progression)
        eta: float   = 0.140,   # N → R (asymptomatic recovery)
        gamma: float = 0.315,   # I → R (symptomatic recovery)
        xi: float    = 0.300,   # R → S (waning immunity)
        kappa: float = 1.0,     # N-state transmission fraction:
                                # a(v,t) = pI(v,t) + κ · pN(v,t)
                                # κ=0: N nodes are fully latent (no transmission)
                                # κ=1: N and I nodes transmit equally (default)
    ):
        self.alpha = alpha
        self.beta  = beta
        self.delta = delta
        self.eta   = eta
        self.gamma = gamma
        self.xi    = xi
        self.kappa = kappa


def compute_influence_range(
    G: nx.DiGraph,
    initial_S: Set,
    initial_N: Set,
    initial_I: Set,
    initial_R: Set,
    params: SNIRParams,
    T: int = 10,
    default_weight: float = 0.05,
) -> Tuple[float, Dict]:
    """
    Algorithm 1: Compute I(G, Ω) = Σ_{t=0}^{T} Σ_u [ pN(u,t) + pI(u,t) ]

    This is the cumulative expected number of infected-state node-steps.
    Lower H(e) ⟹ better containment.

    Parameters
    ----------
    G              : Directed graph. Edge (v→u) means v can infect u.
                     Optional edge attribute 'weight' = p_{v,u}.
    initial_S/N/I/R: Node sets for initial state partition.
    params         : SNIRParams instance.
    T              : Simulation horizon (time steps).
    default_weight : Transmission probability when edge has no 'weight'.

    Returns
    -------
    H        : Scalar I(G,Ω) value.
    history  : {t: {node: {S,N,I,R probs}}} for all t in [0..T].
    """
    nodes = list(G.nodes())
    α, β = params.alpha, params.beta
    δ, η = params.delta, params.eta
    γ, ξ = params.gamma, params.xi
    κ    = params.kappa  # N-state transmission fraction

    # --- Initialise state probability vectors ---
    pS: Dict = {}
    pN: Dict = {}
    pI: Dict = {}
    pR: Dict = {}
    a:  Dict = {}   # a(u,t) = pI + pN  (actively infectious)

    for u in nodes:
        if u in initial_I:
            pS[u], pN[u], pI[u], pR[u] = 0.0, 0.0, 1.0, 0.0
        elif u in initial_N:
            pS[u], pN[u], pI[u], pR[u] = 0.0, 1.0, 0.0, 0.0
        elif u in initial_R:
            pS[u], pN[u], pI[u], pR[u] = 0.0, 0.0, 0.0, 1.0
        else:
            pS[u], pN[u], pI[u], pR[u] = 1.0, 0.0, 0.0, 0.0
        a[u] = pI[u] + κ * pN[u]  # a(u,t) = pI + κ·pN

    H = sum(pI[u] + pN[u] for u in nodes)
    history = {0: {u: {"S": pS[u], "N": pN[u], "I": pI[u], "R": pR[u]} for u in nodes}}

    # Active frontier: only nodes reachable from sources need updates
    active_union = set(initial_N) | set(initial_I)

    for t in range(1, T + 1):
        # Expand frontier by one hop
        new_neighbors: Set = set()
        for v in active_union:
            new_neighbors.update(G.successors(v))
        active_union = active_union | new_neighbors

        new_pS: Dict = {}
        new_pN: Dict = {}
        new_pI: Dict = {}
        new_pR: Dict = {}
        new_a:  Dict = {}

        for u in nodes:
            if u not in active_union:
                new_pS[u] = pS[u]
                new_pN[u] = pN[u]
                new_pI[u] = pI[u]
                new_pR[u] = pR[u]
                new_a[u]  = a[u]
                continue

            # r(u,t+1) = 1 - ∏_{v:(v→u)∈E} [1 - p_{v,u}·a(v,t)]
            r_u = 1.0
            for v in G.predecessors(u):
                p_vu = G[v][u].get("weight", default_weight)
                r_u *= (1.0 - p_vu * a[v])
            r_u = 1.0 - r_u

            # State update equations (from paper Section II.A)
            new_pI_u = (1 - γ) * pI[u] + r_u * β * pS[u] + δ * pN[u]
            new_pR_u = (1 - ξ) * pR[u] + γ * pI[u] + η * pN[u]
            new_pN_u = (1 - δ - η) * pN[u] + α * r_u * pS[u]
            new_pS_u = 1.0 - new_pN_u - new_pI_u - new_pR_u
            new_pS_u = max(0.0, min(1.0, new_pS_u))

            new_pI[u] = new_pI_u
            new_pR[u] = new_pR_u
            new_pN[u] = new_pN_u
            new_pS[u] = new_pS_u
            new_a[u]  = new_pI_u + κ * new_pN_u  # a(u,t) = pI + κ·pN

        pS, pN, pI, pR, a = new_pS, new_pN, new_pI, new_pR, new_a
        H += sum(pI[u] + pN[u] for u in nodes)
        history[t] = {u: {"S": pS[u], "N": pN[u], "I": pI[u], "R": pR[u]} for u in nodes}

    return H, history
