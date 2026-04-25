"""
experiments/exp7_edge_analysis.py
====================================
Experiment 7 — KSCC Edge Selection Analysis

Claim: SG-SNIR selects structurally critical edges inside or entering the KSCC,
while MaxExpectedH wastes budget on peripheral edges.

For each selected edge, classify as:
  - "internal" : both u and v in S*
  - "bridge"   : u ∉ S*, v ∈ S*  (bridge into KSCC)
  - "peripheral": neither endpoint in S*

Also record SpectralDrop value for each selected edge.

Datasets : HIV, Wiki-Vote
Output   : results/experiment7/<dataset>_edge_analysis.json

Usage
-----
  python -m experiments.exp7_edge_analysis
"""

import json
import time
from pathlib import Path
from typing import Dict, FrozenSet, List, Optional, Set, Tuple

import numpy as np

from src.data_loader import load_dataset, seed_fixed_nodes, assign_fixed_weights
from src.snir_model import SNIRParams, compute_influence_range
from src.sg_snir import sg_snir_blocking, get_kscc
from src.baselines import max_expected_h_blocking

RESULTS_DIR = Path("results/experiment7")
RESULTS_DIR.mkdir(parents=True, exist_ok=True)

DATASETS    = ["HIV", "Wiki-Vote"]
K_FIXED     = 20
N_SEEDS     = 5
T_HORIZON   = 10
BASE_SEED   = 42
BASE_PARAMS = SNIRParams(alpha=0.033, beta=0.022, delta=0.020,
                          eta=0.140, gamma=0.315, xi=0.300)


def classify_edge(
    u, v,
    S_star: Optional[FrozenSet],
) -> str:
    if S_star is None:
        return "peripheral"
    u_in = u in S_star
    v_in = v in S_star
    if u_in and v_in:
        return "internal"
    elif (not u_in) and v_in:
        return "bridge"
    else:
        return "peripheral"


def compute_spectral_drop(
    u, v,
    S_star,
    intra: Dict,
    vol_star: float,
    deg_prod: float,
    rho_star: float,
) -> float:
    """
    Compute SpectralDrop (Eq. 17) for an edge.

    For internal KSCC edges (u ∈ S*, v ∈ S*): returns the estimated drop in
    spectral radius from removing this edge.
    For bridge edges (u ∉ S*, v ∈ S*) or peripheral edges: returns 0.0,
    which is their true SpectralDrop per the paper (bridge edges do not
    reduce the KSCC degree-product sum directly). Returning 0.0 instead of
    None keeps the plotting histogram clean without special null handling.
    """
    if S_star is None or u not in S_star or v not in S_star:
        return 0.0   # bridge or peripheral: SpectralDrop = 0 by definition
    if v not in intra or vol_star <= 1:
        return 0.0
    d_in_u  = intra[u][0]
    d_out_v = intra[v][1]
    new_num = deg_prod - d_in_u - d_out_v
    new_rho = new_num / (vol_star - 1)
    return rho_star - new_rho


def annotate_edges(
    selected_edges: List[Tuple],
    G_original,
) -> List[Dict]:
    """
    Classify and annotate each selected edge with structural properties.

    KSCC is computed once on the original graph (before any edge removals) as
    a fixed reference for consistent classification. This allows a fair apples-
    to-apples comparison between SG-SNIR and MaxExpectedH selections, since
    both are classified against the same structural baseline. In practice the
    KSCC evolves as edges are removed, but using a fixed reference avoids the
    confound of each method seeing a different KSCC due to its own removals.
    """
    S_star, rho_star, intra, vol_star, deg_prod = get_kscc(G_original)

    annotations = []
    for (u, v) in selected_edges:
        category      = classify_edge(u, v, S_star)
        spectral_drop = compute_spectral_drop(
            u, v, S_star, intra, vol_star, deg_prod, rho_star
        )
        annotations.append({
            "u":             u,
            "v":             v,
            "category":      category,
            "spectral_drop": round(spectral_drop, 6),
        })
    return annotations


def run_experiment7(datasets: List[str] = DATASETS, verbose: bool = False):
    print("\n" + "="*70)
    print("  EXPERIMENT 7 — KSCC Edge Selection Analysis")
    print("="*70)

    for dataset in datasets:
        print(f"\n{'─'*60}\n  Dataset: {dataset}")
        G_orig = load_dataset(dataset)
        assign_fixed_weights(G_orig)

        np.random.seed(BASE_SEED)
        initial_S, initial_N, initial_I, initial_R = seed_fixed_nodes(
            G_orig, n_infected=N_SEEDS, seed=BASE_SEED
        )

        # --- SG-SNIR ---
        print("  Running SG-SNIR...")
        t0 = time.time()
        sg_edges, sg_H, sg_evals = sg_snir_blocking(
            G_orig, initial_S, initial_N, initial_I, initial_R,
            BASE_PARAMS, K_FIXED, T_HORIZON, verbose=verbose,
        )
        sg_time = round(time.time() - t0, 3)

        sg_annotations = annotate_edges(sg_edges, G_orig)

        # --- MaxExpectedH ---
        print("  Running MaxExpectedH...")
        t0 = time.time()
        meh_edges, meh_H, _ = max_expected_h_blocking(
            G_orig, initial_S, initial_N, initial_I, initial_R,
            BASE_PARAMS, K_FIXED, T_HORIZON, verbose=verbose,
        )
        meh_time = round(time.time() - t0, 3)

        meh_annotations = annotate_edges(meh_edges, G_orig)

        # Compute category fractions
        def category_fractions(annotations: List[Dict]) -> Dict:
            total = max(len(annotations), 1)
            counts = {"internal": 0, "bridge": 0, "peripheral": 0}
            for a in annotations:
                counts[a["category"]] = counts.get(a["category"], 0) + 1
            return {k: round(v / total, 4) for k, v in counts.items()}

        sg_fractions  = category_fractions(sg_annotations)
        meh_fractions = category_fractions(meh_annotations)

        result = {
            "dataset": dataset,
            "k": K_FIXED,
            "SG-SNIR": {
                "time_s":          sg_time,
                "selected_edges":  sg_annotations,
                "category_fractions": sg_fractions,
                "H_final":         sg_H[-1] if sg_H else None,
            },
            "MaxExpectedH": {
                "time_s":          meh_time,
                "selected_edges":  meh_annotations,
                "category_fractions": meh_fractions,
                "H_final":         meh_H[-1] if meh_H else None,
            },
        }

        out_path = RESULTS_DIR / f"{dataset}_edge_analysis.json"
        with open(out_path, "w") as f:
            json.dump(result, f, indent=2)

        print(f"\n  SG-SNIR edge categories:   {sg_fractions}")
        print(f"  MaxExpH edge categories:   {meh_fractions}")
        print(f"  ✓ Saved → {out_path}")


if __name__ == "__main__":
    run_experiment7()
