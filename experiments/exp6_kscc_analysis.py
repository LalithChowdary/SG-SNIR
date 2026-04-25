"""
experiments/exp6_kscc_analysis.py
====================================
Experiment 6 — KSCC Structure Analysis

Claim: Efficiency gain correlates with KSCC concentration — networks with
smaller KSCCs benefit more from filtering.

For each dataset:
  - Compute KSCC size and fraction (|S*| / |V|)
  - Compute mean filtering ratio = mean(|C|) / mean(|Γ(W)|)
  - Reference speedup from Experiment 2 results

Output: results/experiment6/kscc_structure.json

Usage
-----
  python -m experiments.exp6_kscc_analysis
"""

import json
import time
from pathlib import Path
from typing import List

import networkx as nx
import numpy as np

from src.data_loader import (
    load_dataset, seed_fixed_nodes, assign_fixed_weights,
    PRIMARY_DATASETS, LARGE_DATASETS,
)
from src.sg_snir import get_candidates, get_kscc, spectral_filter

RESULTS_DIR = Path("results/experiment6")
RESULTS_DIR.mkdir(parents=True, exist_ok=True)
EXP2_DIR    = Path("results/experiment2")

DATASETS  = PRIMARY_DATASETS + LARGE_DATASETS
K_FIXED   = 20
N_SEEDS   = 5
T_HORIZON = 10
EPS       = 0.01
BASE_SEED = 42


def run_experiment6(datasets: List[str] = DATASETS):
    print("\n" + "="*70)
    print("  EXPERIMENT 6 — KSCC Structure Analysis")
    print("="*70)

    results = []

    for dataset in datasets:
        print(f"\n{'─'*60}\n  Dataset: {dataset}")
        G_orig = load_dataset(dataset)
        assign_fixed_weights(G_orig)

        np.random.seed(BASE_SEED)
        initial_S, initial_N, initial_I, initial_R = seed_fixed_nodes(
            G_orig, n_infected=N_SEEDS, seed=BASE_SEED
        )

        # Compute KSCC on original graph
        S_star, rho_star, intra, vol_star, deg_prod = get_kscc(G_orig)
        kscc_size     = len(S_star) if S_star else 0
        n_nodes       = G_orig.number_of_nodes()
        kscc_fraction = kscc_size / max(n_nodes, 1)

        # Simulate k iterations to gather filtering ratios
        G = G_orig.copy()
        W = set(initial_N) | set(initial_I)
        gamma_sizes, C_sizes = [], []

        for _ in range(min(K_FIXED, 10)):   # cap at 10 iters for analysis only
            gamma_W = get_candidates(G, W, initial_S)
            if not gamma_W:
                break

            S_star_i, rho_i, intra_i, vol_i, deg_i = get_kscc(G)

            if S_star_i is not None:
                gamma_prime = [(u, v) for (u, v) in gamma_W if v in S_star_i]
            else:
                gamma_prime = []

            if not gamma_prime:
                C = gamma_W
            else:
                C = spectral_filter(gamma_prime, S_star_i, rho_i, intra_i, vol_i, deg_i, EPS)

            gamma_sizes.append(len(gamma_W))
            C_sizes.append(len(C))

            # Use degree-product proxy to select the best edge rather than running
            # a full SNIR simulation. This avoids expensive SNIR evaluations during
            # what is purely a structural analysis loop. The filtering ratios
            # (|Γ(W)| vs |C|) are structurally determined by the KSCC topology and
            # are not sensitive to which specific edge is removed each iteration,
            # so the approximation is valid for the purpose of Experiment 6.
            if C:
                best = max(C, key=lambda e: G.out_degree(e[0]) * G.in_degree(e[1]))
                if G.has_edge(*best):
                    G.remove_edge(*best)

        mean_gamma = float(np.mean(gamma_sizes)) if gamma_sizes else 0.0
        mean_C     = float(np.mean(C_sizes))     if C_sizes     else 0.0
        filter_ratio = mean_C / max(mean_gamma, 1.0)

        # Load speedup from Experiment 2 if available
        speedup = None
        exp2_path = EXP2_DIR / f"{dataset}_efficiency.json"
        if exp2_path.exists():
            with open(exp2_path) as f:
                exp2 = json.load(f)
            sg_evals  = sum(r["eval_count"] for r in exp2.get("sg_snir_log", []))
            meh_log   = exp2.get("max_expected_h_log") or {}
            meh_evals = (
                sum(meh_log.get("eval_counts", []))
                if not exp2.get("meh_extrapolated", False)
                else meh_log.get("gamma_W_size", 0)
            )
            speedup = round(meh_evals / max(sg_evals, 1), 2) if sg_evals else None

        row = {
            "dataset":       dataset,
            "n_nodes":       n_nodes,
            "n_edges":       G_orig.number_of_edges(),
            "kscc_size":     kscc_size,
            "kscc_fraction": round(kscc_fraction, 4),
            "rho_star":      round(rho_star, 4),
            "mean_gamma_W":  round(mean_gamma, 1),
            "mean_C":        round(mean_C, 1),
            "filter_ratio":  round(filter_ratio, 4),
            "speedup_evals": speedup,
        }
        results.append(row)

        print(f"  KSCC size: {kscc_size} / {n_nodes} ({kscc_fraction:.2%})")
        print(f"  ρ(S*) = {rho_star:.4f}")
        print(f"  Filter ratio |C|/|Γ(W)| = {filter_ratio:.4f}")
        print(f"  Speedup (from Exp2) = {speedup if speedup else 'N/A (run Exp2 first)'}")

    out_path = RESULTS_DIR / "kscc_structure.json"
    with open(out_path, "w") as f:
        json.dump({"results": results}, f, indent=2)
    print(f"\n  ✓ Saved → {out_path}")

    # Print summary table
    print(f"\n  {'Dataset':<18} {'KSCC':>6} {'Frac%':>7} {'Γ(W)':>8} "
          f"{'|C|':>7} {'Ratio':>7} {'Speedup':>8}")
    print(f"  {'─'*18} {'─'*6} {'─'*7} {'─'*8} {'─'*7} {'─'*7} {'─'*8}")
    for r in results:
        spd = f"{r['speedup_evals']}×" if r['speedup_evals'] else "—"
        print(f"  {r['dataset']:<18} {r['kscc_size']:>6,} "
              f"{r['kscc_fraction']*100:>6.1f}% {r['mean_gamma_W']:>8.0f} "
              f"{r['mean_C']:>7.0f} {r['filter_ratio']:>7.3f} {spd:>8}")


if __name__ == "__main__":
    run_experiment6()
