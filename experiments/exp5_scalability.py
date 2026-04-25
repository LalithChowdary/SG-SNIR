"""
experiments/exp5_scalability.py
=================================
Experiment 5 — Scalability Analysis

Claim: SG-SNIR scales to large networks where MaxExpectedH is impractical.

Datasets: All five, ordered by |V|+|E|
Setup   : k=20, T=10, seed=42, fixed transmission probability
OOT     : MaxExpectedH flagged OOT if projected total runtime exceeds
           OOT_SECONDS (default: 3600s = 1 hour for the full k=20 run).
           This is more practical than the 48h Paper-2 convention —
           state the threshold explicitly in the paper.
           In practice: if one iteration takes >3 minutes, k=20 is OOT.
           (We run 1 real iteration to measure, then extrapolate.)

Output  : results/experiment5/scalability_summary.json

Usage
-----
  python -m experiments.exp5_scalability
  python -m experiments.exp5_scalability --maxh-cap-iters 1
"""

import argparse
import json
import time
from pathlib import Path
from typing import List, Optional

import numpy as np

from src.data_loader import (
    load_dataset, seed_fixed_nodes, assign_fixed_weights,
    ALL_DATASETS,
)
from src.snir_model import SNIRParams, compute_influence_range
from src.sg_snir import sg_snir_blocking, get_candidates
from src.baselines import max_expected_h_blocking

RESULTS_DIR = Path("results/experiment5")
RESULTS_DIR.mkdir(parents=True, exist_ok=True)

K_FIXED      = 20
T_HORIZON    = 10
N_SEEDS      = 5
BASE_SEED    = 42
# OOT threshold: flag MaxExpectedH as Out-of-Time if projected total runtime
# for k=20 exceeds this value. 3600s (1 hour) is practical — far more useful
# than 48h since experiments must finish overnight. State this in the paper.
OOT_SECONDS  = 3600
BASE_PARAMS  = SNIRParams(alpha=0.033, beta=0.022, delta=0.020,
                           eta=0.140, gamma=0.315, xi=0.300)

# Datasets ordered by (|V| + |E|) ascending
DATASET_ORDER = ["HIV", "p2p-Gnutella", "Wiki-Vote", "soc-Epinions1", "email-EuAll"]


def estimate_meh_time(G, initial_S, initial_N, initial_I, initial_R,
                      params, n_iters: int = 1) -> float:
    """Run n_iters of MaxExpectedH, return time-per-iter (seconds)."""
    G = G.copy()
    W = set(initial_N) | set(initial_I)
    total = 0.0
    for _ in range(n_iters):
        gamma_W = get_candidates(G, W, initial_S)
        if not gamma_W:
            break
        t0 = time.time()
        best_H, best_e = float("inf"), None
        for (u, v) in gamma_W:
            attrs = G[u][v].copy()
            G.remove_edge(u, v)
            H_e, _ = compute_influence_range(
                G, initial_S, initial_N, initial_I, initial_R, params, T_HORIZON
            )
            G.add_edge(u, v, **attrs)
            if H_e < best_H:
                best_H, best_e = H_e, (u, v)
        total += time.time() - t0
        if best_e:
            G.remove_edge(*best_e)
    return total / max(n_iters, 1)


def run_experiment5(
    datasets: List[str] = DATASET_ORDER,
    maxh_cap_iters: int = 2,
):
    print("\n" + "="*70)
    print("  EXPERIMENT 5 — Scalability Analysis")
    print("="*70)

    summary = []

    for dataset in datasets:
        print(f"\n{'─'*60}\n  Dataset: {dataset}")
        G = load_dataset(dataset)
        n_nodes = G.number_of_nodes()
        n_edges = G.number_of_edges()
        assign_fixed_weights(G)

        np.random.seed(BASE_SEED)
        initial_S, initial_N, initial_I, initial_R = seed_fixed_nodes(
            G, n_infected=N_SEEDS, seed=BASE_SEED
        )

        # --- SG-SNIR ---
        print("  Running SG-SNIR...")
        t0 = time.time()
        _, sg_H, sg_evals = sg_snir_blocking(
            G, initial_S, initial_N, initial_I, initial_R,
            BASE_PARAMS, K_FIXED, T_HORIZON, eps=0.01,
        )
        sg_time = round(time.time() - t0, 3)
        print(f"  SG-SNIR: {sg_time:.1f}s")

        # --- MaxExpectedH ---
        per_iter_s = estimate_meh_time(
            G, initial_S, initial_N, initial_I, initial_R,
            BASE_PARAMS, n_iters=min(maxh_cap_iters, 1)
        )
        projected_total = per_iter_s * K_FIXED
        oot = projected_total > OOT_SECONDS

        if oot:
            meh_time     = None
            meh_time_str = "OOT"
            speedup      = None
        else:
            print("  Running MaxExpectedH (full)...")
            t0 = time.time()
            max_expected_h_blocking(
                G, initial_S, initial_N, initial_I, initial_R,
                BASE_PARAMS, K_FIXED, T_HORIZON,
            )
            meh_time     = round(time.time() - t0, 3)
            meh_time_str = f"{meh_time}s"
            speedup      = round(meh_time / max(sg_time, 0.001), 2)

        row = {
            "dataset":             dataset,
            "n_nodes":             n_nodes,
            "n_edges":             n_edges,
            "size":                n_nodes + n_edges,
            "sg_snir_time_s":      sg_time,
            "meh_time_s":          meh_time,
            "meh_oot":             oot,
            "meh_projected_s":     round(projected_total, 1) if oot else None,
            "speedup":             speedup,
        }
        summary.append(row)
        print(f"  SG-SNIR={sg_time}s  MaxExpH={meh_time_str}  "
              f"Speedup={speedup if speedup else '—'}")

    out_path = RESULTS_DIR / "scalability_summary.json"
    with open(out_path, "w") as f:
        json.dump({"k": K_FIXED, "rows": summary}, f, indent=2)
    print(f"\n  ✓ Saved → {out_path}")

    # Print table
    print(f"\n  {'Dataset':<18} {'V':>8} {'E':>10} "
          f"{'SG-SNIR':>10} {'MaxExpH':>12} {'Speedup':>8}")
    print(f"  {'─'*18} {'─'*8} {'─'*10} {'─'*10} {'─'*12} {'─'*8}")
    for row in summary:
        meh_str = "OOT" if row["meh_oot"] else f"{row['meh_time_s']}s"
        spd_str = "—" if row["speedup"] is None else f"{row['speedup']}×"
        print(f"  {row['dataset']:<18} {row['n_nodes']:>8,} {row['n_edges']:>10,} "
              f"{row['sg_snir_time_s']:>9.1f}s {meh_str:>12} {spd_str:>8}")


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--datasets", nargs="+", default=DATASET_ORDER)
    p.add_argument("--maxh-cap-iters", type=int, default=1,
                   help="Iterations of MaxExpH to run before extrapolating (default: 1)")
    return p.parse_args()


if __name__ == "__main__":
    args = parse_args()
    run_experiment5(datasets=args.datasets, maxh_cap_iters=args.maxh_cap_iters)
