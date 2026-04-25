"""
experiments/exp2_efficiency.py
================================
Experiment 2 — Computational Efficiency

Claim: SG-SNIR is significantly faster than MaxExpectedH while maintaining quality.

Per-iteration logging:
  |Γ(W)|   — candidates before filtering
  |Γ'(W)|  — after KSCC + bridge filter
  |C|       — after spectral threshold
  eval_counts — SNIR simulations performed
  wall-clock time per iteration

email-EuAll: SG-SNIR runs to completion.
             MaxExpectedH: one iteration extrapolated → marked OOT.

Output: results/experiment2/<dataset>_efficiency.json
        results/experiment2/plots/

Usage
-----
  python -m experiments.exp2_efficiency
  python -m experiments.exp2_efficiency --datasets HIV p2p-Gnutella
"""

import argparse
import json
import time
from pathlib import Path
from typing import Dict, List

import networkx as nx
import numpy as np

from src.data_loader import (
    load_dataset, seed_fixed_nodes,
    assign_fixed_weights,
    PRIMARY_DATASETS, LARGE_DATASETS, SCALABILITY_ONLY,
)
from src.snir_model import SNIRParams, compute_influence_range
from src.sg_snir import sg_snir_blocking, get_candidates, get_kscc, spectral_filter
from src.baselines import max_expected_h_blocking

RESULTS_DIR = Path("results/experiment2")
RESULTS_DIR.mkdir(parents=True, exist_ok=True)
(RESULTS_DIR / "plots").mkdir(exist_ok=True)

N_INFECTED_SEEDS = 5
T_HORIZON        = 10
EPS              = 0.01
BASE_SEED        = 42
DEFAULT_K        = {"HIV": 20, "p2p-Gnutella": 20, "Wiki-Vote": 20,
                    "soc-Epinions1": 20, "email-EuAll": 20}

BASE_PARAMS = SNIRParams(alpha=0.033, beta=0.022, delta=0.020,
                          eta=0.140, gamma=0.315, xi=0.300)


def run_sgsnir_instrumented(
    G: nx.DiGraph,
    initial_S, initial_N, initial_I, initial_R,
    params: SNIRParams,
    k: int,
) -> List[Dict]:
    """
    Run SG-SNIR and record per-iteration filtering statistics.
    Returns list of dicts with keys:
      iteration, gamma_W_size, gamma_prime_size, C_size,
      eval_count, wall_time_s, H_after

    TECHNICAL DEBT: This function reimplements the SG-SNIR inner loop rather
    than calling sg_snir_blocking() with an instrumentation flag. If sg_snir.py
    is updated (e.g. the fallback logic), this function must be updated in sync.
    Tracked for future refactor: add an optional `record_sizes` flag to
    sg_snir_blocking() that returns (selected_edges, H_history, eval_counts,
    gamma_sizes, C_sizes) as a single call.
    """

    G = G.copy()
    W = set(initial_N) | set(initial_I)
    log = []

    for iteration in range(k):
        t_iter = time.time()

        # Step 2
        gamma_W = get_candidates(G, W, initial_S)
        if not gamma_W:
            break

        # Step 3
        S_star, rho, intra, vol, deg_prod = get_kscc(G)

        # Build Γ'(W)
        if S_star is not None:
            gamma_prime = [(u, v) for (u, v) in gamma_W if v in S_star]
        else:
            gamma_prime = []

        # Step 4
        if len(gamma_prime) == 0 and len(gamma_W) > 0:
            C = gamma_W
        else:
            C = spectral_filter(gamma_prime, S_star, rho, intra, vol, deg_prod, EPS)

        # Step 5
        best_edge, best_H = None, float("inf")
        for (u, v) in C:
            attrs = G[u][v].copy()
            G.remove_edge(u, v)
            H_e, _ = compute_influence_range(
                G, initial_S, initial_N, initial_I, initial_R, params, T_HORIZON
            )
            G.add_edge(u, v, **attrs)
            if H_e < best_H:
                best_H, best_edge = H_e, (u, v)

        if best_edge is None:
            break

        G.remove_edge(*best_edge)
        wall_time = time.time() - t_iter

        log.append({
            "iteration":      iteration,
            "gamma_W_size":   len(gamma_W),
            "gamma_prime_size": len(gamma_prime),
            "C_size":         len(C),
            "eval_count":     len(C),
            "wall_time_s":    round(wall_time, 4),
            "H_after":        round(best_H, 6),
        })

    return log


def run_experiment2(
    datasets: List[str],
    verbose: bool = False,
):
    print("\n" + "="*70)
    print("  EXPERIMENT 2 — Computational Efficiency")
    print("="*70)

    all_datasets = PRIMARY_DATASETS + LARGE_DATASETS + SCALABILITY_ONLY

    for dataset in datasets:
        print(f"\n{'─'*60}\n  Dataset: {dataset}")
        G_orig = load_dataset(dataset)
        assign_fixed_weights(G_orig)
        k = DEFAULT_K.get(dataset, 20)

        np.random.seed(BASE_SEED)
        initial_S, initial_N, initial_I, initial_R = seed_fixed_nodes(
            G_orig, n_infected=N_INFECTED_SEEDS, seed=BASE_SEED
        )

        # --- SG-SNIR (instrumented) ---
        print("  Running SG-SNIR (instrumented)...")
        sg_log = run_sgsnir_instrumented(
            G_orig, initial_S, initial_N, initial_I, initial_R,
            BASE_PARAMS, k,
        )

        meh_log = None
        meh_extrapolated = False

        if dataset in SCALABILITY_ONLY:
            # Run only one iteration of MaxExpectedH → extrapolate projected runtime.
            # We evaluate ALL edges in gamma_W (not a capped subset) so the per-iteration
            # time estimate is accurate. On email-EuAll gamma_W may be very large, but
            # we need the true cost — a capped sample would give a misleadingly optimistic
            # projection and could suppress a correct OOT flag.
            print("  MaxExpectedH: running 1 full iteration to estimate cost...")
            G_tmp = G_orig.copy()
            gamma_W = get_candidates(G_tmp, set(initial_N) | set(initial_I), initial_S)
            t0 = time.time()
            for (u, v) in gamma_W:
                attrs = G_tmp[u][v].copy()
                G_tmp.remove_edge(u, v)
                H_e, _ = compute_influence_range(
                    G_tmp, initial_S, initial_N, initial_I, initial_R,
                    BASE_PARAMS, T_HORIZON
                )
                G_tmp.add_edge(u, v, **attrs)
            per_iter_s = time.time() - t0
            meh_log = {
                "extrapolated": True,
                "gamma_W_size": len(gamma_W),
                "one_iter_time_s": round(per_iter_s, 3),
                "projected_total_time_s": round(per_iter_s * k, 1),
            }
            meh_extrapolated = True
            print(f"  MaxExpectedH: 1-iter={per_iter_s:.1f}s → projected k={k}: {per_iter_s*k:.0f}s (OOT)")
        else:
            print("  Running MaxExpectedH...")
            t0 = time.time()
            _, meh_H, meh_evals = max_expected_h_blocking(
                G_orig, initial_S, initial_N, initial_I, initial_R,
                BASE_PARAMS, k, T_HORIZON, verbose=verbose,
            )
            total_t = time.time() - t0
            meh_log = {
                "extrapolated": False,
                "H_history": meh_H,
                "eval_counts": meh_evals,
                "total_time_s": round(total_t, 3),
            }

        payload = {
            "dataset": dataset,
            "k": k,
            "sg_snir_log": sg_log,
            "max_expected_h_log": meh_log,
            "meh_extrapolated": meh_extrapolated,
        }
        out_path = RESULTS_DIR / f"{dataset}_efficiency.json"
        with open(out_path, "w") as f:
            json.dump(payload, f, indent=2)

        total_sg_evals = sum(r["eval_count"] for r in sg_log)
        total_meh_evals = (
            sum(meh_log.get("eval_counts", [])) if not meh_extrapolated
            else meh_log.get("gamma_W_size", "?")
        )
        print(f"  ✓ SG-SNIR total evals: {total_sg_evals}")
        print(f"  ✓ MaxExpH total evals: {total_meh_evals}")
        print(f"  ✓ Saved → {out_path}")


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--datasets", nargs="+",
                   default=PRIMARY_DATASETS + LARGE_DATASETS + SCALABILITY_ONLY)
    p.add_argument("--verbose", action="store_true")
    return p.parse_args()


if __name__ == "__main__":
    args = parse_args()
    run_experiment2(datasets=args.datasets, verbose=args.verbose)
