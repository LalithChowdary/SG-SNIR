"""
experiments/exp1_effectiveness.py
==================================
Experiment 1 — Primary Effectiveness: H vs k

Claim: SG-SNIR achieves epidemic containment quality comparable to
MaxExpectedH on directed networks.

Datasets  : HIV, p2p-Gnutella, Wiki-Vote, soc-Epinions1 (if timing OK)
Baselines : SG-SNIR, MaxExpectedH, Random, DegreeProduct
Output    : results/experiment1/<dataset>_<p_mode>_<eta_range>_H_vs_k.json
            results/experiment1/plots/ (PNG panels)

Usage
-----
  python -m experiments.exp1_effectiveness
  python -m experiments.exp1_effectiveness --datasets HIV p2p-Gnutella
  python -m experiments.exp1_effectiveness --skip-maxh      # fast dry run
"""

import argparse
import json
import time
from pathlib import Path
from typing import Dict, List, Optional

import numpy as np

from src.data_loader import (
    load_dataset, seed_fixed_nodes,
    assign_fixed_weights, assign_variable_weights,
    PRIMARY_DATASETS, LARGE_DATASETS,
)
from src.snir_model import SNIRParams, compute_influence_range
from src.sg_snir import sg_snir_blocking
from src.baselines import (
    max_expected_h_blocking, random_blocking, degree_product_blocking,
)

RESULTS_DIR = Path("results/experiment1")
RESULTS_DIR.mkdir(parents=True, exist_ok=True)
(RESULTS_DIR / "plots").mkdir(exist_ok=True)

# ---------------------------------------------------------------------------
# Plan constants
# ---------------------------------------------------------------------------

# Transition probability ranges mid-points (Paper 1 Table)
ETA_CONFIGS = {
    "eta_low":  0.140,   # mid of [0.1, 0.18]
    "eta_high": 0.215,   # mid of [0.18, 0.25]
}
# Note: SNIRParams constructed directly inside the eta loop — see run_experiment1.

BUDGET_K = {
    "HIV":           [5, 10, 20, 30, 50],
    "p2p-Gnutella":  [5, 10, 20],
    "Wiki-Vote":     [5, 10, 20],
    "soc-Epinions1": [10, 20, 50],
}

TRIAL_COUNT = {
    "HIV":           10,
    "p2p-Gnutella":  10,
    "Wiki-Vote":     5,
    "soc-Epinions1": 3,
}

N_INFECTED_SEEDS = 5
T_HORIZON       = 10
EPS             = 0.01
BASE_SEED       = 42


# ---------------------------------------------------------------------------
# Single trial runner
# ---------------------------------------------------------------------------

def run_trial(
    G_orig,
    dataset: str,
    p_mode: str,
    eta_label: str,
    params: SNIRParams,
    k_values: List[int],
    trial_idx: int,
    seed: int,
    run_maxh: bool,
    verbose: bool,
) -> Dict:
    """Run one trial; returns H_history keyed by method name at each k value."""

    import copy
    G = copy.deepcopy(G_orig)

    # Assign transmission probabilities
    if p_mode == "fixed":
        assign_fixed_weights(G)
    else:
        assign_variable_weights(G)

    # Seed initial states
    initial_S, initial_N, initial_I, initial_R = seed_fixed_nodes(
        G, n_infected=N_INFECTED_SEEDS, seed=seed
    )

    k_max = max(k_values)
    result = {"trial": trial_idx, "seed": seed, "methods": {}}

    # --- SG-SNIR ---
    t0 = time.time()
    _, sg_H, sg_evals = sg_snir_blocking(
        G, initial_S, initial_N, initial_I, initial_R,
        params, k_max, T_HORIZON, EPS, verbose=verbose,
    )
    result["methods"]["SG-SNIR"] = {
        "H_history": sg_H,
        "eval_counts": sg_evals,
        "time_s": round(time.time() - t0, 3),
    }

    # --- MaxExpectedH ---
    if run_maxh:
        t0 = time.time()
        _, meh_H, meh_evals = max_expected_h_blocking(
            G, initial_S, initial_N, initial_I, initial_R,
            params, k_max, T_HORIZON, verbose=verbose,
        )
        result["methods"]["MaxExpectedH"] = {
            "H_history": meh_H,
            "eval_counts": meh_evals,
            "time_s": round(time.time() - t0, 3),
        }

    # --- Random ---
    t0 = time.time()
    _, rand_H, _ = random_blocking(
        G, initial_S, initial_N, initial_I, initial_R,
        params, k_max, T_HORIZON, seed=seed,
    )
    result["methods"]["Random"] = {
        "H_history": rand_H,
        "time_s": round(time.time() - t0, 3),
    }

    # --- DegreeProduct ---
    # degree_product_blocking selects edges greedily by d_out(u)*d_in(v) with
    # no SNIR simulation.  We must run compute_influence_range after EVERY
    # removal so H_history has the same length as SG-SNIR and MaxExpectedH,
    # making all curves directly comparable on the same x-axis.
    t0 = time.time()
    deg_edges, _, _ = degree_product_blocking(
        G, initial_S, initial_N, initial_I, initial_R, k_max,
    )
    deg_H = []
    G_deg = G.copy()
    for edge in deg_edges:
        if G_deg.has_edge(*edge):
            G_deg.remove_edge(*edge)
        H_k, _ = compute_influence_range(
            G_deg, initial_S, initial_N, initial_I, initial_R, params, T_HORIZON
        )
        deg_H.append(H_k)
    result["methods"]["DegreeProduct"] = {
        "H_history": deg_H,
        "time_s": round(time.time() - t0, 3),
    }

    return result


# ---------------------------------------------------------------------------
# Main experiment loop
# ---------------------------------------------------------------------------

def run_experiment1(
    datasets: List[str],
    run_maxh: bool = True,
    verbose: bool = False,
):
    print("\n" + "="*70)
    print("  EXPERIMENT 1 — Primary Effectiveness: H vs k")
    print("="*70)

    for dataset in datasets:
        print(f"\n{'─'*60}")
        print(f"  Dataset: {dataset}")

        G_orig = load_dataset(dataset)

        for eta_label, eta_val in ETA_CONFIGS.items():
            # Construct SNIRParams directly — consistent with all other experiment files.
            # kappa=1.0 default: N-state nodes transmit equally to I-state nodes.
            params = SNIRParams(
                alpha=0.033, beta=0.022, delta=0.020,
                eta=eta_val,  gamma=0.315, xi=0.300,
            )

            for p_mode in ["fixed", "variable"]:
                k_values = BUDGET_K.get(dataset, [5, 10, 20])
                n_trials  = TRIAL_COUNT.get(dataset, 5)

                print(f"\n  [{dataset}] p_mode={p_mode}, {eta_label}")
                all_trials = []

                for trial in range(n_trials):
                    seed = BASE_SEED + trial * 7   # deterministic per-trial seeds
                    np.random.seed(seed)
                    print(f"    Trial {trial + 1}/{n_trials} (seed={seed})")

                    res = run_trial(
                        G_orig, dataset, p_mode, eta_label, params,
                        k_values, trial, seed, run_maxh, verbose,
                    )
                    all_trials.append(res)

                # Save raw results
                out_key = f"{dataset}_{p_mode}_{eta_label}"
                out_path = RESULTS_DIR / f"{out_key}.json"
                payload = {
                    "dataset": dataset,
                    "p_mode": p_mode,
                    "eta_label": eta_label,
                    "eta_value": eta_val,
                    "k_values": k_values,
                    "n_trials": n_trials,
                    "trials": all_trials,
                }
                with open(out_path, "w") as f:
                    json.dump(payload, f, indent=2)
                print(f"    ✓ Saved → {out_path}")


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--datasets", nargs="+",
                   default=PRIMARY_DATASETS + LARGE_DATASETS,
                   help="Datasets to run (default: all qualifying)")
    p.add_argument("--skip-maxh", action="store_true",
                   help="Skip MaxExpectedH baseline (fast dry run)")
    p.add_argument("--verbose", action="store_true")
    return p.parse_args()


if __name__ == "__main__":
    args = parse_args()
    run_experiment1(
        datasets=args.datasets,
        run_maxh=not args.skip_maxh,
        verbose=args.verbose,
    )
