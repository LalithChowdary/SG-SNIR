"""
experiments/exp4_epsilon.py
=============================
Experiment 4 — Sensitivity to ε

Claim: ε = 0.01 is a robust default balancing quality and efficiency.

Datasets : HIV, Wiki-Vote
Setup    : Fix k=10, seed=42, vary ε ∈ {0.001, 0.005, 0.01, 0.02, 0.05,
           0.1, 0.2, 0.5}
Output   : results/experiment4/<dataset>_epsilon_sensitivity.json

Visualization (dual-axis):
  Left y : H_final
  Right y: total eval_count
  x-axis : ε on log scale

Usage
-----
  python -m experiments.exp4_epsilon
"""

import argparse
import json
import time
from pathlib import Path
from typing import List

import numpy as np

from src.data_loader import load_dataset, seed_fixed_nodes, assign_fixed_weights
from src.snir_model import SNIRParams
from src.sg_snir import sg_snir_blocking

RESULTS_DIR = Path("results/experiment4")
RESULTS_DIR.mkdir(parents=True, exist_ok=True)

DATASETS     = ["HIV", "Wiki-Vote"]
K_FIXED      = 10
N_SEEDS      = 5
T_HORIZON    = 10
BASE_SEED    = 42
BASE_PARAMS  = SNIRParams(alpha=0.033, beta=0.022, delta=0.020,
                           eta=0.140, gamma=0.315, xi=0.300)

EPS_VALUES = [0.001, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5]


def run_experiment4(datasets: List[str]):
    print("\n" + "="*70)
    print("  EXPERIMENT 4 — Sensitivity to ε")
    print("="*70)

    for dataset in datasets:
        print(f"\n{'─'*60}\n  Dataset: {dataset}")
        G = load_dataset(dataset)
        assign_fixed_weights(G)

        np.random.seed(BASE_SEED)
        initial_S, initial_N, initial_I, initial_R = seed_fixed_nodes(
            G, n_infected=N_SEEDS, seed=BASE_SEED
        )

        records = []
        print(f"  {'ε':>8}  {'H_final':>10}  {'eval_counts':>12}  {'time_s':>8}")
        print(f"  {'─'*8}  {'─'*10}  {'─'*12}  {'─'*8}")

        for eps in EPS_VALUES:
            t0 = time.time()
            _, H_history, eval_counts = sg_snir_blocking(
                G, initial_S, initial_N, initial_I, initial_R,
                BASE_PARAMS, K_FIXED, T_HORIZON, eps=eps,
            )
            elapsed = round(time.time() - t0, 3)
            H_final     = H_history[-1] if H_history else float("nan")
            total_evals = sum(eval_counts)

            records.append({
                "eps":         eps,
                "H_final":     round(H_final, 6),
                "H_history":   H_history,
                "eval_counts": eval_counts,
                "total_evals": total_evals,
                "time_s":      elapsed,
            })
            print(f"  {eps:>8.3f}  {H_final:>10.4f}  {total_evals:>12}  {elapsed:>8.3f}s")

        out_path = RESULTS_DIR / f"{dataset}_epsilon_sensitivity.json"
        with open(out_path, "w") as f:
            json.dump({"dataset": dataset, "k": K_FIXED, "records": records}, f, indent=2)
        print(f"  ✓ Saved → {out_path}")


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--datasets", nargs="+", default=DATASETS)
    return p.parse_args()


if __name__ == "__main__":
    args = parse_args()
    run_experiment4(datasets=args.datasets)
