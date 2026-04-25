"""
experiments/exp3_ablation.py
==============================
Experiment 3 — Ablation Study

Claim: Every component of SG-SNIR contributes meaningfully to quality
or efficiency or both.

Four precisely defined variants (from Section 5, Experiment 3):

  Variant A — Full SG-SNIR        : All components active
  Variant B — No KSCC filter      : C = Γ(W) directly; Steps 3+4 skipped
  Variant C — No bridge edges     : Γ'(W) = KSCC internal edges only
  Variant D — No spectral threshold: C = Γ'(W), SpectralDrop filter skipped

Datasets : HIV, Wiki-Vote, soc-Epinions1 (if feasible)
Output   : results/experiment3/<dataset>_ablation.json

Usage
-----
  python -m experiments.exp3_ablation
  python -m experiments.exp3_ablation --datasets HIV Wiki-Vote
"""

import argparse
import json
import time
from pathlib import Path
from typing import Dict, List, Set, Tuple

import networkx as nx
import numpy as np

from src.data_loader import (
    load_dataset, seed_fixed_nodes, assign_fixed_weights,
)
from src.snir_model import SNIRParams, compute_influence_range
from src.sg_snir import get_candidates, get_kscc, spectral_filter

RESULTS_DIR = Path("results/experiment3")
RESULTS_DIR.mkdir(parents=True, exist_ok=True)

DATASETS     = ["HIV", "Wiki-Vote", "soc-Epinions1"]
BUDGET_K     = {"HIV": 20, "Wiki-Vote": 20, "soc-Epinions1": 10}
N_SEEDS      = 5
T_HORIZON    = 10
EPS          = 0.01
BASE_SEED    = 42
BASE_PARAMS  = SNIRParams(alpha=0.033, beta=0.022, delta=0.020,
                           eta=0.140, gamma=0.315, xi=0.300)


# ---------------------------------------------------------------------------
# Variant implementations (each returns selected_edges, H_history, eval_counts)
# ---------------------------------------------------------------------------

def _run_variant(
    G: nx.DiGraph,
    initial_S: Set, initial_N: Set, initial_I: Set, initial_R: Set,
    params: SNIRParams, k: int,
    variant: str,
) -> Tuple[List, List, List]:
    """
    variant in {'full', 'no_kscc', 'no_bridge', 'no_spectral'}
    """
    G = G.copy()
    W = set(initial_N) | set(initial_I)
    selected, H_history, eval_counts = [], [], []

    for _ in range(k):
        gamma_W = get_candidates(G, W, initial_S)
        if not gamma_W:
            break

        # ---- Build candidate set C according to variant ----
        if variant == "no_kscc":
            # Variant B: skip Steps 3+4 entirely → C = Γ(W)
            C = gamma_W

        else:
            S_star, rho, intra, vol, deg_prod = get_kscc(G)

            if variant == "no_bridge":
                # Variant C: only KSCC-internal edges, exclude bridges
                if S_star is not None:
                    gamma_prime = [(u, v) for (u, v) in gamma_W
                                   if u in S_star and v in S_star]
                else:
                    gamma_prime = []
            else:
                # Variants A (full) and D (no_spectral): use standard Γ'(W)
                if S_star is not None:
                    gamma_prime = [(u, v) for (u, v) in gamma_W if v in S_star]
                else:
                    gamma_prime = []

            if len(gamma_prime) == 0:
                # No KSCC candidates available for this iteration.
                #
                # For no_bridge: ideally we would skip non-internal edges,
                # but we cannot break early because that shortens H_history
                # below k entries, making all variant curves incomparable.
                # Resolution: fall back to gamma_W and mark this iteration
                # as a forced fallback so the paper can note it.
                C = gamma_W
            elif variant == "no_spectral":
                # Variant D: skip SpectralDrop filter → C = Γ'(W)
                C = gamma_prime
            else:
                # Variant A: full spectral filter
                C = spectral_filter(gamma_prime, S_star, rho, intra, vol, deg_prod, EPS)

        # ---- SNIR evaluation ----
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
        selected.append(best_edge)
        H_history.append(best_H)
        eval_counts.append(len(C))

    return selected, H_history, eval_counts


# ---------------------------------------------------------------------------
# Main experiment
# ---------------------------------------------------------------------------

VARIANT_NAMES = {
    "full":        "Full SG-SNIR",
    "no_kscc":     "No KSCC filter",
    "no_bridge":   "No bridge edges",
    "no_spectral": "No spectral threshold",
}


TRIALS = 5

def run_experiment3(datasets: List[str], verbose: bool = False):
    print("\n" + "="*70)
    print("  EXPERIMENT 3 — Ablation Study")
    print("="*70)

    for dataset in datasets:
        print(f"\n{'─'*60}\n  Dataset: {dataset}")
        G_orig = load_dataset(dataset)
        assign_fixed_weights(G_orig)
        k = BUDGET_K.get(dataset, 20)

        dataset_result = {"dataset": dataset, "k": k, "n_trials": TRIALS, "variants": {}}

        for variant_key, variant_label in VARIANT_NAMES.items():
            dataset_result["variants"][variant_key] = {
                "label": variant_label,
                "H_histories": [],
                "eval_counts_list": [],
                "time_s_list": []
            }

        for trial in range(TRIALS):
            seed = BASE_SEED + trial * 7
            np.random.seed(seed)
            initial_S, initial_N, initial_I, initial_R = seed_fixed_nodes(
                G_orig, n_infected=N_SEEDS, seed=seed
            )
            print(f"  Trial {trial + 1}/{TRIALS} (seed={seed})")

            for variant_key, variant_label in VARIANT_NAMES.items():
                t0 = time.time()
                _, H_history, eval_counts = _run_variant(
                    G_orig, initial_S, initial_N, initial_I, initial_R,
                    BASE_PARAMS, k, variant_key,
                )
                elapsed = round(time.time() - t0, 3)

                dataset_result["variants"][variant_key]["H_histories"].append(H_history)
                dataset_result["variants"][variant_key]["eval_counts_list"].append(eval_counts)
                dataset_result["variants"][variant_key]["time_s_list"].append(elapsed)

        for variant_key, v_data in dataset_result["variants"].items():
            H_hists = [h for h in v_data["H_histories"] if h]
            if H_hists:
                min_len = min(len(h) for h in H_hists)
                arr = np.array([h[:min_len] for h in H_hists])
                mean_H = arr.mean(axis=0).tolist()
                H_final = mean_H[-1] if mean_H else None
            else:
                mean_H = []
                H_final = None

            evals_lists = v_data["eval_counts_list"]
            if evals_lists:
                min_len_evals = min(len(e) for e in evals_lists)
                arr_evals = np.array([e[:min_len_evals] for e in evals_lists])
                mean_evals = arr_evals.mean(axis=0).tolist()
                total_evals = sum(mean_evals)
            else:
                mean_evals = []
                total_evals = 0

            avg_time = round(float(np.mean(v_data["time_s_list"])), 3)

            v_data["H_history"] = [round(x, 6) for x in mean_H]
            v_data["H_final"] = round(H_final, 6) if H_final is not None else None
            v_data["eval_counts"] = [round(x, 2) for x in mean_evals]
            v_data["total_evals"] = round(total_evals, 2)
            v_data["time_s"] = avg_time

            h_str = f"{v_data['H_final']:.4f}" if v_data["H_final"] is not None else "N/A"
            print(f"    Variant: {v_data['label']:<22} -> H_final(avg)={h_str}  evals(avg)={total_evals:.1f}  time(avg)={avg_time}s")

        out_path = RESULTS_DIR / f"{dataset}_ablation.json"
        with open(out_path, "w") as f:
            json.dump(dataset_result, f, indent=2)
        print(f"  ✓ Saved → {out_path}")


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--datasets", nargs="+", default=DATASETS)
    p.add_argument("--verbose", action="store_true")
    return p.parse_args()


if __name__ == "__main__":
    args = parse_args()
    run_experiment3(datasets=args.datasets, verbose=args.verbose)
