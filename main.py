#!/usr/bin/env python3
"""
main.py
=======
SG-SNIR Experiment Runner

Usage examples:
  # Run on HIV network with budget k=20
  python main.py --dataset hiv --k 20

  # Run on all ER synthetic networks
  python main.py --dataset er --k 10

  # Run epsilon sensitivity study on HIV
  python main.py --dataset hiv --k 15 --eps-sweep

  # Disable baselines for speed
  python main.py --dataset hiv --k 10 --no-maxh --no-random

Configurable parameters:
  --k          : Edge blocking budget (default: 10)
  --T          : SNIR simulation horizon (default: 10)
  --eps        : Spectral threshold ε for SG-SNIR (default: 0.01)
  --alpha      : SNIR α — S→N probability (default: 0.033)
  --beta       : SNIR β — S→I probability (default: 0.022)
  --delta      : SNIR δ — N→I probability (default: 0.020)
  --eta        : SNIR η — N→R probability (default: 0.140)
  --gamma      : SNIR γ — I→R probability (default: 0.315)
  --xi         : SNIR ξ — R→S probability (default: 0.300)
  --kappa      : N-state transmission fraction κ (default: 1.0)
                 a(v,t) = pI(v,t) + κ·pN(v,t)
                 κ=0 → N nodes are fully latent; κ=1 → equal to I nodes
  --frac-I     : Fraction of nodes initially in I state (default: 0.02)
  --frac-N     : Fraction of nodes initially in N state (default: 0.02)
  --weight     : Default edge transmission probability (default: 0.05)
  --seed       : Random seed for reproducibility (default: 42)
  --er-n       : Override ER node count (single custom run)
  --er-p       : Override ER edge probability (single custom run)
  --eps-sweep  : Run sensitivity test over ε ∈ {0.01, 0.05, 0.10}
  --no-maxh    : Skip MaxExpectedH baseline (expensive)
  --no-random  : Skip Random baseline
  --no-degprod : Skip DegreeProduct baseline
  --no-plot    : Skip matplotlib plots (print table only)
"""

import argparse
import time
from typing import Dict, List, Optional, Tuple

import networkx as nx
import numpy as np

from src.data_loader import (
    ER_CONFIGS,
    generate_er_network,
    load_hiv_network,
    seed_initial_states,
)
from src.snir_model import SNIRParams, compute_influence_range
from src.sg_snir import sg_snir_blocking
from src.baselines import (
    degree_product_blocking,
    max_expected_h_blocking,
    random_blocking,
)


# ---------------------------------------------------------------------------
# Argument parsing
# ---------------------------------------------------------------------------

def parse_args():
    p = argparse.ArgumentParser(description="SG-SNIR Experiment Runner")

    p.add_argument("--dataset",  type=str, default="hiv",
                   choices=["hiv", "er"],
                   help="Dataset to use: 'hiv' or 'er' (all three ER sizes)")

    # Core algorithm parameters
    p.add_argument("--k",        type=int,   default=10,    help="Edge blocking budget")
    p.add_argument("--T",        type=int,   default=10,    help="SNIR simulation horizon")
    p.add_argument("--eps",      type=float, default=0.01,  help="Spectral threshold ε")
    p.add_argument("--weight",   type=float, default=0.05,  help="Default edge transmission probability")
    p.add_argument("--seed",     type=int,   default=42,    help="Random seed")

    # SNIR parameters (with defaults from published ranges)
    p.add_argument("--alpha",    type=float, default=0.033, help="S→N probability α")
    p.add_argument("--beta",     type=float, default=0.022, help="S→I probability β")
    p.add_argument("--delta",    type=float, default=0.020, help="N→I probability δ")
    p.add_argument("--eta",      type=float, default=0.140, help="N→R probability η")
    p.add_argument("--gamma",    type=float, default=0.315, help="I→R probability γ")
    p.add_argument("--xi",       type=float, default=0.300, help="R→S probability ξ")
    p.add_argument("--kappa",    type=float, default=1.0,
                   help="N-state transmission fraction κ: a(v,t)=pI+κ·pN (0=latent, 1=equal to I)")

    # Initial state seeding
    p.add_argument("--frac-I",   type=float, default=0.02,  help="Fraction of nodes in I state at t=0")
    p.add_argument("--frac-N",   type=float, default=0.02,  help="Fraction of nodes in N state at t=0")

    # Custom ER override
    p.add_argument("--er-n",     type=int,   default=None,  help="Override ER node count")
    p.add_argument("--er-p",     type=float, default=None,  help="Override ER edge probability")

    # Experiment flags
    p.add_argument("--eps-sweep",   action="store_true", help="Run ε sensitivity sweep")
    p.add_argument("--no-maxh",     action="store_true", help="Skip MaxExpectedH baseline")
    p.add_argument("--no-random",   action="store_true", help="Skip Random baseline")
    p.add_argument("--no-degprod",  action="store_true", help="Skip DegreeProduct baseline")
    p.add_argument("--no-plot",     action="store_true", help="Skip plots")
    p.add_argument("--verbose",     action="store_true", help="Print per-iteration details")

    return p.parse_args()


# ---------------------------------------------------------------------------
# Core experiment runner
# ---------------------------------------------------------------------------

def run_experiment(
    G: nx.DiGraph,
    network_name: str,
    params: SNIRParams,
    initial_S, initial_N, initial_I, initial_R,
    k: int,
    T: int,
    eps: float,
    weight: float,
    seed: int,
    run_maxh: bool,
    run_random: bool,
    run_degprod: bool,
    verbose: bool,
) -> Dict:
    """Run all algorithms on a single network and collect results."""

    results = {"network": network_name, "k": k}

    # Baseline H with no blocking
    H_baseline, _ = compute_influence_range(
        G, initial_S, initial_N, initial_I, initial_R, params, T, weight
    )
    results["H_baseline"] = H_baseline
    print(f"\n{'='*60}")
    print(f"Network: {network_name} | Nodes: {G.number_of_nodes()} | Edges: {G.number_of_edges()}")
    print(f"H_baseline (no blocking): {H_baseline:.4f}")
    print(f"{'='*60}")

    # --- SG-SNIR ---
    t0 = time.time()
    sg_edges, sg_H, sg_evals = sg_snir_blocking(
        G, initial_S, initial_N, initial_I, initial_R,
        params, k, T, eps, weight, verbose
    )
    sg_time = time.time() - t0
    results["sg_snir"] = {
        "H_final": sg_H[-1] if sg_H else H_baseline,
        "H_history": sg_H,
        "edges": sg_edges,
        "eval_counts": sg_evals,
        "time_s": sg_time,
        "total_evals": sum(sg_evals),
    }
    print(f"\n[SG-SNIR]       H_final={results['sg_snir']['H_final']:.4f}  "
          f"total_evals={results['sg_snir']['total_evals']}  time={sg_time:.2f}s")

    # --- MaxExpectedH ---
    if run_maxh:
        t0 = time.time()
        maxh_edges, maxh_H, maxh_evals = max_expected_h_blocking(
            G, initial_S, initial_N, initial_I, initial_R,
            params, k, T, weight, verbose
        )
        maxh_time = time.time() - t0
        results["max_expected_h"] = {
            "H_final": maxh_H[-1] if maxh_H else H_baseline,
            "H_history": maxh_H,
            "edges": maxh_edges,
            "eval_counts": maxh_evals,
            "time_s": maxh_time,
            "total_evals": sum(maxh_evals),
        }
        print(f"[MaxExpectedH]  H_final={results['max_expected_h']['H_final']:.4f}  "
              f"total_evals={results['max_expected_h']['total_evals']}  time={maxh_time:.2f}s")

    # --- Random ---
    if run_random:
        t0 = time.time()
        rand_edges, rand_H, _ = random_blocking(
            G, initial_S, initial_N, initial_I, initial_R,
            params, k, T, weight, seed, verbose
        )
        rand_time = time.time() - t0
        results["random"] = {
            "H_final": rand_H[-1] if rand_H else H_baseline,
            "H_history": rand_H,
            "edges": rand_edges,
            "time_s": rand_time,
        }
        print(f"[Random]        H_final={results['random']['H_final']:.4f}  time={rand_time:.2f}s")

    # --- DegreeProduct ---
    if run_degprod:
        t0 = time.time()
        deg_edges, _, _ = degree_product_blocking(
            G, initial_S, initial_N, initial_I, k, verbose
        )
        # Compute H after all removals for fair comparison
        G_deg = G.copy()
        for e in deg_edges:
            if G_deg.has_edge(*e):
                G_deg.remove_edge(*e)
        H_deg, _ = compute_influence_range(
            G_deg, initial_S, initial_N, initial_I, initial_R, params, T, weight
        )
        deg_time = time.time() - t0
        results["degree_product"] = {
            "H_final": H_deg,
            "edges": deg_edges,
            "time_s": deg_time,
        }
        print(f"[DegreeProduct] H_final={H_deg:.4f}  time={deg_time:.2f}s")

    return results


# ---------------------------------------------------------------------------
# Plotting
# ---------------------------------------------------------------------------

def plot_results(all_results: List[Dict], show: bool = True):
    try:
        import matplotlib.pyplot as plt
        import matplotlib.cm as cm
    except ImportError:
        print("matplotlib not installed — skipping plots.")
        return

    fig, axes = plt.subplots(1, len(all_results), figsize=(6 * len(all_results), 5), squeeze=False)

    for col, res in enumerate(all_results):
        ax = axes[0][col]
        k = res["k"]
        x = list(range(1, k + 1))

        if "sg_snir" in res and res["sg_snir"]["H_history"]:
            ax.plot(x[:len(res["sg_snir"]["H_history"])],
                    res["sg_snir"]["H_history"], "b-o", label="SG-SNIR", linewidth=2)

        if "max_expected_h" in res and res["max_expected_h"]["H_history"]:
            ax.plot(x[:len(res["max_expected_h"]["H_history"])],
                    res["max_expected_h"]["H_history"], "g--s", label="MaxExpectedH", linewidth=2)

        if "random" in res and res["random"]["H_history"]:
            ax.plot(x[:len(res["random"]["H_history"])],
                    res["random"]["H_history"], "r:^", label="Random", linewidth=2)

        ax.axhline(res["H_baseline"], color="gray", linestyle="--", alpha=0.5, label="No blocking")
        ax.set_xlabel("Budget k (edges blocked)")
        ax.set_ylabel("H = I(G, Ω)")
        ax.set_title(res["network"])
        ax.legend()
        ax.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig("results_H_vs_k.png", dpi=150, bbox_inches="tight")
    print("\nPlot saved to results_H_vs_k.png")
    if show:
        plt.show()


def plot_efficiency(all_results: List[Dict], show: bool = True):
    """Plot eval_counts per iteration: SG-SNIR vs MaxExpectedH."""
    try:
        import matplotlib.pyplot as plt
    except ImportError:
        return

    fig, axes = plt.subplots(1, len(all_results), figsize=(6 * len(all_results), 4), squeeze=False)

    for col, res in enumerate(all_results):
        ax = axes[0][col]
        if "sg_snir" in res:
            ax.bar([i - 0.2 for i in range(len(res["sg_snir"]["eval_counts"]))],
                   res["sg_snir"]["eval_counts"], width=0.4, label="SG-SNIR", color="steelblue")
        if "max_expected_h" in res:
            ax.bar([i + 0.2 for i in range(len(res["max_expected_h"]["eval_counts"]))],
                   res["max_expected_h"]["eval_counts"], width=0.4, label="MaxExpectedH", color="seagreen", alpha=0.7)

        ax.set_xlabel("Iteration")
        ax.set_ylabel("SNIR evaluations")
        ax.set_title(f"{res['network']} — Evaluation count")
        ax.legend()
        ax.grid(True, alpha=0.3, axis="y")

    plt.tight_layout()
    plt.savefig("results_eval_counts.png", dpi=150, bbox_inches="tight")
    print("Efficiency plot saved to results_eval_counts.png")
    if show:
        plt.show()


# ---------------------------------------------------------------------------
# Epsilon sensitivity sweep
# ---------------------------------------------------------------------------

def run_eps_sweep(G, network_name, params, initial_S, initial_N, initial_I, initial_R,
                  k, T, weight, verbose):
    eps_values = [0.01, 0.05, 0.10]
    print(f"\n{'='*60}")
    print(f"  ε sensitivity sweep — {network_name}")
    print(f"{'='*60}")
    for eps in eps_values:
        _, H_list, evals = sg_snir_blocking(
            G, initial_S, initial_N, initial_I, initial_R,
            params, k, T, eps, weight, verbose=False
        )
        H_final = H_list[-1] if H_list else float("nan")
        total_evals = sum(evals)
        print(f"  ε={eps:.2f}  H_final={H_final:.4f}  total_evals={total_evals}")


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

def main():
    args = parse_args()

    params = SNIRParams(
        alpha=args.alpha, beta=args.beta, delta=args.delta,
        eta=args.eta, gamma=args.gamma, xi=args.xi,
        kappa=args.kappa,
    )

    networks_to_run: List[Tuple[str, nx.DiGraph]] = []

    if args.dataset == "hiv":
        G = load_hiv_network(default_weight=args.weight)
        networks_to_run.append(("HIV-Transmission", G))

    elif args.dataset == "er":
        if args.er_n is not None and args.er_p is not None:
            # Custom override
            G = generate_er_network(args.er_n, args.er_p, args.seed, args.weight)
            networks_to_run.append((f"ER-Custom(n={args.er_n},p={args.er_p})", G))
        else:
            for name, cfg in ER_CONFIGS.items():
                G = generate_er_network(cfg["n"], cfg["p"], args.seed, args.weight)
                networks_to_run.append((name, G))

    all_results = []

    for (name, G) in networks_to_run:
        initial_S, initial_N, initial_I, initial_R = seed_initial_states(
            G, frac_I=args.frac_I, frac_N=args.frac_N, seed=args.seed
        )

        # Epsilon sensitivity sweep
        if args.eps_sweep:
            run_eps_sweep(G, name, params, initial_S, initial_N, initial_I, initial_R,
                          args.k, args.T, args.weight, args.verbose)

        res = run_experiment(
            G=G,
            network_name=name,
            params=params,
            initial_S=initial_S,
            initial_N=initial_N,
            initial_I=initial_I,
            initial_R=initial_R,
            k=args.k,
            T=args.T,
            eps=args.eps,
            weight=args.weight,
            seed=args.seed,
            run_maxh=not args.no_maxh,
            run_random=not args.no_random,
            run_degprod=not args.no_degprod,
            verbose=args.verbose,
        )
        all_results.append(res)

    if not args.no_plot:
        plot_results(all_results, show=False)
        if not args.no_maxh:
            plot_efficiency(all_results, show=False)


if __name__ == "__main__":
    main()
