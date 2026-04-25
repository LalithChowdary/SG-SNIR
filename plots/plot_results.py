"""
plots/plot_results.py
======================
Plotting script for all 8 SG-SNIR experiments.

Reads JSON results from results/experiment<N>/ and produces publication-ready
figures saved to plots/<exp_name>/<figure>.png.

Usage
-----
  # Plot all experiments
  python plots/plot_results.py

  # Plot specific experiments
  python plots/plot_results.py --experiments 1 2 3 4
"""

import argparse
import json
from pathlib import Path
from typing import Dict, List, Optional

import numpy as np

try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import matplotlib.ticker as ticker
    HAS_MPL = True
except ImportError:
    HAS_MPL = False
    print("WARNING: matplotlib not installed — cannot generate plots.")

RESULTS_ROOT = Path("results")
PLOTS_ROOT   = Path("plots")

METHOD_STYLES = {
    "SG-SNIR":      {"color": "#2563EB", "linestyle": "-",  "marker": "o", "lw": 2.5},
    "MaxExpectedH": {"color": "#16A34A", "linestyle": "--", "marker": "s", "lw": 2.0},
    "Random":       {"color": "#DC2626", "linestyle": ":",  "marker": "^", "lw": 1.8},
    "DegreeProduct":{"color": "#D97706", "linestyle": "-.", "marker": "D", "lw": 1.8},
}


# ---------------------------------------------------------------------------
# Experiment 1 — H vs k
# ---------------------------------------------------------------------------

def plot_exp1():
    exp_dir  = RESULTS_ROOT / "experiment1"
    out_dir  = PLOTS_ROOT / "experiment1"
    out_dir.mkdir(parents=True, exist_ok=True)

    files = sorted(exp_dir.glob("*.json"))
    if not files:
        print("[Exp1] No result files found — run exp1 first.")
        return

    # Group by dataset
    datasets: Dict[str, List] = {}
    for f in files:
        data = json.loads(f.read_text())
        ds = data["dataset"]
        datasets.setdefault(ds, []).append(data)

    for ds, configs in datasets.items():
        fig, axes = plt.subplots(2, 2, figsize=(14, 10))
        fig.suptitle(f"Experiment 1 — H vs k: {ds}", fontsize=14, fontweight="bold")

        for config in configs:
            p_mode    = config["p_mode"]
            eta_label = config["eta_label"]
            k_values  = config["k_values"]
            trials    = config["trials"]
            row = 0 if eta_label == "eta_low" else 1
            col = 0 if p_mode    == "fixed"   else 1
            ax  = axes[row][col]

            # Aggregate across trials
            for method_key, style in METHOD_STYLES.items():
                H_per_trial = []
                for trial in trials:
                    h_list = trial["methods"].get(method_key, {}).get("H_history", [])
                    if h_list:
                        H_per_trial.append(h_list)

                if not H_per_trial:
                    continue

                min_len = min(len(h) for h in H_per_trial)
                arr     = np.array([h[:min_len] for h in H_per_trial])
                means   = arr.mean(axis=0)
                stds    = arr.std(axis=0)
                xs      = k_values[:min_len]

                ax.plot(xs, means, label=method_key, **style)
                ax.fill_between(xs, means - stds, means + stds,
                                alpha=0.12, color=style["color"])

            p_label = "Fixed p" if p_mode == "fixed" else "Variable p"
            η_label = "η ∈ [0.1, 0.18]" if eta_label == "eta_low" else "η ∈ [0.18, 0.25]"
            ax.set_title(f"{p_label} | {η_label}", fontsize=10)
            ax.set_xlabel("Budget k (edges blocked)")
            ax.set_ylabel("H = I(G, Ω)")
            ax.legend(fontsize=8)
            ax.grid(True, alpha=0.3)

        plt.tight_layout()
        out_path = out_dir / f"{ds}_H_vs_k.png"
        plt.savefig(out_path, dpi=150, bbox_inches="tight")
        plt.close()
        print(f"  [Exp1] Saved → {out_path}")


# ---------------------------------------------------------------------------
# Experiment 2 — Efficiency
# ---------------------------------------------------------------------------

def plot_exp2():
    exp_dir = RESULTS_ROOT / "experiment2"
    out_dir = PLOTS_ROOT / "experiment2"
    out_dir.mkdir(parents=True, exist_ok=True)

    files = sorted(exp_dir.glob("*_efficiency.json"))
    if not files:
        print("[Exp2] No result files found."); return

    all_data = [json.loads(f.read_text()) for f in files]

    # --- Visualization 1: Filtering reduction (stacked bar) ---
    for data in all_data:
        ds  = data["dataset"]
        log = data.get("sg_snir_log", [])
        if not log:
            continue

        iters    = [r["iteration"]        for r in log]
        gamma_W  = [r["gamma_W_size"]     for r in log]
        gamma_p  = [r["gamma_prime_size"] for r in log]
        C_sizes  = [r["C_size"]           for r in log]

        fig, ax = plt.subplots(figsize=(10, 5))
        ax.bar(iters, gamma_W,  label="|Γ(W)|",    color="#93C5FD", alpha=0.9)
        ax.bar(iters, gamma_p,  label="|Γ'(W)|",   color="#3B82F6", alpha=0.9)
        ax.bar(iters, C_sizes,  label="|C|",        color="#1D4ED8", alpha=0.9)
        ax.set_title(f"{ds} — Candidate Set Reduction per Iteration")
        ax.set_xlabel("Iteration"); ax.set_ylabel("Candidate set size")
        ax.legend(); ax.grid(True, alpha=0.3, axis="y")
        plt.tight_layout()
        p = out_dir / f"{ds}_filtering_reduction.png"
        plt.savefig(p, dpi=150, bbox_inches="tight"); plt.close()
        print(f"  [Exp2] Saved → {p}")

    # --- Visualization 3: Speedup bar chart ---
    speedups, ds_labels = [], []
    for data in all_data:
        sg_evals  = sum(r["eval_count"] for r in data.get("sg_snir_log", []))
        meh_log   = data.get("max_expected_h_log", {})
        extrapolated = data.get("meh_extrapolated", False)
        if extrapolated:
            meh_evals = meh_log.get("gamma_W_size", 0)
        else:
            meh_evals = sum(meh_log.get("eval_counts", []))
        if sg_evals > 0:
            speedups.append(meh_evals / sg_evals)
            ds_labels.append(data["dataset"])

    if speedups:
        fig, ax = plt.subplots(figsize=(8, 5))
        bars = ax.bar(ds_labels, speedups, color="#2563EB", alpha=0.85)
        ax.bar_label(bars, fmt="%.1f×", fontsize=9)
        ax.set_title("Speedup: MaxExpH Evaluations / SG-SNIR Evaluations")
        ax.set_ylabel("Speedup ratio"); ax.grid(True, alpha=0.3, axis="y")
        plt.tight_layout()
        p = out_dir / "speedup_bar.png"
        plt.savefig(p, dpi=150, bbox_inches="tight"); plt.close()
        print(f"  [Exp2] Saved → {p}")


# ---------------------------------------------------------------------------
# Experiment 3 — Ablation
# ---------------------------------------------------------------------------

def plot_exp3():
    exp_dir = RESULTS_ROOT / "experiment3"
    out_dir = PLOTS_ROOT / "experiment3"
    out_dir.mkdir(parents=True, exist_ok=True)

    for f in sorted(exp_dir.glob("*_ablation.json")):
        data = json.loads(f.read_text())
        ds   = data["dataset"]
        variants = data["variants"]

        # H vs k
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 5))
        fig.suptitle(f"Experiment 3 — Ablation: {ds}", fontsize=13, fontweight="bold")
        colors = ["#2563EB", "#DC2626", "#16A34A", "#D97706"]

        for idx, (vkey, vdata) in enumerate(variants.items()):
            H = vdata.get("H_history", [])
            xs = list(range(1, len(H) + 1))
            if H:
                ax1.plot(xs, H, label=vdata["label"], color=colors[idx], linewidth=2,
                         linestyle=["-", "--", "-.", ":"][idx], marker=["o","s","^","D"][idx])

            evals = vdata.get("eval_counts", [])
            if evals:
                ax2.plot(list(range(1, len(evals)+1)), evals, label=vdata["label"],
                         color=colors[idx], linewidth=2,
                         linestyle=["-", "--", "-.", ":"][idx])

        ax1.set_title("H vs k"); ax1.set_xlabel("k"); ax1.set_ylabel("H")
        ax1.legend(fontsize=8); ax1.grid(True, alpha=0.3)
        ax2.set_title("Eval count vs k"); ax2.set_xlabel("k"); ax2.set_ylabel("Evaluations")
        ax2.legend(fontsize=8); ax2.grid(True, alpha=0.3)

        plt.tight_layout()
        p = out_dir / f"{ds}_ablation.png"
        plt.savefig(p, dpi=150, bbox_inches="tight"); plt.close()
        print(f"  [Exp3] Saved → {p}")


# ---------------------------------------------------------------------------
# Experiment 4 — ε sensitivity
# ---------------------------------------------------------------------------

def plot_exp4():
    exp_dir = RESULTS_ROOT / "experiment4"
    out_dir = PLOTS_ROOT / "experiment4"
    out_dir.mkdir(parents=True, exist_ok=True)

    for f in sorted(exp_dir.glob("*_epsilon_sensitivity.json")):
        data    = json.loads(f.read_text())
        ds      = data["dataset"]
        records = data["records"]

        eps_vals    = [r["eps"]         for r in records]
        H_finals    = [r["H_final"]     for r in records]
        eval_totals = [r["total_evals"] for r in records]

        fig, ax1 = plt.subplots(figsize=(8, 5))
        color_H, color_e = "#2563EB", "#DC2626"

        l1, = ax1.plot(eps_vals, H_finals, color=color_H, marker="o", lw=2, label="H_final")
        ax1.axvline(0.01, color="gray", linestyle="--", alpha=0.6, label="ε=0.01 (default)")
        ax1.set_xscale("log"); ax1.set_xlabel("ε (log scale)")
        ax1.set_ylabel("H_final", color=color_H); ax1.tick_params(axis="y", labelcolor=color_H)

        ax2 = ax1.twinx()
        l2, = ax2.plot(eps_vals, eval_totals, color=color_e, marker="s", lw=2,
                       linestyle="--", label="Total evals")
        ax2.set_ylabel("Total SNIR evaluations", color=color_e)
        ax2.tick_params(axis="y", labelcolor=color_e)

        lines = [l1, l2]
        ax1.legend(lines, [l.get_label() for l in lines], fontsize=9)
        ax1.set_title(f"{ds} — ε Sensitivity")
        ax1.grid(True, alpha=0.3)
        plt.tight_layout()
        p = out_dir / f"{ds}_epsilon.png"
        plt.savefig(p, dpi=150, bbox_inches="tight"); plt.close()
        print(f"  [Exp4] Saved → {p}")


# ---------------------------------------------------------------------------
# Experiment 5 — Scalability
# ---------------------------------------------------------------------------

def plot_exp5():
    exp_dir = RESULTS_ROOT / "experiment5"
    out_dir = PLOTS_ROOT / "experiment5"
    out_dir.mkdir(parents=True, exist_ok=True)

    summary_path = exp_dir / "scalability_summary.json"
    if not summary_path.exists():
        print("[Exp5] No summary file found."); return

    data = json.loads(summary_path.read_text())
    rows = data["rows"]

    sizes      = [r["size"]           for r in rows]
    sg_times   = [r["sg_snir_time_s"] for r in rows]
    meh_times  = [r.get("meh_time_s") for r in rows]  # None if OOT
    labels     = [r["dataset"]        for r in rows]

    fig, ax = plt.subplots(figsize=(10, 6))
    ax.plot(sizes, sg_times, "-o", color="#2563EB", lw=2.5, label="SG-SNIR")

    # MaxExpH — draw available points, mark OOT
    avail_x = [sizes[i] for i, t in enumerate(meh_times) if t is not None]
    avail_y = [t for t in meh_times if t is not None]
    oot_x   = [sizes[i] for i, t in enumerate(meh_times) if t is None]

    if avail_x:
        ax.plot(avail_x, avail_y, "--s", color="#16A34A", lw=2, label="MaxExpectedH")
    for x in oot_x:
        ax.axvline(x, color="#16A34A", linestyle=":", alpha=0.5)
        ax.text(x, ax.get_ylim()[1] * 0.9, "OOT", ha="center",
                color="#16A34A", fontsize=9, fontstyle="italic")

    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlabel("|V| + |E| (log scale)"); ax.set_ylabel("Total runtime (s, log scale)")
    ax.set_title("Experiment 5 — Scalability: Runtime vs Network Size")
    ax.legend(); ax.grid(True, alpha=0.3)

    for i, (x, y, label) in enumerate(zip(sizes, sg_times, labels)):
        ax.annotate(label, (x, y), textcoords="offset points",
                    xytext=(5, 5), fontsize=8, color="#2563EB")

    plt.tight_layout()
    p = out_dir / "scalability.png"
    plt.savefig(p, dpi=150, bbox_inches="tight"); plt.close()
    print(f"  [Exp5] Saved → {p}")


# ---------------------------------------------------------------------------
# Experiment 6 — KSCC Structure scatter
# ---------------------------------------------------------------------------

def plot_exp6():
    exp_dir = RESULTS_ROOT / "experiment6"
    out_dir = PLOTS_ROOT / "experiment6"
    out_dir.mkdir(parents=True, exist_ok=True)

    struct_path = exp_dir / "kscc_structure.json"
    if not struct_path.exists():
        print("[Exp6] No result file found."); return

    data    = json.loads(struct_path.read_text())["results"]
    fracs   = [r["kscc_fraction"]  for r in data]
    speedup = [r.get("speedup_evals") or 1.0 for r in data]
    labels  = [r["dataset"]        for r in data]

    fig, ax = plt.subplots(figsize=(7, 5))
    ax.scatter(fracs, speedup, s=120, color="#2563EB", zorder=5)
    for x, y, label in zip(fracs, speedup, labels):
        ax.annotate(label, (x, y), textcoords="offset points",
                    xytext=(6, 4), fontsize=9)

    # Trend line
    if len(fracs) > 1:
        z = np.polyfit(fracs, speedup, 1)
        p_line = np.poly1d(z)
        xs = np.linspace(min(fracs), max(fracs), 100)
        ax.plot(xs, p_line(xs), "--", color="#93C5FD", lw=1.5, label="Trend")

    ax.set_xlabel("KSCC Size Fraction (|S*| / |V|)")
    ax.set_ylabel("Eval Speedup (MEH evals / SG evals)")
    ax.set_title("Experiment 6 — KSCC Concentration vs Efficiency Gain")
    ax.legend(); ax.grid(True, alpha=0.3)
    plt.tight_layout()
    p = out_dir / "kscc_scatter.png"
    plt.savefig(p, dpi=150, bbox_inches="tight"); plt.close()
    print(f"  [Exp6] Saved → {p}")


# ---------------------------------------------------------------------------
# Experiment 7 — Edge category breakdown
# ---------------------------------------------------------------------------

def plot_exp7():
    exp_dir = RESULTS_ROOT / "experiment7"
    out_dir = PLOTS_ROOT / "experiment7"
    out_dir.mkdir(parents=True, exist_ok=True)

    for f in sorted(exp_dir.glob("*_edge_analysis.json")):
        data = json.loads(f.read_text())
        ds   = data["dataset"]

        categories = ["internal", "bridge", "peripheral"]
        colors     = ["#2563EB", "#16A34A", "#DC2626"]
        methods    = ["SG-SNIR", "MaxExpectedH"]

        sg_fracs  = data["SG-SNIR"]["category_fractions"]
        meh_fracs = data["MaxExpectedH"]["category_fractions"]

        fig, ax = plt.subplots(figsize=(8, 5))
        x = np.arange(len(methods))
        width = 0.25
        bottom_sg, bottom_meh = 0, 0

        for i, (cat, col) in enumerate(zip(categories, colors)):
            sg_val  = sg_fracs.get(cat, 0)
            meh_val = meh_fracs.get(cat, 0)
            ax.bar(x[0], sg_val,  width, bottom=bottom_sg,  color=col,
                   label=cat, alpha=0.85)
            ax.bar(x[1], meh_val, width, bottom=bottom_meh, color=col, alpha=0.55)
            bottom_sg  += sg_val
            bottom_meh += meh_val

        ax.set_xticks(x); ax.set_xticklabels(methods)
        ax.set_ylabel("Fraction of budget")
        ax.set_title(f"{ds} — Edge Category Breakdown")
        ax.legend(fontsize=9); ax.grid(True, alpha=0.3, axis="y")
        plt.tight_layout()
        p = out_dir / f"{ds}_edge_categories.png"
        plt.savefig(p, dpi=150, bbox_inches="tight"); plt.close()
        print(f"  [Exp7] Saved → {p}")


# ---------------------------------------------------------------------------
# Experiment 8 — Robustness box plots
# ---------------------------------------------------------------------------

def plot_exp8():
    exp_dir = RESULTS_ROOT / "experiment8"
    out_dir = PLOTS_ROOT / "experiment8"
    out_dir.mkdir(parents=True, exist_ok=True)

    for f in sorted(exp_dir.glob("*_robustness.json")):
        data = json.loads(f.read_text())
        ds   = data["dataset"]
        seed_sizes = data["seed_sizes"]

        sg_boxes  = [data["results"][str(s)]["sg_snir"]["H_finals"]      for s in seed_sizes]
        meh_boxes = [data["results"][str(s)]["max_expected_h"]["H_finals"] for s in seed_sizes]

        fig, axes = plt.subplots(1, 2, figsize=(13, 5), sharey=True)
        fig.suptitle(f"Experiment 8 — Robustness to Seed Selection: {ds}",
                     fontsize=13, fontweight="bold")

        for ax, boxes, method, color in [
            (axes[0], sg_boxes,  "SG-SNIR",      "#2563EB"),
            (axes[1], meh_boxes, "MaxExpectedH",  "#16A34A"),
        ]:
            valid = [b for b in boxes if b]
            if not valid:
                ax.text(0.5, 0.5, "No data", ha="center", va="center",
                        transform=ax.transAxes)
                continue
            bp = ax.boxplot(valid, labels=[str(s) for s in seed_sizes[:len(valid)]],
                            patch_artist=True,
                            boxprops=dict(facecolor=color, alpha=0.5),
                            medianprops=dict(color="black", lw=2))
            ax.set_title(method, fontsize=11, color=color)
            ax.set_xlabel("Initial infected seed count")
            ax.set_ylabel("H after k removals")
            ax.grid(True, alpha=0.3, axis="y")

        plt.tight_layout()
        p = out_dir / f"{ds}_robustness_boxplot.png"
        plt.savefig(p, dpi=150, bbox_inches="tight"); plt.close()
        print(f"  [Exp8] Saved → {p}")


# ---------------------------------------------------------------------------
# Main dispatcher
# ---------------------------------------------------------------------------

PLOTTERS = {
    1: plot_exp1,
    2: plot_exp2,
    3: plot_exp3,
    4: plot_exp4,
    5: plot_exp5,
    6: plot_exp6,
    7: plot_exp7,
    8: plot_exp8,
}


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--experiments", nargs="+", type=int,
                   default=list(PLOTTERS.keys()),
                   help="Which experiment numbers to plot (default: all)")
    args = p.parse_args()

    if not HAS_MPL:
        print("matplotlib is required for plotting. Install with: pip install matplotlib")
        return

    for exp_num in args.experiments:
        fn = PLOTTERS.get(exp_num)
        if fn:
            print(f"\n{'─'*60}")
            print(f"  Plotting Experiment {exp_num}...")
            fn()
        else:
            print(f"  Unknown experiment number: {exp_num}")

    print("\n✓ Done.")


if __name__ == "__main__":
    main()
