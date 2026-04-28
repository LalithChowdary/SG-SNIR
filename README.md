# SG-SNIR

**Spectral-Guided Edge Blocking for Epidemic Containment in Directed Networks Using the SNIR Model**


---


## What is SG-SNIR?

SG-SNIR is a hybrid epidemic containment algorithm for directed networks. Given a budget of `k` edge removals, it finds the best edges to block in order to minimise total epidemic spread — measured by the SNIR epidemic model.

The core idea: instead of evaluating every possible edge for blocking (expensive), SG-SNIR uses the network's **spectral structure** to prune the candidate set first, then runs full epidemic simulations only on the most promising edges.

---

## Project Structure

```
combo/
├── src/                        # Core algorithm and model source code
│   ├── sg_snir.py              # Main SG-SNIR algorithm
│   ├── snir_model.py           # SNIR epidemic model
│   ├── baselines.py            # Baseline algorithms (MaxExpH, DegProd, Random)
│   ├── data_loader.py          # Dataset loading utilities
│   └── tests.py                # Unit tests
│
├── experiments/                # One script per experiment (exp1–exp9)
│   ├── exp1_effectiveness.py
│   ├── exp2_efficiency.py
│   ├── exp3_ablation.py
│   ├── exp4_epsilon.py
│   ├── exp5_scalability.py
│   ├── exp6_kscc_analysis.py
│   ├── exp7_edge_analysis.py
│   ├── exp8_robustness.py
│   └── exp9_baseline_comparison.py
│
├── results/                    # JSON result files + run logs, per experiment
│   ├── experiment1/
│   │   ├── logs/               # Terminal logs from experiment runs
│   │   └── *.json              # Result data files
│   ├── experiment2/ ... experiment9/
│   └── (same structure for each)
│
├── plots/                      # Plotting scripts and generated figures
│   ├── plot_results.py
│   ├── generate_diagrams.py
│   └── diagrams/               # Generated PDF/image outputs
│
├── data/                       # Raw network datasets
│   ├── snap/                   # SNAP datasets (p2p-Gnutella, Wiki-Vote, etc.)
│   └── ICPSR_22140/            # HIV Transmission Network dataset
│
├── docs/                       # Documentation
│   ├── EXPERIMENTS.md          # ← Detailed experiment guide (how to run, results, logs)
│   └── Spectral-Guided Edge Blocking for Epidemic Containment in Directed Networks Using the SNIR Model.pdf  # Paper PDF
│
├── latex/                      # LaTeX source for the paper
│   ├── conference_101719.tex   # Main paper source
│   └── conference_101719.pdf   # Compiled PDF
│
└── main.py                     # Standalone runner (quick single-dataset test)
```

---

## Installation

**Requirements:** Python 3.9+

```bash
pip install networkx numpy scipy matplotlib
```

No additional setup is needed. Datasets are downloaded automatically on first run via `data_loader.py`.

---

## Quick Start

Run SG-SNIR on the HIV network with a budget of 20 edge removals:

```bash
python main.py --dataset hiv --k 20
```

Run on p2p-Gnutella and compare against all baselines:

```bash
python main.py --dataset gnutella --k 20
```

Run an epsilon sensitivity sweep:

```bash
python main.py --dataset hiv --k 15 --eps-sweep
```

See `main.py` header for the full list of configurable flags (`--alpha`, `--beta`, `--seed`, `--no-maxh`, etc.).

---

## Running Experiments

Each of the 9 experiments has a dedicated script in the `experiments/` folder.  

📄 **For full details on what each experiment tests, how to run it, expected output, and where to find logs, see:**

👉 **[docs/EXPERIMENTS.md](docs/EXPERIMENTS.md)**

### Quick Reference

| Experiment | Script | What it proves |
|-----------|--------|---------------|
| Exp 1 | `exp1_effectiveness.py` | SG-SNIR beats DegProd and Random baselines |
| Exp 2 | `exp2_efficiency.py` | SG-SNIR uses 2.3–2.4× fewer SNIR evaluations |
| Exp 3 | `exp3_ablation.py` | Every component (KSCC, bridges, SpectralDrop) contributes |
| Exp 4 | `exp4_epsilon.py` | Results are robust across ε ∈ [0.001, 0.5] |
| Exp 5 | `exp5_scalability.py` | Wall-clock speedup on large networks |
| Exp 6 | `exp6_kscc_analysis.py` | KSCC structure analysis across datasets |
| Exp 7 | `exp7_edge_analysis.py` | Which edges are selected and why |
| Exp 8 | `exp8_robustness.py` | Stability across different random seeds |
| Exp 9 | `exp9_baseline_comparison.py` | Comparison against DINO node-immunisation |

---

## Key Results Summary

| Dataset | KSCC% | SG-SNIR vs MaxExpH (quality gap) | Speedup vs MaxExpH | SG-SNIR vs DINO (same disruption budget) |
|---------|:-----:|:---------------------------------:|:------------------:|:-----------------------------------------:|
| HIV | 0.18% | 0.000% (exact match) | 1.00× (graceful) | — |
| p2p-Gnutella | 32.8% | +0.106% | **2.41×** | **17% fewer infections** (300-seed outbreak) |
| Wiki-Vote | 18.3% | +2.256% (3-trial upper bound) | **2.29×** | — |
| soc-Epinions1 | 42.5% | — | 1.31× | — |

---

## Paper

The full paper PDF is available at:
- `docs/Spectral-Guided Edge Blocking for Epidemic Containment in Directed Networks Using the SNIR Model.pdf`
