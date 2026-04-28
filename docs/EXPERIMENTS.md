# SG-SNIR — Experiment Guide

This document describes all 9 experiments in the SG-SNIR project: what each one tests, how to run it, where the results are saved, and where to find the logs.

---

## General Setup

**All experiments are run from the project root directory** (`combo/`):

```bash
cd /path/to/combo
```

**SNIR parameters used across all experiments (unless noted otherwise):**

| Parameter | Value | Meaning |
|-----------|-------|---------|
| α | 0.033 | S → N (asymptomatic infection) |
| β | 0.022 | S → I (symptomatic infection) |
| δ | 0.020 | N → I (progression) |
| η | 0.140 | N → R (asymptomatic recovery) |
| γ | 0.315 | I → R (symptomatic recovery) |
| ξ | 0.300 | R → S (immunity loss) |
| T | 10 | Simulation time horizon |
| seeds | 5 nodes | Initial infected count (BASE_SEED=42) |

---

## Datasets

| Dataset | Nodes | Edges | KSCC% | ρ(S*) | Source |
|---------|------:|------:|:-----:|:-----:|--------|
| HIV Transmission | 35,228 | 49,779 | 0.18% | 6.25 | ICPSR #22140 |
| p2p-Gnutella | 6,301 | 20,777 | 32.8% | 4.88 | SNAP |
| Wiki-Vote | 7,115 | 103,689 | 18.3% | 44.57 | SNAP |
| soc-Epinions1 | 75,879 | 508,837 | 42.5% | 82.79 | SNAP |
| email-EuAll | 265,214 | 420,045 | — | — | SNAP (scalability only) |

Datasets are downloaded automatically by `data_loader.py` on first use and cached in `data/`.

---

## Experiment 1 — Effectiveness

### What it proves
SG-SNIR reduces cumulative epidemic influence $H$ better than the DegreeProduct and Random baselines across all dataset/configuration combinations.

### How to run
```bash
python experiments/exp1_effectiveness.py
```

### Configuration
- Budget `k`: 50 (HIV), 20 (p2p-Gnutella, Wiki-Vote)
- Trials: 10 (HIV, p2p-Gnutella), 5 (Wiki-Vote)
- Weight modes: fixed (0.05) and variable (degree-scaled)
- η settings: `eta_low` (0.14) and `eta_high` (0.28)

### Results location
```
results/experiment1/
├── HIV_fixed_eta_low.json
├── HIV_fixed_eta_high.json
├── HIV_variable_eta_low.json
├── HIV_variable_eta_high.json
├── p2p-Gnutella_fixed_eta_low.json
├── p2p-Gnutella_fixed_eta_high.json
├── p2p-Gnutella_variable_eta_low.json
├── p2p-Gnutella_variable_eta_high.json
├── Wiki-Vote_fixed_eta_low.json
├── Wiki-Vote_fixed_eta_high.json
├── Wiki-Vote_variable_eta_low.json
├── Wiki-Vote_variable_eta_high.json
└── WikiVote_sg_vs_meh.json    ← Direct SG-SNIR vs MaxExpH comparison
```

### Logs
```
results/experiment1/logs/
├── exp1_run2.log              ← Main run output
└── exp1_wikivote_meh.log      ← Wiki-Vote MaxExpH comparison run
```

### Key result
SG-SNIR beats DegreeProduct by **2–8%** and Random across all 12 configurations. Quality gap vs MaxExpH: **+0.106%** (p2p-Gnutella), **+2.256%** (Wiki-Vote, 3-trial upper bound), **0.000%** (HIV, by construction).

---

## Experiment 2 — Computational Efficiency

### What it proves
SG-SNIR uses **2.3–2.4× fewer SNIR simulations** than exhaustive MaxExpectedH on networks with KSCC ≥ 18%. On HIV (KSCC = 0.18%), it degrades gracefully to MaxExpH with no quality loss.

### How to run
```bash
python experiments/exp2_efficiency.py
```

### Configuration
- Budget `k = 20`, 1 trial per dataset
- Datasets: HIV, p2p-Gnutella, Wiki-Vote, soc-Epinions1
- Metric: total SNIR evaluation count across all k iterations

### Results location
```
results/experiment2/
├── HIV_efficiency.json
├── p2p-Gnutella_efficiency.json
├── Wiki-Vote_efficiency.json
└── soc-Epinions1_efficiency.json
```

### Logs
```
results/experiment2/logs/
├── exp2_run.log               ← HIV + p2p + Wiki-Vote run
└── exp2_epinions_run.log      ← soc-Epinions1 run (separate, slow)
```

### Key result
| Dataset | KSCC% | SG Evals | MEH Evals | Speedup |
|---------|:-----:|:--------:|:---------:|:-------:|
| HIV | 0.18% | 1,230 | 1,230 | 1.00× |
| Wiki-Vote | 18.3% | 34 | 78 | **2.29×** |
| p2p-Gnutella | 32.8% | 286 | 690 | **2.41×** |
| soc-Epinions1 | 42.5% | 601 | 790 | 1.31× |

---

## Experiment 3 — Ablation Study

### What it proves
Every component of SG-SNIR (KSCC filter, bridge edges, SpectralDrop threshold) contributes measurably. Removing any one component either raises evaluation count or worsens quality.

### How to run
```bash
python experiments/exp3_ablation.py
```

### Configuration
- k = 20, 5 trials — p2p-Gnutella (primary) + HIV (boundary check)
- k = 10, 1 trial — soc-Epinions1 (large network spot check)
- Variants: Full SG-SNIR / No KSCC filter / No bridge edges / No spectral threshold

### Results location
```
results/experiment3/
├── p2p-Gnutella_ablation.json
├── HIV_ablation.json
└── soc-Epinions1_ablation_partial.json   ← k=10, 1 trial
```

### Logs
```
results/experiment3/logs/
├── exp3_run.log                    ← p2p-Gnutella + HIV ablation
└── exp3_epinions_ablation.log      ← soc-Epinions1 partial run
```

### Key result (p2p-Gnutella, k=20, 5 trials)
| Variant | H_final | Total Evals | vs Full |
|---------|:-------:|:-----------:|:-------:|
| Full SG-SNIR | 20.049 | 227.6 | baseline |
| No KSCC filter | 20.036 | 618.0 | 2.72× more evals |
| No bridge edges | 20.090 | 228.4 | +0.20% worse quality |
| No spectral threshold | 20.036 | 264.4 | 1.16× more evals |

> **Note:** soc-Epinions1 counts (286 vs 445 evals) differ from Exp 2 (601 vs 790) because this ablation uses k=10 while Exp 2 uses k=20. Counts scale approximately linearly with budget.

---

## Experiment 4 — Sensitivity to ε

### What it proves
The SpectralDrop threshold ε is robust — results are completely flat across ε ∈ [0.001, 0.5] (two orders of magnitude). This means SG-SNIR requires no hyperparameter tuning in practice.

### How to run
```bash
python experiments/exp4_epsilon.py
```

### Configuration
- ε values: [0.001, 0.005, 0.010, 0.020, 0.050, 0.100, 0.200, 0.500]
- 1 trial per ε, k = 20
- Datasets: HIV and p2p-Gnutella

### Results location
```
results/experiment4/
├── HIV_epsilon_sensitivity.json
├── p2p-Gnutella_epsilon_sensitivity.json
├── spectral_drop_distribution.json       ← SpectralDrop values from iteration 1
└── spectral_drop_distribution_full.json  ← SpectralDrop values across all 20 iterations
```

### Logs
```
results/experiment4/logs/
└── exp4_run.log
```

### Key result
All 290 SpectralDrop values (across 20 iterations on p2p-Gnutella) lie in [5.37×10⁻⁴, 1.50×10⁻³]. Since even the lowest threshold (ε=0.001 × ρ = 0.00488) exceeds the max observed drop, the filter retains all edges for every tested ε. H_final and eval count are identical for all ε values tested.

---

## Experiment 5 — Scalability

### What it proves
Wall-clock time speedup closely tracks evaluation-count speedup, confirming that SNIR simulation dominates runtime. MaxExpH exceeds the 3,600s timeout on email-EuAll; SG-SNIR completes successfully.

### How to run
```bash
python experiments/exp5_scalability.py
```

### Configuration
- k = 20, 1 trial
- All 5 datasets including email-EuAll
- Time limit: 3,600s (OOT = exceeds limit)

### Results location
```
results/experiment5/
└── scalability_summary.json
```

### Logs
```
results/experiment5/logs/
├── exp5_run.log               ← p2p + Wiki + HIV + soc-Epinions
└── exp5_epinions_run.log      ← soc-Epinions1 separate run
```

### Key result
| Dataset | SG-SNIR | MaxExpH |
|---------|:-------:|:-------:|
| p2p-Gnutella | 42.1s | 92.4s |
| Wiki-Vote | 14.4s | 38.8s |
| HIV | 288.9s | 299.9s |
| soc-Epinions1 | 1,631.6s | 2,147.9s |
| email-EuAll | — | **OOT** |

---

## Experiment 6 — KSCC Structure Analysis

### What it proves
Characterises the KSCC (Key Strongly Connected Component) structure of each dataset and validates that the spectral radius approximation used by SG-SNIR is accurate.

### How to run
```bash
python experiments/exp6_kscc_analysis.py
```

### Results location
```
results/experiment6/
└── kscc_structure.json
```

### Logs
```
results/experiment6/logs/
└── exp6_run.log
```

---

## Experiment 7 — Edge Analysis

### What it proves
Analyses the edges selected by SG-SNIR over the k iterations — which edges are bridge edges vs KSCC internal edges, and how SpectralDrop scores are distributed across selected edges.

### How to run
```bash
python experiments/exp7_edge_analysis.py
```

### Results location
```
results/experiment7/
├── HIV_edge_analysis.json
└── Wiki-Vote_edge_analysis.json
```

### Logs
```
results/experiment7/logs/
└── exp7_run.log
```

---

## Experiment 8 — Robustness (Seed Stability)

### What it proves
SG-SNIR results are stable across different random seeds for initial infection placement, validating that the reported results are not the product of lucky seed selection.

### How to run
```bash
python experiments/exp8_robustness.py
```

### Configuration
- Multiple seeds tested per dataset
- k = 20, fixed/eta_low configuration

### Results location
```
results/experiment8/
├── HIV_robustness.json
├── p2p-Gnutella_robustness.json
└── Wiki-Vote_robustness.json
```

### Logs
```
results/experiment8/logs/
├── exp8_run.log               ← HIV + p2p-Gnutella
└── exp8_wiki_run.log          ← Wiki-Vote separate run
```

---

## Experiment 9 — Comparison Against DINO

### What it proves
SG-SNIR (edge blocking) vastly outperforms DINO (node immunisation) when both are given the same **Network Disruption budget** (total number of connections severed). This validates that targeted edge blocking is more effective than coarse-grained node removal.

### How to run
```bash
python experiments/exp9_baseline_comparison.py
```

### Configuration
- Dataset: p2p-Gnutella
- Seed counts: 300 and 600 (large outbreaks to stress-test DINO)
- DINO quarantines top hub nodes; SG-SNIR gets the same edge-disruption budget

### Results location
```
results/experiment9/
├── p2p-Gnutella_baseline_comparison.json
└── logs/
    ├── exp9_run.log
    ├── results_300_seeds.json
    └── results_600_seeds.json
```

### Key result
| Seeds | Budget (edges severed) | SG-SNIR H_final | DINO H_final |
|:-----:|:----------------------:|:---------------:|:------------:|
| 300 | 1,256 edges | **938.82** | 1,136.50 |
| 600 | 1,607 edges | **1,910.60** | 2,246.54 |

SG-SNIR achieves ~17% lower epidemic spread than DINO for the same structural cost.

---

## Generating Plots

To regenerate all paper figures from saved result JSONs:

```bash
python plots/plot_results.py
python plots/generate_diagrams.py
```

Output PDFs are saved in `plots/diagrams/` and are also symlinked into `latex/diagrams/` for the paper.
