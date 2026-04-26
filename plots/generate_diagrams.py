import os
import matplotlib.pyplot as plt
import numpy as np

# Ensure directory exists
out_dir = '/Users/lalith/snu/sem6/sin/project/code/combo/plots/diagrams'
os.makedirs(out_dir, exist_ok=True)

# ---------------------------------------------------------
# Plot 1: Efficiency vs. Effectiveness Trade-off Scatter Plot
# ---------------------------------------------------------
# Data for p2p-Gnutella (k=20)
methods = ['MaxExpectedH\n(Paper 1)', 'SG-SNIR\n(Ours)', 'DegreeProduct\n(Heuristic)']
times = [92.4, 42.1, 1.5]  # Wall-clock time in seconds
h_final = [20.775, 20.797, 21.465]  # Expected Influence
colors = ['#e74c3c', '#2ecc71', '#95a5a6']

plt.figure(figsize=(7, 5))
plt.scatter(times, h_final, c=colors, s=200, edgecolor='black', zorder=3)

# Add labels to points with better placement
plt.annotate(methods[0], (times[0], h_final[0]), xytext=(-10, -35), 
             textcoords='offset points', ha='center', fontsize=10)
plt.annotate(methods[1], (times[1], h_final[1]), xytext=(0, -35), 
             textcoords='offset points', ha='center', fontweight='bold', fontsize=11)
plt.annotate(methods[2], (times[2], h_final[2]), xytext=(15, -5), 
             textcoords='offset points', fontsize=10)

# Formatting
# Remove log scale, use linear from 0 to 100
plt.xlim(0, 100) # Give more horizontal space
plt.ylim(20.65, 21.55) # Give more vertical space

plt.grid(True, linestyle='--', alpha=0.6, zorder=0)
plt.title('Quality vs. Efficiency Trade-off (p2p-Gnutella, k=20)', fontsize=12, pad=15)
plt.xlabel('$\\leftarrow$ Faster          Wall-clock Runtime (seconds)          Slower $\\rightarrow$', fontsize=11)
plt.ylabel('$\\leftarrow$ Better          Final Epidemic Influence ($H_{final}$)          Worse $\\rightarrow$', fontsize=11)

# Annotate the trade-off
plt.annotate('Optimal\nTrade-off', xy=(times[1], h_final[1] - 0.02), xytext=(10, 20.9),
             arrowprops=dict(facecolor='#27ae60', shrink=0.05, width=1.5, headwidth=6),
             fontsize=10, ha='center', color='#27ae60', fontweight='bold')

# Draw Pareto frontier line between SG-SNIR and MaxExpH
plt.plot([times[1], times[0]], [h_final[1], h_final[0]], 'k--', alpha=0.5, zorder=2)
plt.plot([times[2], times[1]], [h_final[2], h_final[1]], 'k--', alpha=0.5, zorder=2)

plt.tight_layout()
plt.savefig(os.path.join(out_dir, 'tradeoff_scatter.pdf'), format='pdf', dpi=300)
plt.close()

# ---------------------------------------------------------
# Plot 2: Equivalent Contact Budget Line Charts (300 & 600 seeds)
# ---------------------------------------------------------
import json

# Plot 300 Seeds (Grouped Bar Chart)
try:
    with open('../results/experiment9/logs/results_300_seeds.json', 'r') as f:
        exp9_300_data = json.load(f)['data']

    budgets_300 = [str(d['budget_edges']) + ' Edges\n(' + str(d['k_nodes']) + ' Nodes)' for d in exp9_300_data]
    dino_H_300 = [d['dino_H'] for d in exp9_300_data]
    sg_H_300 = [d['sg_H'] for d in exp9_300_data]

    x = np.arange(len(budgets_300))
    width = 0.35

    plt.figure(figsize=(8, 5))
    bars1 = plt.bar(x - width/2, dino_H_300, width, label='DINO (Node Removal)', color='#e74c3c', edgecolor='black', zorder=3)
    bars2 = plt.bar(x + width/2, sg_H_300, width, label='SG-SNIR (Edge Removal)', color='#3498db', edgecolor='black', zorder=3)

    plt.ylabel('Final Epidemic Influence ($H_{final}$) $\\leftarrow$ Better', fontsize=11)
    plt.title('Equivalent Contact Budget (300 Seeds - Massive Outbreak)', fontsize=12, pad=15)
    plt.xticks(x, budgets_300)
    plt.legend()
    plt.grid(axis='y', linestyle='--', alpha=0.6, zorder=0)

    # Add text labels
    for bar in bars1:
        plt.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 5, f'{bar.get_height():.1f}', ha='center', va='bottom', fontsize=9, fontweight='bold')
    for bar in bars2:
        plt.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 5, f'{bar.get_height():.1f}', ha='center', va='bottom', fontsize=9, fontweight='bold')

    # Keep y-axis scalable to show gap clearly, but ensure 0 isn't required if differences are small.
    # We will let matplotlib auto-scale the y-axis, but add some padding at top
    plt.ylim(min(sg_H_300)*0.9, max(dino_H_300)*1.05)

    plt.tight_layout()
    plt.savefig(os.path.join(out_dir, 'budget_comparison_300_seeds.pdf'), format='pdf', dpi=300)
    plt.close()
except FileNotFoundError:
    print("300 seeds data not found.")

# Plot 600 Seeds (Grouped Bar Chart)
try:
    with open('../results/experiment9/logs/results_600_seeds.json', 'r') as f:
        exp9_600_data = json.load(f)['data']

    budgets_600 = [str(d['budget_edges']) + ' Edges\n(' + str(d['k_nodes']) + ' Nodes)' for d in exp9_600_data]
    dino_H_600 = [d['dino_H'] for d in exp9_600_data]
    sg_H_600 = [d['sg_H'] for d in exp9_600_data]

    x = np.arange(len(budgets_600))
    width = 0.35

    plt.figure(figsize=(8, 5))
    bars1 = plt.bar(x - width/2, dino_H_600, width, label='DINO (Node Removal)', color='#e74c3c', edgecolor='black', zorder=3)
    bars2 = plt.bar(x + width/2, sg_H_600, width, label='SG-SNIR (Edge Removal)', color='#3498db', edgecolor='black', zorder=3)

    plt.ylabel('Final Epidemic Influence ($H_{final}$) $\\leftarrow$ Better', fontsize=11)
    plt.title('Equivalent Contact Budget (600 Seeds - Extreme Outbreak)', fontsize=12, pad=15)
    plt.xticks(x, budgets_600)
    plt.legend()
    plt.grid(axis='y', linestyle='--', alpha=0.6, zorder=0)

    # Add text labels
    for bar in bars1:
        plt.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 5, f'{bar.get_height():.1f}', ha='center', va='bottom', fontsize=9, fontweight='bold')
    for bar in bars2:
        plt.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 5, f'{bar.get_height():.1f}', ha='center', va='bottom', fontsize=9, fontweight='bold')

    plt.ylim(min(sg_H_600)*0.9, max(dino_H_600)*1.05)

    plt.tight_layout()
    plt.savefig(os.path.join(out_dir, 'budget_comparison_600_seeds.pdf'), format='pdf', dpi=300)
    plt.close()
except FileNotFoundError:
    print("600 seeds data not found.")

print("Successfully generated plots in plots/diagrams/")
