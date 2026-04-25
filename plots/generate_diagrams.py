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
# Plot 2: Equivalent Contact Budget Line Chart
# ---------------------------------------------------------
import json

with open('../results/experiment9/p2p-Gnutella_baseline_comparison.json', 'r') as f:
    exp9_data = json.load(f)['data']

budgets = [d['budget_edges'] for d in exp9_data]
dino_H = [d['dino_H'] for d in exp9_data]
sg_H = [d['sg_H'] for d in exp9_data]

plt.figure(figsize=(7, 5))

# Plot lines
plt.plot(budgets, dino_H, marker='s', markersize=8, linewidth=2, color='#e74c3c', label='DINO (Paper 2) [Node Removal]')
plt.plot(budgets, sg_H, marker='o', markersize=8, linewidth=2, color='#3498db', label='SG-SNIR (Ours) [Edge Removal]')

# Formatting
plt.grid(True, linestyle='--', alpha=0.6)
plt.title('Equivalent Contact Budget (p2p-Gnutella)', fontsize=12, pad=15)
plt.xlabel('Equivalent Contact Budget (Edges Removed)', fontsize=11)
plt.ylabel('Final Epidemic Influence ($H_{final}$) $\\leftarrow$ Better', fontsize=11)

# Annotate the gap
plt.annotate('21.1% Reduction\nin Epidemic Spread!', 
             xy=(budgets[-1], sg_H[-1]), xytext=(budgets[-1]-100, 17.5),
             arrowprops=dict(facecolor='#2c3e50', shrink=0.05, width=2, headwidth=8),
             fontsize=11, ha='center', fontweight='bold', color='#2c3e50',
             bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="black", lw=1))

plt.ylim(15, 20.5)
plt.legend(loc='upper right')
plt.tight_layout()
plt.savefig(os.path.join(out_dir, 'budget_comparison_line.pdf'), format='pdf', dpi=300)
plt.close()

print("Successfully generated plots in plots/diagrams/")
