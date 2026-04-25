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

# Add labels to points
for i, method in enumerate(methods):
    if method == 'SG-SNIR\n(Ours)':
        plt.annotate(method, (times[i], h_final[i]), xytext=(10, -20), 
                     textcoords='offset points', fontweight='bold')
    else:
        plt.annotate(method, (times[i], h_final[i]), xytext=(10, 5), 
                     textcoords='offset points')

# Formatting
plt.xscale('log')
plt.grid(True, linestyle='--', alpha=0.6, zorder=0)
plt.title('Quality vs. Efficiency Trade-off (p2p-Gnutella, k=20)', fontsize=12, pad=15)
plt.xlabel('Wall-clock Runtime (seconds) [Log Scale] $\\rightarrow$ Faster', fontsize=11)
plt.ylabel('Final Epidemic Influence ($H_{final}$) $\\rightarrow$ Better', fontsize=11)
# Invert x-axis so faster (smaller time) is to the right? No, standard is left-to-right.
# Actually, lower H is better (bottom), lower time is better (left).
# So bottom-left is the optimal corner.
plt.annotate('Optimal Region\n(Fast & High Quality)', xy=(5, 20.75), xytext=(2, 20.9),
             arrowprops=dict(facecolor='black', shrink=0.05, width=1.5, headwidth=6),
             fontsize=10, ha='center', color='#27ae60')

plt.tight_layout()
plt.savefig(os.path.join(out_dir, 'tradeoff_scatter.pdf'), format='pdf', dpi=300)
plt.close()

# ---------------------------------------------------------
# Plot 2: Equivalent Contact Budget Bar Chart
# ---------------------------------------------------------
# Data for Exp 9: 459 edges removed
methods_bar = ['DINO (Paper 2)\n[Node Removal]', 'SG-SNIR (Ours)\n[Edge Removal]']
h_bar = [19.80, 15.63]
colors_bar = ['#e74c3c', '#3498db']

plt.figure(figsize=(6, 5))
bars = plt.bar(methods_bar, h_bar, color=colors_bar, width=0.5, edgecolor='black', zorder=3)

# Add value labels on top of bars
for bar in bars:
    yval = bar.get_height()
    plt.text(bar.get_x() + bar.get_width()/2, yval + 0.2, f'{yval:.2f}', 
             ha='center', va='bottom', fontweight='bold', fontsize=11)

plt.grid(axis='y', linestyle='--', alpha=0.6, zorder=0)
plt.title('Equivalent Contact Budget (p2p-Gnutella)', fontsize=12, pad=15)
plt.ylabel('Final Epidemic Influence ($H_{final}$) $\\rightarrow$ Lower is Better', fontsize=11)

# Annotate the improvement
plt.annotate('21.1% Reduction\nin Epidemic Spread!', 
             xy=(1, 15.63), xytext=(0.5, 12),
             arrowprops=dict(facecolor='#2c3e50', shrink=0.05, width=2, headwidth=8),
             fontsize=11, ha='center', fontweight='bold', color='#2c3e50',
             bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="black", lw=1))

plt.ylim(0, 22)
plt.tight_layout()
plt.savefig(os.path.join(out_dir, 'budget_comparison_bar.pdf'), format='pdf', dpi=300)
plt.close()

print("Successfully generated plots in plots/diagrams/")
