#!/usr/bin/env python3
"""Cost function vs. M1QN3 iteration for the full-year hybrid Darwin
(forward)/BLING(adjoint) surrogate-gradient run (12-month SOCAT cost,
DIC+ALK only, numiter=1). Values also available in
../results/cost_trajectory_results.csv."""
import os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

OUT_DIR = os.environ.get("FIGURES_OUT_DIR", "../figures")

cycles = [0, 1]
bling_fc = [212372.7, 164280.1]
darwin_J = [340152.9, 304597.1]

fig, ax = plt.subplots(figsize=(6, 5))
ax.plot(cycles, bling_fc, 'o-', color='#D55E00', label='BLING (adjoint source)', linewidth=2, markersize=9)
ax.plot(cycles, darwin_J, 's-', color='#0072B2', label='Darwin (true objective)', linewidth=2, markersize=9)
ax.set_xlabel('M1QN3 iteration')
ax.set_ylabel('SOCAT full-year (12-month) pCO2 cost function')
ax.set_title('Full-year hybrid pilot: cost trajectory')
ax.legend(frameon=False)
ax.spines[['top', 'right']].set_visible(False)
ax.set_xticks(cycles)
darwin_pct = 100 * (darwin_J[-1] - darwin_J[0]) / darwin_J[0]
ax.annotate(f'{darwin_pct:+.1f}%', (cycles[-1], darwin_J[-1]),
            textcoords='offset points', xytext=(-45, 8), fontsize=10, color='#0072B2')
fig.tight_layout()

os.makedirs(OUT_DIR, exist_ok=True)
OUT = f'{OUT_DIR}/cost_trajectory_1yr.png'
fig.savefig(OUT, dpi=150)
print(f"saved to {OUT}")
print(f"Darwin J: {darwin_J[0]:.1f} -> {darwin_J[-1]:.1f}  ({darwin_pct:+.1f}%)")
