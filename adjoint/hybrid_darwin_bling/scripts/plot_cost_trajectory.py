#!/usr/bin/env python3
"""Cost function vs. M1QN3 iteration for the 1-month hybrid Darwin
(forward)/BLING(adjoint) surrogate-gradient pilot: both BLING's own cost
(what the adjoint gradient was actually computed from) and Darwin's own
cost (the true objective the surrogate gradient is being used to reduce).
Values also available in ../results/cost_trajectory_results.csv."""
import os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

OUT_DIR = os.environ.get("FIGURES_OUT_DIR", "../figures")

cycles = [0, 1, 2, 3, 4, 5]
bling_fc = [35041.3, 28365.3, 15988.2, 12451.4, 7874.6, 3725.4]
darwin_J = [40908.9, 34686.1, 22281.0, 18438.7, 13096.4, 6321.4]

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

ax1.plot(cycles, bling_fc, 'o-', color='#D55E00', label='BLING (adjoint source)', linewidth=2, markersize=7)
ax1.plot(cycles, darwin_J, 's-', color='#0072B2', label='Darwin (true objective)', linewidth=2, markersize=7)
ax1.set_xlabel('M1QN3 iteration')
ax1.set_ylabel('SOCAT month-01 pCO2 cost function')
ax1.set_title('Cost trajectory (linear)')
ax1.legend(frameon=False)
ax1.spines[['top', 'right']].set_visible(False)
ax1.set_xticks(cycles)
darwin_pct = 100 * (darwin_J[-1] - darwin_J[0]) / darwin_J[0]
ax1.annotate(f'{darwin_pct:+.1f}%', (cycles[-1], darwin_J[-1]),
             textcoords='offset points', xytext=(10, -10), fontsize=10, color='#0072B2')

ax2.semilogy(cycles, bling_fc, 'o-', color='#D55E00', label='BLING (adjoint source)', linewidth=2, markersize=7)
ax2.semilogy(cycles, darwin_J, 's-', color='#0072B2', label='Darwin (true objective)', linewidth=2, markersize=7)
ax2.set_xlabel('M1QN3 iteration')
ax2.set_ylabel('SOCAT month-01 pCO2 cost function (log scale)')
ax2.set_title('Cost trajectory (log)')
ax2.legend(frameon=False)
ax2.spines[['top', 'right']].set_visible(False)
ax2.set_xticks(cycles)

fig.suptitle('Hybrid Darwin(forward)/BLING(adjoint) surrogate-gradient pilot\n'
             'DIC+ALK only, 1-month SOCAT cost, numiter=5 (stopped via M1QN3 iteration cap, not error)',
             fontsize=11)
fig.tight_layout(rect=[0, 0, 1, 0.92])

os.makedirs(OUT_DIR, exist_ok=True)
OUT = f'{OUT_DIR}/cost_trajectory.png'
fig.savefig(OUT, dpi=150)
print(f"saved to {OUT}")
print(f"Darwin J: {darwin_J[0]:.1f} -> {darwin_J[-1]:.1f}  ({darwin_pct:+.1f}%)")
