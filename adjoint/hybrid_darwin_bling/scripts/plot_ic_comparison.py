#!/usr/bin/env python3
"""Initial vs. optimized DIC/ALK surface maps from the 1-month hybrid
pilot: prior (cycle 0, unperturbed) vs. final (cycle 5, last accepted
M1QN3 guess) Darwin IC, plus the difference map.

Configure via environment variables:
  HYBRID_DIR    -- hybrid_darwin_bling/ root (default: parent of scripts/)
  BLING_PRIOR_DIR -- global_oce_biogeo_bling_SOCAT/input_ad/ (prior IC)
  FIGURES_OUT_DIR -- where to write the output PNG
  FINAL_CYCLE   -- last cycle's IC to compare against (default 5, matches
                    the 1-month pilot's numiter=5 -> 6 cycles, 0..5)
"""
import os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

NZ, NY, NX = 15, 64, 128
HYBRID = os.environ.get("HYBRID_DIR", "..")
BLING_PRIOR = os.environ.get("BLING_PRIOR_DIR", "../../global_oce_biogeo_bling_SOCAT/input_ad")
OUT_DIR = os.environ.get("FIGURES_OUT_DIR", "../figures")
FINAL_CYCLE = os.environ.get("FINAL_CYCLE", "5")

FIELDS = {
    'DIC': dict(prior=f'{BLING_PRIOR}/dic_init.bin', units='mol m$^{-3}$'),
    'ALK': dict(prior=f'{BLING_PRIOR}/alk_init.bin', units='eq m$^{-3}$'),
}
FINAL = {
    'DIC': f'{HYBRID}/run_darwin/dic_darwin_init.bin',
    'ALK': f'{HYBRID}/run_darwin/alk_darwin_init.bin',
}


def load(path):
    return np.fromfile(path, dtype='>f4').reshape(NZ, NY, NX)[0]


fig, axes = plt.subplots(2, 3, figsize=(15, 7.5))

for row, (name, cfg) in enumerate(FIELDS.items()):
    prior = load(cfg['prior'])  # BLING units (mol/m3 or eq/m3)
    final_darwin_units = load(FINAL[name])  # mmol/m3
    final = final_darwin_units / 1000.0  # -> BLING units for apples-to-apples

    prior_m = np.ma.masked_where(prior == 0, prior)
    final_m = np.ma.masked_where(final == 0, final)
    delta_m = np.ma.masked_where(prior == 0, final - prior)

    vmin, vmax = prior_m.min(), prior_m.max()
    im0 = axes[row, 0].imshow(prior_m, origin='lower', extent=[0, 360, -90, 90],
                               vmin=vmin, vmax=vmax, cmap='viridis', aspect='auto')
    axes[row, 0].set_title(f'{name}: prior (cycle 0)', fontsize=10)
    fig.colorbar(im0, ax=axes[row, 0], label=cfg['units'], shrink=0.85)

    im1 = axes[row, 1].imshow(final_m, origin='lower', extent=[0, 360, -90, 90],
                               vmin=vmin, vmax=vmax, cmap='viridis', aspect='auto')
    axes[row, 1].set_title(f'{name}: optimized (cycle {FINAL_CYCLE})', fontsize=10)
    fig.colorbar(im1, ax=axes[row, 1], label=cfg['units'], shrink=0.85)

    dmax = np.abs(delta_m).max()
    im2 = axes[row, 2].imshow(delta_m, origin='lower', extent=[0, 360, -90, 90],
                               vmin=-dmax, vmax=dmax, cmap='RdBu_r', aspect='auto')
    axes[row, 2].set_title(f'{name}: optimized - prior', fontsize=10)
    fig.colorbar(im2, ax=axes[row, 2], label=cfg['units'], shrink=0.85)

    for ax in axes[row]:
        ax.set_ylabel('latitude')
    if row == 1:
        for ax in axes[row]:
            ax.set_xlabel('longitude')

fig.suptitle(f'Hybrid Darwin/BLING pilot: surface IC before vs. after {FINAL_CYCLE} M1QN3 iterations\n'
             '(BLING-adjoint-driven optimization of Darwin\'s own DIC/ALK ICs against SOCAT)',
             fontsize=12)
fig.tight_layout(rect=[0, 0, 1, 0.93])

os.makedirs(OUT_DIR, exist_ok=True)
OUT = f'{OUT_DIR}/ic_comparison.png'
fig.savefig(OUT, dpi=140)
print(f"saved to {OUT}")

for name in FIELDS:
    prior = load(FIELDS[name]['prior'])
    final = load(FINAL[name]) / 1000.0
    wet = prior != 0
    print(f"{name}: mean prior={prior[wet].mean():.4f}  mean final={final[wet].mean():.4f}  "
          f"mean |delta|={np.abs(final[wet]-prior[wet]).mean():.4f}  "
          f"max |delta|={np.abs(final[wet]-prior[wet]).max():.4f}")
