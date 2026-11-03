#!/usr/bin/env python3
"""Kelp blade optics for radtrans -> ../input/data.kelpoptics (kelp_optics_mode = 1).

Default (since 2026-10-07): MEASURED spectra, 13 radtrans bands from
../work/blade_optics/macrocystis_blade_ART_radtrans13.csv:
  A(l): mean relative absorbance D(l)/D(675) of 9 M. pyrifera blade spectra
        (0-18 m, spring/summer; Colombo-Pallotta 2006, PhD thesis, CICESE,
        Figs. 8, 15-16, pp. 57, 89-91; digitised from the PDF vector paths),
        A = 1 - 10^(-k D), k set so the PAR-band mean A = 0.80 (same thesis)
  R(l): mean of 1,770 measured M. pyrifera canopy-blade reflectance spectra
        (SBC LTER knb-lter-sbc.157.1, Bell & Siegel,
        doi:10.6073/pasta/3b0c41b7b20e87c5d2ecd8e9030607ac); PAR mean 0.038
  T(l) = 1 - A - R, set to 0 where A + R > 1 (406-465 nm, at most -0.016);
        no measured M. pyrifera blade T spectrum was found.
  PAR means: A 0.80, R 0.04, T 0.16. (The thesis PAR-mean R = 0.10 is not
  consistent with the 1,770 SBC spectra: mean 0.038, max 0.075.)
--assumed: the earlier set below.

Earlier (--assumed): blade reflectance R is spectrally flat, R = R_PAR; absorptance
A(l) = (1 - R) (1 - exp(-tau s(l))), with s(l) the spectral shape of the
diatom absorption spectrum in OPTICS_COEFF2/optics_plankton.txt (same pigments
as brown algae: chl a, chl c, fucoxanthin; ASSUMED shape, no absolute
Macrocystis blade spectrum found) normalised to a mean of 1, tau set so that
the PAR-band mean absorptance equals A_PAR; transmittance T = 1 - A - R.

  A_PAR = 0.80       : measured PAR-mean absorptance of M. pyrifera blades,
                       0-18 m, Baja California, no depth trend (Colombo-Pallotta
                       2006, PhD thesis, CICESE, Fig. 9B: "La absortancia en la
                       region PAR en promedio fue de 80%")
  R_PAR = 0.10       : measured PAR-mean blade reflectance, 10% (+/-0.7), spectrally
                       and with depth uniform (same thesis, text with Fig. 9)
                       -> T_PAR = 0.10. For comparison, Fucus vesiculosus: A 0.92,
                       R ~0.03, T 0.05 (T 0.10-0.13 near 570 nm, ~0 in blue/red;
                       Lichtenberg et al. 2017, Plant Physiol. 175:721, Fig. 4A);
                       14 tropical macroalgae: A 0.64-0.92, R 0.07-0.23 (Beach et
                       al. 2006, Phycologia 45:450)
  sfa  = 1/58 m2/gDW : blade area per dry weight, Macrocystis dry_sa = 58 gDW/m2
                       in G-MACMODS (mag_species.py; Arzeno-Soltero et al. 2023,
                       Commun. Earth Environ. 4:185, Suppl. Table 1, range
                       43.5-72.5); measured 44 gDW/m2 in the Santa Barbara Channel
                       (Rodriguez 2014, UCSB PhD thesis, Fig. 1D)
  proj = 0.5         : ASSUMED mean projection of randomly oriented blades

Also writes the particle backscatter slope [K6] and the river CDOM fraction
[K7] (namelists DARWIN_BBP_OPTICS, DARWIN_CDOM_RIVER):
  bbp_eta = 1.0      : median PACE OCI bbp_s over Santa Monica Bay, Jan 2025
                       (IQR 0.76-1.13; our PACE data); typical of the literature
                       (Sbp median 1, 80% within 0.5-1.5: Werdell et al. 2018,
                       Prog. Oceanogr. 160:186, citing Reynolds et al. 2016);
                       bbp_lambda0 = 550 nm (ASSUMED reference, mid-visible)
  river CDOM         : a_CDOM(440) = 2 m-1 in Malibu Creek (PLACEHOLDER: no
                       measured CDOM for Malibu Creek or other southern California
                       creeks found; Harringmeyer et al. 2021, Front. Environ. Sci.
                       9:647966, report higher a_g in north Santa Monica Bay after
                       storms but no creek end-member),
                       as CDOM [mmol P/m3] = a440 / (CDOMcoeff exp(-Sdom (440-450)))
                       divided by the creek DOP (make_rivers.py)

Usage: python3 make_kelp_optics.py [--assumed [A_PAR] [R_PAR]]
"""
import os
import sys

import numpy as np
from scipy.optimize import brentq

HERE = os.path.dirname(os.path.abspath(__file__))
INP = os.path.join(HERE, '..', 'input')
WL = np.arange(400, 701, 25)          # radtrans RT_wbRefWLs (13 bands)
BBP_ETA, BBP_LAMBDA0 = 1.0, 550.0
TSS_BETA = round(float(np.log(1.025) / np.log(555.0 / 412.0)), 3)   # Babin et al. 2003 Table 5
RIVER_A440 = 2.0                      # m-1, PLACEHOLDER
CDOM_COEFF, CDOM_SDOM = 100.0, 0.014  # darwin3 defaults (a at 450 nm)


def diatom_shape():
    lines = open(os.path.join(INP, 'OPTICS_COEFF2', 'optics_plankton.txt')).read().splitlines()
    i = next(n for n, l in enumerate(lines) if l.startswith('*** Diatom'))
    tab = {}
    for l in lines[i + 2:i + 2 + len(WL)]:
        w, a = l.split()[:2]
        tab[int(w)] = float(a)
    a = np.array([tab[w] for w in WL])
    return a / a.mean()


def measured():
    import csv
    rows = list(csv.DictReader(open(os.path.join(HERE, '..', 'work', 'blade_optics',
                                                 'macrocystis_blade_ART_radtrans13.csv'))))
    assert [int(float(r['wavelength_nm'])) for r in rows] == list(WL)
    return (np.array([float(r[k]) for r in rows]) for k in ('A', 'R', 'T'))


def main():
    assumed = '--assumed' in sys.argv
    args = [x for x in sys.argv[1:] if x != '--assumed']
    if assumed:
        a_par = float(args[0]) if len(args) > 0 else 0.80
        r_par = float(args[1]) if len(args) > 1 else 0.10
        s = diatom_shape()
        tau = brentq(lambda t: np.mean((1 - r_par) * (1 - np.exp(-t * s))) - a_par, 1e-3, 50)
        A = (1 - r_par) * (1 - np.exp(-tau * s))
        R = np.full_like(A, r_par)
        Tr = 1 - A - R
        head = (f"# A(l) = (1-R)(1 - exp(-tau s(l))), s = diatom absorption shape (ASSUMED shape),\n"
                f"# PAR-mean A = {a_par}, R = {r_par} flat (M. pyrifera, Colombo-Pallotta 2006 PhD thesis,\n"
                f"# CICESE), T = 1 - A - R")
    else:
        A, R, Tr = measured()
        tau = float('nan')
        head = ("# MEASURED M. pyrifera blades: A shape Colombo-Pallotta 2006 (PhD thesis, CICESE), PAR mean\n"
                "# 0.80; R mean of 1,770 SBC LTER blade spectra (knb-lter-sbc.157.1); T = 1 - A - R (>= 0)")
    fmt = lambda v: ', '.join(f'{x:.4f}' for x in v)
    from make_rivers import RIVERS
    dop = RIVERS[0]['conc']['DOP']                       # mmol P/m3 (Malibu Creek)
    cdom_riv = RIVER_A440 / (CDOM_COEFF * np.exp(-CDOM_SDOM * (440 - 450)))
    frac = cdom_riv / dop
    txt = f"""# Kelp blade optics for radtrans (python/make_kelp_optics.py; malibu_v01 [K1,K2])
# bands: {', '.join(map(str, WL))} nm
{head}
# sfa = 1/58 m2/gDW (G-MACMODS dry_sa), proj = 0.5 (ASSUMED)
 &DARWIN_KELP_OPTICS
 kelp_optics_mode = 1,
 kelp_proj = 0.5,
 kelp_sfa = {1 / 58:.6f},
 kelp_absBlade = {fmt(A)},
 kelp_refBlade = {fmt(R)},
 kelp_trnBlade = {fmt(Tr)},
 /
# TSS (mineral particles, g/m3) [K4]: a*(l) = a443 exp(-S (l - 443)), a443 = 0.05 m2/g
# (within 0.03-0.1 m2/g at 443 nm for most natural mineral assemblages, Babin & Stramski
# 2004, L&O 49:756, p.765), S = 0.012 nm-1 (typical S_NAP, Babin et al. 2003, L&O 48:843,
# Eq. 14); b*(l) = 0.5 (l/555)^-{TSS_BETA} m2/g (case-2 b_pm(555) 0.51 m2/g, Table 4;
# spectral exponent from bp(412):bp(555) = 1.025, all stations, Table 5, same paper);
# bb ratio 0.025 (low end of 0.025-0.044 for minerals, Mie base model, Wozniak & Stramski
# 2004, Appl. Opt. 43:3489, p.3494)
 &DARWIN_TSS_OPTICS
 tss_optics_mode = 1,
 tss_aStar = {fmt(0.05 * np.exp(-0.012 * (WL - 443)))},
 tss_bStar = {fmt(0.5 * (WL / 555.0) ** -TSS_BETA)},
 tss_bbRatio = 0.025,
 /
# particle backscatter slope [K6]: bbp(l) = bbp_darwin (l / bbp_lambda0)^-bbp_eta
# bbp_eta = median PACE OCI bbp_s, Santa Monica Bay Jan 2025 (IQR 0.76-1.13)
 &DARWIN_BBP_OPTICS
 bbp_eta = {BBP_ETA},
 bbp_lambda0 = {BBP_LAMBDA0},
 /
# river CDOM [K7] = cdom_riv_frac x DOP runoff: a_CDOM(440) = {RIVER_A440} m-1 (PLACEHOLDER)
# -> CDOM {cdom_riv:.4f} mmol P/m3 / creek DOP {dop:.4f} mmol P/m3
 &DARWIN_CDOM_RIVER
 cdom_riv_frac = {frac:.5f},
 /
"""
    open(os.path.join(INP, 'data.kelpoptics'), 'w').write(txt)
    print(txt)
    print(f'tau = {tau:.3f}; A {A.min():.2f}-{A.max():.2f} (PAR mean {A.mean():.2f}), '
          f'T {Tr.min():.2f}-{Tr.max():.2f} (PAR mean {Tr.mean():.2f}), R PAR mean {R.mean():.3f}; '
          f'TSS b* exponent {TSS_BETA}; Lk at B = 100 gDW/m3: {0.5 / 58 * 100:.2f} m-1')


if __name__ == '__main__':
    main()
