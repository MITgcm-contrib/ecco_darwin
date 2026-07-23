#!/usr/bin/env python3
"""Derive a river carbonate boundary (ALK, DIC, pH) from USGS discrete water-quality
samples, via the Water Quality Portal (waterqualitydata.us).

Motivation: the shipped boundaries used median ALK + median pH solved SEPARATELY, which
understates the flux-relevant pCO2 -- pH and ALK covary, and the median pH is biased by
high-pH low-flow samples while the CO2-rich moments drive the air-sea flux. This tool
instead computes pCO2 PER SAMPLE from co-located (same-activity) ALK + pH + T, takes the
median, and back-solves DIC from the well-constrained median ALK and that median pCO2.
The result reproduces the observed CO2 supersaturation these rivers actually carry.

Alkalinity is reported by USGS as mg/L CaCO3; 50.04 mg/L CaCO3 = 1 meq/L, and at S=0
(rho~1) meq/L = mmol/m^3 = the model unit. Carbonate constants are fun_module's (Weiss
1974 / Millero 1995) at S=0, pb=0 (a shallow river; the pb pressure term is ~1).

Usage:
    python tools/usgs_carbonate_boundary.py siteid USGS-15896000     # Kuparuk gauge
    python tools/usgs_carbonate_boundary.py bbox -146.5,69.0,-143.5,70.6   # Canning/ANWR pool

Re-run after any change to fun_module's carbonate constants (e.g. the _pbar_rho fix).
"""
import sys
import csv
import io
import math
import urllib.request
from collections import defaultdict
from statistics import median

BASE = ("https://www.waterqualitydata.us/data/Result/search?{sel}"
        "&characteristicName=Alkalinity;pH;Temperature,%20water"
        "&mimeType=csv&zip=no&dataProfile=resultPhysChem")


def _k0(T):
    tk = 273.15 + T; t1 = tk / 100
    return math.exp(93.4517 / t1 - 60.2409 + 23.3585 * math.log(t1))

def _k1(T):
    tk = 273.15 + T
    return 10 ** (-(3670.7 / tk - 62.008 + 9.7944 * math.log(tk)))

def _k2(T):
    tk = 273.15 + T
    return 10 ** (-(1394.7 / tk + 4.777))


def _pco2_from_alk_ph(alk_mmol, pH, T):
    """Water pCO2 [uatm] from carbonate alkalinity [mmol/m^3], pH, T [C] at S=0."""
    alk = alk_mmol * 1e-6
    H = 10 ** (-pH)
    hco3 = alk / (1 + 2 * _k2(T) / H)
    co2s = hco3 * H / _k1(T)
    return co2s / _k0(T) * 1e6


def _dic_from_alk_pco2(alk_mmol, pco2, T):
    """Back-solve DIC [mmol/m^3] and pH from ALK [mmol/m^3] and target pCO2 [uatm]."""
    alk = alk_mmol * 1e-6
    co2s = pco2 * 1e-6 * _k0(T)
    k1, k2 = _k1(T), _k2(T)
    H = 1e-8
    for _ in range(200):                       # relax H so HCO3(1+2K2/H) matches ALK
        hco3 = co2s * k1 / H
        H *= (hco3 * (1 + 2 * k2 / H)) / alk
    hco3 = co2s * k1 / H
    co3 = hco3 * k2 / H
    return (co2s + hco3 + co3) * 1e6, -math.log10(H)


def fetch(sel):
    with urllib.request.urlopen(BASE.format(sel=sel), timeout=120) as r:
        return list(csv.DictReader(io.StringIO(r.read().decode("utf-8", "replace"))))


def main(mode, arg):
    sel = f"siteid={arg}" if mode == "siteid" else f"bBox={arg}"
    rows = fetch(sel)
    act = defaultdict(dict)                     # activity -> {characteristic: value}
    alks, temps, nsta = [], [], set()
    for r in rows:
        ch = r["CharacteristicName"]; d = r["ActivityStartDate"]
        try:
            v = float(r["ResultMeasureValue"])
        except (ValueError, KeyError):
            continue
        if not (d and 6 <= int(d[5:7]) <= 9):   # open water Jun-Sep
            continue
        act[r["ActivityIdentifier"]][ch] = v
        if ch == "Alkalinity":
            alks.append(v); nsta.add(r["MonitoringLocationIdentifier"])
        if ch == "Temperature, water":
            temps.append(v)
    pcs = []
    for vals in act.values():
        if "Alkalinity" in vals and "pH" in vals:
            alk = vals["Alkalinity"] / 50.04 * 1000.0
            T = vals.get("Temperature, water", 7.0)
            p = _pco2_from_alk_ph(alk, vals["pH"], T)
            if 0 < p < 1e5:
                pcs.append(p)
    if not alks or not pcs:
        sys.exit(f"{sel}: insufficient co-located Alkalinity+pH samples.")
    alk_med = median(alks) / 50.04 * 1000.0     # mmol/m^3
    T_med = median(temps) if temps else 7.0
    pco2_med = median(pcs)
    dic, pH = _dic_from_alk_pco2(alk_med, pco2_med, T_med)
    print(f"{sel}  open-water (Jun-Sep):")
    print(f"  n(ALK)={len(alks)} from {len(nsta)} station(s)   n(paired ALK+pH)={len(pcs)}")
    print(f"  ALK        = {alk_med:.1f} mmol/m^3   T = {T_med:.1f} C")
    print(f"  median pCO2 = {pco2_med:.0f} uatm   ({'SUPERSATURATED -> outgas' if pco2_med > 415 else 'undersaturated'})")
    print(f"  -> DIC = {dic:.1f} mmol/m^3, pH = {pH:.2f}  (reproduce median pCO2 at median ALK)")
    print()
    print(f'  BOUNDARIES: ("pH", {pH:.2f}), ("ALK", {alk_med:.1f}), ("DIC", {dic:.1f})')


if __name__ == "__main__":
    if len(sys.argv) == 3 and sys.argv[1] in ("siteid", "bbox"):
        main(sys.argv[1], sys.argv[2])
    else:
        # back-compat: a bare gauge number
        main("siteid", f"USGS-{sys.argv[1]}" if len(sys.argv) > 1 else "USGS-15896000")
