"""
Before/after for adopting the multi-channel width-role separation, all four rivers.

  before  runs/singlechannel_archive/<site>     legacy: B is the delta SUM at the mouth
                                                and the PER-CHANNEL width upstream
  after   runs/multichannel_preicefix/<site>    B is the TOTAL conveyance width; the
                                                per-thread width B/n_chan feeds dispersion

NOTE the "after" side is deliberately NOT runs/definitive. Both sides here carry the
PRE-FIX bottom-fast ice cap, so this stays a controlled one-variable comparison of the
width change alone. The ice fix landed afterwards and is isolated by
tools/compare_icefix.py (multichannel_preicefix -> definitive).

The expected signature, from config.MULTICHANNEL: DEPTH is untouched (ZZ = B*depth), so
the per-AREA carbon flux should barely move while the BASIN-INTEGRATED flux tracks the
corrected water-surface area. Kuparuk and Canning are near-single-thread controls -- their
upstream width barely changes, so only their mouth per-thread width (dispersion) does.

Usage: python tools/compare_adoption.py
"""
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "tools"))
from compare_multichannel import load, summarize  # noqa: E402

SITES = ["colville", "kuparuk", "sagavanirktok", "canning"]
BEFORE = ROOT / "runs/singlechannel_archive"
AFTER = ROOT / "runs/multichannel_preicefix"


def main():
    rows = []
    for s in SITES:
        b_p, a_p = BEFORE / s / "output.nc", AFTER / s / "output.nc"
        if not (b_p.exists() and a_p.exists()):
            print(f"skip {s}: missing {'before' if not b_p.exists() else 'after'} run")
            continue
        rows.append((s, summarize(s, load(b_p)), summarize(s, load(a_p))))

    hdr = (f"{'site':15s}{'surface km2':>22s}{'gC/m2/yr (per-area)':>26s}"
           f"{'basin tC/yr':>24s}{'S summer med':>20s}")
    sub = (f"{'':15s}{'before   after  ratio':>22s}{'before   after  ratio':>26s}"
           f"{'before   after  ratio':>24s}{'before  after':>20s}")
    print("\nAdopting multi-channel geometry — year 2, all four rivers\n" + "=" * len(hdr))
    print(hdr); print(sub); print("-" * len(hdr))
    for s, b, a in rows:
        print(f"{s:15s}"
              f"{b['surf_km2']:8.1f}{a['surf_km2']:8.1f}{a['surf_km2']/b['surf_km2']:6.2f}x"
              f"{b['FCO2_ann']:10.1f}{a['FCO2_ann']:8.1f}{a['FCO2_ann']/b['FCO2_ann']:7.2f}x"
              f"{b['FCO2_total']:9.0f}{a['FCO2_total']:8.0f}{a['FCO2_total']/b['FCO2_total']:6.2f}x"
              f"{b['S_summer_med']:11.2f}{a['S_summer_med']:8.2f}")
    print("\nPer-area flux should stay ~flat (DEPTH unchanged); the basin total is the"
          "\nreal correction. Salinity is NOT constrained -- see docs/multichannel_test.md.")


if __name__ == "__main__":
    main()
