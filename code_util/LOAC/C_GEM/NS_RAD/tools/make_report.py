"""
Combine the six NS-RAD section PDFs into one document, docs/ns_rad_report.pdf,
in reading order (headline summary -> schematic -> geometry -> river networks ->
diagnostics -> validation). This combined report is the kept deliverable; the six
section PDFs are throwaway intermediates that build_all.sh deletes after this stitch.

A missing input is skipped with a warning rather than aborting, so it still works
after a partial build (e.g. no runs/regression_bnd -> no validation PDF). Prefers the
`pypdf` package if installed; otherwise falls back to poppler's `pdfunite` on PATH.

Runs as the final figure step of tools/build_all.sh. Because the section PDFs are
deleted after each build, rebuild the report with `tools/build_all.sh --figures-only`
(which recreates them first) rather than calling this script on its own.
"""
import shutil
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
DOCS = ROOT / "docs"
OUT = DOCS / "ns_rad_report.pdf"

# The individual report PDFs, in the order they should appear in the combined document.
PARTS = [
    "ns_rad_model_summary.pdf",     # headline summary
    "ns_rad_model_schematic.pdf",   # framework + data flow + call graph
    "ns_rad_geometry.pdf",          # channel geometry per river
    "ns_rad_river_networks.pdf",    # river-network maps
    "ns_rad_diagnostics.pdf",       # modelled fields, all four rivers
    "ns_rad_validation.pdf",        # model vs observations
]


def main():
    present = [DOCS / n for n in PARTS if (DOCS / n).exists()]
    for n in PARTS:
        if not (DOCS / n).exists():
            print(f"  [make_report] skipping missing {n}")
    if not present:
        print("  [make_report] no input PDFs found -- build the figures first; nothing written")
        return 1

    try:
        from pypdf import PdfWriter
        w = PdfWriter()
        for p in present:
            w.append(str(p))
        with open(OUT, "wb") as f:
            w.write(f)
    except ImportError:
        exe = shutil.which("pdfunite")
        if not exe:
            print("  [make_report] need either the `pypdf` package or poppler's "
                  "`pdfunite` on PATH", file=sys.stderr)
            return 2
        subprocess.run([exe, *(str(p) for p in present), str(OUT)], check=True)

    print(f"  [make_report] wrote {OUT.relative_to(ROOT)} "
          f"({len(present)} of {len(PARTS)} parts, {OUT.stat().st_size / 1e6:.1f} MB)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
