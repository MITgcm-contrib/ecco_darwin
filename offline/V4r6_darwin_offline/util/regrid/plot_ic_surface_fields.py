"""Surface-field QC plots for every PTRACERS initial condition -- ALL 36
tracers, one page each, into a single PDF.

Scope: initial conditions only (data.ptracers), not forcing (wind/iron/
apCO2/oasim/SIarea are time-varying boundary conditions, a different
category -- not plotted here).

Two kinds of page, one per tracer (a third, "zero", existed while c01-c06
were briefly cold-started at exactly zero; reverted, see
data.ptracers's own header for why -- growth in darwin_plankton.F is
purely multiplicative, so a group starting at exactly zero can never
grow, confirmed empirically after 120 hours of real integration):
  - "file": read directly from regridded_llc90/ (already on the LLC90
    grid -- no regridding happens in this script), surface (level 0),
    plotted via ecco_v4_py.plot_proj_to_latlon_grid (nearest-neighbor
    reprojection onto a Robinson map with real coastlines via cartopy).
    Now includes all 10 plankton tracers (c01-c10), not just c07-c10.
  - "derived_unknown": Chl01-06 are not file-based at all -- darwin_init_chl.F
    (pkg/darwin) unconditionally recomputes them at cold start as
    chltmp=MAX(carbon*chl2cmin,...) from the c01-c06 biomass, which is now
    real/nonzero -- genuine "instant acclimation" (Steph's own term), not
    a guaranteed-zero value. The real cold-start value isn't computable
    here without reimplementing that formula, so this still plots as a
    placeholder zero field -- read as "not shown", not "confirmed zero".

Land (hFacC==0 at the surface) is masked to NaN on every page so it reads
as blank, not a spurious 0.
"""
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import ecco_v4_py as ecco
import regrid_to_llc90 as r

IC_DIR = "regridded_llc90"
OUT_PDF = "figures/IC_surface_fields.pdf"

# (index, name, units, kind, filename-or-None, source label, note)
# Cross-checked against input_darwin_offline/data.ptracers as of.
TRACERS = [
    (1, "DIC", "mmol C/m^3", "file", "DIC_llc270_to_llc90.bin",
     "LLC270-native (nearest-wet-neighbor regrid)", "steph_input/v06/ICs/DIC.R4_270x3510x50.bin"),
    (2, "NO3", "mmol N/m^3", "file", "lev01_nitrate_ann-3d_llc90.bin",
     "1-degree Levitus (lat-lon regrid)", "lev01_nitrate_ann-3d.bin"),
    (3, "NO2", "mmol N/m^3", "file", "ptr03_run33_4_28800_llc90.bin",
     "1-degree DARWIN_ECCO_run33_4 restart (lat-lon regrid)", "ptr03_run33_4_28800.bin"),
    (4, "NH4", "mmol N/m^3", "file", "ptr04_run33_4_28800_llc90.bin",
     "1-degree DARWIN_ECCO_run33_4 restart (lat-lon regrid)", "ptr04_run33_4_28800.bin"),
    (5, "PO4", "mmol P/m^3", "file", "lev01_phosphate_ann-3d_llc90.bin",
     "1-degree Levitus (lat-lon regrid)", "lev01_phosphate_ann-3d.bin"),
    (6, "FeT", "mmol Fe/m^3", "file", "ptr06_run33_4_28800_llc90.bin",
     "1-degree DARWIN_ECCO_run33_4 restart (lat-lon regrid)", "ptr06_run33_4_28800.bin"),
    (7, "SiO2", "mmol Si/m^3", "file", "lev01_silicate_ann-3d_llc90.bin",
     "1-degree Levitus (lat-lon regrid)", "lev01_silicate_ann-3d.bin"),
    (8, "DOC", "mmol C/m^3", "file", "ptr08_run33_4_28800_llc90.bin",
     "1-degree DARWIN_ECCO_run33_4 restart (lat-lon regrid)", "ptr08_run33_4_28800.bin"),
    (9, "DON", "mmol N/m^3", "file", "ptr09_run33_4_28800_llc90.bin",
     "1-degree DARWIN_ECCO_run33_4 restart (lat-lon regrid)", "ptr09_run33_4_28800.bin"),
    (10, "DOP", "mmol P/m^3", "file", "ptr10_run33_4_28800_llc90.bin",
     "1-degree DARWIN_ECCO_run33_4 restart (lat-lon regrid)", "ptr10_run33_4_28800.bin"),
    (11, "DOFe", "mmol Fe/m^3", "file", "ptr11_run33_4_28800_llc90.bin",
     "1-degree DARWIN_ECCO_run33_4 restart (lat-lon regrid)", "ptr11_run33_4_28800.bin"),
    (12, "POC", "mmol C/m^3", "file", "ptr12_run33_4_28800_llc90.bin",
     "1-degree DARWIN_ECCO_run33_4 restart (lat-lon regrid)", "ptr12_run33_4_28800.bin"),
    (13, "PON", "mmol N/m^3", "file", "ptr13_run33_4_28800_llc90.bin",
     "1-degree DARWIN_ECCO_run33_4 restart (lat-lon regrid)", "ptr13_run33_4_28800.bin"),
    (14, "POP", "mmol P/m^3", "file", "ptr14_run33_4_28800_llc90.bin",
     "1-degree DARWIN_ECCO_run33_4 restart (lat-lon regrid)", "ptr14_run33_4_28800.bin"),
    (15, "POFe", "mmol Fe/m^3", "file", "ptr15_run33_4_28800_llc90.bin",
     "1-degree DARWIN_ECCO_run33_4 restart (lat-lon regrid)", "ptr15_run33_4_28800.bin"),
    (16, "POSi", "mmol Si/m^3", "file", "ptr16_run33_4_28800_llc90.bin",
     "1-degree DARWIN_ECCO_run33_4 restart (lat-lon regrid)", "ptr16_run33_4_28800.bin"),
    (17, "PIC", "mmol C/m^3", "file", "ptr17_run33_4_28800_llc90.bin",
     "1-degree DARWIN_ECCO_run33_4 restart (lat-lon regrid)", "ptr17_run33_4_28800.bin"),
    (18, "ALK", "mmol eq/m^3", "file", "ALK_llc270_to_llc90.bin",
     "LLC270-native (nearest-wet-neighbor regrid)", "steph_input/v06/ICs/ALK.R4_270x3510x50.bin"),
    (19, "O2", "mmol O/m^3", "file", "lev01_oxygen_micromolar_ann-3d_llc90.bin",
     "1-degree Levitus (lat-lon regrid)", "lev01_oxygen_micromolar_ann-3d.bin"),
    (20, "CDOM", "mmol C/m^3", "file", "ptr20_run33_4_28800_llc90.bin",
     "1-degree DARWIN_ECCO_run33_4 restart (lat-lon regrid)", "ptr20_run33_4_28800.bin"),
]

# c01-c10 (6 phytoplankton + 4 Zoo): all 10 read the SAME restart file
# (READ_FLD_XYZ_RL always reads record 1, and PTRACERS_initialFile
# repeats one filename), matching Steph's own original data.ptracers
# exactly. REVERTED: c01-c06 briefly
# had no PTRACERS_initialFile at all, cold-starting at exactly zero --
# caught on the first real run attempt: darwin_plankton.F's growth/uptake
# is purely multiplicative (every term scales with the group's own
# existing biomass), so a group cold-started at exactly zero can never
# grow (confirmed empirically -- still bit-exact zero after 120 real
# hours of integration) and would have stayed zero for the entire 34-year
# run, not just at t=0. Reverted to match Steph's real practice.
for i, cname in zip(range(21, 31), ["c01", "c02", "c03", "c04", "c05", "c06",
                                     "c07", "c08", "c09", "c10"]):
    grp = (["PicoCyano", "PicoEuk", "Cocco", "Diazo", "Diatom", "Dino"][i - 21]
           if i <= 26 else "Zoo")
    TRACERS.append((
        i, f"{cname} ({grp})", "mmol C/m^3", "file", "biomass_run33_4_28800_llc90.bin",
        "1-degree DARWIN_ECCO_run33_4 restart (lat-lon regrid)",
        "biomass_run33_4_28800.bin -- same file/record read for all 10 "
        "plankton tracers (PTRACERS_initialFile(21:30) repeats one "
        "filename), so c01-c10 are identical fields at t=0 by Steph's own "
        "original design.",
    ))

# Chl01-06: not file-based -- darwin_init_chl.F unconditionally computes
# chltmp=MAX(carbon*chl2cmin,...) at cold start from the REAL c01-c06
# biomass above (now nonzero), i.e. genuine "instant acclimation" (Steph's
# own term) -- NOT guaranteed zero (that was only true -
#, when the underlying carbon was itself zero). The real
# cold-start value isn't computable here without reimplementing that
# formula, so this still plots as a placeholder zero field -- read as
# "not shown", not "confirmed zero".
for i, cname in zip(range(31, 37), ["Chl01", "Chl02", "Chl03", "Chl04", "Chl05", "Chl06"]):
    TRACERS.append((
        i, cname, "mg Chl a/m^3", "derived_unknown", None,
        "Derived at cold start (darwin_init_chl.F) from real biomass, not read from a file",
        "darwin_init_chl.F (pkg/darwin) computes chltmp=MAX(carbon*chl2cmin,...) "
        "at nIter0 from the real (now nonzero) c01-c06 biomass above -- "
        "genuine acclimation, not a guaranteed-zero value. Plotted as a "
        "placeholder zero field below only because the real value isn't "
        "computable without running the model; do not read this page as "
        "confirming Chl is actually zero.",
    ))

assert len(TRACERS) == 36, len(TRACERS)
assert [t[0] for t in TRACERS] == list(range(1, 37)), [t[0] for t in TRACERS]


def load_surface(filename, wet0):
    path = os.path.join(IC_DIR, filename)
    raw = np.fromfile(path, dtype=">f4").reshape(50, 1170, 90)
    surf = ecco.llc_compact_to_tiles(raw[0], less_output=True)  # (13,90,90)
    return np.where(wet0, surf, np.nan)


def zero_surface(wet0):
    return np.where(wet0, 0.0, np.nan)


def plot_field_page(pdf, xc, yc, field, title, subtitle, cmap="viridis"):
    fig = plt.figure(figsize=(11, 7))
    finite = field[np.isfinite(field)]
    if finite.size and finite.max() > finite.min():
        cmin, cmax = finite.min(), finite.max()
    else:
        # uniform field (e.g. all-zero): give pcolormesh a nonzero span
        center = finite[0] if finite.size else 0.0
        cmin, cmax = center - 1.0, center + 1.0
    # dx/dy=0.5 (not the 0.25 default) and rasterizing the data mesh below
    # are both needed -- without them, the first attempt at this script
    # produced a >500MB, still-growing PDF after 8+ minutes: pcolormesh at
    # 0.25 degrees on a Robinson projection embeds each of ~1e6 grid cells
    # as a separate vector polygon, per page, times 36 pages. Rasterizing
    # only the data mesh keeps coastlines/gridlines/colorbar as clean
    # vectors while collapsing the data layer to a single embedded image.
    result = ecco.plot_proj_to_latlon_grid(
        xc, yc, field, projection_type="robin", show_colorbar=True,
        cmap=cmap, cmin=cmin, cmax=cmax, dx=0.5, dy=0.5, less_output=True,
    )
    mesh = result[2]
    if mesh is not None:
        mesh.set_rasterized(True)
    fig.suptitle(title, fontsize=13, y=0.98)
    plt.figtext(0.5, 0.02, subtitle, ha="center", fontsize=8, wrap=True)
    pdf.savefig(fig, dpi=150)
    plt.close(fig)


def main():
    os.makedirs(os.path.dirname(OUT_PDF), exist_ok=True)
    xc, yc, hfacc = r.target_grid()
    wet0 = hfacc[0] > 0
    zero_field = zero_surface(wet0)

    with PdfPages(OUT_PDF) as pdf:
        # Title / legend page
        fig = plt.figure(figsize=(11, 7))
        legend = (
            "6+4+0 offline Darwin on ECCOv4r6/LLC90 -- ALL 36 PTRACERS "
            "initial conditions, surface (level 0),.\n\n"
            "Source categories:\n"
            "  - LLC270-native: real cube-sphere LLC270 field, regridded to "
            "LLC90 by nearest-wet-neighbor (regrid_llc270_to_llc90). "
            "[DIC, ALK]\n"
            "  - 1-degree (lat-lon regrid): Steph's original 1x1-degree "
            "field, regridded to LLC90 via RegularGridInterpolator, with "
            "source-side land-fill and polar-clamp fixes applied. "
            "[18 tracers + c01-c10]\n"
            "  - Derived (unknown value here): computed internally by "
            "darwin_init_chl.F at cold start from the real c01-c06 biomass "
            "-- genuine acclimation, plotted as a placeholder zero field "
            "only because the real value isn't computable without running "
            "the model. [Chl01-06]\n\n"
            "Land (hFacC==0 at the surface) is masked to blank on every map "
            "-- it is not a literal zero value in the underlying file.\n\n"
            "NOTE: c01-c06 were briefly cold-started at exactly zero; "
            "reverted to the real biomass restart "
            "(matching c07-c10 and Steph's own original config) after "
            "confirming a group starting at exactly zero can never grow "
            "under darwin_plankton.F's purely multiplicative growth "
            "kinetics -- see data.ptracers's own header for the full story."
        )
        fig.text(0.5, 0.9, "PTRACERS Initial Conditions -- Surface Fields (all 36)",
                  ha="center", va="top", fontsize=16, weight="bold")
        fig.text(0.1, 0.75, legend, ha="left", va="top", fontsize=11, wrap=True)
        pdf.savefig(fig)
        plt.close(fig)

        for idx, name, units, kind, filename, source, note in TRACERS:
            if kind == "file":
                field = load_surface(filename, wet0)
            else:
                field = zero_field
            title = f"Tracer {idx}: {name}  [{units}]"
            subtitle = f"Source: {source}\n{note}"
            plot_field_page(pdf, xc, yc, field, title, subtitle)
            finite = field[np.isfinite(field)]
            print(f"Plotted tracer {idx} ({name}, {kind}): "
                  f"min/max {finite.min():.4g}/{finite.max():.4g}"
                  if finite.size else f"Plotted tracer {idx} ({name}, {kind}): no wet points")

    print(f"\nDone: {OUT_PDF} -- {len(TRACERS)} tracer pages + 1 legend page")


if __name__ == "__main__":
    main()
