#!/usr/bin/env python3
"""plot_nares_grid.py -- review figure for the Nares Strait grid + bathymetry
(same layout as west_AO/grid/plot_westAO_grid.py).

(a) model bathymetry on a map, with the wet open-boundary cells of each edge;
(b) the sNx x sNy tiling in index space, active vs all-land (blank) tiles;
(c) wet cells across the grid at each along-axis position (channel width);
(d) Robeson Channel / Hall Basin at native resolution, every cell drawn.

Usage: python3 plot_nares_grid.py <grid dir> <out.png> [sNx sNy]
"""
import os
import sys
import json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm, ListedColormap
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
import cartopy.crs as ccrs
import cartopy.feature as cfeature

sys.path.insert(0, "/Users/carrolld/Documents/research/FORTE/ECCO_FORTE/code")
import agu_style as S

OB_COLORS = {"south": "#0072B2", "north": "#D55E00", "west side": "#009E73", "east side": "#CC79A7"}


def main(gdir, out, snx=40, sny=40):
    S.use()
    m = json.load(open(os.path.join(gdir, "Nares_grid_meta.json")))
    nx, ny = m["nx"], m["ny"]
    g = np.fromfile(os.path.join(gdir, "Nares.mitgrid"), ">f8").reshape(16, ny + 1, nx + 1)
    xC, yC = g[0][:-1, :-1], g[1][:-1, :-1]
    xG, yG = g[5], g[6]
    b = np.fromfile(os.path.join(gdir, "Nares_bathymetry.bin"), ">f8").reshape(ny, nx)
    wet = b < 0
    depth = np.where(wet, -b, np.nan)

    fig = plt.figure(figsize=(S.DOUBLE, 7.6))
    proj = ccrs.NorthPolarStereo(central_longitude=-62)
    pc = ccrs.PlateCarree()

    # (a) map
    ax = fig.add_axes([0.06, 0.43, 0.62, 0.52], projection=proj)
    ax.set_extent([-86, -38, 73.8, 85.2], crs=pc)
    ax.add_feature(cfeature.LAND.with_scale("50m"), facecolor=S.LAND, zorder=0)
    pm = ax.pcolormesh(xG, yG, depth, transform=pc, cmap="Blues", norm=LogNorm(10, 2500),
                       rasterized=True, zorder=1)
    ax.add_feature(cfeature.COASTLINE.with_scale("50m"), linewidth=0.3, edgecolor="0.35", zorder=2)
    # i = 0 is the south (North Water) end, i = nx-1 the north (Lincoln Sea) end;
    # j = 0 / ny-1 are the two side edges (named by which shore they parallel)
    side0 = "west side" if np.mean(xC[0, :]) < np.mean(xC[-1, :]) else "east side"
    side1 = "east side" if side0 == "west side" else "west side"
    ob = {"south": (xC[:, 0], yC[:, 0], wet[:, 0]), "north": (xC[:, -1], yC[:, -1], wet[:, -1]),
          side0: (xC[0, :], yC[0, :], wet[0, :]), side1: (xC[-1, :], yC[-1, :], wet[-1, :])}
    handles = []
    for k in ("south", "north", "west side", "east side"):
        x, y, msk = ob[k]
        ax.scatter(x[msk], y[msk], s=1.0, color=OB_COLORS[k], transform=pc, zorder=4, linewidths=0)
        handles.append(Line2D([], [], marker="s", ls="", ms=4, color=OB_COLORS[k],
                              label=f"{k.capitalize()} OB ({int(msk.sum())} wet cells)"))
    ex = np.r_[xG[0, :], xG[:, -1], xG[-1, ::-1], xG[::-1, 0]]
    ey = np.r_[yG[0, :], yG[:, -1], yG[-1, ::-1], yG[::-1, 0]]
    ax.plot(ex, ey, transform=pc, color="k", lw=0.6, zorder=3)
    handles.append(Line2D([], [], color="k", lw=0.6, label="Domain edge"))
    gl = ax.gridlines(draw_labels=True, x_inline=False, y_inline=False, linewidth=0.3, color="0.6",
                      xlocs=np.arange(-90, -30, 10), ylocs=np.arange(75, 86, 2.5))
    gl.top_labels = gl.right_labels = False
    gl.xlabel_style = gl.ylabel_style = {"size": 6, "rotation": 0}
    # --- simulated glaciers -------------------------------------------------
    ice_f = os.path.join(gdir, "Nares_shelficeTopo.bin")
    if os.path.exists(ice_f):
        ice = np.fromfile(ice_f, ">f8").reshape(ny, nx)
        ax.contourf(xC, yC, (ice < 0).astype(float), levels=[0.5, 1.5], colors="none", hatches=["////"],
                    transform=pc, zorder=3.5)
        ax.contour(xC, yC, (ice < 0).astype(float), levels=[0.5], colors="#E69F00", linewidths=0.8,
                   transform=pc, zorder=3.6)
        handles.append(Patch(fc="none", ec="#E69F00", hatch="////", lw=0.8,
                             label=f"Petermann ice-shelf cavity ({int((ice < 0).sum())} cells)"))
    outl = "/Users/carrolld/Documents/research/downscaling/raw_data/mankoff_freshwater_v8/ice/outlets.csv"
    if os.path.exists(outl):
        import pandas as pd
        o = pd.read_csv(outl)
        # 'Sermersuaq' names more than one Greenland glacier: keep Humboldt's, i.e. inside this grid
        o = o[(o.B2015_name == "Sermersuaq") & (o.elev < 0)]
        sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
        import gen_nares_mitgrid as GN
        e1, e2, nn = (np.array(m[k]) for k in ("e1", "e2", "pole"))
        lamO, phiO = GN.to_rot(o.lon.values, o.lat.values, e1, e2, nn)
        lamG, _ = GN.to_rot(xG[0, :], yG[0, :], e1, e2, nn)
        _, phiG = GN.to_rot(xG[:, 0], yG[:, 0], e1, e2, nn)
        o = o[(lamO > lamG[0]) & (lamO < lamG[-1]) & (phiO > phiG[0]) & (phiO < phiG[-1])]
        ax.scatter(o.lon, o.lat, s=6, marker="^", color="#CC0000", transform=pc, zorder=5, linewidths=0)
        handles.append(Line2D([], [], marker="^", ls="", ms=4, color="#CC0000",
                              label=f"Humboldt plume outlets ({len(o)}, Mankoff v8)"))
    # --- regions, glaciers, arches -----------------------------------------
    halo = dict(fontsize=5.5, transform=pc, ha="center", va="center", zorder=6,
                bbox=dict(fc="white", ec="none", alpha=0.75, pad=0.6))

    def channel_centre(i):
        """lon, lat of the centre of the widest contiguous wet run across column i."""
        col = wet[:, i]
        runs, j = [], 0
        while j < ny:
            if col[j]:
                k = j
                while k < ny and col[k]:
                    k += 1
                runs.append((k - j, j, k)); j = k
            else:
                j += 1
        _, j0, j1 = max(runs)
        jc = (j0 + j1 - 1) // 2
        return float(xC[jc, i]), float(yC[jc, i])

    # water bodies placed from the model geometry: i = along-strait column chosen
    # from the channel-width profile (narrows: Smith ~78.3N, Kennedy ~80.5N,
    # Robeson ~82.0N; basins between them). Narrow channels get a leader line.
    ncr = wet.sum(axis=0)
    def narrowest(i0, i1):
        return i0 + int(np.argmin(ncr[i0:i1]))
    i_smith, i_kenn, i_rob = narrowest(370, 430), narrowest(630, 700), narrowest(800, 880)
    wide = [("Baffin Bay", 60), ("North Water", 230), ("Kane Basin", (i_smith + i_kenn) // 2),
            ("Hall Basin", (i_kenn + i_rob) // 2), ("Lincoln Sea", 1000)]
    for name, i in wide:
        lo, la = channel_centre(i)
        ax.text(lo, la, name, **halo)
    for name, i, dlon in (("Smith Sound", i_smith, -6.0), ("Kennedy Ch.", i_kenn, -6.5),
                          ("Robeson Ch.", i_rob, -6.5)):
        lo, la = channel_centre(i)
        ax.annotate(name, xy=(lo, la), xycoords=pc._as_mpl_transform(ax),
                    xytext=(lo + dlon, la), textcoords=pc._as_mpl_transform(ax),
                    fontsize=5.5, ha="right", va="center", zorder=6,
                    bbox=dict(fc="white", ec="none", alpha=0.75, pad=0.6),
                    arrowprops=dict(arrowstyle="-", lw=0.5, color="0.2"))
    for name, lat, lon in (("Ellesmere I.", 79.6, -79.0), ("Devon I.", 75.3, -84.0),
                           ("Greenland", 78.6, -58.0)):
        ax.text(lon, lat, name, style="italic", **halo)
    ax.annotate("Jones Sound", xy=(-81.0, 76.1), xycoords=pc._as_mpl_transform(ax),
                xytext=(-86.5, 77.0), textcoords=pc._as_mpl_transform(ax), fontsize=5.5,
                ha="center", va="center", zorder=6, bbox=dict(fc="white", ec="none", alpha=0.75, pad=0.6),
                arrowprops=dict(arrowstyle="-", lw=0.5, color="0.2"))
    # glaciers: label beside the cavity / the plume outlets
    if os.path.exists(ice_f):
        jj, ii = np.where(ice < 0)
        plo, pla = float(xC[jj, ii].mean()), float(yC[jj, ii].mean())
        ax.annotate("Petermann Gl.", xy=(plo, pla), xycoords=pc._as_mpl_transform(ax),
                    xytext=(plo + 6.0, pla + 0.1), textcoords=pc._as_mpl_transform(ax),
                    fontsize=5.5, fontweight="bold", ha="left", va="center", zorder=6,
                    bbox=dict(fc="white", ec="none", alpha=0.75, pad=0.6),
                    arrowprops=dict(arrowstyle="-", lw=0.5, color="0.2"))
    if os.path.exists(outl) and len(o):
        hlo, hla = float(o.lon.mean()), float(o.lat.mean())
        ax.annotate("Humboldt Gl.", xy=(hlo, hla), xycoords=pc._as_mpl_transform(ax),
                    xytext=(hlo + 6.0, hla), textcoords=pc._as_mpl_transform(ax),
                    fontsize=5.5, fontweight="bold", ha="left", va="center", zorder=6,
                    bbox=dict(fc="white", ec="none", alpha=0.75, pad=0.6),
                    arrowprops=dict(arrowstyle="-", lw=0.5, color="0.2"))
    # typical arch sites: northern arch at the Lincoln Sea end of Robeson Channel,
    # southern arch across Smith Sound / southern Kane Basin
    i_narch = i_rob + int(np.argmax(ncr[i_rob:i_rob + 120] > 150))
    for name, i in (("N arch", i_narch), ("S arch", i_smith)):
        lo, la = channel_centre(i)
        ax.plot(lo, la, marker="D", ms=4, mfc="white", mec="k", mew=0.8, transform=pc, zorder=6)
        ax.text(lo + 3.0, la, name, **{**halo, "ha": "left"})
    handles.append(Line2D([], [], marker="D", ls="", ms=4, mfc="white", mec="k", label="Ice arch site (typical)"))
    ax.legend(handles=handles, loc="lower right", fontsize=5.5, frameon=True, framealpha=0.95, edgecolor="0.6")
    ax.text(0.01, 0.99, f"{nx} x {ny} cells, {m['dx_km']:.0f} km, {int(wet.sum())} wet ({100 * wet.mean():.0f}%)",
            transform=ax.transAxes, ha="left", va="top", fontsize=6, bbox=dict(fc="white", ec="none", alpha=0.85))
    fig.canvas.draw()
    pa = ax.get_position()
    cax = fig.add_axes([pa.x0 + 0.02, pa.y0 - 0.045, pa.width - 0.04, 0.013])
    cb = fig.colorbar(pm, cax=cax, orientation="horizontal", extend="both")
    cb.set_label("Depth (m)")

    # (b) tiling
    axb = fig.add_axes([0.76, 0.62, 0.22, 0.22])
    tiles = wet.reshape(ny // sny, sny, nx // snx, snx).any(axis=(1, 3))
    axb.imshow(np.where(wet, 1.0, 0.0), origin="lower", cmap=ListedColormap([S.LAND, "#9ec5f4"]),
               extent=[0, nx, 0, ny], interpolation="nearest", aspect="auto")
    for j in range(tiles.shape[0]):
        for i in range(tiles.shape[1]):
            if not tiles[j, i]:
                axb.add_patch(plt.Rectangle((i * snx, j * sny), snx, sny, fc="none", ec="0.3", hatch="////", lw=0.0))
    for i in range(0, nx + 1, snx):
        axb.axvline(i, color="0.3", lw=0.3)
    for j in range(0, ny + 1, sny):
        axb.axhline(j, color="0.3", lw=0.3)
    axb.set_xlim(0, nx); axb.set_ylim(0, ny)
    axb.set_xlabel("i (along strait)"); axb.set_ylabel("j")
    axb.legend(handles=[Patch(fc="#9ec5f4", ec="0.3", lw=0.3, label="Wet"),
                        Patch(fc=S.LAND, ec="0.3", lw=0.3, label="Land"),
                        Patch(fc="none", ec="0.3", hatch="////", lw=0.3, label="Blank tile")],
               loc="upper center", bbox_to_anchor=(0.5, -0.32), ncol=3, fontsize=6, frameon=False)
    axb.text(0.5, 1.04, f"sNx={snx}, sNy={sny}\n{int(tiles.sum())} of {tiles.size} tiles active",
             transform=axb.transAxes, ha="center", va="bottom", fontsize=6.5)

    # (c) channel width: wet cells across the grid at each along-axis position
    axc = fig.add_axes([0.08, 0.07, 0.38, 0.24])
    xkm = (np.arange(nx) + 0.5) * m["dx_km"]
    ncross = wet.sum(axis=0)
    axc.plot(xkm, ncross, color="#0072B2", lw=1.0, label="Wet cells across grid")
    axc.axhline(10, color=S.NEUTRAL, lw=0.8, ls="--", label="10 cells")
    axc.set_xlabel("Distance along strait axis from S edge (km)")
    axc.set_ylabel("Wet cells across (= km)")
    axc.set_xlim(0, xkm[-1]); axc.set_ylim(0, ny)
    axc.grid(lw=0.3, color="0.85")
    # mark the narrowest point within the strait proper (middle 60% of the axis)
    mid = slice(int(0.2 * nx), int(0.8 * nx))
    k = int(np.argmin(ncross[mid])) + mid.start
    axc.plot(xkm[k], ncross[k], "o", ms=4, color="#D55E00",
             label=f"Narrowest: {int(ncross[k])} cells, {yC[ny // 2, k]:.1f}N")
    axc.legend(fontsize=6, loc="upper center")

    # (d) Robeson Channel / Hall Basin zoom, cell edges drawn
    axd = fig.add_axes([0.56, 0.06, 0.42, 0.27], projection=proj)
    # centre on the narrowest section found for panel (c)
    jw = np.where(wet[:, k])[0]
    lat0, lon0 = float(yC[jw, k].mean()), float(xC[jw, k].mean())
    axd.set_extent([lon0 - 2.2, lon0 + 2.2, lat0 - 0.30, lat0 + 0.30], crs=pc)
    sel = (np.abs(yC - lat0) < 0.5) & (np.abs(xC - lon0) < 3.5)
    jj, ii = np.where(sel)
    j0, j1, i0, i1 = jj.min(), jj.max() + 1, ii.min(), ii.max() + 1
    XG, YG = xG[j0:j1 + 1, i0:i1 + 1], yG[j0:j1 + 1, i0:i1 + 1]
    sub = depth[j0:j1, i0:i1]
    axd.pcolormesh(XG, YG, np.where(np.isnan(sub), 0, np.nan), transform=pc,
                   cmap=ListedColormap([S.LAND]), edgecolors="0.6", linewidth=0.15, rasterized=True, zorder=0)
    axd.pcolormesh(XG, YG, sub, transform=pc, cmap="Blues", norm=LogNorm(10, 2500),
                   edgecolors="0.6", linewidth=0.15, rasterized=True, zorder=1)
    axd.add_feature(cfeature.COASTLINE.with_scale("10m"), linewidth=0.5, edgecolor="k", zorder=2)
    gl = axd.gridlines(draw_labels=True, x_inline=False, y_inline=False, linewidth=0.3, color="0.6",
                       xlocs=np.arange(-70, -50, 1), ylocs=np.arange(80, 84, 0.25))
    gl.top_labels = gl.right_labels = False
    gl.xlabel_style = gl.ylabel_style = {"size": 6, "rotation": 0}
    axd.legend(handles=[Patch(fc="#6da7ec", ec="0.55", lw=0.3, label="Ocean cell (1 km)"),
                        Patch(fc=S.LAND, ec="0.55", lw=0.3, label="Land cell"),
                        Line2D([], [], color="k", lw=0.5, label="Coastline (NE 10m)")],
               loc="lower left", fontsize=6, frameon=True, framealpha=0.95, edgecolor="0.6")

    # panel letters in FIGURE coordinates from each panel's drawn box (the map axes
    # shrink to keep their aspect, so axes-fraction offsets land in the wrong place)
    fig.canvas.draw()
    for a, l in zip([ax, axb, axc, axd], "abcd"):
        bb = a.get_position()
        dx = {"a": 0.0, "b": -0.045, "c": -0.05, "d": 0.0}[l]
        dy = {"a": 0.008, "b": 0.045, "c": 0.008, "d": 0.008}[l]   # (b) clears its 2-line title
        fig.text(bb.x0 + dx, bb.y1 + dy, f"({l})", fontsize=9, fontweight="bold", ha="left", va="bottom")
    fig.savefig(out, dpi=S.DPI)
    print("wrote", out, "| narrowest", int(ncross[k]), "cells at", round(float(yC[ny // 2, k]), 2), "N")


if __name__ == "__main__":
    a = sys.argv[1:]
    main(a[0], a[1], *(int(v) for v in a[2:4]))
