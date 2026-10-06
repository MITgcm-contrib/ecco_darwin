# Reading and plotting MITgcm / ECCO output

## Formats

- **MDS**: `name.<iter10>.data` + `.meta`. Big-endian, `>f4` (float32) or `>f8`; the
  `.meta` gives `dimList`, `dataprec`, `nrecords`, `fldList` (diagnostics) and timeInterval.
  With `useSingleCpuIO=.FALSE.` there is one file per tile (`.001.001.data`); readers
  stitch them.
- **MNC/NetCDF**: per-tile `*.t001.nc` files unless gluing; ECCO PO.DAAC products are
  already NetCDF on the native LLC grid or regridded 0.5°.
- Grid: `XC, YC, XG, YG, RC, RF, DRF, hFacC/W/S, rA, DXG, DYG, Depth` written at start
  (or in the ECCO grid NetCDF). Use hFac-weighted volumes for averages and budgets.

## Python (preferred)

- `from MITgcmutils import rdmds` — `rdmds('dir/T', itrs=np.nan)` reads all iterations;
  `rdmds('diags/state_3d', 36000, rec=[0,1])` for selected records; returns (nrec,nz,ny,nx).
- `xmitgcm.open_mdsdataset(run_dir, grid_dir=..., geometry='llc', delta_t=..., ref_date=...)`
  for xarray datasets with grid metrics; `geometry='llc'` handles faces.
- `ecco_v4_py`: `llc_compact_to_tiles(arr)` → (13, 90, 90) tiles;
  `plot_proj_to_latlon_grid(XC, YC, field, ...)` for maps; budget/transport helpers.
- Minimal reader in the user's style: `ECCO/BBL/figures/scripts/mitio.py` (parses `.meta`, reads
  global files or a single `.001.001` tile, no multi-tile stitching; shared plotting style +
  categorical palette). Reuse it for global output; use rdmds/xmitgcm for multi-tile output.
- LLC regridding: nearest-wet-neighbour with `scipy.spatial.cKDTree` on XC/YC (see
  `ECCO/offline/scripts/regrid_to_llc90.py`); mask land with hFacC.
- EOS: compute density/sigma with the model's own EOS (e.g. JMD95Z for ECCO v4), not TEOS-10.

LLC compact layout (LLC90): faces 1–2 are 90×270 each, face 3 (Arctic cap) 90×90, faces
4–5 are 270×90 and rotated. Compact array shape is (…, 1170, 90). Don't plot it raw; convert
to tiles or regrid.

## MATLAB

R2025b is at `/Applications/MATLAB_R2025b.app/bin/matlab` (not on PATH). Run headless
as `matlab -nodisplay -batch "run('script.m')"`, one process per script. Headless render
loops crash intermittently, so retry and check the output file exists. `rdmds.m` lives in
`MITgcm/utils/matlab`.

## Analysis conventions the user wants

- Drop the partial last-period averages that `dumpAtLast` produces before computing means.
- Check budgets close (tendency = advection + diffusion + forcing) before trusting a signal.
- Compare tiling/restart invariance on fields, not on `%MON` lines (global-sum order differs).
- Figures: every panel with more than one element gets a legend; minimal titles; real
  geometry and coastlines, never invented ones; prefer physics-native data over
  interpolated stand-ins. Make a movie (matplotlib animation + ffmpeg) for each test case.
- darwin v05→v06 diagnostic renames bite: `cDIC`→`C_DIC`, `cDIC_PIC`→`C_DICPIC`, `freeFe`→`freeFeLs`
  (a different quantity). v05/v06 ptracer indices diverge from index 20 onward — match by name.
- When code changes, regenerate the existing figures and decks in place rather than
  creating new copies (decks use the `dustin-talk-slides` skill).
- Large outputs on Pleiades: analyse there, copy back small results.

## PO.DAAC / NetCDF publishing

UDUNITS-parseable SI units (`mmol m-3`, not `meq` or `uM C`; `1` for dimensionless),
official CF `standard_name` or none, GCMD keywords looked up rather than copied, spelled-out
`long_name`. Details in `ECCO/notes/ian_metadata_instructions.docx` and `ECCO/CLAUDE.md`.
