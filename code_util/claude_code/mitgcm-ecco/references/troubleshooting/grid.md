# Troubleshooting: grids, bathymetry, hFac, vertical coordinates, r*
Distilled from answered mitgcm-support threads (2003-2026) and MITgcm GitHub issues/PRs. Names verified against the current-source index and `origin/master` (fetched 2026-10-05); "Era" says when a bug/message is old or fixed. Thread URLs are the month index; search the subject within.

## Input files, namelists, rebuilds

### Custom `CPP_OPTIONS.h` (or any header in `code/`) ignored; CONFIG_CHECK says `#undef NONLIN_FRSURF`
- Cause: genmake2 symlinks were stale; headers moved into `code/` after the first genmake2/`make depend`.
- Fix: `make CLEAN` (not `make clean`) then `make depend`; or `make makefile` to rerun genmake2 with the same options. Check `ls -l CPP_OPTIONS.h` -> `../code/CPP_OPTIONS.h`.
- Era: any (2021 report; Martin Losch).
- Src: mitgcm-support 2021-July 'Custom CPP_OPTIONS.h ignored' http://mailman.mitgcm.org/pipermail/mitgcm-support/2021-July/thread.html

### `S/R LOAD_GRID_SPACING: No value for delX at i = N` (or delY) after changing resolution
- Cause: SIZE.h is compile-time. Changed Nx/Ny without rebuilding (model still wants the old count), or delX list shorter than Nx. Also `delX = 200.` is NOT `100*200.` (gives one value + unset entries).
- Fix: rebuild after any SIZE.h edit (only workaround for varying domains is one big executable with masked land, wasteful). Use `dxSpacing = 200., dySpacing = ...` (units m or deg) for constant spacing, independent of sNx.
- Era: all versions; message text unchanged in `model/src/load_grid_spacing.F`.
- Src: 2017-January 'Digest Vol 163 Issue 1' (Ed Doddridge) http://mailman.mitgcm.org/pipermail/mitgcm-support/2017-January/thread.html ; 2006-October 'Spinup in a 2-D channel surface waves?' (Martin: dxSpacing) http://mailman.mitgcm.org/pipermail/mitgcm-support/2006-October/thread.html ; 2011-November 'frequent grid modifications' http://mailman.mitgcm.org/pipermail/mitgcm-support/2011-November/thread.html

### `S/R LOAD_GRID_SPACING: delX(i=k)= 0.0 : MUST BE >0` for a block of i
- Cause: binary file written in a different precision than `readBinaryPrec` (float32 file, 64-bit read -> zeros/garbage; "found 30 invalid delX" for half the points).
- Fix: write `>f8` and keep `readBinaryPrec=64`, or write `>f4` and set `readBinaryPrec=32`. Big-endian always. Verify with `cksum` against a clean verification copy.
- Era: 2021; message still in `load_grid_spacing.F`.
- Src: 2021-April 'internal wave verification: problem with output for delX' http://mailman.mitgcm.org/pipermail/mitgcm-support/2021-April/thread.html

### `S/R INI_VERTICAL_GRID: neither delR nor delRc are defined` / `S/R INI_PARMS: No value for delZ/delP/delR at k = N`
- Cause: almost never delR itself; an earlier namelist error corrupts parsing (too many delX values, stray space/typo/missing comma in the delZ list). Users also "fixed" it by renaming delZ<->delR or writing `60*3.3333` explicitly.
- Fix: check the &PARM04 list for typos; print delX/delY/delR right after the read in `ini_parms.F`; use dxSpacing/dySpacing; count entries == Nr in SIZE.h.
- Era: 2015-2023 (an empty space in delZ in 2023). Messages still in `ini_vertical_grid.F`, `ini_parms.F`.
- Src: 2015-October 'Problem with vertical grid.' http://mailman.mitgcm.org/pipermail/mitgcm-support/2015-October/thread.html ; 2016-January 'error about deltR' http://mailman.mitgcm.org/pipermail/mitgcm-support/2016-January/thread.html ; 2023-June 'How to customize vertical grid' http://mailman.mitgcm.org/pipermail/mitgcm-support/2023-June/thread.html

### `namelist read: read unexpected character` / `invalid reference to variable in NAMELIST input` / `Cannot match namelist object name`
- Cause: parameter in the wrong namelist (e.g. hFacMin goes in PARM01, not PARM04), misspelled, `useMDSIO` (not a namelist variable), `debugMode` in `data` (moved to `eedata`), or a `#` comment not in column 1.
- Fix: look the parameter up in `model/src/ini_parms.F`; put `#` in column 1; `debugMode` -> eedata. MDS output is the default when `useMNC=.FALSE.`; `outputTypesInclusive=.TRUE.` for both.
- Era: any; `debugMode` error since ~c66.
- Src: 2006-June 'How to comment out the lines "implicSurfPress"?' http://mailman.mitgcm.org/pipermail/mitgcm-support/2006-June/thread.html ; 2005-September 'exp0 with non-constant bathymetry' http://mailman.mitgcm.org/pipermail/mitgcm-support/2005-September/thread.html ; 2018-October 'U, V, W, T, S, and Eta files on the curvilinear coordinates' http://mailman.mitgcm.org/pipermail/mitgcm-support/2018-October/thread.html ; 2025-June '(no subject)' http://mailman.mitgcm.org/pipermail/mitgcm-support/2025-June/thread.html

### Input files opened read-write; fails on read-only filesystem
- Cause: older code opened input files with read/write access.
- Fix: since PR #596 (2022-04-27) inputs open with `ACTION='read'` via the `_READONLY_ACTION` macro (`eesupp/inc/CPP_EEMACROS.h`); defining `EXCLUDE_OPEN_ACTION` restores read-write. Namelists: PR #621.
- Era: fixed upstream in PR #596 / #621 (c68-era).
- Src: 2023-April 'Open all input files in Read only mode' http://mailman.mitgcm.org/pipermail/mitgcm-support/2023-April/thread.html ; https://github.com/MITgcm/MITgcm/pull/596

### Immediate segfault / `Killed` at startup of a large setup
- Cause: static memory footprint exceeds RAM, or stack limit; changing nSx/nPx alone does not reduce memory per process.
- Fix: `size mitgcmuv` (or `size -A mitgcmuv`) to see per-process static memory; `ulimit -s unlimited` (`limit stacksize unlimited`); raise nPx/nPy if MPI; rule of thumb ~50 double-precision 3-D fields.
- Era: 2007-2008, still valid.
- Src: 2007-March 'Killed' http://mailman.mitgcm.org/pipermail/mitgcm-support/2007-March/thread.html ; 2008-October 'Memory estimates' http://mailman.mitgcm.org/pipermail/mitgcm-support/2008-October/thread.html

## Bathymetry, masks, binary I/O

### Depth.data all zeros / no ocean; or model depth = input depth + constant
- Cause: bathyFile must be negative in ocean (positive/zero = land); or a leftover `Ro_SeaLevel` (retired) offset the depths; or NaNs in the file.
- Fix: write negative depths; delete `Ro_SeaLevel` from `data` (use `top_Pres` for p-coord or `seaLev_Z`); replace NaN by 0.
- Era: 2008-2016; Ro_SeaLevel is retired (ini_parms.F comment: "replaced by top_Pres or seaLev_Z").
- Src: 2016-April 'problem reading binary grid file ...' http://mailman.mitgcm.org/pipermail/mitgcm-support/2016-April/thread.html ; 2011-November 'curvilinear grid' (Ro_SeaLevel=1.E5) http://mailman.mitgcm.org/pipermail/mitgcm-support/2011-November/thread.html ; 2017-January 'bathymetry file in a .bin format' http://mailman.mitgcm.org/pipermail/mitgcm-support/2017-January/thread.html

### Segfault, `forrtl: error (65): floating invalid` in `ini_masks_etc`, or run dies at once, after editing bathymetry / T,S init
- Cause: NaN in bathyFile (land coded as NaN) or a NaN in a wet cell of hydrogThetaFile/hydrogSaltFile. Land-point values in init files are ignored (masked), but zeros in wet cells make INI_THETA stop.
- Fix: land = 0 in bathy; fill T/S with valid numbers (extrapolate into land) so no NaN/0 where hFacC>0.
- Era: 2008-2018, still valid.
- Src: 2014-February 'Segmentation Fault' http://mailman.mitgcm.org/pipermail/mitgcm-support/2014-February/thread.html ; 2017-January 'Error while creating bathymery : ini_masks_etc' http://mailman.mitgcm.org/pipermail/mitgcm-support/2017-January/thread.html ; 2008-October 'mask_values' http://mailman.mitgcm.org/pipermail/mitgcm-support/2008-October/thread.html ; 2018-August '(no subject)' http://mailman.mitgcm.org/pipermail/mitgcm-support/2018-August/thread.html

### Bathymetry/forcing map appears transposed/rotated, wrong number of points
- Cause: array order. MITgcm files are Fortran order, x fastest, big-endian, no record markers. Matlab `(nx,ny)` + `fwrite(...,'ieee-be')` is already right (transpose only for plotting); numpy needs `(ny,nx)` C-order or `order='F'`, `astype('>f8')`.
- Fix: check `size = Nx*Ny*prec`; `np.fromfile(f,'>f4').reshape(ny,nx)`; compare with a `verification/*/input/gendata.{m,py}` output via `cksum`.
- Era: any.
- Src: 2017-March 'Bathymetry file' http://mailman.mitgcm.org/pipermail/mitgcm-support/2017-March/thread.html ; 2018-August 'Questions Setting Up Bathymetry File' http://mailman.mitgcm.org/pipermail/mitgcm-support/2018-August/thread.html ; 2012-August 'topography' http://mailman.mitgcm.org/pipermail/mitgcm-support/2012-August/thread.html ; 2021-April (above)

### Initial/forcing fields interpolated to a finer land-sea mask have zero/garbage at coast
- Cause: coarse-grid land points hold fill values; interpolation smears them into new ocean points. Also OBC files (see obcs.md) with land = 0.
- Fix: extrapolate fields into land on the coarse grid first (Dimitris: `xpolate.m`, nearest neighbour), then interpolate and apply the new mask. Superimposing several bathymetries = pre-processing only.
- Era: any.
- Src: 2013-November 'regridding question' http://mailman.mitgcm.org/pipermail/mitgcm-support/2013-November/thread.html ; 2012-August 'Superimposing topography files' http://mailman.mitgcm.org/pipermail/mitgcm-support/2012-August/thread.html

### Single-row (2-D x-z) domain is periodic in y; need varying cross-section width
- Cause: with no walls/OBCS the domain is doubly periodic regardless of `no_slip_sides`; one y-row has no wall.
- Fix: add a second y-row with bathy=0, or use thin walls (`addSwallFile`/`addWwallFile` in PARM05, e.g. one row of 1s; `add_walls2masks.F`). dyF varying with x: hardwire dyG/dyF in `ini_cartesian_grid.F` or use a curvilinear grid (+ explicit f via `selectCoriMap=3`).
- Era: 2016-2020.
- Src: 2020-November '2D model - Varying cell width (delY)' http://mailman.mitgcm.org/pipermail/mitgcm-support/2020-November/thread.html ; 2016-May 'Default BC's at the Sides for Curvilinear Grid' http://mailman.mitgcm.org/pipermail/mitgcm-support/2016-May/thread.html

### T,S,U,V = 0 at the bottom level / "bottom" output is zero
- Cause: level k=Nr is entirely land (hFacC(:,:,Nr)=0) because bathymetry never reaches rF(Nr); W at Nr+1 does not exist.
- Fix: check hFacC/hFacW/hFacS at k=Nr in the grid output; bathymetry equal to sum(delR) is fine (solid ground below domain).
- Era: any.
- Src: 2017-December 'Zeroing quantities in the bottom' http://mailman.mitgcm.org/pipermail/mitgcm-support/2017-December/thread.html ; 2015-January 'Lev spacing and bottom depth' http://mailman.mitgcm.org/pipermail/mitgcm-support/2015-January/thread.html

## Horizontal grids: Cartesian, spherical, rotated, curvilinear

### Wind stress / forcing on a rotated or curvilinear grid points the wrong way; UVEL is not "eastward"
- Cause: model only knows +i/+j. `zonalWindFile`/ustress, uVel, vVel are along grid i/j (positive UVEL = toward increasing i).
- Fix: rotate lat/lon vectors into grid direction before writing (or use EXF, which rotates with angleCosC/angleSinC); rotate output back: `uE = AngleCS*uc - AngleSN*vc`, `vN = AngleSN*uc + AngleCS*vc` with u,v interpolated to C-points first (Martin never remembers the sign; check). AngleCS/SN are 1/0 if the grid files do not supply them (hydrostatic runs do not need them).
- Era: any.
- Src: 2012-August 'cartesian grid' http://mailman.mitgcm.org/pipermail/mitgcm-support/2012-August/thread.html ; 2011-April 'definition of UVEL' http://mailman.mitgcm.org/pipermail/mitgcm-support/2011-April/thread.html ; 2014-January 'Local rotation of co-ordinates' http://mailman.mitgcm.org/pipermail/mitgcm-support/2014-January/thread.html ; 2008-June 'Velocity cubed sphere + pickup seaice' http://mailman.mitgcm.org/pipermail/mitgcm-support/2008-June/thread.html

### `usingCurvilinearGrid=.TRUE.`: `file not found tile002.mitgrid`, every proc reads tile001, or `attempt to access non-existent record`
- Cause: new-style grid input (`horizGridFile`, `tileNNN.mitgrid`) is built for cubed-sphere topology: one file per tile/face, each field (Nx+1)x(Ny+1) (extra row/col only for the two "missing corners" of a cube face). Without pkg/exch2 (plain exch1) it still wants a file per tile; "more bug than feature" (JM, 2020).
- Fix: compile pkg/exch2 (defaults need no `data.exch2` for one facet) and supply one file with all fields; or `#define OLD_GRID_IO` in CPP_OPTIONS.h (off by default) and supply 16 individual (Nx,Ny) ieee-be files. OLD_GRID_IO has no angleCosC/angleSinC.
- Era: 2007-2020; OLD_GRID_IO still in `model/inc/CPP_OPTIONS.h` and `ini_curvilinear_grid.F`. JM planned a PR (2020) to relax the per-tile file requirement; not verified.
- Src: 2020-August 'Curvilinear grid file formats' http://mailman.mitgcm.org/pipermail/mitgcm-support/2020-August/thread.html ; 2016-June 'Internal-Wave: ... only one tile001.mitgrid' http://mailman.mitgcm.org/pipermail/mitgcm-support/2016-June/thread.html ; 2016-April 'How does MITgcm readin grid from tile001.mitgrid file' http://mailman.mitgcm.org/pipermail/mitgcm-support/2016-April/thread.html ; 2013-July 'Format of curvilinear grid files' http://mailman.mitgcm.org/pipermail/mitgcm-support/2013-July/thread.html

### OLD_GRID_IO file names; dummy/copy metrics (DXF=DXG=DXC...) give noise or blow-up
- Cause: with OLD_GRID_IO the model reads (2-D, ieee-be) `DXC DXF DXG DXV DYC DYF DYG DYU LATC LATG LONC LONG RA RAS RAW RAZ`.bin. Setting DXF=DXG=DXV, RA=RAZ=RAW=RAS etc. "to get running" breaks all conservation properties (noise near ice shelf in 2024 report; "EXTREME Pot.Temp" blow-up in 2011).
- Fix: compute the metrics properly (average the available dx/dy as `ini_spherical_polar_grid.F` / `ini_cartesian_grid.F` do; DXF/DYF/DXV/DYU are simple averages). Examples: MITgcm-contrib/arctic (cs_36km readme).
- Era: 2008-2024.
- Src: 2011-October 'curvilinear grid' http://mailman.mitgcm.org/pipermail/mitgcm-support/2011-October/thread.html ; 2024-April 'Curvilinear grid generation' http://mailman.mitgcm.org/pipermail/mitgcm-support/2024-April/thread.html ; 2016-April (above)

### `delX/delY/ygOrigin/f0` have no effect once `usingCurvilinearGrid=.TRUE.`
- Cause: all geometry comes from the grid files; spacing/origin parameters mean something only for Cartesian or spherical-polar grids; `f0` only for Cartesian.
- Fix: put the geometry in the grid files, or avoid the problem: rotate a spherical-polar grid with `rotateGrid=.TRUE., phiEuler, thetaEuler, psiEuler` (`rotate_spherical_polar_grid.F`, called from `ini_spherical_polar_grid.F`). Check spelling: a misspelled Euler-angle name silently leaves the grid unrotated.
- Era: 2017-2026.
- Src: 2018-October 'to give delX=360*0.25 ... in curvilinear coordinates' http://mailman.mitgcm.org/pipermail/mitgcm-support/2018-October/thread.html ; 2026-February 'Resources for generating curvilinear grid for MITgcm (displaced-pole Arctic)' http://mailman.mitgcm.org/pipermail/mitgcm-support/2026-February/thread.html ; 2017-October 'Invalid rotateGrid setting' http://mailman.mitgcm.org/pipermail/mitgcm-support/2017-October/thread.html

### Where do curvilinear/cs/llc grid files come from? (SPGrid dead; mitgrid with `rAz=0`, `dxC~22 m`)
- Cause: MITgcm itself never generates grid files. Hand-built mitgrids frequently have bad metrics (zero rAz, constant dyC) -> immediate cg2d divergence.
- Fix: use existing sets: ASTE (github crios-ut/aste), MITgcm-contrib/arctic, SASSIE-ECCO, cs-derived Arctic faces (Menemenlis); cs grids `MITgcm_contrib/dyncore_ASP/csgrids` (incl. ref_96), JM's cs96 `https://stuff.mit.edu/~jm_c/bin_files/cs96_dxC3_dXYa.tar.gz`; generator `MITgcm_contrib/high_res_cube/matlab-grid-generator` (gengrids.m/convertMITgrid.m: fiddly; per-face vs combined file mismatch); llc grids via ECCO/llc_hires. exch2 topology = SIZE.h + `data.exch2` only (matlab topology driver.m is obsolete).
- Era: 2007-2026; SPGrid (WildMagic/boost) unmaintained since ~2010.
- Src: 2026-February (above); 2025-February 'cubed sphere question - higher resolution grid files?' http://mailman.mitgcm.org/pipermail/mitgcm-support/2025-February/thread.html ; 2015-October 'exch2 and autogenerated grids (llc, cs, etc.)' http://mailman.mitgcm.org/pipermail/mitgcm-support/2015-October/thread.html ; 2012-February 'CS Grid Generator and Convertor Matlab Scripts' http://mailman.mitgcm.org/pipermail/mitgcm-support/2012-February/thread.html ; 2018-August 'creating grids for curvilinear orthogonal coordinates' http://mailman.mitgcm.org/pipermail/mitgcm-support/2018-August/thread.html

### AngleCS/AngleSN (angleCosC/angleSinC) wrong, jump to 1/0, in last row/column (i=Nx, j=Ny)
- Cause: `CALC_GRID_ANGLES` needs yG, dxG, dyG at i+1/j+1; only periodic/cs/llc domains fill them via exchange. Non-periodic regional grids get bad values at Nx, Ny. The angles are used by EXF_SET_UV, GMREDI, NH Coriolis, rotate_uv2en.
- Fix: ignore the last row/col (should be wall or open boundary anyway); or compute angles offline (`utils/matlab/cs_grid/cubeCalcAngle.m`) and supply in grid files (not OLD_GRID_IO). Martin's proposed `yG(i+1,j)=2*yG(i,j)-yG(i-1,j)` extrapolation after `EXCH_Z_3D_RS(yG,...)` in `ini_curvilinear_grid.F` / `ini_spherical_polar_grid.F` is not in `origin/master` (grep).
- Era: 2019-2024 (c65z, c67x); issue #828 closed Aug 2024 without a code fix; see also issue #362.
- Src: https://github.com/MITgcm/MITgcm/issues/828 ; 2019-October 'some problems about curvilinear grid' http://mailman.mitgcm.org/pipermail/mitgcm-support/2019-October/thread.html

### Locations of U/V points on curvilinear grid; no XU/YU/XV/YV output
- Cause: XG,YG are the south-west corner; U sits mid-way between (XG(i,j),YG(i,j)) and (XG(i,j+1),YG(i,j+1)); the write_grid lines for XU/XV are commented out.
- Fix: `XU(i,j)=0.5*(XG(i,j)+XG(i,j+1))`, `XV(i,j)=0.5*(XG(i,j)+XG(i+1,j))` (same for Y); T/S sit at XC,YC; W at XC,YC,RF(1:Nr). For divergences use model dx/dy, never rebuilt from coordinates.
- Era: any.
- Src: 2012-November 'XU, YU and XV, YV coordinates in the curvilinear grid...' http://mailman.mitgcm.org/pipermail/mitgcm-support/2012-November/thread.html ; 2009-May 'Location of U/V in grid' http://mailman.mitgcm.org/pipermail/mitgcm-support/2009-May/thread.html ; 2007-February 'vertical grid points and velocities' http://mailman.mitgcm.org/pipermail/mitgcm-support/2007-February/thread.html

### MNC grid output: last row/column of XG/YG equals the first; gluemnc puts curvilinear tiles side by side
- Cause: mnc pads u/v-point variables with the neighbour's first row/col (periodic); for non-exch2 curvilinear grids the X/Y coordinate variables are lat/lon not indices, so gluemnc cannot place tiles. MDS output has no padding.
- Fix: ignore the extra column; use MDS + rdmds/xmitgcm; Martin's patch used grid indices for X/Y when `usingCurvilinearGrid` in `pkg/mnc/mnc_cw_cvars.F`.
- Era: 2008-2013; mnc is legacy.
- Src: 2013-February 'MNC grid boundaries with curvilinear coordinates' http://mailman.mitgcm.org/pipermail/mitgcm-support/2013-February/thread.html ; 2008-January 'netcdf with curvilinear grid' http://mailman.mitgcm.org/pipermail/mitgcm-support/2008-January/thread.html

## Cubed-sphere and LLC specifics

### `forrtl: severe (36): attempt to access non-existent record` reading `grid_csNN.faceNNN.bin`
- Cause: file too short for what `ini_curvilinear_grid.F` reads: 18 records of (N+1)x(N+1) real*8, in order xC yC dxF dyF rA xG yG dxV dyU rAz dxC dyC rAw rAs dxG dyG angleCosC angleSinC. Grid sets lacking the two angle fields (the 16-record ones) are too short.
- Fix: use sets that include angles (dyncore_ASP/csgrids, JM's cs96 tar); `rA` is the 5th record if you need `RA.bin`; extra row/col are the 2 "missing corners" of each face.
- Era: 2013-2020; current code still reads 17,18.
- Src: 2020-February 'Problem with CS96 model set-up' http://mailman.mitgcm.org/pipermail/mitgcm-support/2020-February/thread.html ; 2016-November 'CS64-128 grid files' http://mailman.mitgcm.org/pipermail/mitgcm-support/2016-November/thread.html ; 2013-July (above)

### `EXCH1_RS_CUBE: Wrong Tiling ... works only with sNx=sNy & nSx=6 & nSy=nPx=nPy=1`
- Cause: `useCubedSphereExchange=.TRUE.` left in `eedata` for a non-cube (e.g. curvilinear regional) domain.
- Fix: remove it from eedata. Otherwise the only cube-specific part is in `ini_curvilinear_grid.F` (guarded by useCubedSphereExchange).
- Era: 2004-2011; routine is now templated, message reads `EXCH1_RX_CUBE: Wrong Tiling`.
- Src: 2011-October 'curvilinear grid' http://mailman.mitgcm.org/pipermail/mitgcm-support/2011-October/thread.html ; 2004-April 'Curvilinear grid' http://mailman.mitgcm.org/pipermail/mitgcm-support/2004-April/thread.html

### Blank tiles / `data.exch2` / "cpus on land points"
- Cause: confusion over when pkg/exch2 + data.exch2 are needed.
- Fix: data.exch2 only needed for blank tiles or non-regular facet sizes; two examples: `adjustment.cs-32x32x1` (curvilinear cs32) and `global_ocean.90x40x15` (lat-lon), compare `SIZE.h` vs `code/SIZE.h_mpi` and `input/data.exch2.mpi`; obsolete W2_EXCH2_SIZE.h / matlab driver.m must not be linked.
- Era: 2012-2023.
- Src: https://github.com/MITgcm/MITgcm/issues/791 (JM comment 2023-12-20); 2012-May 'cubed sphere: init format, subdomains and tiles' http://mailman.mitgcm.org/pipermail/mitgcm-support/2012-May/thread.html

### Cubed-sphere velocity / vector output looks broken across tile edges; initializing winds
- Cause: u,v are on the C-grid along face axes; quiver/plot of raw U,V is wrong; initial uVel/vVel from lat-lon winds need rotation + C-grid placement.
- Fix: interpolate to C-points, rotate with AngleCS/SN; use `utils/matlab/cs_grid` (`uvLatLon2cube.m`, `rotate_uv2uvEN.m`, `cubeCalcAngle.m`), MITgcmutils (`rdmds`/`wrmds`, `cs` plotting), xmitgcm cs support (PR #98) + xgcm, MeshArrays.jl/MITgcm.jl. Non-divergent init winds: evaluate a streamfunction at XG,YG (see `verification/advect_cs/code/ini_vel.F`).
- Era: 2004-2024.
- Src: 2024-October 'Processing and visualising cubed sphere output' http://mailman.mitgcm.org/pipermail/mitgcm-support/2024-October/thread.html ; 2012-May 'Initial Wind on Cubed Sphere' http://mailman.mitgcm.org/pipermail/mitgcm-support/2012-May/thread.html

### Cubed-sphere/LLC unified `.meta/.data` per time step; tiled output merging
- Cause: default MDS output is one file per tile.
- Fix: `useSingleCpuIO=.TRUE.` (faster, more robust than `globalFiles`, if the global field fits on a node); `globalFiles` not safe in multi-processor runs; to merge existing tiles read with `rdmds` and write with `wrmds` (MITgcmutils/mds.py). `useSingleCpuIO=.FALSE.` is still needed on machines where per-node local disk or memory dictates (xmitgcm assumes single-file).
- Era: 2016-2025.
- Src: 2025-July 'Question about generating unified .meta and .data files in Cubed-Sphere simulations' http://mailman.mitgcm.org/pipermail/mitgcm-support/2025-July/thread.html ; 2016-October 'xmitgcm: python package for reading mitgcm mds files' http://mailman.mitgcm.org/pipermail/mitgcm-support/2016-October/thread.html

### MOC / overturning on cs or llc grid gives nonsense; broken-line file mismatch
- Cause: (1) interpolating u,v to a lat-lon grid breaks non-divergence (large artificial error); (2) broken-line file made by the new `utils/matlab/cs_grid/mk_isoLat_bkl.m` (PR #729, fixes PR #792, 2023-24) has a different storage convention than the old `cs_grid/bk_line` files.
- Fix: integrate transports on the model grid along a broken line near each latitude; with new-format files follow `use_isoLat_bkl.m` (reshape at lines 92-96); old files keep the old recipe.
- Era: 2023-2025.
- Src: 2023-August 'AMOC estimation in cs grid' http://mailman.mitgcm.org/pipermail/mitgcm-support/2023-August/thread.html ; 2025-May 'MITgcmutil - cubed sphere - creating a broken line file' http://mailman.mitgcm.org/pipermail/mitgcm-support/2025-May/thread.html ; https://github.com/MITgcm/MITgcm/pull/729

### Non-Earth planet on cubed sphere (grid files assume Earth radius)
- Cause: shipped cs grid files were generated for rSphere=6370e3 m; `rSphere`/`omega` in `data` do not rescale curvilinear metrics.
- Fix: scale dx*/dy* by (r/rEarth) and rA* by (r/rEarth)^2 in the grid files, or use the run-time `radius_fromHorizGrid` (`ini_curvilinear_grid.F` rescales when it differs from rSphere; added ~2011).
- Era: 2008-2012; radius_fromHorizGrid present today.
- Src: 2008-June 'Velocity cubed sphere + pickup seaice' (Chris Hill, Martin) http://mailman.mitgcm.org/pipermail/mitgcm-support/2008-June/thread.html ; 2012-April 'Best practice for cubed-sphere grid generation?' http://mailman.mitgcm.org/pipermail/mitgcm-support/2012-April/thread.html

### `deepAtmosphere`: where is `rSphere`? p* + deepAtmosphere?
- Cause: scaling is horizontally uniform, a function of k only.
- Fix: Z-coord: distance from centre = rSphere + rF/rC, rF(1)=seaLev_Z (default 0) so rSphere = model top. Ocean in P-coord: also top. Atmosphere in P-coord: rSphere = bottom. Works the same in p* as p (approximation); little tested; gravity variation with height is NOT part of it.
- Era: 2013 / 2023 (Jean-Michel Campin).
- Src: 2023-February 'rSphere definition when using deepAtmosphere option for ocean modelling' http://mailman.mitgcm.org/pipermail/mitgcm-support/2023-February/thread.html ; 2013-January 'DeepAtmosphere and rstar coordinates' http://mailman.mitgcm.org/pipermail/mitgcm-support/2013-January/thread.html

## Vertical grid, partial cells (hFac), cavities

### Blow-up (`SOLUTION IS HEADING OUT OF BOUNDS`, `EXTREME Pot.Temp`, NaN, cg2d) after enabling partial cells or refining dz
- Cause: partial cells thinner than hFacMin*delR break the vertical CFL; abrupt dz jumps; aspect ratio. Min cell = `MAX(hFacMin, MIN(hFacMinDr*recip_drF(k),1))` (ini_masks_etc.F).
- Fix: raise hFacMin (0.9 almost disables PC; 0.1-0.3 typical), use hFacMinDr (metres) to bound thin cells; keep dz(k+1)/dz(k) < 1.4; `monitorFreq` = deltaT and look at `advcfl_W_hf_max` (includes hFac) not only `advcfl_wvel_max`; run with `debugLevel`>=1; temporarily set hFacMin=1 to test; use `viscAhGrid` not `viscAh`.
- Era: 2004-2023; Naughten 2018 (hFacMinDr=5 + GM, flux limiters) never resolved on-list (settled on hFacMinDr=20).
- Src: 2012-September 'Problem using partial cells' http://mailman.mitgcm.org/pipermail/mitgcm-support/2012-September/thread.html ; 2018-December 'Tracer instabilities when I reduce hFacMinDr' http://mailman.mitgcm.org/pipermail/mitgcm-support/2018-December/thread.html ; 2023-November 'Calculation overflow problem' http://mailman.mitgcm.org/pipermail/mitgcm-support/2023-November/thread.html ; 2004-January 'consult' http://mailman.mitgcm.org/pipermail/mitgcm-support/2004-January/thread.html

### High-resolution horizontal stripes / grid noise; how big is viscAhGrid?
- Cause: viscAh too small for resolution; `viscAh_eff = viscAhGrid*0.25*L^2/deltaT` (`mom_calc_visc.F`); e.g. L=3.7 km, deltaT=90 s: viscAh=0.7 needs viscAhGrid=1.8e-5.
- Fix: viscAhGrid~0.01-0.001, `viscAhGridMax` just below 1 (numerical stability bound; contributions add), optional `viscA4Grid` (~0.01); prefer 2-D Leith to Smagorinsky; dz ratio < 1.4; no_slip_bottom plus bottom drag double-counts drag.
- Era: 2020-2023.
- Src: 2023-June 'Query regarding horizontal strips in currents' http://mailman.mitgcm.org/pipermail/mitgcm-support/2023-June/thread.html ; 2020-April 'Regional high-res model configuration' http://mailman.mitgcm.org/pipermail/mitgcm-support/2020-April/thread.html

### `OBCS_CHECK: Inside Mask and OB locations disagree` with pkg/shelfice (ice draft meets bathymetry at hFacMin*delR)
- Cause: bug in `ini_masks_etc.F`: bathymetry exactly at rF(k)-drF(k)*hFacMin/2 with a deeper ice draft gave nonzero hFacW/S outside the wet domain.
- Fix: upgrade (fixed in PR #312, merged Feb 2020); workaround for old code: hFacMin 0.11, or comment out the duplicated lines 84-85, 94-95 of ini_masks_etc.F; `useMin4hFacEdges=.TRUE.` changes the symptom.
- Era: before checkpoint67-ish; fixed upstream by PR #312 (2020-02).
- Src: https://github.com/MITgcm/MITgcm/issues/324 ; 2019-December 'OBCS_CHECK' http://mailman.mitgcm.org/pipermail/mitgcm-support/2019-December/thread.html

### Thin ice shelf (hFacC(k=1)>0): thsice/seaice forms under shelf, short-wave warms cavity
- Cause: shelf thinner than the top level.
- Fix: PR #187 (Dec 2018) prevents pkg/thsice ice growth where a shelf exists and removes short-wave below shelf; pkg/seaice still needed separate fix (issue #99). Static-ice setups usually "dig" the seabed (not the draft) so each wet column keeps >=2 open cells at velocity points (Paul Holland, Kaitlin Naughten code).
- Era: fixed upstream in PR #187 (2018-12) for thsice.
- Src: https://github.com/MITgcm/MITgcm/pull/187 ; 2020-November 'Treatment of Antarctic grounding zone regions' http://mailman.mitgcm.org/pipermail/mitgcm-support/2020-November/thread.html

### Isolated coastal cells get salinity -200 ... +90 (surface fresh-water/EmPmR forcing)
- Cause: cell connected to the rest by a single face, linear free surface applies virtual salt flux with no lateral flux to balance it; partial cells only shrink the volume.
- Fix: `nonlinFreeSurf`+`useRealFreshWaterFlux` (volume changes instead of salt), remove the isolated points from bathymetry, or open a connection.
- Era: 2019 (Martin Losch).
- Src: 2019-November 'Bulk forcing and partial cells' http://mailman.mitgcm.org/pipermail/mitgcm-support/2019-November/thread.html

### Coastal currents veer / flow follows shelf wrongly in regional setup
- Cause: hFacMin=0.9 (partial cells essentially off) on 5 m dz; no OBC inflow; `OBCSfixTopo=.FALSE.`; high aspect ratio.
- Fix: hFacMin=0.3 or 0.1, `OBCSfixTopo=.TRUE.`, prescribe U,V,T,S at open boundaries (see obcs.md), grid aspect ratio ~1, no stress-forcing misalignment.
- Era: 2020.
- Src: 2020-July 'Ocean currents veering due to bathymetry' http://mailman.mitgcm.org/pipermail/mitgcm-support/2020-July/thread.html

## r*, nonlinear free surface, SSH

### `STOP in CALC_R_STAR : too SMALL rStarFac[C,W,S] !` (`WARNING: r*FacC < hFacInf`)
- Cause: eta driven to ~ -(water depth): ice loading with `useRealFreshWaterFlux` (200 m thick ice, 1-cell fjord where ice cannot flow away), tidal/forcing blow-up, near-dry shallow cells, or plain explosion/CFL. The i,j printed in `fail at i,j=` is the tile-local index; the "bi,bj,Thid" line only counts.
- Fix: find the STDERR.NNNN holding the message; `monitorFreq=1`, `grep _eta_min STDOUT.0000`; fill 1-cell bays / raise minimum depth; smooth dz; reduce deltaT; hFacInf/hFacSup defaults 0.2/2.0; do not rely on capping ice thickness: `MAX_HEFF` is now a retired parameter (error if set in data.seaice).
- Era: 2004-2020; message still in `calc_r_star.F`.
- Src: 2016-June 'too SMALL rStarFac' http://mailman.mitgcm.org/pipermail/mitgcm-support/2016-June/thread.html ; 2015-March 'small rStarFac on land....' http://mailman.mitgcm.org/pipermail/mitgcm-support/2015-March/thread.html ; 2019-March 'Simulation Breaking with r*' http://mailman.mitgcm.org/pipermail/mitgcm-support/2019-March/thread.html ; 2008-July 'ice layer modeling' http://mailman.mitgcm.org/pipermail/mitgcm-support/2008-July/thread.html

### Which free-surface settings are valid? `select_rStar`, `nonlinFreeSurf`, `DISABLE_RSTAR_CODE`
- Cause: r* is only implemented with the full nonlinear free surface; users want to know the CPP switches.
- Fix: `nonlinFreeSurf=4`, `select_rStar=0` (z) or `2` (z*); `config_check.F` enforces `select_rStar>=1 -> nonlinFreeSurf>0`, `=2 -> nonlinFreeSurf=4`, and `exactConserv=.TRUE.`. `DISABLE_RSTAR_CODE`/`DISABLE_SIGMA_CODE` (commented lines in CPP_OPTIONS.h) are adjoint-only: for TAF/ECCO set `#define NONLIN_FRSURF`, `#undef DISABLE_RSTAR_CODE`, `#define DISABLE_SIGMA_CODE` (removes many recomputations). Without r*, `SURF_ADJUSTMENT` in STDOUT means water added where hFac<hFacInf.
- Era: 2012-2020 (Martin, Patrick Heimbach).
- Src: 2020-May 'SSH problems with nonlinFreeSurf and Star' http://mailman.mitgcm.org/pipermail/mitgcm-support/2020-May/thread.html ; 2014-May 'online calculation of energy flux (Remi Tailleux)' (Patrick/Gael on NLFS+TAF) http://mailman.mitgcm.org/pipermail/mitgcm-support/2014-February/thread.html ; 2012-September 'Non-linear free-surface and vertical resolution' http://mailman.mitgcm.org/pipermail/mitgcm-support/2012-September/thread.html

### Negative / depressed SSH where sea ice is thick; SSH correlates with sea-ice thickness
- Cause: `useRealFreshWaterFlux=.TRUE.` makes ice load the surface (Campin et al 2008 ice loading), not "levitating".
- Fix: expected physics; set `useRealFreshWaterFlux=.FALSE.` for levitating ice (ice thickness results similar in Kara Sea test); Europa-type thick ice needs eta >> bottom to be avoided.
- Era: 2008-2020.
- Src: 2020-May (above); 2008-July 'ice layer modeling' (above)

### Mean SSH drifts continuously (regional model)
- Cause: unbalanced OBC net transport or imbalanced EmPmR (evaporation without return).
- Fix: `balanceEmPmR=.TRUE.` (and `balanceSaltClimRelax`) in PARM01; OBC: net inflow adjusted when `useOBCSbalance` with OBCS_balanceFac*=1 (default; this also suppresses SSH response to a prescribed barotropic forcing); alternatively remove regional mean ETAN offline, or tune boundary normal velocity to the true freshwater input (not with balance on).
- Era: 2012-2014.
- Src: 2014-August 'Continuous decrease of sea surface height in model results' http://mailman.mitgcm.org/pipermail/mitgcm-support/2014-August/thread.html ; 2012-October 'Modeling of the surface waves' http://mailman.mitgcm.org/pipermail/mitgcm-support/2012-October/thread.html

### No wetting/drying; SSH lower than the top cell or the whole depth
- Cause: linear free surface lets eta exceed top-level thickness (unphysical); there is no wet/dry.
- Fix: nonlinear free surface with z* (cells rescale; will not go dry but fails if scaled depth < ~10-20% of initial); NH + z* only if supported by your checkpoint (2006 answer: r* needs nonlinFreeSurf).
- Era: 2006-2018.
- Src: 2018-January 'Simple problems about the output of salinity and eta' http://mailman.mitgcm.org/pipermail/mitgcm-support/2018-January/thread.html ; 2006-October 'rstar coordinates' http://mailman.mitgcm.org/pipermail/mitgcm-support/2006-October/thread.html

### r* with ice-shelf cavities / thick sea ice and JMD95Z EOS
- Cause: tilted r* surfaces produce spurious density gradient with `eosType='JMD95Z'` (depth-based pressure); r* as terrain-following cavity coordinate gives pressure-gradient error.
- Fix: `selectP_inEOS_Zc=2` (dynamic hydrostatic pressure) for JMD95Z; `MDJWF` uses it by default; r* with pkg/shelfice otherwise OK.
- Era: 2018.
- Src: 2018-October 'rStar with cavities' http://mailman.mitgcm.org/pipermail/mitgcm-support/2018-October/thread.html

### r* grad-P ("sloping term") wrong when model top is not at r=0
- Cause: 2004 formulation assumed Ro_surf=0 (z) / R_low=0 (p).
- Fix: update; PR #655 (Sep 2022) generalizes for P or Z coordinate with non-zero top.
- Era: fixed upstream in PR #655 (2022-09); older checkpoints affected only for r_top != 0 (e.g. ocean in p-coord with top_Pres, ice-shelf/atm loading setups).
- Src: https://github.com/MITgcm/MITgcm/pull/655

### Which tRef/hydrogThetaFile? Kelvin-wave test gives V != 0 at step 1
- Cause: `hydrogThetaFile` is the full initial theta; `tRef` is only the reference profile (init if no file); tRef matters only for `eosType='LINEAR'` (`rho = rhoNil*(1-tAlpha*(T-tRef)+...)`), ignored for non-linear EOS. For non-zero V at t1: ICs live at their own C-grid points (eta at centre, u at west face, v at south face) and GRID.h shows the layout; fields are Nx*Ny (never Nx+1), extra points come from periodicity.
- Fix: put full 3-D T in hydrogThetaFile; prescribe ICs at the right staggered points; uniform U over topography is divergent on a C-grid (initial W stripes), so give a non-divergent field.
- Era: 2006-2020.
- Src: 2020-May 'Vertical temperature profile tRef' http://mailman.mitgcm.org/pipermail/mitgcm-support/2020-May/thread.html ; 2017-November 'Non-Zero Meridional Velocity Kelvin Wave' http://mailman.mitgcm.org/pipermail/mitgcm-support/2017-November/thread.html ; 2009-January 'initial W caused by bathymetry?' http://mailman.mitgcm.org/pipermail/mitgcm-support/2009-January/thread.html

## Budgets, pressure, diagnostics tied to the grid

### Volume/heat/tracer budget does not close with grid-file hFacW/hFacS
- Cause: with `nonlinFreeSurf>0` hFacW/S vary in time; grid-file hFacs are eta=0. With r*, div(u) is not zero; closing requires exactConserv (eta tendency).
- Fix: use diagnostics UVELMASS/VVELMASS (include hFac) or flux diagnostics ADVx_TH/ADVy_TH/ADVr_TH + DFxE_TH etc (already x dy*drF*hFac, units degC m^3/s); heat transport across latitude with VTHMASS*dxG*drF (no ETAN) or ADVy_TH/DIFy_TH; tendency = -(1/vol)*(sum flux div - TH*div(U)). div(u)=0 exactly in finite volumes: use rA, dxG, dyG, hFac, drF not dxF/dyF.
- Era: 2006-2021.
- Src: 2017-June 'On the calculation of the volume conservation from the output of a regional model' http://mailman.mitgcm.org/pipermail/mitgcm-support/2017-June/thread.html ; 2016-April 'calculation of divu' http://mailman.mitgcm.org/pipermail/mitgcm-support/2016-April/thread.html ; 2006-May 'questions re computation of meridional heat flux' http://mailman.mitgcm.org/pipermail/mitgcm-support/2006-May/thread.html ; 2017-September 'On the exact calculation of ADVx/y/r_TH and DFrx/y/rE_TH' http://mailman.mitgcm.org/pipermail/mitgcm-support/2017-September/thread.html

### PHIHYD / PHL / bottom pressure: values huge, topography imprint, PH != pressure
- Cause: phiHyd is a potential anomaly relative to rhoConst (default 999.8): `P/rhoConst = -g*rC + PH + PNH`; `PHL`/`PHIBOT` (phiHydLow) = bottom pressure/rhoConst - g*D; contains eta and atmospheric loading; mean offset and rho'*D term imprint topography. `totPhiHyd` includes eta, `PHIHYD` does not.
- Fix: remove horizontal mean (`REMOVE_MEAN_RL` in `model/src/remove_mean.F`; its call in dynamics.F is commented out), multiply by rhoConst; for non-hydrostatic add phi_nh at bottom (approximate); `PH` is at centres, `PHL` at the bottom; with z* need care (see next).
- Era: 2004-2021 (Dimitris, Martin, Jean-Michel).
- Src: 2013-August 'difference between PH and PHL' http://mailman.mitgcm.org/pipermail/mitgcm-support/2013-August/thread.html ; 2011-September 'phihydlow without topography signal' http://mailman.mitgcm.org/pipermail/mitgcm-support/2011-September/thread.html ; 2004-August 'phiHydLow' http://mailman.mitgcm.org/pipermail/mitgcm-support/2004-August/thread.html ; 2019-December 'Bottom pressure, non hydrostatic run, in oceanic set up' http://mailman.mitgcm.org/pipermail/mitgcm-support/2019-December/thread.html ; 2021-October (PHL question, Song) http://mailman.mitgcm.org/pipermail/mitgcm-support/2021-October/thread.html

### Closing the momentum budget in LLC4320 with self-computed hydrostatic pressure (noisy stripes)
- Cause: rebuilding pressure from T,S depends on eosType, selectP_inEOS_Zc, integr_GeoPot, z* choices; not straightforward.
- Fix: save `PHIHYD` and `Um_dPhiX`, `Vm_dPhiY` (momentum tendency from pressure gradient); for z* (`select_rStar>0`) ask for the relation of each diagnostic to the momentum equation. Re-running is needed if not saved; cutouts: MITgcm-contrib/llc_hires/llc_4320/regions/Box56.
- Era: 2024.
- Src: 2024-November 'Hydrostatic pressure in LLC runs' http://mailman.mitgcm.org/pipermail/mitgcm-support/2024-November/thread.html

### GM_Psi (bolus streamfunction) weighting with hFac and deepAtmosphere
- Cause: GM_PsiX/Y are already vertically integrated from the bottom and include hFac.
- Fix: do NOT multiply by hFacS; when integrating along longitude multiply by `dxG*deepFacF(k)` (matches `gmredi_residual_flow.F`, deepFacC only for volume).
- Era: 2025 (Jean-Michel Campin).
- Src: 2025-August 'GM_Psi diagnostic when using HFacC and DeepFacC' http://mailman.mitgcm.org/pipermail/mitgcm-support/2025-August/thread.html

### pkg/layers: overturning in density coordinates in p-coordinates / wrong reference pressure; budget with LinFSConserveTr=F
- Cause: `layers_krho` (formerly layers_kref) defaults to 1 for any coordinate and gives the k-level of the reference pressure for potential density; the reference pressure must coincide with a tracer-level pressure; `LAYERS_THERMODYNAMICS` code cannot be trusted in p-coordinates (issue #825 open) and budget not closed if `linFSConserveTr=F` (PR #988).
- Fix: z-coord: e.g. `layers_krho(1)=37` for ECCO 50 levels (~1934 dbar); p-coord equivalent level 14; need `LAYERS_UFLUX`, `LAYERS_VFLUX`, `LAYERS_THICKNESS` (default) for MOC; test with `global_ocean.cs32x15/input.in_p`; diagnostics names depend on layers_name order (LTto1RHO...).
- Era: 2024-2026; open issues #825, #987.
- Src: https://github.com/MITgcm/MITgcm/issues/825 ; https://github.com/MITgcm/MITgcm/issues/987

### SEAICE LSR solver: ice does not move at fine grid spacing (<~500 m), depends on tile edges
- Cause: poor LSR convergence (restricted additive Schwarz per tile); free-drift floe touching a tile edge is held back.
- Fix: keep LSR with `#define SEAICE_ALLOW_FREEDRIFT` and `LSR_mixIniGuess=2`, or use mEVP/aEVP (not original EVP); avoid JFNK/Krylov (cost).
- Era: 2020 (Martin Losch).
- Src: 2020-January 'Ice not Moving when Resolution is Very High' http://mailman.mitgcm.org/pipermail/mitgcm-support/2020-January/thread.html

### Tides on a prescribed grid: forcing via ATMOSPHERIC_LOADING instead of editing external_forcing.F
- Cause: users add tidal potential gradients in external_forcing.F and blow r*.
- Fix: `#define ATMOSPHERIC_LOADING` (CPP_OPTIONS.h) and give the potential via `pLoadFile`/`apressurefile` (EXF), only bi-linear interpolation to tracer points; use Crank-Nicolson (implicSurfPress=implicDiv2DFlow=0.5) for tides; strong shallow-water energy needs bottom drag (no_slip_bottom + bottomDragQuadratic).
- Era: 2008-2011.
- Src: 2008-April 'calc_r_star' http://mailman.mitgcm.org/pipermail/mitgcm-support/2008-April/thread.html ; 2011-November 'very large velocity near the land boundary' http://mailman.mitgcm.org/pipermail/mitgcm-support/2011-November/thread.html

### Byte-swapped input on little-endian/Pathscale: `-D_BYTESWAPIO`
- Cause: compiler default endianness vs big-endian input.
- Fix: use either `DEFINES='-D_BYTESWAPIO'` in the optfile or the compiler flag (`-byteswapio` PGI, `-convert big_endian` ifort/pathf90), never both (double swap).
- Era: 2008; the supported optfiles now set this.
- Src: 2008-October 'Error reading dx file with pathscale' http://mailman.mitgcm.org/pipermail/mitgcm-support/2008-October/thread.html
