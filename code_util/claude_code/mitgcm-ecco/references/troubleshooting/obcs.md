# Troubleshooting: pkg/obcs (open boundaries, Orlanski, Stevens, sponge, tides, balance, OB input files)
Distilled from answered mitgcm-support threads (2003-2026) and MITgcm issues/PRs. Names checked against `origin/master` (fetched 2026-10-05, includes checkpoint69o); "Era" says when a bug/message is old or fixed. Month URLs are the thread index; search the subject within. Grid/hFac/bathymetry-at-boundary items that overlap are in grid.md.

## Enabling the package and CPP/namelist mismatches

### genmake2 "In ../code/CPP_OPTIONS.h there is an illegal line: #define ALLOW_OBCS" / runtime "Run-time control flag useOBCS was used when CPP flag ALLOW_OBCS was unset"
- Cause: obcs is a package, not a CPP_OPTIONS.h switch (old tutorial text was wrong).
- Fix: put `obcs` in `code/packages.conf` (or genmake2 `-enable=obcs`), `useOBCS=.TRUE.` in `data.pkg`, and no `#define ALLOW_OBCS` anywhere; select boundaries and options in a copy of `OBCS_OPTIONS.h` in `code/`.
- Era: since 2004 (c52+); doc fixed later.
- Src: 2004-October 'problem with defining ALLOW_OBCS flag' http://mailman.mitgcm.org/pipermail/mitgcm-support/2004-October/thread.html ; 2008-May 'obcs in CPP_OPTIONS.h?' http://mailman.mitgcm.org/pipermail/mitgcm-support/2008-May/thread.html ; 2007-June 'OBCS CPP flag' http://mailman.mitgcm.org/pipermail/mitgcm-support/2007-June/thread.html

### `OBCS_CHECK: Cannot set useOBCSsponge=.TRUE. ... with #undef ALLOW_OBCS_SPONGE` (same for prescribe, tides, Stevens)
- Cause: runtime flag needs the matching CPP flag in `OBCS_OPTIONS.h`. In current master `ALLOW_OBCS_PRESCRIBE`, `ALLOW_OBCS_BALANCE`, `ALLOW_ORLANSKI`, `ALLOW_OBCS_N/S/E/W` are #define, but `ALLOW_OBCS_SPONGE`, `ALLOW_OBCS_STEVENS`, `ALLOW_OBCS_TIDES`, `ALLOW_OBCS_SEAICE_SPONGE` are `#undef`. Misspelled flag names (cpp does not complain) or a stale link of OBCS_OPTIONS.h in the build dir give the same error.
- Fix: `#define` the flag in `code/OBCS_OPTIONS.h`; `ls -l OBCS_OPTIONS.h` in build; `make makefile && make CLEAN && make depend && make`; to be sure compare `obcs_check.F` vs `obcs_check.f`. Undefining `ALLOW_OBCS_NORTH/SOUTH` does not make a wall (domain stays periodic).
- Era: 2008-2026, unchanged.
- Src: 2008-June 'Obcs problems' http://mailman.mitgcm.org/pipermail/mitgcm-support/2008-June/thread.html ; 2015-November 'Problem with ALLOW_OBCS_SPONGE.' http://mailman.mitgcm.org/pipermail/mitgcm-support/2015-November/thread.html ; 2021-February 'OBCS add_tide does not change open boundary velocities' http://mailman.mitgcm.org/pipermail/mitgcm-support/2021-February/thread.html

### `Fortran runtime error: End of file` in obcs_readparms / `Cannot match namelist object name` / `"tidalPeriod" no longer allowed in file "data.obcs"`
- Cause: data.obcs namelist malformed (missing `&` terminator, extra/misplaced namelist; OB_* given as reals like `450*-0.1`, which must be integers), or an old parameter name.
- Fix: terminate each namelist; OB_I*/OB_J* are integers; `tidalPeriod` was replaced by `OBCS_tidalPeriod`; old tide files OBNamFile/OBSphFile etc. are "retired-and-unset" -> `OB[N,S,E,W]_[u,v]Tid[Am,Ph]File` (checkpoint68y, 2024-06).
- Era: tide renames fixed upstream in checkpoint68y; the old names stop the run.
- Src: 2019-June 'Error when introducing OBCS-Sponge in my config' http://mailman.mitgcm.org/pipermail/mitgcm-support/2019-June/thread.html ; 2014-January 'OBCS (open boundary forcing) problem' http://mailman.mitgcm.org/pipermail/mitgcm-support/2014-January/thread.html ; 2021-May 'Open boundary conditions for Internal Waves.' http://mailman.mitgcm.org/pipermail/mitgcm-support/2021-May/thread.html

## Where the boundary is: OB_Ieast/Iwest/Jnorth/Jsouth, OBCS_CHECK, corners

### Meaning and length of OB_Ieast, OB_Iwest, OB_Jnorth, OB_Jsouth
- Cause: vectors indexed along the boundary, entries are the *tracer-cell* index of the OB along the normal: OB_Jnorth/OB_Jsouth have Nx entries (value = J index), OB_Ieast/OB_Iwest have Ny entries (value = I index). 0 = no OB; `-1` is shorthand for Nx (east) / Ny (north); values must be integers. Mixing Nx/Ny is a classic error (`OB_Jsouth=801*1` for Nx=400).
- Fix: `OB_Jsouth=21*1,104*0` (not `OB_Jsouth(1:21)=21*1`, compiler-dependent); partial OBs: `OB_Jnorth=549*0,2*-100,...`; irregular boundaries (OB_Ieast varying with j) work in newer code (changes Oct 2012). The same N-vector means the OB input files are always full-length (see file entry below).
- Era: 2003-2021 (doc example in Fig. 8.6 had an index typo, flagged 2025).
- Src: 2008-May 'obcs in CPP_OPTIONS.h?' http://mailman.mitgcm.org/pipermail/mitgcm-support/2008-May/thread.html ; 2009-June 'cal, exf, obcs: problem at the tile boundaries' (Patrick Heimbach) http://mailman.mitgcm.org/pipermail/mitgcm-support/2009-June/thread.html ; 2007-January 'obcs' (-1 == nx) http://mailman.mitgcm.org/pipermail/mitgcm-support/2007-January/thread.html ; 2016-June 'Quick Question about Setting OBCS' http://mailman.mitgcm.org/pipermail/mitgcm-support/2016-June/thread.html

### `OBCS_CHECK: Inside Mask and OB locations disagree` / `N errors in OB location vs Mask` / `N errors in tile OB set-up`
- Cause: the OB index list does not match the wet/dry mask built from the bathymetry: OB points over land or next to land, transposed bathymetry (NY*NX vs NX*NY), changed SIZE.h/tiling, OBCSfixTopo modifying depths, partial-cell effects, or CPP flags (`ALLOW_OBCS_NORTH` etc.) contradicting data.obcs. "N errors in tile OB set-up" is a CPP-vs-namelist conflict whose exact line (e.g. `tile bi,bj has Northern OB`) is in STDERR.
- Fix: read the individual messages in STDERR.*/STDOUT.* (per-tile `OB_Jn/OB_Js/OB_Ie/OB_Iw` local indices are printed after `OBCS_CHECK: start summary`); to derive the right list run one step with `useOBCS=.FALSE.` and run-length-encode `hFacW(2,:,1)`, `hFacW(Nx,:,1)`, `hFacS(:,2,1)`, `hFacS(:,Ny,1)` (1s are -1 on N/E); try one OB point first, then add; with exch2 an `insideOBmaskFile` can define the OB region; make bathymetry flat/consistent across the boundary or use OBCSfixTopo=.TRUE. Check input orientation (bathymetry is Fortran-ordered).
- Era: 2012-2021; pkg/shelfice + OBCS version of the message was a real bug fixed in PR #312 (see grid.md).
- Src: 2014-September 'OBCS errors in OB location vs Mask' http://mailman.mitgcm.org/pipermail/mitgcm-support/2014-September/thread.html ; 2020-February 'S/R OBCS_CHECK: Inside Mask and OB locations disagree' http://mailman.mitgcm.org/pipermail/mitgcm-support/2020-February/thread.html ; 2017-April 'Error Mask' (OBCSfixTopo default .TRUE.) http://mailman.mitgcm.org/pipermail/mitgcm-support/2017-April/thread.html ; 2019-December 'OBCS_CHECK' http://mailman.mitgcm.org/pipermail/mitgcm-support/2019-December/thread.html

### OB point in a domain corner or inside a bay: OBCS_CHECK fails, large w, or SSH drift with four open sides
- Cause: at a corner the code sees two boundaries (OB_Iwest=0 still counts); prescribed values at corners are applied twice (once per side) so u/v are not divergence-free and w is forced; a boundary row with land on both sides is not supported.
- Fix: put a land cell at corners (depth 0) or start OBs at 2..N-1; specify the neighbouring row explicitly, e.g. `OB_Iwest=1*0,1*2`, `OB_Jsouth=1*0,1*2` for a 1-cell step (Estanislao, Martin: logically there should be flux through that interface); make corner u/v consistent (obwu(1)=obsu(1), obwv(1)=obsv(1)).
- Era: 2005-2019.
- Src: 2019-October 'Unstructured open boundaries' http://mailman.mitgcm.org/pipermail/mitgcm-support/2019-October/thread.html ; 2005-October 'SSH drift and obcs again' http://mailman.mitgcm.org/pipermail/mitgcm-support/2005-November/thread.html ; 2014-July 'Logarithmic decrease in Eta without OBCS/EXF' (Matt: corners land) http://mailman.mitgcm.org/pipermail/mitgcm-support/2014-July/thread.html

### Cannot prescribe a river/"inner" OB inside a domain whose perimeter is all open, or a second OB on the same row
- Cause: one OB index per row/column and one file per side; OB points with water on both sides are rejected.
- Fix: use surface runoff via EXF (maybe at a cell a few rows in), or Orlanski with `insideOBmaskFile` (still cannot read OB values from files); for a second western OB at a different j you need zeros elsewhere in the single full-length file. OBCS for deep sewage/point sources: put the source in ptracers_forcing/external_forcing instead.
- Era: 2004-2025 (Jean-Michel Campin confirmed 2025-03).
- Src: 2025-March 'Rivers imposed as OBCs - problem with fully-open boundary domains' http://mailman.mitgcm.org/pipermail/mitgcm-support/2025-March/thread.html ; 2021-August 'Backtrace error while implementing NOBCS' http://mailman.mitgcm.org/pipermail/mitgcm-support/2021-August/thread.html ; 2004-December 'OBCS for deep sewage water source' http://mailman.mitgcm.org/pipermail/mitgcm-support/2004-December/thread.html

### Inclined / rotated open boundaries
- Cause: OB_I*/OB_J* are grid directions.
- Fix: don't rotate the domain if exch2 can drop land tiles; otherwise use a curvilinear grid (OLD_GRID_IO, e.g. MITgcm-contrib/arctic/cs_36km), rotate boundary velocities into grid direction with AngleCS/SN (see grid.md).
- Era: 2019-2024.
- Src: 2024-December 'How to set inclined open boundaries' http://mailman.mitgcm.org/pipermail/mitgcm-support/2024-December/thread.html ; 2019-August 'On the set up of the OBCS volume conservation in a regional simulation with curvilinear grid' http://mailman.mitgcm.org/pipermail/mitgcm-support/2019-August/thread.html

## OB input files: dimensions, records, time, precision, defaults

### `forrtl: severe (36)/ Fortran runtime error: Non-existing record number` in mdsio_read_section (file OBNu...), `MDSREADFIELD_XZ_GL: File does not exist`
- Cause: OB file shorter than what the run reads. Each OB file must be the full-length section times Nr times nrec: north/south (Nx, Nr, nrec), east/west (Ny, Nr, nrec), even for a single OB point (zeros elsewhere; 3-D vars need Nr, 2-D (eta, ice) don't); `nrec >= nTimeSteps*deltaT/externForcingPeriod + 1` (linear interpolation reads one record beyond the end); file precision != `readBinaryPrec` (real*4 file read as real*8 gives half the records); with EXF a record before and after the run window; boundary data that simply end before your run (Arctic cs_36km OB ends Jul 2002).
- Fix: rebuild files with the above shape; check byte size = Nx*Nr*nrec*prec; extend/loop the OB series (monthly climatology repeated N*12 records).
- Era: 2003-2021 all; message text from `mdsio_read_section.F`.
- Src: 2021-August 'Backtrace error while implementing NOBCS' http://mailman.mitgcm.org/pipermail/mitgcm-support/2021-August/thread.html ; 2019-June 'Error while implementing OBCS' http://mailman.mitgcm.org/pipermail/mitgcm-support/2019-June/thread.html ; 2014-June 'OBCS Error with 36km Arctic config' http://mailman.mitgcm.org/pipermail/mitgcm-support/2014-June/thread.html ; 2013-December 'format of input files wrong in regional model with open boundary conditions?' http://mailman.mitgcm.org/pipermail/mitgcm-support/2013-December/thread.html

### OB velocities/tracers appear as ±1e32, zeros, or "model never reads OB files"
- Cause: wrong precision/endianness (OB files use `readBinaryPrec`; with pkg/exf the precision is `exf_iprec_obcs`, defaults to `exf_iprec`, so `exf_iprec=32` with real*8 files gives 1e30 garbage); OB*File names missing (no uFile -> OB u,v = 0; STDOUT monitor `obc_N_vVel_*` all zero); `useOBCSprescribe` not set or `ALLOW_OBCS_PRESCRIBE` undef.
- Fix: set `exf_iprec_obcs`/`exf_iprec` to match the files; provide u,v files; `useOBCSprescribe=.TRUE.`; check the `MDS_READ_SEC_XZ: opening global file` lines in STDOUT. Issue #309 (open): no obcs-own precision flag; tidal OB files are read in obcs_init_variables with readBinaryPrec even when EXF is on.
- Era: 2005-2025; #309 open.
- Src: 2005-March 'obcs_prescribe and *.h conflict?' (exf_iprec 32->64) http://mailman.mitgcm.org/pipermail/mitgcm-support/2005-March/thread.html ; 2006-September 'OBCS settings' http://mailman.mitgcm.org/pipermail/mitgcm-support/2006-September/thread.html ; 2017-August 'Open Boundary' http://mailman.mitgcm.org/pipermail/mitgcm-support/2017-August/thread.html ; https://github.com/MITgcm/MITgcm/issues/309

### What are the OB values when you give no files? (excess cooling / relaxation to tRef at the boundary)
- Cause: defaults are T,S = tRef,sRef profiles and u,v = 0 (not zero T); a T or S far from tRef at the boundary then drives strong adjustment or cooling.
- Fix: set tRef/sRef close to the interior, or prescribe T,S files; to hard-wire values edit a copy of `obcs_calc.F` in code/ (check the values aren't overwritten later). Constant forcing without files: `periodicExternalForcing=.FALSE.` (one record) or EXF `obcs?period=0`.
- Era: 2008-2021.
- Src: 2008-November 'How to prescribe constant OBCS' http://mailman.mitgcm.org/pipermail/mitgcm-support/2008-November/thread.html ; 2021-August 'Alternative flag instead of OBXyFile in data.obcs' http://mailman.mitgcm.org/pipermail/mitgcm-support/2021-August/thread.html ; 2011-November 'Excess cooling at the surface layers' http://mailman.mitgcm.org/pipermail/mitgcm-support/2011-November/thread.html

### OB time stamps: how records map to time (no EXF vs EXF), monthly/yearly series
- Cause: without EXF OB records follow `periodicExternalForcing/externForcingPeriod/externForcingCycle` in PARM03, stamped at the *middle* of each interval and linearly interpolated (so t=0 is the mean of first and last record); with EXF the timing comes from `&EXF_NML_OBCS` (`obcs{N,S,E,W}startdate1/2`, `period`, `repCycle`), one per boundary.
- Fix: non-EXF: constant -> `periodicExternalForcing=.FALSE.`; set cycle/period so records = cycle/period. EXF: constant -> `obcs?period=0`; place data mid-interval (`obcsSstartdate1=19790116`) or `obcs?period=-12` for monthly climatology (no year files); with `useOBCSYearlyFields` the year of startdate is replaced by the file year, so use a Jan-1 startdate (a Dec-31 start froze the OB at the first record); Gregorian monthly OB: constant period 2629800 s (=0.25*(366+3*365)*86400/12); interannual monthly with cal: `period=-1` (added c68e; EXF bug at c68d `CAL_CHECKDATE: Invalid month in date(1)= 0` for undefined boundaries fixed by PR #562, Nov 2021).
- Era: 2004-2021; period=-1 and PR #562 in checkpoint68e.
- Src: 2016-February 'variables (.)_F1 in obcs_calc.F' (Jody: records at mid-interval) http://mailman.mitgcm.org/pipermail/mitgcm-support/2016-February/thread.html ; 2014-October 'OBCS_PRESCRIBE_READ' http://mailman.mitgcm.org/pipermail/mitgcm-support/2014-October/thread.html ; 2019-October 'Monthly averaged OBCS with useOBCSYearlyFields' http://mailman.mitgcm.org/pipermail/mitgcm-support/2019-October/thread.html ; 2011-April 'useOBCSYearlyFields' http://mailman.mitgcm.org/pipermail/mitgcm-support/2011-April/thread.html ; 2021-January 'issues with SST/SSS restoring and OBCS in EXF package' http://mailman.mitgcm.org/pipermail/mitgcm-support/2021-January/thread.html ; 2021-November 'OBCS/EXF in Checkpoint68d' http://mailman.mitgcm.org/pipermail/mitgcm-support/2021-November/thread.html ; https://github.com/MITgcm/MITgcm/pull/562

### Which grid point/hFac do OB values refer to? (staggering, inward point, OBCS_uvApplyFac)
- Cause: tracer OB at index OB; normal velocity is imposed at the OB and at the point *inward* (u(OB+1) for a western OB at OB, u(OB) for eastern; v same for N/S); "true" OB velocity is the inward point, the outer one is for historical reasons (`OBCS_uvApplyFac=0` gives identical results). eta (cg2d) is zero at OBs; hFacW/S on the outer face is a periodic copy.
- Fix: interpolate boundary data to the velocity point inside the tracer OB; for fluxes use `hFacW(2,:,:)` (west), `hFacW(Nx,:,:)` (east), `hFacS(:,2,:)` (south), `hFacS(:,Ny,:)` (north) with dyG/dxG at the same index; keep depth flat across the OB (Depth(1,:)=Depth(2,:)); fields in IC/RBCS/OBCS go on their own grid points (U on U, V on V).
- Era: 2005-2017 (Martin Losch, Matt Mazloff, Jean-Michel Campin).
- Src: 2005-September 'obcs' (Martin: what is really applied) http://mailman.mitgcm.org/pipermail/mitgcm-support/2005-September/thread.html ; 2016-April 'stagger grid and obcs' http://mailman.mitgcm.org/pipermail/mitgcm-support/2016-April/thread.html ; 2014-April 'how to control obcs balance flow ?' http://mailman.mitgcm.org/pipermail/mitgcm-support/2014-April/thread.html ; 2017-March 'Aligning hFacW, UVEL, dyG' http://mailman.mitgcm.org/pipermail/mitgcm-support/2017-March/thread.html

### Bad T/S/u/v at land-ocean transitions along the boundary (e.g. salinity ~0 next to coast)
- Cause: OB file values over land points are masked, not extrapolated; mismatch of OB mask and file; also a sponge (40 pts) restoring to zeros (2020 thread unresolved on-list; the thread ended with RBCS replacing the OBCS sponge).
- Fix: extrapolate OB fields into land before writing; check OB location vs masks; test with `OBCSsponge_Salt=.FALSE., OBCSsponge_Theta=.FALSE.`; make bathymetry flat through the sponge; Martin/Dimitris: toggle sponge/balance one at a time. See also the 2026 sponge fix below.
- Era: 2020; sponge bug fixed in checkpoint69o (2026-07).
- Src: 2020-June 'topography with OBCS' http://mailman.mitgcm.org/pipermail/mitgcm-support/2020-June/thread.html

### OB inputs on an irregular vertical grid; OB from ECCO/other model output
- Cause: pkg/obcs expects inputs on the model's own vertical levels; ECCO "mass-weighted" velocities include hFac.
- Fix: write OB fields on the model's delR grid; ECCO: use UVELMASS/hFacW to back out UVEL (ignore bolus velocity), and balance after interpolation.
- Era: 2017-2019.
- Src: 2019-March 'Irregular Vertical Grid and OBCS' http://mailman.mitgcm.org/pipermail/mitgcm-support/2019-March/thread.html ; 2017-April 'Open boundary conditions using ECCO data' http://mailman.mitgcm.org/pipermail/mitgcm-support/2017-April/thread.html

### OB eta files: `OB*etaFile(s) only allowed with nonlinFreeSurf`; SSH at the boundary does not drive the interior; OBNhFile
- Cause: eta OB values only set layer thickness (volume flux) with the nonlinear free surface; the pressure gradient across the OB never enters the equations, and eta at the OB point is zero in cg2d. OB[N,S,E,W]hFile is ice thickness, not SSH.
- Fix: nonlinFreeSurf=4 (+select_rStar=2), or convert SSH to a barotropic normal flow (u = c*eta/H, c=sqrt(gH), or via a streamfunction) and prescribe velocity (+Stevens); relaxing eta does not change interior SSH.
- Era: 2004-2017.
- Src: 2016-February 'only eta as open boundary ?' http://mailman.mitgcm.org/pipermail/mitgcm-support/2016-February/thread.html ; 2005-October 'eta at OBCS' http://mailman.mitgcm.org/pipermail/mitgcm-support/2005-October/thread.html ; 2017-July 'Ask for boundary conditions (OB[N/E/W/S]hFile)' http://mailman.mitgcm.org/pipermail/mitgcm-support/2017-July/thread.html ; 2016-November 'obcs with non linear free surface' http://mailman.mitgcm.org/pipermail/mitgcm-support/2016-November/thread.html

### Using the internal_wave verification code: constant OB forcing is ignored
- Cause: that experiment's `code/obcs_calc.F` has periodic forcing hard-coded.
- Fix: remove `obcs_calc.F` from your code dir and rebuild (`make makefile && make CLEAN && make depend && make`) to use file-driven OBs. Editing obcs_calc.F directly is supported (Martin: "that's what this file is for"). Prescribing W: `OB[N,S,E,W]wFile` now exist (they didn't before ~2009); w at an OB is reset by obcs_apply_w.
- Era: 2005-2020.
- Src: 2020-June 'consatant open boundary flow but interior current changing periodicaly' http://mailman.mitgcm.org/pipermail/mitgcm-support/2020-June/thread.html ; 2020-June 'Tidal input at open boundary' http://mailman.mitgcm.org/pipermail/mitgcm-support/2020-June/thread.html ; 2008-January 'Specify W velocity with OBCS?' http://mailman.mitgcm.org/pipermail/mitgcm-support/2008-January/thread.html

## Volume balance and SSH drift

### Sea level rises/falls by metres to tens of metres in a regional run with prescribed OBs (or Orlanski)
- Cause: net volume inflow != 0 (also tiny offsets, or single-precision OB files); Orlanski cannot control net flow; in addition unbalanced E-P/runoff.
- Fix: balance the boundary transports offline (preferred; use hFacW/hFacS at the inward index, dyG/dxG, drF, double precision `readBinaryPrec=64`, balance over the period you want, e.g. a year) or `useOBCSbalance=.TRUE.` with `OBCS_balanceFacN/S/E/W` (balances every time step); `exactConserv=.TRUE.`; for freshwater: `balanceEmPmR` (needs `ALLOW_BALANCE_FLUXES`, which is `#undef` in CPP_OPTIONS.h) or regional-mean removal of ETAN offline. Diagnose with OBCS monitor (`OBCS_monitorFreq`, `obc_*_Int`).
- Era: 2003-2022.
- Src: 2014-April 'how to control obcs balance flow ?' http://mailman.mitgcm.org/pipermail/mitgcm-support/2014-April/thread.html ; 2014-October 'Open boundary conditions' (Martin) http://mailman.mitgcm.org/pipermail/mitgcm-support/2014-October/thread.html ; 2008-December 'obcs balance' http://mailman.mitgcm.org/pipermail/mitgcm-support/2008-December/thread.html ; 2012-September 'Anamalous Decrease in Sea Surface Level and SPONGE Layer' http://mailman.mitgcm.org/pipermail/mitgcm-support/2012-September/thread.html ; 2022-August 'model unstability caused by increasing sea surface elevation(eta)' http://mailman.mitgcm.org/pipermail/mitgcm-support/2022-August/thread.html

### What do the OBCS_balanceFac values mean? (`-1`, `0`, `1`, `2`)
- Cause: the factor selects where the net-flow correction is applied.
- Fix: `0.` no correction on that side; `>0` relative size of the share of the correction taken by that side; `-1` forces zero net flow through *that* boundary (removes a real through-flow). Typical: 0 everywhere and 1 on one side (outflow side) or split between two. With a through-current (Kuroshio/ACC) never put -1 on inflow/outflow sides. Orlanski + balance: `OBCS_balanceFacN=1., OBCS_balanceFacS=0.`. Balance is applied before tides are added, so tidal volume oscillation of domain-mean ETA remains.
- Era: 2012-2021.
- Src: 2014-April 'how to control obcs balance flow ?' http://mailman.mitgcm.org/pipermail/mitgcm-support/2014-April/thread.html ; 2021-November '[EXTERNAL] Warm start and cold start yielding same results' http://mailman.mitgcm.org/pipermail/mitgcm-support/2021-November/thread.html ; 2019-September 'Orlanski boundary condition' http://mailman.mitgcm.org/pipermail/mitgcm-support/2019-September/thread.html ; 2019-December 'Volume conservation with open boundaries and tides' http://mailman.mitgcm.org/pipermail/mitgcm-support/2019-December/thread.html

### Balance does not correct non-OBCS volume sources (ice-sheet/shelf volume change, E-P)
- Cause: obcs_balance_flow only removes the OB-prescribed net flux; `balanceEmPmR` is global and not passed to obcs.
- Fix: no upstream mechanism (issue #145 open); idea: add the mean EmPmR/volume defect to `inFlow` in `obcs_balance_flow.F` (Dan Goldberg's branch); `OBCSbalanceSurf=.TRUE.` adds the surface mass flux into the balance.
- Era: 2018-2019 (issue #145), `OBCSbalanceSurf` present in master.
- Src: https://github.com/MITgcm/MITgcm/issues/145 ; 2019-April 'SSH drift with sponge-OBCS & OBCSbalance' http://mailman.mitgcm.org/pipermail/mitgcm-support/2019-April/thread.html

### `useOBCSbalance=.TRUE.` crash (0/0) with nSx>1 and OBs at domain edge
- Cause: tiles without boundary points had area=0 -> division by zero (compiler-dependent).
- Fix: fixed in pkg/obcs 2007 (guard Ar>0).
- Era: fixed upstream 2007-05.
- Src: 2007-April 'OBS_CALC bug: divide by zero at internal boundaries' http://mailman.mitgcm.org/pipermail/mitgcm-support/2007-April/thread.html

### dynstat_eta_max = 0 / eta exactly 0 on open boundary
- Cause: pressure/eta is not solved at OB points (cg2d_b, cg2d_x zeroed); max SSH reads 0 if the whole basin is negative.
- Fix: not a bug; look at the basin's volume budget (above).
- Era: any.
- Src: 2008-August 'problem with obcs package' http://mailman.mitgcm.org/pipermail/mitgcm-support/2008-August/thread.html ; 2009-January 'data precision' http://mailman.mitgcm.org/pipermail/mitgcm-support/2009-January/thread.html

## Orlanski radiation

### `OBCS_CHECK: useOrlanski* OBC not yet implemented for pTracers` / `nonlinFreeSurf not yet implemented in Orlanski OBC` / with seaice / `OBCS not yet implemented in CD-Scheme`
- Cause: these combinations are still hard errors in `obcs_check.F` (Orlanski+pTracers, +nonlinFreeSurf, +useSEAICE; any OBCS with `useCDscheme`).
- Fix: use Dirichlet(+sponge) or Stevens instead; for pTracers either comment out the stop in `obcs_check.F`/`obcs_calc.F` (ptracers then take a Neumann-like OB from the interior) and "keep fingers crossed", per-tracer `OBCS_u1_adv_Tr(iTr)=1` (1st-order upwind at outflow; `verification/so_box_biogeo`), or RBCS to damp tracers at the edge; for CD scheme replace by small biharmonic viscosity. Orlanski can be switched per side.
- Era: all versions up to master.
- Src: 2012-October 'Way around ptracers & Orlanski OBCS?' http://mailman.mitgcm.org/pipermail/mitgcm-support/2012-October/thread.html ; 2019-April 'Behavior of Ptracers at Orlanski boundary' http://mailman.mitgcm.org/pipermail/mitgcm-support/2019-April/thread.html ; 2012-April 'nonlinear free surface and orlanski boundary' http://mailman.mitgcm.org/pipermail/mitgcm-support/2012-April/thread.html ; 2012-January 'CD scheme/OBCS' http://mailman.mitgcm.org/pipermail/mitgcm-support/2012-January/thread.html

### Orlanski: reflection of internal waves, no net flow, SSH drift, tides ignored
- Cause: single (adaptive) phase speed per variable, barotropic mode not handled; "Sommerfeld-type" so volume is not conserved; barotropic tidal velocity added on top is not filtered out of the signal Orlanski diagnoses (issue #810).
- Fix: use a sponge (Dirichlet+sponge; wide, gentle) for internal waves/tides; Orlanski only for mode-1/outflow cases; set `CMAX` ~0.45 and `cVelTimeScale` ~ deltaT if no intrinsic timescale (2004 advice); use `useOBCSbalance` + `OBCS_balanceFac*`; sensitivity to viscosity: SSH blow-up with small deltaT cured by switching viscAhGrid to `viscC2smag`.
- Era: 2004-2024; issue #810 closed Apr 2024 with the sponge recommendation.
- Src: https://github.com/MITgcm/MITgcm/issues/810 ; 2017-November 'Orlanski vs sponge obcs in ocean modelling' http://mailman.mitgcm.org/pipermail/mitgcm-support/2017-November/thread.html ; 2022-April 'Orlanski radiation condition parameters' http://mailman.mitgcm.org/pipermail/mitgcm-support/2022-April/thread.html ; 2004-April 'FW: obcs for lock-exchange problem (fwd)' http://mailman.mitgcm.org/pipermail/mitgcm-support/2004-April/thread.html ; 2019-July 'SSH blow-up at Orlanski BC with small deltaT' http://mailman.mitgcm.org/pipermail/mitgcm-support/2019-July/thread.html

### NaN in one variable (e.g. T) along an Orlanski boundary after long runs (divide by zero)
- Cause: `CL = (phi - phi_old)/(ab1*S2+ab2*S3)` with denominator exactly 0 for uniform fields.
- Fix: upstream guards `denom .NE. 0` and sets CL=0 (PR #826/#824, issue #822); whether CL should be CMAX instead is undecided (changing it breaks `dome` and `tutorial_plume_on_slope`).
- Era: fixed upstream in checkpoint68y (2024-06).
- Src: https://github.com/MITgcm/MITgcm/issues/822 ; https://github.com/MITgcm/MITgcm/pull/826

### Orlanski restart and OB pickup/section files written per tile
- Cause: old code did not restart Orlanski properly (phase speeds), NH+Orlanski lost wVel at the OB; OBCS slice (xz/yz) files and pickups are per-tile because `useSingleCpuIO` does not apply to sections/vectors.
- Fix: current master writes/reads `pickup_orlanski{N,S,E,W}`; for one global file use `globalFiles=.TRUE.` together with `useSingleCpuIO=.TRUE.` (warning is benign); ECCO OBCS controls also need it.
- Era: Orlanski restart fixed 2008-2009; per-tile sections still current.
- Src: 2006-January 'mnc obcs pickup files?' http://mailman.mitgcm.org/pipermail/mitgcm-support/2006-January/thread.html ; 2008-February 'Non-hydrostatic Orlanski pickup and W velocity' http://mailman.mitgcm.org/pipermail/mitgcm-support/2008-February/thread.html ; 2014-March 'Orlanksi Boundary Condition Files' http://mailman.mitgcm.org/pipermail/mitgcm-support/2014-March/thread.html ; 2008-May 'is this bug in mitgcm code?' http://mailman.mitgcm.org/pipermail/mitgcm-support/2008-May/thread.html

## Sponge layer

### How to choose spongeThickness, Urelaxobcsinner/bound; sponge restores to OB values, not to "nothing"
- Cause: `obcs_sponge.F` relaxes u,v,T,S in `spongeThickness` cells inside each OB toward the OB value with a timescale varying linearly from `*relaxobcsbound` at the OB to `*relaxobcsinner` at the inner edge. bound=0 is Dirichlet (reflects at sponge edge), inner=inf has no effect; defaults 0 thickness. By default all four sides and all of U,V,T,S are sponged (switches `OBCSsponge_N/S/E/W`, `_UatNS`, `_UatEW`, `_VatNS`, `_VatEW`, `_Theta`, `_Salt`, `useLinearSponge`).
- Fix: inner long (e.g. 10 d), bound short (e.g. 12 h; for internal waves hour-scale at the edge); width of order one mode-1 wavelength (telescope dx to get that with few cells; stretch <~1-2%/cell); OB values must equal the interior state you want (tRef profile for T/S, zero for v) or the sponge itself generates motion; flat bathymetry through the sponge; zero tangential velocity; for time-propagating waves inside the sponge use pkg/rbcs.
- Era: 2006-2026.
- Src: 2017-November 'parameters in OBCS sponge layer in ocean modeling' http://mailman.mitgcm.org/pipermail/mitgcm-support/2017-November/thread.html ; 2018-March 'parameters in obcs sponge layer' http://mailman.mitgcm.org/pipermail/mitgcm-support/2018-March/thread.html ; 2012-September 'Anamalous Decrease in Sea Surface Level and SPONGE Layer' (Jody) http://mailman.mitgcm.org/pipermail/mitgcm-support/2012-September/thread.html ; 2012-February 'instability issue with sponge' http://mailman.mitgcm.org/pipermail/mitgcm-support/2012-February/thread.html ; https://github.com/MITgcm/MITgcm/issues/1031

### Sponge bugs fixed upstream: South salinity timescale, T/S relaxation skipped where OB value is 0, east OB `isl-1`
- Cause: (1) `lambda_obcs_s` at the Southern OB used `float(spongeThickness)` instead of `float(spongeThickness-jsl)` (issue #989); (2) T/S relaxation was not applied if the OB value was zero, inconsistent with `obcs_apply_ts.F` (PR #990); (3) in `obcs_sponge_v` (eastern OB) `float(isl-1)` instead of `isl` (2010).
- Fix: update; (1)+(2) fixed in checkpoint69o (2026-07-06), (3) fixed 2010. Results change only where those branches are used (verification `shelfice_2d_remesh`, `obcs_ctrl` unchanged).
- Era: bug present from the introduction of pkg/obcs sponge until 2026-07 (c69n and older).
- Src: https://github.com/MITgcm/MITgcm/issues/989 ; https://github.com/MITgcm/MITgcm/pull/990 ; 2010-April 'obcs_sponge.F, line 330' http://mailman.mitgcm.org/pipermail/mitgcm-support/2010-April/thread.html

### Sponge thicker than a tile starts "suddenly" / sponge appears at every tile edge (parallel)
- Cause: sponge cells are counted in tile-local indices; a sponge wider than sNx/sNy is truncated to the boundary tile.
- Fix: use fewer, larger tiles or coarser dx in the sponge, or pkg/rbcs (tile-independent, 3-D relaxation-strength map; relaxation toward spatially varying targets; can relax U,V,T,S, ptracers).
- Era: 2016 (not re-verified in current source).
- Src: 2016-August 'sponge limited to one process?' http://mailman.mitgcm.org/pipermail/mitgcm-support/2016-August/thread.html ; 2004-August 'question on sponge-layer/mpi' http://mailman.mitgcm.org/pipermail/mitgcm-support/2004-August/thread.html

### Sponge for passive tracers, T-only sponge, "closed" wall + sponge, or sea-ice sponge
- Cause: OBCS sponge relaxes only u,v,T,S (and optionally sea-ice via `useSeaiceSponge`, `seaiceSpongeThickness` in OBCS_PARM05, needs `ALLOW_OBCS_SEAICE_SPONGE`).
- Fix: ptracers -> pkg/rbcs; T-only: `OBCSsponge_Theta` plus switch off the velocity sponge flags; with a bathymetry wall plus a T sponge you still sponge all variables; surface-current/SST nudging -> rbcs (pkg/ctrl cannot restore).
- Era: 2015-2026.
- Src: 2018-July 'OBCS sponge for ptracers' http://mailman.mitgcm.org/pipermail/mitgcm-support/2018-July/thread.html ; 2015-August 'Sponge BCs and passive tracers?' http://mailman.mitgcm.org/pipermail/mitgcm-support/2015-August/thread.html ; 2015-February 'Using a sponge layer to relax temperature' http://mailman.mitgcm.org/pipermail/mitgcm-support/2015-February/thread.html ; 2019-July 'solid boundary with sponge layer' http://mailman.mitgcm.org/pipermail/mitgcm-support/2019-July/thread.html ; 2026-April 'How to restore the surface currents ...' http://mailman.mitgcm.org/pipermail/mitgcm-support/2026-April/thread.html

### Spurious vertical velocity / waveguide / upwelling along an open boundary
- Cause: prescribed U,V,T,S not dynamically consistent (not in geostrophic balance); nonzero tangential velocity whose divergence forces w; bathymetry gradient normal to the OB (prescribed normal flow implies w); corner effects; also one-sided Dirichlet with a narrow sponge.
- Fix (Matt Mazloff): restoring "sponge" region ~10 cells with no normal bathymetry gradient; tangential velocity along the OB = 0 (omit tangential OB files, optionally restore it); corners land; Stevens BCs reduce the problem (vertical shear consistent with dynamics) but are less-tested; prefer pkg/rbcs for more flexible relaxation timescales.
- Era: 2011-2020.
- Src: 2020-July 'Need assistance for tuning currents in the physical model' http://mailman.mitgcm.org/pipermail/mitgcm-support/2020-July/thread.html ; 2014-July 'Logarithmic decrease in Eta without OBCS/EXF' http://mailman.mitgcm.org/pipermail/mitgcm-support/2014-July/thread.html ; 2016-June 'Questions about the OBCS Package' http://mailman.mitgcm.org/pipermail/mitgcm-support/2016-June/thread.html ; 2009-August 'Re: MITgcm-support Digest, Vol 74, Issue 10' http://mailman.mitgcm.org/pipermail/mitgcm-support/2009-August/thread.html

## Stevens (1990) boundary conditions

### `useStevens*`: `stops due to EXTREME Pot.Temp` near the boundary; warnings about nonlinFreeSurf/pTracers
- Cause: normal velocity is vertically averaged and the rest computed from the interior; tangential OB velocities are not reset to zero (Martin has never run Stevens with nonzero tangential velocity); Stevens is only partly supported with `nonlinFreeSurf>0` ("not yet implemented", error text only), pTracers ("expect the unexpected"), sea ice.
- Fix: need `#define ALLOW_OBCS_STEVENS` in OBCS_OPTIONS.h (default #undef), `useStevensWest/...`, `TrelaxStevens/SrelaxStevens` in OBCS_PARM04; do not give tangential velocity files; for SSH-only data convert SSH to barotropic normal velocity (stream function) + Stevens; Martin uses it routinely for climate regional runs, with an (d eta/dn=0) hack for NLFS. Can be mixed per boundary with Orlanski.
- Era: 2009-2021; PR #265 (2019, still open in digest) would restructure the code and change results.
- Src: 2021-May 'Stevens BC' http://mailman.mitgcm.org/pipermail/mitgcm-support/2021-May/thread.html ; 2018-May 'Stevens boundary conditions with sea ice' http://mailman.mitgcm.org/pipermail/mitgcm-support/2018-May/thread.html ; 2016-November 'MITgcm-support Digest, Vol 161, Issue 21' http://mailman.mitgcm.org/pipermail/mitgcm-support/2016-November/thread.html ; https://github.com/MITgcm/MITgcm/pull/265 ; https://github.com/MITgcm/MITgcm/pull/103

## Tides

### Setting up barotropic tides at OBs (units, files, nothing happens)
- Cause: needs `#define ALLOW_OBCS_TIDES` in OBCS_OPTIONS.h and `useOBCStides=.TRUE.` in data.obcs (else tides are silently absent), `OBCS_tidalPeriod(1:n)` in seconds, amplitude files in m/s (normal-velocity amplitude, not elevation) and phase files in seconds (not radians); `obcs_add_tides.F` adds the barotropic tide after balance/sponge targets are set. 1-hour/0.1 m/s "small response" reports were setup issues.
- Fix: TPXO `u,v` are cm/s (divide by 100) or use transports/depth; phase_seconds = phase_cycles * period (phase corrected to start date); files are Nx(or Ny) x nTidalComp (old `OB?amFile/OB?phFile` names are retired: use `OB[N,S,E,W]_[u,v]Tid[Am,Ph]File` since checkpoint68y; `verification/seaice_obcs/input.tides`, MITgcm_contrib/tides). Division by zero when fewer than 10 components given (`tidalPeriod=0`) fixed in checkpoint68g (PR #602). Tides raise CFL: check `advcfl_W`, reduce deltaT.
- Era: 2013-2026.
- Src: 2019-April 'Tidal forcing at the open boundaries' http://mailman.mitgcm.org/pipermail/mitgcm-support/2019-April/thread.html ; 2017-June 'Problem with tides in the OBCS' http://mailman.mitgcm.org/pipermail/mitgcm-support/2017-June/thread.html ; 2016-October 'open boundary prescribe' http://mailman.mitgcm.org/pipermail/mitgcm-support/2016-October/thread.html ; 2025-March 'Rivers imposed as OBCs ...' (PR #752, issue #617 new tide naming) http://mailman.mitgcm.org/pipermail/mitgcm-support/2025-March/thread.html ; https://github.com/MITgcm/MITgcm/pull/602 ; 2017-September 'Tide in the open boundary' http://mailman.mitgcm.org/pipermail/mitgcm-support/2017-September/thread.html

### Tides + sponge/daily boundary data; tidal currents too strong in a narrow channel; no SSH forcing
- Cause: `obcs_add_tides` adds the barotropic tide to the OB velocity *before* sponge relaxation, so a sponge relaxes toward (daily-mean + tide); the tide is depth-independent (no baroclinic tide). Model has no SSH nudging: tidal elevation at the OB is not connected to elevation inside, only volume flux matters. In narrow shallow channels strong currents are volume conservation + too little friction.
- Fix: use `obcs_add_tides` with sponge; for strong channels: `no_slip_sides=.TRUE.`, bottomDragQuadratic >=1e-3 or `zRoughBot`, `viscAhGridMax=0.5`, `viscAhGrid=0.01`, `viscA4Grid=0.01`, widen/deepen channel; GGL90 mixing; for SSH-driven forcing convert to velocity (u = c eta/H) or add an SSH nudging term (Ken Hughes' patch, not in a package).
- Era: 2017-2026 (Dimitris Menemenlis, Martin Losch).
- Src: 2026-September 'Open boundaries and sponge layers' http://mailman.mitgcm.org/pipermail/mitgcm-support/2026-September/thread.html ; 2024-May 'Problem: Extreme Tidal Currents in the Simulation at Narrow Channels.' http://mailman.mitgcm.org/pipermail/mitgcm-support/2024-May/thread.html ; 2022-August 'model unstability caused by increasing sea surface elevation(eta)' http://mailman.mitgcm.org/pipermail/mitgcm-support/2022-August/thread.html ; 2020-May 'Open boundary tidal amplitude' http://mailman.mitgcm.org/pipermail/mitgcm-support/2020-May/thread.html

### Internal-tide sponge gives artificial U/T disturbances; OB internal-wave reflection
- Cause: T/S sponge toward a constant Tref (not the stratification profile) or wrong boundary file; sponge too short/fast/slow; barotropic OB only.
- Fix: sponge toward the time-mean stratification (T/S files), width ~ one mode-1 wavelength, hour-scale timescales at the outer edge; telescoped grid; or pkg/rbcs with a propagating signal. (`Trelax*/Srelax*` are not OBCS-sponge parameters: use `Urelaxobcs*, Vrelaxobcs*` and the `OBCSsponge_Theta/_Salt` switches.)
- Era: 2026-08 (issue #1031 open).
- Src: https://github.com/MITgcm/MITgcm/issues/1031 ; 2021-May 'Open boundary conditions for Internal Waves.' (Jody) http://mailman.mitgcm.org/pipermail/mitgcm-support/2021-May/thread.html ; 2023-May 'Problem with Obcs when simulating tide' http://mailman.mitgcm.org/pipermail/mitgcm-support/2023-May/thread.html

## Sea ice at open boundaries

### Sea ice piles up/gets stuck at OBs; `Orlanski ... useSEAICE`; periodic ice with ocean OBs
- Cause: pkg/seaice + OBCS is incomplete (hacks): ice strength (`press`) and viscosities are not set at OBs; prescribed ice at low frequency; ice convergence at edges.
- Fix: choose among: `useSeaiceNeumann=.TRUE.` (Neumann BC for sea ice variables, checkpoint68z, 2024-07); CPP `OBCS_SEAICE_AVOID_CONVERGENCE`, `OBCS_SEAICE_COMPUTE_UVICE`, `OBCS_SEAICE_SMOOTH_*` (all #undef by default); ice sponge `useSeaiceSponge=.TRUE.` (+ `ALLOW_OBCS_SEAICE_SPONGE`, `seaiceSpongeThickness`; test `seaice_obcs.seaiceSponge`); Orlanski + seaice is an error; prescribed ice in obcs files (OBN/S/E/W: ice thickness `h`, area `a`, snow `sn`, `uice`,`vice`); free slip on ice at OBs would need zero shear/bulk viscosity (not coded); to keep sea ice periodic while the ocean has OBs, make `obcs_adjust_uvice.F`, `obcs_apply_seaice.F`, `obcs_apply_uvice.F`, `obcs_seaice_sponge.F` return immediately. RBCS is not extended to sea ice (would need coding).
- Era: 2013-2024.
- Src: 2022-November 'RBCS for sea ice boundary conditions' (Martin Losch) http://mailman.mitgcm.org/pipermail/mitgcm-support/2022-November/thread.html ; 2018-February 'OBCS for ocean, but periodic for sea ice?' http://mailman.mitgcm.org/pipermail/mitgcm-support/2018-February/thread.html ; 2013-June 'OBCS for sea ice!' http://mailman.mitgcm.org/pipermail/mitgcm-support/2013-June/thread.html ; 2005-June 'obcs and seaice' http://mailman.mitgcm.org/pipermail/mitgcm-support/2005-June/thread.html

## Passive tracers at open boundaries

### Passive tracer concentration explodes (1e38, negative) near an OB; what is the default tracer OB?
- Cause: default OB for pTracers is Neumann-like (value copied from the interior), not zero; inflow then carries interior values; no tracer OB file; unlimited advection scheme.
- Fix: prescribe `OB[N,S,E,W]ptrFile(iTr)` (see `verification/exp4`, `so_box_biogeo`) to get Dirichlet (zeros for no inflow); flux-limited scheme (advScheme 33) to cut undershoots; `OBCS_u1_adv_Tr(iTr)=1`; RBCS for sources/sinks; with Orlanski see the combination error above; offline tracer pkgs (e.g. cfc, dic offline) do not support OBCS.
- Era: 2004-2019.
- Src: 2019-January 'OBCS Tracer Errors' http://mailman.mitgcm.org/pipermail/mitgcm-support/2019-January/thread.html ; 2017-April 'ptracer and default boundary conditions' http://mailman.mitgcm.org/pipermail/mitgcm-support/2017-April/thread.html ; 2016-January 'OBCS for PTracers' http://mailman.mitgcm.org/pipermail/mitgcm-support/2016-January/thread.html ; 2008-April 'obcs for ptracers' http://mailman.mitgcm.org/pipermail/mitgcm-support/2008-April/thread.html

## Domain, parallel layout, misc

### 2-D slab (sNy=1) with OBCS: T/S go to zero or the run blows up after the first step
- Cause: overlap too small for the advection scheme (scheme 33 needs OLy>=2-3, scheme 7 / 77 needs 4) when sNy=1; or only E/W OBs with periodic N/S (model is doubly periodic by default: put a wall row or rely on Ny=1 exchange rules); OLy=1 "worked" only for centered schemes.
- Fix: `OLx=OLy>=2` in general, `>=4` for 7-point schemes; Jean-Michel: Ny=1 is supported by the default exchange routines (no need for sNy>OLy); if you only want E/W OBs, leave OB_Jnorth/Jsouth=0 (no need to undef ALLOW_OBCS_NORTH/SOUTH). For zero meridional velocity set f0=0 or wall row.
- Era: 2010-2017.
- Src: 2010-July '2D setup' http://mailman.mitgcm.org/pipermail/mitgcm-support/2010-July/thread.html ; 2017-May 'problems with OBCS' http://mailman.mitgcm.org/pipermail/mitgcm-support/2017-May/thread.html

### Default boundaries are doubly periodic; OB values "appear on both sides"; extra Nx+1/Ny+1 column in netCDF output
- Cause: with `useOBCS=.FALSE.` (or OB_*=0 and no wall) a point on one side is a copy of the opposite side; NetCDF/mnc u,v output includes the periodic halo (U(1) copied into U(Nx+1)). With OBs, OB values override periodicity; the extra column is harmless.
- Fix: wall = depth 0 at i=1 (and j=1); ignore Nx+1/Ny+1 in mnc output; MDS output has no padding; to damp a closed box use rbcs or a sponge with OB=1,Nx.
- Era: any.
- Src: 2012-July 'doubt regarding OBCs' http://mailman.mitgcm.org/pipermail/mitgcm-support/2012-July/thread.html ; 2014-January 'OBCS (open boundary forcing) problem' http://mailman.mitgcm.org/pipermail/mitgcm-support/2014-January/thread.html ; 2017-March 'Strange U-velocity and V-velocity on North, East' http://mailman.mitgcm.org/pipermail/mitgcm-support/2017-March/thread.html

### Freshwater (river) through an OB: negative salinity, cg2d diverges at step 1
- Cause: vertical CFL (thin surface cells + inflow) and advection undershoot at a zero-salinity front, not a package limit; freshwater OBs work in glacial-fjord setups.
- Fix: reduce deltaT (3 s fixed a 10 s case), check `advcfl_*` (< ~0.5), ramp the inflow, flux-limited scheme, extra mixing (MY82/GGL90).
- Era: 2012.
- Src: 2012-December 'Freshwater OBCS problem' http://mailman.mitgcm.org/pipermail/mitgcm-support/2012-December/thread.html

### OBCS + EXF: `EXF_CHECK_RANGE` / hflux out of range after enabling OBs; OBCS not combined with CD scheme; time/NaN near boundary
- Cause: boundary-driven spurious surface temperatures/fluxes (strong convergence near OBs), not a coding error in EXF; check range is a diagnostic stop.
- Fix: inspect snapshots near the OB (not time means); `useExfCheckRange=.FALSE.` to continue; reduce deltaT; keep OB data consistent with interior; viscAhGrid/stability checks; see also Orlanski error list for CD scheme.
- Era: 2009-2025.
- Src: 2025-June 'Questions about EXF and OBCS' http://mailman.mitgcm.org/pipermail/mitgcm-support/2025-June/thread.html ; 2009-June 'cal, exf, obcs: problem at the tile boundaries' http://mailman.mitgcm.org/pipermail/mitgcm-support/2009-June/thread.html

### Non-hydrostatic + OBCS stop in CONFIG_CHECK / w at the OB
- Cause: (2009) `config_check` stopped nonHydrostatic + nonlinFreeSurf; w at an OB is reset each step by `obcs_apply_w` (needed for advective w tendencies); an ALLOW_OBCS_SOUTH guard around the western OB in `obcs_apply_w.F` (2008 typo) disabled the western w-BC for E-W channels when N/S flags were undef.
- Fix: look at STDERR for the real message; the typo was fixed 2008-11; no relaxation of w in the sponge is possible (w from continuity).
- Era: 2008-2010.
- Src: 2008-November 'Possible OBCS bug' http://mailman.mitgcm.org/pipermail/mitgcm-support/2008-November/thread.html ; 2010-April 'nonhydrostatic w boundary condition' http://mailman.mitgcm.org/pipermail/mitgcm-support/2010-April/thread.html ; 2009-July 'NonHydrostatic and OBCs' http://mailman.mitgcm.org/pipermail/mitgcm-support/2009-July/thread.html

### One-way nesting / regional ECCO-style downscaling: what pkg/obcs can do
- Cause: expectations of two-way nesting.
- Fix: no two-way nesting; run coarse model, save U,V,T,S (+eta), interpolate to boundary (inward velocity points), balance net flow, then pkg/obcs prescribe (+sponge); stretched curvilinear grids limit deltaT by smallest cell. `verification/seaice_obcs` (EXF+cal+obcs+ice), `obcs_ctrl` (adjoint controls), `exp4`, `dome`, `tutorial_plume_on_slope`, `so_box_biogeo` are the examples.
- Era: 2008-2013.
- Src: 2008-June 'grid nesting help' http://mailman.mitgcm.org/pipermail/mitgcm-support/2008-June/thread.html ; 2013-September 'regional zooming method' http://mailman.mitgcm.org/pipermail/mitgcm-support/2013-September/thread.html ; 2010-August 'Problems about OBCS and EXF' http://mailman.mitgcm.org/pipermail/mitgcm-support/2010-August/thread.html
