# pkg/iceplume (as compiled for Nares)

Subglacial-discharge plume parameterization (Cowton et al. 2015) for MITgcm,
packaged by Mike Wood: https://github.com/mhwood/iceplume, commit e4d4ac9
(https://github.com/mhwood/iceplume/commit/e4d4ac9), with 13 fixes made for
Nares (local branch `nares_fixes`, HEAD c2e6dbd, 2026-10-06: volume, salt and
every ptracer close to <= 1.5e-14 of inventory for open and shelf plumes; 1 vs
4 ranks agree to round-off; restart bit-for-bit; heat closure not
independently verified). The 15 files here are byte-identical to
`darwin3/pkg/iceplume` in the Nares build tree (checked 2026-11-16).

Install: copy this directory to `<darwin3>/pkg/iceplume` before running
genmake2. `iceplume` is listed in `../../code/packages.conf`. The core files
that the package needs patched are in `../../code/` (see `../../PROVENANCE.md`
section 1): `do_oceanic_phys.F`, `exf_getffields.F`, `external_fields_load.F`,
`gchem_forcing_sep.F`, `obcs_balance_flow.F`, `packages_boot.F`,
`packages_check.F`, `packages_init_fixed.F`, `packages_init_variables.F`,
`packages_readparms.F`, `ptracers_apply_forcing.F`, plus the iceplume blocks
in `apply_forcing.F` and `PARAMS.h`. `ALLOW_ADDFLUID` must be defined and
`selectAddFluid = 1` set in `data`.

Runtime: `data.iceplume` (plume mask, plume length, `runoffQsg` yearly files;
see `../../PIPELINE.md`). `useICEPLUME = .FALSE.` in every Nares run so far.

License: `LICENSE` (BSD 3-Clause, Mike Wood, 2024) covers the iceplume files.
`opkda1.F`, `opkda2.F` and `opkdmain.F` are ODEPACK (Alan Hindmarsh, LLNL,
https://computing.llnl.gov/projects/odepack); they carry an author line but no
license header. Not verified: the ODEPACK redistribution terms.
