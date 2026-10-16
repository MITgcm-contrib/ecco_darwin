# Build the iteration-2 closure GF control dirs (sourced by job_gf_cl2 + job_gf_cl2_test; 2026-10-15, D118/D121 = A).
# Per-site control = v11_nn regional set (gyre: HOT BATS HydroS; highlat: PAPA PAP SOTS K2 MUNIDA OWSM; CARIACO, DYFAMED:
# v05 baseline) + closure (scripts/make_closure_runs.py --inplace). Mac dry run (2026-10-15) of exactly these steps
# reproduces runs/closure/<SITE> + runs/diazo/s2_{HOT,BATS} byte for byte (namelists, traits, ICs, relax, diffkr).
# Needs: $R (runs dir), $SITES. v11 sets use the 19 v05 controls, so GF_CONTROLS must be unset for gf_make_opt.
unsetenv GF_CONTROLS
rm -rf $R
python3 scripts/make_runs.py sites input $R/base build/mitgcmuv $SITES || exit 1
set G = (); set H = (); set B = ()
foreach s ($SITES)
  switch ($s)
    case HOT:
    case BATS:
    case HydroS:
      set G = ($G $s); breaksw
    case CARIACO:
    case DYFAMED:
      set B = ($B $s); breaksw
    default:
      set H = ($H $s)
  endsw
end
if ($#G > 0) then
  python3 scripts/gf_make_opt.py results/optimized/gyre_step1.csv $R/base $R/opt gyre $G
  if ($status != 0) exit 1
endif
if ($#H > 0) then
  python3 scripts/gf_make_opt.py results/optimized/highlat_step1.csv $R/base $R/opt highlat $H
  if ($status != 0) exit 1
endif
mkdir -p $R/ctrl
foreach s ($G)
  mv $R/opt/gyre_$s $R/ctrl/$s
end
foreach s ($H)
  mv $R/opt/highlat_$s $R/ctrl/$s
end
foreach s ($B)
  mv $R/base/$s $R/ctrl/$s
end
foreach s ($SITES)
  python3 scripts/make_closure_runs.py --inplace $R/ctrl/$s || exit 1
end
setenv GF_CONTROLS closure
