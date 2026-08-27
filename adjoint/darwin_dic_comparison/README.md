# darwin_dic_comparison

Companion to `darwin_bling_comparison/`: the same 38-point centered-FD
methodology, comparing Darwin's own forward finite-difference sensitivity
of the SOCAT-pCO2 cost function against **pkg/dic's** adjoint (instead of
BLING's). pkg/dic carries 5 prognostic tracers (DIC, ALK, PO4, DOP, O2 --
no NO3 or Fe), all 5 compared here.

Darwin's own FD values for DIC, ALK, PO4, and O2 are reused verbatim from
`darwin_bling_comparison`'s `multipoint_comparison_results.csv` -- Darwin's
side of that computation never depended on which adjoint model it's being
compared against. Only DOP required fresh Darwin perturbation runs
(`global_oce_biogeo_darwin/run_multi_DOP_plus`/`_minus`), since BLING has
no DOP control to have motivated computing it earlier.

## Retraction: there is no PO4/DOP adjoint corruption

**An earlier version of this README reported that pkg/dic's adjoint had a
"real, field-specific numerical instability" corrupting PO4 and DOP, and
excluded the affected points. That claim is withdrawn -- it was an
artifact of the diagnostic, not a property of the model.** Both the
original numbers and the "corrected" ones in that version are wrong; the
corrected ones were worse.

The flag was `abs(adxx) > 1e6`, applied in
`compute_multipoint_comparison_dic.py` to `a`, the **raw control-space**
adjoint. Control-space scale differs enormously between tracers. Surface
(k=1) median `|adxx|` over all wet points:

| field | tracer | median &#124;adxx&#124; | 1e6 sits at |
|---|---|---|---|
| ptr1 | DIC | 1.1e3 | far above any point |
| ptr2 | ALK | 8.9e2 | far above any point |
| ptr3 | PO4 | 1.3e4 | ~99th pct |
| ptr4 | DOP | **4.1e5** | **~60th pct** |

So one absolute cut, implicitly calibrated on DIC/ALK, landed in the
*middle* of DOP's ordinary distribution and removed 1028 of 4447 wet
surface points globally, 15 of 38 in the sample. Excluding points based
on the magnitude of the dependent variable is selection bias, and it did
what selection bias does: it inverted DOP's regression (slope -5.57 ->
+27.95, r -0.932 -> +0.652).

### Evidence that the flagged points are ordinary

1. **Same distribution shape for all four tracers.** Surface
   `log10|adxx|` is unimodal and lognormal for every tracer, with a
   near-identical interpercentile width -- the 25th-to-99th percentile
   spans 2.52 decades for DIC, 2.53 for ALK, 2.85 for PO4, 2.53 for DOP.
   The fields differ only by a *shift* along the log axis (medians 1e3.04,
   1e2.95, 1e4.13, 1e5.61). There is no second population.

2. **Smooth vertical structure.** `|adxx|` decays monotonically by
   ~20-30x per level from the surface to 1e-19 at k=15, for both PO4 and
   DOP. A numerical instability does not produce fifteen levels of clean
   geometric decay.

3. **Spatial coherence.** Ratio of a point's `|adxx|` to the median of its
   wet 8-neighbours: field-wide median 0.99, p95 6.8. For the flagged DOP
   points that ratio is **1.35** -- i.e. they sit smoothly inside their
   neighbourhood. Flagged PO4 points give 3.50, still well below the
   field's own p95. Corruption produces isolated speckle; this is not
   speckle.

4. **Darwin independently confirms them.** At every flagged point,
   Darwin's forward FD -- computed with no reference whatsoever to
   pkg/dic -- reproduces the sign and rough magnitude, with a consistent
   ratio (~29-139, clustered near 40, for the 14 clean flagged DOP
   points; 1.39 for the single flagged PO4 point). Two independent models
   do not agree on garbage.

The `run_dic_C`/`run_dic_D` single-control reruns were correctly executed
and correctly reported to be bit-for-bit identical to the 7-control
extraction. That result was *misread* as confirming a bug. It is equally
consistent with -- and, given the above, actually evidence for -- the
values simply being correct and reproducible. It remains true that the
slot the field occupies matters: see the slot-order constraint in
`hybrid_darwin_dic/run_dic/data.ctrl`. That is a separate issue from this
retraction, and slots 5 and 6 are not affected by it.

## Results (all 38 points, no exclusions)

Because OLS on these samples can be dominated by a single high-leverage
point, the leverage-immune **Theil-Sen** slope and **Spearman rank r** are
reported alongside. Where the two columns agree, the relationship is real;
where they diverge, OLS is being carried by one or two points.

| Tracer | OLS slope | OLS r | Theil-Sen slope | rank r | reading |
|---|---|---|---|---|---|
| DIC | 1.157 | 0.999 | 1.128 | 0.954 | **genuine agreement** |
| ALK | 1.176 | 0.998 | 1.142 | 0.954 | **genuine agreement** |
| PO4 | 1.394 | 0.998 | -1.094 | **-0.075** | OLS is pure leverage; **no real correlation** |
| DOP | -5.575 | -0.932 | **43.4** | **0.702** | OLS sign is a 1-point artifact; **consistent positive ~40x** |
| O2 | -- | -- | -- | -- | both models exactly 0 |

**DIC and ALK are the internal control, and they pass.** Both estimators
agree closely for exactly the two tracers that reach the cost through
carbonate chemistry with no dependence on which biogeochemistry package
computed the tracer state. This is what a real agreement looks like in
this table, and it validates the control-space -> physical conversion
(`dJ/dc = adxx * sqrt(w)`): if that direction were inverted, DIC's slope
would land near 58, not 1.16.

**PO4 does not actually agree.** The reported r=0.998 came entirely from
one point (Darwin -2.00e5 vs pkg/dic -2.78e5) that is ~300x the sample
median and dominates the fit. That point is *not* corrupt -- the two
models agree there to within 39%, the best single agreement in the PO4
sample -- but across the other 37 points there is **no rank correlation at
all** (Spearman -0.075). The previous README's instinct that PO4's
apparent agreement was an artifact was right; its diagnosis (corruption)
and its replacement number (0.65/0.322, itself computed on a
selection-biased sample) were both wrong.

**DOP is consistently proportional, not sign-flipped.** The OLS slope of
-5.575 was likewise a single-point artifact: one point has Darwin
dJ = -5.00e4, about 60x larger than any other Darwin DOP value in the
sample, and sign-flipped relative to pkg/dic. Robustly, pkg/dic's DOP
adjoint is **positive, well rank-correlated with Darwin's (0.702), and
larger by a consistent factor of ~40** (Theil-Sen 43.4, median ratio
39.8). This is visible directly in `figures/darwin_vs_dic_multipoint.png`:
the DOP points sit ~1.5 decades off the 1:1 line in *both* quadrants,
which is a scale offset, not a disagreement about direction.

**O2 is exactly zero on both sides, for two different structural
reasons.** Darwin's FD is zero because its carbon-chemistry solver has no
O2 term (same finding as the BLING comparison). pkg/dic's adjoint is zero
because, per `dic_biotic_forcing.F:315`, `GO2 = R_OP*GPO4` -- O2 is
computed *from* `GPO4` but never feeds back into `GDIC`/`GALK`. It is a
strict one-way diagnostic tracer with no adjoint pathway at all, unlike
BLING, whose smooth Michaelis-Menten remineralization gives O2 a tiny but
nonzero sensitivity. This finding is unchanged and was verified directly
against the source.

## Why the singularity hypothesis was checked and rejected

The natural suspect for a PO4-specific blow-up is `bio_export.F:150`:

```fortran
nutlimit = tmppo4/(tmppo4+KPO4)     ! d/dPO4 = KPO4/(PO4+KPO4)**2
```

The adjoint build uses pkg/dic's stock `DIC_OPTIONS.h` (`build_ad/DIC_OPTIONS.h`
is a symlink to it) with `DIC_NO_NEG`, `ALLOW_FE` and `DIC_AD_SAFE` all
**undefined**, so `tmppo4` is the raw tracer with no clamp, and the
derivative is singular at `PO4 = -KPO4 = -5.0e-4 mol/m3` (KPO4 is the
default; `data.dic` does not override it). Given `PTRACERS_ref(3) =
5.438e-4`, PO4 lives right at that scale, so this looked promising.

**It does not happen.** Checking the forward field directly
(`PTRACER03` at iterations 0, 4 and 2880) gives `min(PO4) = 0.0` exactly
at every point and every level, so `PO4 + KPO4 >= 5e-4` always and the
derivative is bounded by `1/KPO4 = 2000`. No singularity is reached. This
is recorded because it is the obvious hypothesis and it is worth not
re-deriving it.

## Open item

The **~40x DOP scale offset is real and unexplained.** It is not a
bookkeeping error: `run_ad/data.ctrl` registers slot 6 as
`xx_ptr4` <-> `wt_DOP.bin`, matching what the analysis script reads, and
the weight fields are uniform constants per tracer (DOP 1.5e-5, so
`sqrt(w) = 3.873e-3`), so the conversion cannot introduce spatial
structure. `DOPfraction = 0.67` does not account for a factor of 40, and
neither does `KDOPRemin = 1/(180 d)`. Identifying its source -- most
plausibly a stoichiometric or units mismatch between pkg/dic's
single-pool DOP and Darwin's multi-pool organic phosphorus -- is the
natural next step.

The secondary open item is **PO4's lack of correlation** (rank r =
-0.075), which the corruption narrative previously obscured. That is now
the more interesting of the two disagreements and has not been
investigated.

## Directory layout

```
scripts/    compute_multipoint_comparison_dic.py  (reports OLS + Theil-Sen + Spearman)
            plot_multipoint_comparison_dic.py, plot_multipoint_nutrients_dic.py
figures/    darwin_vs_dic_multipoint.png (main 5-tracer comparison)
            darwin_vs_dic_nutrients.png (PO4/DOP only, linear axes -- pulled out
            separately since their range is hard to read against DIC/ALK's)
results/    multipoint_comparison_results_dic.csv
            pilot_adxx_fields_dic.npz -- pkg/dic's raw adxx_ptr1-5 fields
            (whole-domain), extracted once via MITgcmutils.rdmds from
            global_oce_biogeo_dic_SOCAT/run_ad's cycle-0 run. Archived here
            for provenance; compute_multipoint_comparison_dic.py reads its
            own copy from the live experiment tree, not from results/.
```

The `corrupted` column has been removed from the results CSV along with
the flag that produced it.
