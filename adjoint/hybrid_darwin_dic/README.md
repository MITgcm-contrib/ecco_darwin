# hybrid_darwin_dic

Hybrid surrogate-gradient optimization using **pkg/dic**'s adjoint against
a **Darwin** cost function -- the pkg/dic counterpart to
`hybrid_darwin_bling/`.

Darwin has no adjoint. The hybrid loop therefore takes the *cost* from a
Darwin forward run and the *gradient* from a cheaper package that does
have one, and asks whether the mismatched pair still converges.
`hybrid_darwin_bling/` answers that with BLING; this directory answers it
with pkg/dic.

Everything here is the **1-month** (month01) SOCAT pCO2 cost, matching the
`pilot_1month` rows in `hybrid_darwin_bling/results/cost_trajectory_results.csv`,
not the full-year cost. Controls are DIC and ALK initial conditions.

## Two attempts

### 1. M1QN3-driven (`scripts/run_hybrid_loop_dic.sh`) -- stalled

The first attempt reused the BLING loop's mechanism: patch Darwin's `J`
into `ecco_cost_MIT_CE_000.optNNNN` leaving pkg/dic's gradient payload
intact, then step with `optim.x`. It did not move. Cycles 0 and 1 returned
Darwin `J = 4.090890e+04` identically, and the control guess was unchanged
between them (`xx_ptr1` surface min/mean/max `0/1081.24/2143.28` in both).

The cause is that M1QN3 packs *every registered control* into a single
vector. With all 7 fields registered, slot 2 (salt) returns surface values
around 1e12 against a normal 1e2-1e5 (see **Slot order** below), and that
one field dominates the two-norm and destroys the step sizing -- even
though the fields actually being optimized, DIC and ALK, are clean.

### 2. Manual steepest descent (`scripts/run_hybrid_loop_dic_manual.sh`) -- converges noisily

The working loop bypasses `optim.x` entirely and reads the raw `adxx_ptr*`
output directly (`scripts/manual_gradient_step.py`):

1. run pkg/dic forward+adjoint and take `adxx_ptr1` (DIC), `adxx_ptr2` (ALK)
2. convert to physical space: `dJ/dphysical = adxx * sqrt(ctrl_weight)`
3. steepest-descent step, clipped per point to +/-1 prior sigma (the
   `wt_*.bin` sigma), so no single grid point can move more than one prior
   standard deviation per cycle
4. write the updated fields back to `dic_init.bin` / `alk_init.bin`
5. run Darwin forward and evaluate its own SOCAT cost

## Results

`results/cost_trajectory_dic.csv`. Darwin's own cost `J` is the objective
that matters; `dic_fc` is pkg/dic's internal cost, shown for reference.

| steps taken | pkg/dic `fc` | Darwin `J` |
|---|---|---|
| 0 | 35214.3 | 40908.9 |
| 1 | 35214.3 | 11808.0 |
| 2 | 9117.2 | 9885.9 |
| 3 | 9653.0 | 6662.8 |
| 4 | 4172.8 | 10290.5 |
| 5 | 10317.9 | 7200.9 |

**Read the cycle numbering carefully.** The manual loop takes its gradient
step *before* the Darwin run within the same cycle, so its "cycle 0" is
already one step in. The CSV keeps the loop's own cycle labels; the table
above is re-indexed by steps taken so it lines up with the BLING pilot,
whose loop applies the control at the *start* of a cycle.

Against the BLING pilot on the same cost (40909 -> 34686 -> 22281 -> 18439
-> 13096 -> 6321), pkg/dic drops much faster initially -- 71% in a single
step, where BLING needs five cycles to reach comparable ground -- but then
oscillates: `J` rises at step 4, and the final value (7201) is worse than
its own best (6663 at step 3). BLING ends slightly lower and monotonically.

The gradient direction is clearly carrying real information; what is
missing is any damping on the step length. Bypassing `optim.x` removed the
line search along with the bad packed vector, and the oscillation is the
cost of that. A line search on the manual step, or restoring M1QN3 with
only clean fields registered, is the obvious next step.

## Slot order in `data.ctrl`

**The gradient returned for `genarr3d` slot 2 is not usable in this
configuration.** Whichever field occupies slot 2 comes back with surface
values around 1e12 against a normal 1e2-1e5; the corruption is confined to
level k=1 and follows the slot, not the field.

Practical consequences:

- Slot 2 must hold a field whose gradient you do not need. In `run_dic/`
  that is salt, and `adxx_salt` is simply ignored.
- Do **not** close the gap by leaving slot 2 blank. An interior blank slot
  aborts the run in `MDS_READ_FIELD` on an empty filename.
- Slots 1, 3 and 4 are verified clean: DIC in slot 3 and ALK in slot 4
  agree **bit-for-bit** with the same fields run as the sole registered
  control in slot 1.

That last point means **one 7-control run is sufficient**. The
`run_dic_A` + `run_dic_B` two-run arrangement that produced the committed
trajectory registers DIC and ALK separately in slot 1, at roughly twice
the adjoint cost (measured 23.3 min each, versus 22.9 min for a single
7-control run at the same resolution and length) for bit-identical
gradients. New work should read `adxx_ptr1` and `adxx_ptr2` from a single
`run_dic` run instead.

## Layout

```
scripts/    run_hybrid_loop_dic.sh          M1QN3-driven loop (stalled; kept for the record)
            run_hybrid_loop_dic_manual.sh   manual steepest-descent loop (the working one)
            manual_gradient_step.py         adxx -> physical gradient, clipped step
            apply_control_to_darwin_ic_dic.py   xx_ptr1/xx_ptr2 -> Darwin's physical IC
            compute_darwin_cost_dic.py      Darwin's own month01 SOCAT pCO2 cost
            patch_ecco_cost_dic.py          splice Darwin's J into ecco_cost, gradient untouched
results/    cost_trajectory_dic.csv         both attempts, 0-based cycle labels as logged
run_dic/    7-control config (DIC slot 3, ALK slot 4, salt sacrificial in slot 2)
run_dic_A/  DIC only, slot 1     run_dic_C/  PO4 only, slot 1
run_dic_B/  ALK only, slot 1     run_dic_D/  DOP only, slot 1
```

Only namelists are committed. Run output, adjoint tapes and binary IC
files are not: a cycle writes ~10 GB of tape scratch that is deleted at
the end of each cycle.

## Caveats

- `run_dic_A` and `run_dic_B` never advance `optimcycle`, so each cycle
  **overwrites** `adxx_*.0000000000` in place. The files left in those
  directories after a run are the *last* cycle's gradients, not cycle 0's.
  Check file mtimes against the loop log before comparing them to anything.
- pkg/dic's O2 adjoint is identically zero by construction: per
  `dic_biotic_forcing.F`, `GO2 = R_OP*GPO4` -- O2 is computed *from* `GPO4`
  and never feeds back into `GDIC`/`GALK`, so it has no adjoint pathway.
- pkg/dic optimizes all 7 registered fields internally, so its own
  trajectory sees theta/salt/PO4/DOP/O2 perturbations. The DIC+ALK-only
  scope is enforced downstream, in `apply_control_to_darwin_ic_dic.py`,
  which reads only `xx_ptr1` and `xx_ptr2`.
