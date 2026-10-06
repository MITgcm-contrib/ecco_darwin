# pkg/admtlm

Adjoint/tangent-linear combined (ADM-TLM) driver for singular-vector / Hessian-type calculations (legacy).

**runtime switch:** `useADMTLM`-style flag in `data.pkg` (check exact name in packages_boot.F)

## Headers
- `ADMTLM_OPTIONS.h` — BOP
- `arpack_debug.h` — \SCCS Information: @(#) FILE: debug.h   SID: 2.3   DATE OF SID: 11/16/95   RELEASE: 2 %---------------------------------% See debug.doc for documentat

## Routines (8)
`admtlm_bypassad.F`, `admtlm_driver.F`, `admtlm_dsvd.F`, `admtlm_dsvd2model.F`, `admtlm_init_fixed.F`, `admtlm_map.F`, `admtlm_metric.F`, `admtlm_model2dsvd.F`

## Called from outside the package
- `ADMTLM_DSVD` ← `eesupp/src/main.F:185`
