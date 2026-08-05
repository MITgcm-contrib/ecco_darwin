C------------------------------------------------------------------------------|
C                           DIAGNOSTICS_VEC_SIZE.h
C------------------------------------------------------------------------------|

C     VEC_points : max number of points storable in one diagnostics_vec mask.
C                  Must be >= the largest per-boundary point count printed by
C                  gen_dvmasks.py -v (see STEP1, Section V.c).
C     nVEC_mask  : number of lateral masks, i.e. nml_vecFiles entries in
C                  data.diagnostics_vec. The accompanying data.diagnostics_vec
C                  declares 20 (4 boundaries x 5 slots).
C     nSURF_mask : number of surface masks, i.e. nml_sfFiles entries.
C                  The accompanying data.diagnostics_vec declares none.

      INTEGER, PARAMETER :: VEC_points = 2000
      INTEGER, PARAMETER :: nVEC_mask = 20
      INTEGER, PARAMETER :: nSURF_mask = 0
