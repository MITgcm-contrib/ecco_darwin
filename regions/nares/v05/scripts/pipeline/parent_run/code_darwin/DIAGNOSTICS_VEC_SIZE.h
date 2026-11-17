C------------------------------------------------------------------------------|
C                           DIAGNOSTICS_VEC_SIZE.h
C------------------------------------------------------------------------------|
C pan_AO downscale. Resized 2026-10-01, third and final time, after the masks
C were regenerated along the DISC RIM by gen_dvmasks_disc.py (job 25215198):
C
C   boundary  2026-09-17   2026-09-29   2026-10-01 (disc, CURRENT)
C             (permuted)   (square)
C   east         215          265            141
C   west         284            0            736   <-- now the LARGEST
C   north        661          358            256
C   south        622          324            685
C   total       1782          947           1818
C
C THE WEST BOUNDARY IS BACK. Dropping it on 2026-09-29 was correct for the
C config as it then stood and wrong about the domain: this is a polar
C stereographic DISC, so the model array's square perimeter is not the open
C boundary (only 174 of its 6720 cells are wet; its east edge is 0 wet). The
C square's west edge happens to run across continental Siberia, which is why
C that mask came out empty -- not because there is no western open boundary.
C On the rim circle west carries more points than any other side. The
C RECL=0 / `forrtl: severe (37)` crash that forced the drop (total_valid = 0
C makes INQUIRE(IOLENGTH=LEN) return 0) cannot arise now: no mask is empty.
C
C VEC_points = 600, UNCHANGED, and this is the number to be careful about,
C because the global totals above are NOT what it bounds. The check in
C diagnostics_vec_check.F is PER-PROCESS:
C   vec_numPnts_allproc(vec_id, mpiMyId+1) .GT. VEC_points
C i.e. per mask entry, per rank. The parent is LLC270 at sNx=sNy=30, nSx=nSy=1,
C so one 30x30 tile per rank, 900 cells, and a single mask can at most put 900
C points on one rank. Measured directly off the four .bin files on 2026-10-01
C by binning the compact (13,270,270) arrays into 30x30 tiles:
C
C   boundary   global   max points on any one tile   tiles touched
C   west        736                 101                   26
C   east        141                  30                    7
C   south       685                 107                   19
C   north       256                  30                   13
C
C Worst case is 107 (south), so 400 leaves ~3.7x headroom on the quantity that
C is actually bounded. It was 700 originally, sized against the old bogus
C global north count of 661 -- sizing VEC_points against a GLOBAL count is the
C conservative error, not the dangerous one, but it is still the wrong
C quantity. The genuinely global buffer is the gather in
C diagnostics_vec_write_bin.F, bufferR4_2D((VEC_points)*(nPx*nPy)) =
C 400*767 = 306,800 slots, against a largest total_valid of 736. Ample.
C
C nVEC_mask = 40 boundaries x 5 mask-entries per boundary (ETAN+PTRACE01-08,
C UICE/VICE+PTRACE09-16, AREA+PTRACE17-24, HSNOW+PTRACE25-31,
C HEFF+THETA/SALT/UVEL/VVEL); no surface (SF) fields. Mirrors kerguelen's
C structure exactly. Must accompany
C pan_AO/parent_run/run/data.diagnostics_vec (nml_vecFiles dimensioned
C nVEC_mask = 40 there too) -- these two files are a matched pair, not
C independently editable. Changing nVEC_mask REQUIRES A REBUILD.

      INTEGER, PARAMETER :: VEC_points = 600
      INTEGER, PARAMETER :: nVEC_mask = 40
      INTEGER, PARAMETER :: nSURF_mask = 0
