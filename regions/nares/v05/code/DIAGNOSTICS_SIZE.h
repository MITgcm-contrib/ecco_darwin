C     Diagnostics Array Dimension
C     ---------------------------
C     ndiagMax   :: maximum total number of available diagnostics
C     numlists   :: maximum number of diagnostics list (in data.diagnostics)
C     numperlist :: maximum number of active diagnostics per list (data.diagnostics)
C     numLevels  :: maximum number of levels to write    (data.diagnostics)
C     numDiags   :: maximum size of the storage array for active 2D/3D diagnostics
C     nRegions   :: maximum number of regions (statistics-diagnostics)
C     sizRegMsk  :: maximum size of the regional-mask (statistics-diagnostics)
C     nStats     :: maximum number of statistics (e.g.: aver,min,max ...)
C     diagSt_size:: maximum size of the storage array for statistics-diagnostics
C Note : may need to increase "numDiags" when using several 2D/3D diagnostics,
C  and "diagSt_size" (statistics-diags) since values here are deliberately small.
      INTEGER    ndiagMax
      INTEGER    numlists, numperlist, numLevels
      INTEGER    numDiags
      INTEGER    nRegions, sizRegMsk, nStats
      INTEGER    diagSt_size
C 2026-09-30: 1500 -> 5000, matching ecco_darwin v06/1deg
C code_darwin_v4r5_v2_AO. ndiagMax counts every AVAILABLE diagnostic
C registered by DIAGNOSTICS_ADDTOLIST, not just the active ones, and adding
C pkg radtrans + pkg oasim + DARWIN_DIAG_IOP registers a block of them PER
C WAVEBAND (nlam=13: Ed/Es/Eu, and absorption/scattering/backscatter), on top
C of the already-enabled DARWIN_DIAG_PERTYPE and DARWIN_DIAG_TENDENCIES.
C Overflow is a loud init-time stop, not silent, but there is no reason to
C discover it in the queue.
      PARAMETER( ndiagMax = 5000 )
      PARAMETER( numlists = 1000, numperlist = 1000, numLevels=2*Nr )
      PARAMETER( numDiags = 25+500*Nr )
      PARAMETER( nRegions = 0 , sizRegMsk = 1 , nStats = 4 )
      PARAMETER( diagSt_size = 100*Nr )

CEH3 ;;; Local Variables: ***
CEH3 ;;; mode:fortran ***
CEH3 ;;; End: ***
