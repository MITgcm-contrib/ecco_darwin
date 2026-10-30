C     LR17 river point sources (data.rivers, namelist RIVERS_PARM):
C     rivN          :: number of rivers
C     rivI, rivJ    :: global indices of the source cell of each river
C     rivNrec       :: number of discharge records
C     rivT0, rivDt  :: model time of record 1, record interval [s]
C     rivQ(r,n)     :: discharge of river r at record n [m^3/s]
C                      (linear in time between records, held at the ends)
C     rivHalf       :: half-width of the block of cells that shares each
C                      river's discharge (0: the source cell only; 1: 3x3)
C     rivCells(r)   :: number of wet cells in river r's block (computed)
C     The discharge enters level kSurfC of the block's wet cells as
C     addMass (selectAddFluid=1), shared equally, at temp_addMass,
C     salt_addMass (PARM01).
      INTEGER rivMax, rivRecMax
      PARAMETER ( rivMax = 10, rivRecMax = 400 )
      INTEGER rivN, rivNrec, rivI(rivMax), rivJ(rivMax), rivHalf
      _RL rivT0, rivDt, rivQ(rivMax,rivRecMax), rivCells(rivMax)
      COMMON /LR17_RIV_I/ rivN, rivNrec, rivI, rivJ, rivHalf
      COMMON /LR17_RIV_R/ rivT0, rivDt, rivQ, rivCells
