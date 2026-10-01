      INTEGER PTRACERS_num
      PARAMETER(PTRACERS_num = 4 )
#ifdef ALLOW_AUTODIFF_TAMC
      INTEGER    maxpass
      PARAMETER( maxpass     = PTRACERS_num + 2 )
#endif
