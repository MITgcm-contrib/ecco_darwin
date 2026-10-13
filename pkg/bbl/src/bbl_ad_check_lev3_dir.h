C     BBL state carried between time steps (pkg/bbl)
#ifdef ALLOW_BBL
CADJ STORE bbl_eta   = tapelev3, key = ilev_3
CADJ STORE bbl_theta = tapelev3, key = ilev_3
CADJ STORE bbl_salt  = tapelev3, key = ilev_3
# ifdef ALLOW_PTRACERS
CADJ STORE bbl_ptr   = tapelev3, key = ilev_3
# endif
#endif /* ALLOW_BBL */
