C     BBL state carried between time steps (pkg/bbl)
#ifdef ALLOW_BBL
CADJ STORE bbl_eta   = tapelev2, key = ilev_2
CADJ STORE bbl_theta = tapelev2, key = ilev_2
CADJ STORE bbl_salt  = tapelev2, key = ilev_2
# ifdef ALLOW_PTRACERS
CADJ STORE bbl_ptr   = tapelev2, key = ilev_2
# endif
#endif /* ALLOW_BBL */
