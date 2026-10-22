C     BBL state carried between time steps (pkg/bbl)
#ifdef ALLOW_BBL
CADJ STORE bbl_eta   = tapelev4, key = ilev_4
CADJ STORE bbl_theta = tapelev4, key = ilev_4
CADJ STORE bbl_salt  = tapelev4, key = ilev_4
CADJ STORE bbl_rho_nr = tapelev4, key = ilev_4
# ifdef ALLOW_PTRACERS
CADJ STORE bbl_ptr   = tapelev4, key = ilev_4
# endif
#endif /* ALLOW_BBL */
