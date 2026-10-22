C     BBL state carried between time steps (pkg/bbl)
#ifdef ALLOW_BBL
CADJ STORE bbl_eta   = comlev1, key = ikey_dynamics, kind = isbyte
CADJ STORE bbl_theta = comlev1, key = ikey_dynamics, kind = isbyte
CADJ STORE bbl_salt  = comlev1, key = ikey_dynamics, kind = isbyte
CADJ STORE bbl_rho_nr = comlev1, key = ikey_dynamics, kind = isbyte
# ifdef ALLOW_PTRACERS
CADJ STORE bbl_ptr   = comlev1, key = ikey_dynamics, kind = isbyte
# endif
#endif /* ALLOW_BBL */
