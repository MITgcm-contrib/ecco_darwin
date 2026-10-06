# tools: genmake2 and testreport

## genmake2 -help
```
GENMAKE :

A program for GENerating MAKEfiles for the MITgcm project.
   For a quick list of options, use "genmake2 -h"
or for more detail see the documentation, section "Building the model"
   (under "Getting Started") at:  https://mitgcm.readthedocs.io/

===  Processing options files and arguments  ===

Usage: "tools/genmake2" [OPTIONS]
  where [OPTIONS] can be:

    -help | --help | -h | --h
	  Print this help message and exit.

    -tap | --tap
	  Generate a Makefile for a Tapenade build
    -tap_extra NAME | --tap_extra NAME | -tap_extra=NAME | --tap_extra=NAME
	  Here, "NAME" specifies a list of extra arguments to pass to the
	  Tapenade command.

    -oad | --oad
	  Generate a Makefile for an OpenAD built
    -oadsingularity NAME | --oadsingularity NAME | -oadsngl NAME | --oadsngl NAME
      -oadsingularity=NAME | --oadsingularity=NAME | -oadsngl=NAME | --oadsngl=NAME
	  Here, "NAME" specifies the singularity file
	  that contains the OpenAD execulable.

    -nocat4ad | -dog4ad | -ncad | -dad
	  do not concatenate (cat) source code sent to TAF
	  resulting in compilation of multiple files

    -adoptfile NAME | --adoptfile NAME | -adof NAME | --adof NAME
      -adoptfile=NAME | --adoptfile=NAME | -adof=NAME | --adof=NAME
	  Use "NAME" as the adoptfile.  By default, use file
	  "tools/adjoint_options/adjoint_default" or
	  "tools/adjoint_options/adjoint_tap" (for Tapenade built) or
	  "tools/adjoint_options/adjoint_oad" (for OpenAD built).

    -optfile NAME | --optfile NAME | -of NAME | --of NAME
      -optfile=NAME | --optfile=NAME | -of=NAME | --of=NAME
	  Use "NAME" as the optfile.  By default, an attempt will be
	  made to find an appropriate "standard" optfile in the
	  tools/build_options/ directory.

    -pdepend NAME | --pdepend NAME
      -pdepend=NAME | --pdepend=NAME
	  Get package dependency information from "NAME".

    -pgroups NAME | --pgroups NAME
      -pgroups=NAME | --pgroups=NAME
	  Get the package groups information from "NAME".

    -bash NAME
	  Explicitly specify the Bourne or BASH shell to use

    -make NAME | -m NAME
      --make=NAME | -m=NAME
	  Use "NAME" for the MAKE program. The default is "make" but
	  many platforms, "gmake" is the preferred choice.

    -makefile NAME | -mf NAME
      --makefile=NAME | -mf=NAME
	  Call the makefile "NAME".  The default is "Makefile".

    -makedepend NAME | -md NAME
      --makedepend=NAME | -md=NAME
	  Use "NAME" for the MAKEDEPEND program.

    -rootdir NAME | --rootdir NAME | -rd NAME | --rd NAME
      -rootdir=NAME | --rootdir=NAME | -rd=NAME | --rd=NAME
	  Specify the location of the MITgcm root directory as 'NAME'.
	  By default, genmake will try to find the location by
	  looking in parent directories (up to the 5th parent).

    -mods NAME | --mods NAME | -mo NAME | --mo NAME
      -mods=NAME | --mods=NAME | -mo=NAME | --mo=NAME
	  Here, "NAME" specifies a list of directories that are used for
	  additional source code.  Files found in the "mods list" are given
	  preference over files of the same name found elsewhere.

    -disable NAME | --disable NAME
      -disable=NAME | --disable=NAME
	  Here "NAME" specifies a list of packages that we don't
	  want to use.  If this violates package dependencies,
	  genmake will exit with an error message.

    -enable NAME | --enable NAME
      -enable=NAME | --enable=NAME
	  Here "NAME" specifies a list of packages that we wish
	  to specifically enable.  If this violates package
	  dependencies, genmake will exit with an error message.

    -standarddirs NAME | --standarddirs NAME
      -standarddirs=NAME | --standarddirs=NAME
	  Here, "NAME" specifies a list of directories to be
	  used as the "standard" code.

    -fortran NAME | --fortran NAME | -fc NAME | --fc NAME
      -fc=NAME | --fc=NAME
	  Use "NAME" as the fortran compiler.  By default, genmake
	  will search for a working compiler by trying a list of
	  "usual suspects" such as g77, f77, etc.

    -cc NAME | --cc NAME | -cc=NAME | --cc=NAME
	  Use "NAME" as the C compiler.  By default, genmake
	  will search for a working compiler by trying a list of
	  "usual suspects" such as gcc, c89, cc, etc.

    -use_real4 | -use_r4 | -ur4 | --use_real4 | --use_r4 | --ur4
	  Use "real*4" type for _RS variable (#undef REAL4_IS_SLOW)
	  *only* works if CPP_EEOPTIONS.h allows this.

    -ignoretime | -ignore_time | --ignoretime | --ignore_time
	  Ignore all the "wall clock" routines entirely.  This will
	  not in any way hurt the model results -- it simply means
	  that the code that checks how long the model spends in
	  various routines will give junk values.

    -ts | --ts
	  Produce timing information per timestep
    -papis | --papis
	  Produce summary MFlop/s (and IPC) with PAPI per timestep
    -pcls | --pcls
	  Produce summary MFlop/s etc. with PCL per timestep
    -foolad | --foolad
	  Fool the AD code generator
    -papi | --papi
	  Performance analysis with PAPI
    -pcl | --pcl
	  Performance analysis with PCL
    -hpmt | --hpmt
	  Performance analysis with the HPM Toolkit

    -ieee | --ieee
	  use IEEE numerics.  Note that this option *only* works
	  if it is supported by the OPTFILE that is being used.
    -devel | --devel
	  Add additional warning and debugging flags for development
	  (if supported by the OPTFILE); also switch to IEEE numerics.
    -gsl | --gsl
	  Use GSL to control floating point rounding and precision

    -mpi | --mpi
	  Include MPI header files and link to MPI libraries
    -mpi=PATH | --mpi=PATH
	  Include MPI header files and link to MPI libraries using MPI_ROOT
	  set to PATH. i.e. Include files from $PATH/include, link to libraries
	  from $PATH/lib and use binaries from $PATH/bin.

    -omp | --omp
	  Activate OpenMP code
    -omp=OMPFLAG | --omp=OMPFLAG
	  Activate OpenMP code + use Compiler option OMPFLAG

    -es | --es | -embed-source | --embed-source
	  Embed a tarball containing the full source code (including the
	  Makefile, etc.) used to build the executable [off by default]

    -ds | --ds
	  Report genmake internal variables status (DUMPSTATE)
	  to file "genmake_state" (for debug purpose)

  While it is most often a single word, the "NAME" variables specified
  above can in many cases be a space-delimited string such as:

    --enable pkg1   --enable 'pkg1 pkg2'   --enable 'pkg1 pkg2 pkg3'
    -mods=dir1   -mods='dir1'   -mods='dir1 dir2 dir3'
    -foptim='-Mvect=cachesize:512000,transform -xtypemap=real:64,double:64,integer:32'

  which, depending upon your shell, may need to be single-quoted.

  For more detailed genmake documentation, please see section "Building the model"
    (under "Getting Started") at:  https://mitgcm.readthedocs.io/
```

## testreport -help
```
parsing options...  
Usage:  ./testreport [OPTIONS]

where possible OPTIONS are:
  (-help|-h)               print usage
 ---- type of test : ----
  (-tlm)                   perform a Tangent-Linear run (defaut: using TAF)
  (-adm|-ad)               perform an Adjoint run (default: using TAF)
  (-tap)                   use Tapenade for Adjoint or Tangent-Linear run
  (-oad)                   perform an OpenAD adjoint run
  (-oadsingularity|-oadsngl) STRING
                           path to singularity container with OpenAD
  (-mth)                   run multi-threaded (using eedata.mth)
  (-mpi)                   use MPI to compile and run on 2 processors
  (-MPI)  NUMBER           use MPI to compile and run on max NUMBER procs
  (-mfile|-mf) STRING      MPI: file with list of possible machines to run on
  (-command|-c) STRING     command to run (e.g., if non-standard MPI setting)
                            DEF='mitgcmuv' or ='mpirun -np TR_NPROC mitgcmuv'
 ---- testing options : ----
  (-optfile|-of) STRING    optfile to use
  (-fast)                  use optfile default for compiler flags (no '-ieee')
                            DEF=off => use IEEE numerics option (if available)
  (-devel)                 use optfile developement flags (if available)
  (-gsl)                   compile with "-gsl" flag
  (-use_r4|-ur4)           if allowed, use real*4 type for '_RS' variable
  (-tdir|-t) STRING        list of group and/or exp. dirs to test
                             (recognized groups: basic, tutorials)
                             (DEF="" which test all)
                             (if list= 'start_from THIS_EXP' then
                              test THIS_EXP + all the following)
  (-skipdir|-skd) STRING   list of exp. dirs to skip
                             (DEF="" which test all)
  (-ts)                    provide timing information per timestep
  (-papis)                 provide MFlop/s per timestep using PAPI
  (-pcls)                  provide MFlop/s per timestep using PCL
 ---- system options : ----
  (-bash|-b) STRING        preferred location of a "bash" or "sh" shell
                             (DEF="" for "bash")
  (-ef) STRING             used as genmake2 "-extra_flag" argument
  (-ncad)                  use genmake2 option "-nocat4ad" (-ncad)
  (-small_f)               make target small_f before making target all
  (-makedepend|-md) STRING command to use for "makedepend"
  (-make|-m) STRING        command to use for "make"
                             (DEF="make")
  (-j) JOBS                use "make -j JOBS" for parallel builds
 ---- output options : ----
  (-match) NUMBER          Matching Criteria (number of digits)
                             (DEF="10")
  (-pass)                  return non-zero exit code if any exp do not pass
  (-odir) STRING           used to build output directory name
                             (DEF="hostname")
  (-addr|-a) STRING        list of email recipients
                             (DEF="" no email is sent)
  (-send)       STRING     sending command (instead of locally built mpack)
  (-savdir|-sd) STRING     location to save output tar file to send (DEF='.')
  (-mpackdir|-mpd) DIR     location of the mpack utility
                             (DEF='../tools/mpack-1.6')
 -- do only some parts: --
  (-clean)                 *ONLY* run "make CLEAN" & clean run-dir
  (-norun|-nr)             skip the "runmodel" stage (stop after make)
  (-obj)                   only produces objects (=norun & no executable)
  (-src)                   only produces small '*.f' src files (not even obj)
                            + with: '-adm/-tlm', also makes taf outp src code
  (-runonly|-ro)           *ONLY* run stage (="-quick" without make)
  (-quick|-q)              same as "-nogenmake -noclean -nodepend"
  (-nogenmake|-ng)         skip the genmake stage
  (-noclean|-nc)           skip the "make clean" stage
  (-nodepend|-nd)          skip the "make depend" stage
  (-postclean|-pc)         after each exp. test, clean build-dir & run-dir
  (-deloutp|-do)           delete output files after successful run
  (-deldir|-dd)            on success, delete the output directory

and where STRING can be a whitespace-delimited list
such as:

  -t 'exp0 exp2 exp3' 
  -addr='abc@123.com testing@home.org'

provided that the expression is properly quoted within the current
shell (note the use of single quotes to protect white space).
```

## Optfiles shipped (tools/build_options)
`SUPER-UX_SX-6_f90+mpi_caspur` `SUPER-UX_SX-6_sx90_dkrz` `SUPER-UX_SX-ACE_sxf90_awi` `cygwin_ia32_g77` `darwin_absoft_f77` `darwin_amd64_gfortran` `darwin_arm64_gfortran` `darwin_ia32_g95` `darwin_ia32_gfortran` `darwin_ia32_ifort` `darwin_ia32_pgf95_trane` `darwin_ppc_f95` `darwin_ppc_g77` `darwin_ppc_xlf` `darwin_ppc_xlf_panther` `darwin_ppc_xlf_panther+wienders` `darwin_ppc_xlf_panther_baylor` `darwin_ppc_xlf_tiger_baylor` `linux_alpha_g77` `linux_amd64_absoft` `linux_amd64_g77` `linux_amd64_g95` `linux_amd64_gfortran` `linux_amd64_gfortran_albedo` `linux_amd64_gfortran_avx2` `linux_amd64_gfortran_greenplanet` `linux_amd64_ifort` `linux_amd64_ifort+gcc` `linux_amd64_ifort+impi` `linux_amd64_ifort+impi_engaging` `linux_amd64_ifort+impi_stampede2_knl` `linux_amd64_ifort+impi_stampede2_skx` `linux_amd64_ifort+mpi_ice_nas` `linux_amd64_ifort+mpi_sal_oxford` `linux_amd64_ifort+mpi_yellowstone` `linux_amd64_ifort10-` `linux_amd64_ifort11` `linux_amd64_ifort_albedo` `linux_amd64_ifort_beagle` `linux_amd64_ifort_discover` `linux_amd64_ifort_fimm_default` `linux_amd64_ifort_fimm_emic` `linux_amd64_ifort_uv100` `linux_amd64_ifx+impi_avx2` `linux_amd64_open64` `linux_amd64_pathf90` `linux_amd64_pgf77` `linux_amd64_pgf77+mpi_ncar` `linux_amd64_pgf77+mpi_xd1` `linux_amd64_pgf77_ocl` `linux_amd64_pgf90+mpi_greenplanet` `linux_amd64_pgf90+mpi_xd1` `linux_amd64_sunf90` `linux_arm64_gfortran` `linux_ia32_absoft` `linux_ia32_g77` `linux_ia32_g95` `linux_ia32_gfortran` `linux_ia32_gfortran+mpi_fc_lam` `linux_ia32_ifort` `linux_ia32_ifort10.1` `linux_ia32_ifort11` `linux_ia32_lf95` `linux_ia32_open64` `linux_ia32_pathf90` `linux_ia32_pgf77` `linux_ia32_pgf77+mpi_aer` `linux_ia32_sunf90` `linux_ia64_cray_archer2` `linux_ia64_cray_cca` `linux_ia64_cray_ollie` `linux_ia64_efc` `linux_ia64_g77` `linux_ia64_ifort` `linux_ia64_ifort+mpi_altix_nas` `linux_ia64_ifort+mpi_swell` `linux_ia64_ifort9.0+mpt_altix3700_stommel` `linux_ia64_ifort_altix_jpl` `linux_ia64_ifort_ollie` `linux_ia64_open64` `linux_ia64_pgf77+mpi_cray_xt3_jaguar` `linux_ppc64_xlf` `linux_ppc_xlf` `sp5+mpi` `sp5+mpi_nas` `sp6+mpi_iblade` `sp6_ncar` `sunos_amd64_f77_awi` `sunos_i86pc_f95` `sunos_sparc_sunf90` `sunos_sparc_sunf90_m64` `sunos_sun4u_f77` `sunos_sun4u_f90` `sunos_sun4u_g77` `sunos_sun4u_mpf77+mpi_sunfire` `unsupported`
