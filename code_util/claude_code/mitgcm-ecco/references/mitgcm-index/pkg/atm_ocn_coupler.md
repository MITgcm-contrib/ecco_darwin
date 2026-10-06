# pkg/atm_ocn_coupler

Stand-alone coupler component exchanging fields between atmosphere and ocean MITgcm components.

**runtime switch:** `useATM_OCN_COUPLER`-style flag in `data.pkg` (check exact name in packages_boot.F)

## README
```

***************************************************************************************
|                               Initialisation                                        |
***************************************************************************************

MITCOMPONENT_init
=================

      CALL MITCOMPONENT_init( 
     I                  name, 
     O                  comm )
      
      name - Name of component to register e.g.
             'ocean', 'atmos', 'ice', 'land'. Up to
             MAX_COMPONENTS ( see "CPLR_SIG.h" )
             are allowed. Default is MAX_COMPONENTS = 10.

      comm - MPI_Communicator which includes all procs.
             that registered as component type 'name.

      Each process can only register as one component type.

      Every process has to call MITCOMPONENT_init at the
      same point otherwise everything deadlocks or dies.
      Except the coupler process which calls MITCOUPLER_init.

      Initialises the MPI context for a particular component
      model. On return the component is given a communicator
      that can be used for communication based on MPI
      for the processes of the component. Internal 
```

## Namelist parameters
### COUPLER_PARAMS
- `cpl_sequential` — =0/1 : selects Synchronous/Sequential Coupling
- `cpl_exchange_RunOff` — controls exchange of RunOff fields
- `cpl_exchange1W_sIce` — controls 1-way exchange of seaice (step fwd in ATM)
- `cpl_exchange2W_sIce` — controls 2-way exchange of ThSIce variables
- `cpl_exchange_SaltPl` — controls exchange of Salt-Plume fields
- `cpl_exchange_DIC` — controls exchange of DIC variables
- `runOffMapSize` — Nunber of connected grid points in RunOff-Map
- `runOffMapFile` — Input file for setting runoffmap

## Headers
- `ATMIDS.h` — are used to identify this component and the fields it exchanges with other components.
- `ATMVARS.h` — grid. Arrays may need adding or removing different couplings.
- `CPLIDS.h` — /==========================================================\ are used to identify this component and the fields it exchanges with other components. No
- `CPL_MAP2GRIDS.h` — Declare arrays used for mapping coupling fields from one grid (atmos., ocean) to the other grid
- `CPL_PARAMS.h` — - parameter for the Coupler, holds in common block
- `CPP_EEOPTIONS.h` — BOP
- `CPP_OPTIONS.h` — BOP
- `OCNIDS.h` — are used to identify this component and the fields it exchanges with other components.
- `OCNVARS.h` — grid. Arrays may need adding or removing different couplings.

## Routines (25)
`accept_component_registrations.F`, `atm_to_ocn_maprunoff.F`, `atm_to_ocn_mapxyr8.F`, `cpl_check_cplconfig.F`, `cpl_init_atm_vars.F`, `cpl_init_ocn_vars.F`, `cpl_read_params.F`, `cpl_recv_atm_atmconfig.F`, `cpl_recv_atm_fields.F`, `cpl_recv_ocn_fields.F`, `cpl_recv_ocn_ocnconfig.F`, `cpl_register_atm.F`, `cpl_register_ocn.F`, `cpl_send_atm_cplparms.F`, `cpl_send_atm_fields.F`, `cpl_send_atm_ocnconfig.F`, `cpl_send_ocn_atmconfig.F`, `cpl_send_ocn_cplparms.F`, `cpl_send_ocn_fields.F`, `exch_component_configs.F`, `initialise.F`, `mds_byteswap.F`, `ocn_to_atm_mapxyr8.F`, `set_runoffmap.F`

## Called from outside the package
- `CPL_REGISTER_OCN` ← `pkg/atm2d/accept_component_registrations.F:42`
- `ACCEPT_COMPONENT_REGISTRATIONS` ← `pkg/atm2d/atm2d_init_fixed.F:89`
- `EXCH_COMPONENT_CONFIGS` ← `pkg/atm2d/atm2d_init_fixed.F:94`
- `INITIALISE` ← `pkg/atm2d/atm2d_init_fixed.F:86`
- `CPL_RECV_OCN_OCNCONFIG` ← `pkg/atm2d/exch_component_configs.F:42`
- `CPL_SEND_OCN_ATMCONFIG` ← `pkg/atm2d/exch_component_configs.F:50`
- `CPL_RECV_OCN_FIELDS` ← `pkg/atm2d/forward_step_atm2d.F:198`
- `CPL_SEND_OCN_FIELDS` ← `pkg/atm2d/forward_step_atm2d.F:218`
- `MDS_BYTESWAPR8` ← `pkg/dic/dic_read_co2_pickup.F:73`
- `MDS_BYTESWAPR4` ← `pkg/exf/exf_interp_read.F:142`
- `MDS_BYTESWAPR8` ← `pkg/exf/exf_interp_read.F:151`
- `MDS_BYTESWAPR4` ← `pkg/fizhi/fizhi_init_chem.F:171`
- `MDS_BYTESWAPR4` ← `pkg/fizhi/fizhi_init_veg.F:124`
- `MDS_BYTESWAPR8` ← `pkg/fizhi/fizhi_init_vegsurftiles.F:77`
- `MDS_BYTESWAPR8` ← `pkg/fizhi/fizhi_readwrite_vegtiles.F:131,139,147,155`
- `MDS_BYTESWAPR4` ← `pkg/fizhi/update_ocean_exports.F:657,658,721,722`
- `MDS_BYTESWAPR4` ← `pkg/mdsio/mdsio_facef_read.F:86,121`
- `MDS_BYTESWAPR8` ← `pkg/mdsio/mdsio_facef_read.F:97,131`
- `MDS_BYTESWAPR4` ← `pkg/mdsio/mdsio_gl.F:278,358,610,710`
- `MDS_BYTESWAPR8` ← `pkg/mdsio/mdsio_gl.F:295,368,620,727`
- `MDS_BYTESWAPR4` ← `pkg/mdsio/mdsio_gl_slice.F:181,421,655,896`
- `MDS_BYTESWAPR8` ← `pkg/mdsio/mdsio_gl_slice.F:198,438,672,913`
- `MDS_BYTESWAPR4` ← `pkg/mdsio/mdsio_rd_rec_rl.F:57`
- `MDS_BYTESWAPR8` ← `pkg/mdsio/mdsio_rd_rec_rl.F:65`
- `MDS_BYTESWAPR4` ← `pkg/mdsio/mdsio_rd_rec_rs.F:57`
- `MDS_BYTESWAPR8` ← `pkg/mdsio/mdsio_rd_rec_rs.F:65`
- `MDS_BYTESWAPR4` ← `pkg/mdsio/mdsio_read_field.F:296,513`
- `MDS_BYTESWAPR8` ← `pkg/mdsio/mdsio_read_field.F:301,515`
- `MDS_BYTESWAPR4` ← `pkg/mdsio/mdsio_read_section.F:232,517`
- `MDS_BYTESWAPR8` ← `pkg/mdsio/mdsio_read_section.F:249,534`
- `MDS_BYTESWAPR4` ← `pkg/mdsio/mdsio_read_tape.F:153,239`
- `MDS_BYTESWAPR8` ← `pkg/mdsio/mdsio_read_tape.F:158,244`
- `MDS_BYTESWAPR4` ← `pkg/mdsio/mdsio_wr_rec_rl.F:59`
- `MDS_BYTESWAPR8` ← `pkg/mdsio/mdsio_wr_rec_rl.F:67`
- `MDS_BYTESWAPR4` ← `pkg/mdsio/mdsio_wr_rec_rs.F:59`
- `MDS_BYTESWAPR8` ← `pkg/mdsio/mdsio_wr_rec_rs.F:67`
- `MDS_BYTESWAPR4` ← `pkg/mdsio/mdsio_write_field.F:345,417`
- `MDS_BYTESWAPR8` ← `pkg/mdsio/mdsio_write_field.F:350,419`
- `MDS_BYTESWAPR4` ← `pkg/mdsio/mdsio_write_section.F:214,474`
- `MDS_BYTESWAPR8` ← `pkg/mdsio/mdsio_write_section.F:231,491`
- `MDS_BYTESWAPR4` ← `pkg/mdsio/mdsio_write_tape.F:190,268`
- `MDS_BYTESWAPR8` ← `pkg/mdsio/mdsio_write_tape.F:195,273`
- `MDS_BYTESWAPR4` ← `pkg/mdsio/mdsio_writelocal.F:291`
- `MDS_BYTESWAPR8` ← `pkg/mdsio/mdsio_writelocal.F:293`
- `MDS_BYTESWAPR8` ← `pkg/obsfit/obsfit_active_file_control.F:117,122,165,175`
- `MDS_BYTESWAPR8` ← `pkg/obsfit/obsfit_init_equifiles.F:121`
- `MDS_BYTESWAPR8` ← `pkg/profiles/active_file_control_profiles.F:127,135,183,199`
- `MDS_BYTESWAPR8` ← `pkg/profiles/profiles_init_ncfile.F:138`

## Verification experiments compiling it (1)
`cpl_aim+ocn`
