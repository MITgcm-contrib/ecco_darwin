   Run with Freshwater Fluxes
   1) In code dir add the modified: <obcs_exf_load.F>  <obcs_calc.F> <OBCS_PARAMS.h> <OBCS_FIELDS.h> <obcs_readparms.F> <obcs_balance_flow.F>
   2) Compile
   3) In datas.obcs add: OBEetaFile_tr = '/project/k1254/vasoup/Glorys12V1_eta_tr_AG_19930101-20181225.bin',
                         useOBCSbalance=.TRUE.,
                         OBCS_balanceFacE=-1,
                         and use OBEuFile = Glorys12V1_UobcsE_1km_19930101-20171225.bin'
