import os
import argparse
import numpy as np
import xarray as xr
import simplegrid as sg
from MITgcmutils import mds, llc
#### Interpolation
from scipy.spatial import Delaunay
from scipy.ndimage import gaussian_filter
from scipy.interpolate import LinearNDInterpolator, interp1d

############################################################
#               INTERPOLATION FUNCTIONS                    #
############################################################
def get_coast_indexes(WCmask0):
    ### Get parameters
    coast_map = np.zeros(WCmask0.shape)
    llc = WCmask0.shape[-1]
    extN = int(llc/90) #number of coast cell to extrapolate
    for k in range(WCmask0.shape[0]):
        id0 = np.where(WCmask0[k]==0)
        for i in range(len(id0[0])):
            if id0[0][i]==0 and id0[1][i]==0:
                pass
            if id0[0][i] in [0,1,3*llc,3*llc+1,7*llc-1,7*llc-2,10*llc-1,10*llc-2]:
                pass
            else:
                if np.any(WCmask0[k,id0[0][i]-extN:id0[0][i]+extN+1,id0[1][i]-extN:id0[1][i]+extN+1])==1:
                    coast_map[k,id0[0][i],id0[1][i]]=1
    return coast_map

def Hinterp(data0, lmask0, coast0, Dtri, XC1, YC1, sigmaG):
    llc = lmask0.shape[-1]
    if len(data0.shape) == 3:
        data1 = np.zeros((data0.shape[0],XC1.shape[0],XC1.shape[1]))
        for k in range(data0.shape[0]):
            # Get index of coast at the depth k
            id_coast = np.where(coast0[k]==1)
            # Extrapolate data to the coast
            tmp = data0[k].copy()
            tmp[lmask0[k]==0] = np.nan
            for i in range(len(id_coast[0])):
                tmp[id_coast[0][i],id_coast[1][i]] = np.nanmean(tmp[id_coast[0][i]-3:id_coast[0][i]+4,
                                                                id_coast[1][i]-3:id_coast[1][i]+4])
            # Downscal data (barycentric interpolation + Gaussian filter)
            tmp[np.isnan(tmp)] = 0
            tranform = LinearNDInterpolator(Dtri, tmp.ravel())
            data1[k] = tranform(XC1, YC1)
            data1[k] = gaussian_filter(data1[k], sigma=sigmaG)
    else:
        data1 = np.zeros(XC1.shape)
        # Get index of coast at the surface
        id_coast = np.where(coast0[0]==1)
        # Extrapolate data to the coast
        tmp = data0.copy()
        tmp[lmask0[0]==0] = np.nan
        for i in range(len(id_coast[0])):
            tmp[id_coast[0][i],id_coast[1][i]] = np.nanmean(tmp[id_coast[0][i]-3:id_coast[0][i]+4,
                                                                id_coast[1][i]-3:id_coast[1][i]+4])
        # Downscal data (barycentric interpolation + Gaussian filter)
        tmp[np.isnan(tmp)] = 0
        tranform = LinearNDInterpolator(Dtri, tmp.ravel())
        data1 = tranform(XC1, YC1)
        data1 = gaussian_filter(data1, sigma=sigmaG)
    return data1

def Zinterp(data0, drF0, drF1, HFacC1):
    ### Initialization
    depth0 = np.cumsum(drF0)
    depth1 = np.cumsum(drF1)
    lmask1 = HFacC1.copy(); lmask1[lmask1>0] = 1
    data1 = np.zeros(HFacC1.shape)
    ## Interpolate vertical profiles
    for j in range(data0.shape[1]):
        for i in range(data0.shape[2]):
            if HFacC1[0,j,i] == 0:
                pass
            else:
                id_bot1 = np.where(HFacC1[:,j,i]>0)[0][-1]
                id_bot0 = np.where(depth0>=depth1[id_bot1])[0][1]
                interp = interp1d(depth0[:id_bot0+1], data0[:id_bot0+1,j,i], fill_value='extrapolate')
                data1[:id_bot1+1,j,i] = interp(depth1[:id_bot1+1])
    data1 *= lmask1 
    return(data1)

############################################################
#                  GRID WORK FUNCTIONS                     #
############################################################

def transp_tiles(data):
    nx = data.shape[1]
    ny = data.shape[0]
    tmp = data[7*nx:,::-1]
    transpo = np.concatenate((tmp[2::3,:].transpose(),tmp[1::3,:].transpose(),tmp[0::3,:].transpose()))
    data_out = np.concatenate((data[:7*nx],np.flipud(transpo[:,:nx]),np.flipud(transpo[:,nx:])))
    return data_out

def UVrot(uvel, vvel, angle_cos, angle_sin):
    zvel = np.zeros_like(uvel)
    mvel = np.zeros_like(vvel)
    uvelC, vvelC = llc.uv2c(uvel,vvel)
    if len(uvelC.shape)==2:
        zvel = angle_cos * uvelC - angle_sin * vvelC
        mvel = angle_sin * uvelC + angle_cos * vvelC
    else:    
        for k in range(np.shape(uvel)[0]):
            zvel[k,:,:] = angle_cos * uvelC[k,:,:] - angle_sin * vvelC[k,:,:]
            mvel[k,:,:] = angle_sin * uvelC[k,:,:] + angle_cos * vvelC[k,:,:]
    return zvel, mvel

############################################################
#                READING/WRITING FUNCTIONS                 #
############################################################

def read_eccogrid(config_dir):
    grid_dir = os.path.join(config_dir, 'parent/outputs/grid/')
    XC = transp_tiles(mds.rdmds(grid_dir+'XC'))
    YC = transp_tiles(mds.rdmds(grid_dir+'YC'))
    dXC = transp_tiles(mds.rdmds(grid_dir+'DXC'))
    dYC = transp_tiles(mds.rdmds(grid_dir+'DYC'))
    AngleCS = transp_tiles(mds.rdmds(grid_dir+'AngleCS'))
    AngleSN = transp_tiles(mds.rdmds(grid_dir+'AngleSN'))
    tmp = mds.rdmds(grid_dir+'hFacC')
    HFacC = np.zeros(tmp.shape)
    for i in range(len(HFacC)):
        HFacC[i] = transp_tiles(tmp[i])
    delR = np.asanyarray(mds.rdmds(grid_dir+'DRF')[:,0,0], dtype=np.float32)
    return XC, YC, dXC, dYC, AngleCS, AngleSN, HFacC, delR

def read_ncgrid(config_dir, model_name, var_ls):
    ds = xr.open_dataset(os.path.join(config_dir,model_name+'_ncgrid.nc'))
    var_mats = []
    for var in var_ls:
        var_mats.append(ds[var].values)
    ds.close()
    return var_mats

def gen_nc(XC1, YC1, drF1, Nr1, pickup_phy, pickup_sic, pickup_ggl90, pickup_darwin, pickup_ptracers):
    #### Generate xarray with forcings
    dims = dict(X=np.arange(XC1.shape[1]), Y=np.arange(XC1.shape[0]), 
                Z=np.arange(drF1.shape[0]),TRC=np.arange(pickup_ptracers.shape[0]))
    data = dict(XC=(["Y","X"],XC1), YC=(["Y","X"],YC1), drF=(["Z"],drF1),
                U=(["Z","Y","X"],pickup_phy[:Nr1]), V=(["Z","Y","X"],pickup_phy[Nr1:2*Nr1]),
                THETA=(["Z","Y","X"],pickup_phy[2*Nr1:3*Nr1]), SALT=(["Z","Y","X"],pickup_phy[3*Nr1:4*Nr1]),
                ETAN=(["Y","X"],pickup_phy[-1]), siTICE=(["Y","X"],pickup_sic[0]), siAREA=(["Y","X"],pickup_sic[1]),
                siHEFF=(["Y","X"],pickup_sic[2]), siHSNOW=(["Y","X"],pickup_sic[3]), siUICE=(["Y","X"],pickup_sic[4]),
                siVICE=(["Y","X"],pickup_sic[5]), ggl90=(["Z","Y","X"],pickup_ggl90), pH=(["Z","Y","X"],pickup_darwin),
                tracers=(["TRC","Z","Y","X"],pickup_ptracers))
    ds = xr.Dataset(data_vars=data, coords=dims)
    #### Generate netcdf file
    output_dir = os.path.join(config_dir, 'forcings/pickups/')
    ds.to_netcdf(output_dir+'pickups_check.nc', mode='w')
    return

############################################################
#                 INIT CONDITION FUNCTIONS                 #
############################################################

############## PHYSICS ##############
def gen_pickup(config_dir, pickup_itr, Nr0, AngleCS0, AngleSN0, WCmask0, drF0, coast, 
               XC1, YC1, Nr1, drF1, HFacC1, Dtri, sigmaG, print_level):
    #### Read pickup file
    pickup_dir = os.path.join(config_dir, 'parent/outputs/pickups/')
    if print_level>=1:
        print('        > Reading in the ECCO pickup file')
    pickup0, _, metap = mds.rdmds(os.path.join(pickup_dir,'pickup'), itrs=pickup_itr, returnmeta=True)
    for i in range(len(pickup0)):
        pickup0[i] = transp_tiles(pickup0[i])
    #### Center & Rotate velocity fields
    if print_level>=1:
        print('        > Rotating orignal velocity fields into regional coordinates')
    Uid = metap['fldlist'].index('Uvel')
    Vid = metap['fldlist'].index('Vvel')
    Zvel, Mvel = UVrot(pickup0[Uid*Nr0:(Uid+1)*Nr0], pickup0[Vid*Nr0:(Vid+1)*Nr0], AngleCS0, AngleSN0)
    #### Interpolate ICs on the regional grid
    ICvar = ['Zvel', 'Mvel','Theta','Salt']
    pickup1 = np.zeros((Nr1*len(ICvar)+1,XC1.shape[0],XC1.shape[1]))
    if print_level>=1:
        print('        > Downscaling initial conditions')
    for i, var in enumerate(ICvar):
        if var in metap['fldlist']:
            vid = metap['fldlist'].index(var)
            tmp = Hinterp(pickup0[vid*Nr0:(vid+1)*Nr0], WCmask0, coast, Dtri, XC1, YC1, sigmaG)
            pickup1[i*Nr1:(i+1)*Nr1] = Zinterp(tmp, drF0, drF1, HFacC1)
        else:
            tmp = Hinterp(locals()[var], WCmask0, coast, Dtri, XC1, YC1, sigmaG)
            pickup1[i*Nr1:(i+1)*Nr1] = Zinterp(tmp, drF0, drF1, HFacC1)
    vid = metap['fldlist'].index('EtaN')
    pickup1[-1] = Hinterp(pickup0[vid*Nr0], WCmask0, coast, Dtri, XC1, YC1, sigmaG)*HFacC1[0]
    #### Generate new pickup files#### Generate new pickup files
    out_pckf = ['pickup_U','pickup_V','pickup_THETA','pickup_SALT','pickup_ETAN']
    output_dir = os.path.join(config_dir, 'forcings/pickups/')
    if print_level>=1:
        print('        > Generating sea-ice initial conditions file')
    for i in range(len(out_pckf)-1):
        data = pickup1[i*Nr1:(i+1)*Nr1]
        mds.wrmds(output_dir+out_pckf[i], data, ndims=3, dataprec='float64',  nrecords=1, 
                  dimlist=[data.shape[2], data.shape[1], data.shape[0]], itr=pickup_itr)
    data = pickup1[-1]
    mds.wrmds(output_dir+out_pckf[-1], data, ndims=2, dataprec='float64',  nrecords=1, 
              dimlist=[data.shape[1], data.shape[0]], itr=pickup_itr)
    return pickup1

############## SEA ICE ##############
def gen_pickup_seaice(config_dir, pickup_itr, AngleCS0, AngleSN0, WCmask0, coast, 
                      XC1, YC1, HFacC1, Dtri, sigmaG, print_level):
    #### Read pickup file
    pickup_dir = os.path.join(config_dir, 'parent/outputs/pickups/')
    if print_level>=1:
        print('        > Reading in the ECCO pickup_seaice file')
    pickup0, _, metap = mds.rdmds(os.path.join(pickup_dir,'pickup_seaice'), itrs=pickup_itr, returnmeta=True)
    for i in range(len(pickup0)):
        pickup0[i] = transp_tiles(pickup0[i])
    #### Center & Rotate velocity fields
    if print_level>=1:
        print('        > Rotating orignal velocity fields into regional coordinates')
    Uid = metap['fldlist'].index('siUICE')
    Vid = metap['fldlist'].index('siVICE')
    Zvel, Mvel = UVrot(pickup0[Uid], pickup0[Vid], AngleCS0, AngleSN0)
    #### Interpolate ICs on the regional grid
    idsel = [metap['fldlist'].index(nm) for nm in ['siTICE', 'siAREA', 'siHEFF', 'siHSNOW']]
    ICvar = np.array(metap['fldlist'])[idsel].tolist()+['Zvel', 'Mvel']
    pickup1 = np.zeros((len(ICvar),XC1.shape[0],XC1.shape[1]))
    if print_level>=1:
        print('        > Downscaling initial conditions')
    for i, var in enumerate(ICvar):
        if var in metap['fldlist']:
            vid = metap['fldlist'].index(var)
            pickup1[i] = Hinterp(pickup0[vid], WCmask0, coast, Dtri, XC1, YC1, sigmaG)*HFacC1[0]
        else:
            pickup1[i] = Hinterp(locals()[var], WCmask0, coast, Dtri, XC1, YC1, sigmaG)*HFacC1[0]
    #### Generate new pickup files
    output_dir = os.path.join(config_dir, 'forcings/pickups/')
    if print_level>=1:
        print('        > Generating sea-ice initial conditions file')
    mds.wrmds(output_dir+'pickup_seaice', pickup1, ndims=metap['ndims'], dataprec='float64',  nrecords=metap['nrecords'], 
              fields=metap['fldlist'], dimlist=[pickup1[0].shape[1], pickup1[0].shape[0]], itr=pickup_itr)
    return pickup1

############## 1D ##############
def gen_pickup_1D(config_dir, pickup_itr, Nr0, AngleCS0, AngleSN0, WCmask0, drF0, coast,
                  XC1, YC1, Nr1, drF1, HFacC1, Dtri, sigmaG, print_level, init_fnm):
    #### Read pickup file
    pickup_dir = os.path.join(config_dir, 'parent/outputs/pickups/')
    if print_level>=1:
        print(f'        > Reading in the ECCO {init_fnm} file')
    pickup0, _, metap = mds.rdmds(os.path.join(pickup_dir,init_fnm), itrs=pickup_itr, returnmeta=True)
    for i in range(len(pickup0)):
        pickup0[i] = transp_tiles(pickup0[i])
    #### Interpolate ICs on the regional grid
    pickup1 = np.zeros((Nr1,XC1.shape[0],XC1.shape[1]))
    if print_level>=1:
        print('        > Downscaling initial conditions')
    tmp = Hinterp(pickup0, WCmask0, coast, Dtri, XC1, YC1, sigmaG)
    pickup1 = Zinterp(tmp, drF0, drF1, HFacC1)
    #### Generate new pickup files
    output_dir = os.path.join(config_dir, 'forcings/pickups/')
    if print_level>=1:
        print(f'        > Generating {init_fnm[7:]} initial conditions file')
    mds.wrmds(output_dir+init_fnm, pickup1, ndims=3, dataprec='float64',  nrecords=1, 
              dimlist=[pickup1.shape[2], pickup1.shape[1], pickup1.shape[0]], itr=pickup_itr)
    return pickup1

############## PTRACERS ##############
def gen_pickup_ptracers(config_dir, pickup_itr, Nr0, AngleCS0, AngleSN0, WCmask0, drF0, coast,
                        XC1, YC1, Nr1, drF1, HFacC1, Dtri, sigmaG, print_level, init_fnm):
    #### Read pickup file
    pickup_dir = os.path.join(config_dir, 'parent/outputs/pickups/')
    if print_level>=1:
        print(f'        > Reading in the ECCO {init_fnm} file')
    pickup0, _, metap = mds.rdmds(os.path.join(pickup_dir, init_fnm), itrs=pickup_itr, returnmeta=True)
    for i in range(pickup0.shape[0]):
        for j in range(pickup0.shape[1]):
            pickup0[i,j] = transp_tiles(pickup0[i,j])
    pickup1 = np.zeros((len(pickup0),Nr1,XC1.shape[0],XC1.shape[1]))
    if print_level>=1:
        print('        > Downscaling initial conditions')
    for i in range(len(pickup0)):
        tmp = Hinterp(pickup0[i], WCmask0, coast, Dtri, XC1, YC1, sigmaG)
        pickup1[i] = Zinterp(tmp, drF0, drF1, HFacC1)
    #### Generate new pickup files
    output_dir = os.path.join(config_dir, 'forcings/pickups/')
    if print_level>=1:
        print(f'        > Generating {init_fnm[7:]} initial conditions file')
    mds.wrmds(output_dir, pickup1, ndims=3, dataprec='float64',  nrecords=metap['nrecords'][0], 
              dimlist=[pickup1.shape[3], pickup1.shape[2], pickup1.shape[1]], itr=pickup_itr)
    return pickup1

############################################################
#                     MAIN FUNCTION                        #
############################################################

def gen_pickup_files(config_dir, model_name, pickup_itr, sigmaG, bgc, print_level, gennc):

    #################################################
    ############# Create forcing folder #############
    #################################################
    if 'forcings' not in os.listdir(config_dir):
        os.mkdir(os.path.join(config_dir,'forcings'))
    if 'pickups' not in os.listdir(os.path.join(config_dir,'forcings')):
        os.mkdir(os.path.join(os.path.join(config_dir,'forcings'),'pickups'))

    #################################################
    ####### Read regional grid file (Level 1) #######
    #################################################
    if print_level>=1:
        print('    - Reading in the regional model tile geometry')
    var_ls = ['XC', 'YC', 'AngleCS', 'AngleSN', 'HFacC', 'drF']
    data = read_ncgrid(config_dir, model_name, var_ls)
    for i, nm in enumerate(var_ls):
        globals()[nm+'1'] = data[i]
    Nr1 = len(drF1)
    
    ###############################################
    ####### Read parent grid file (Level 0) #######
    ###############################################
    if print_level>=1:
        print('    - Reading in the parent model tile geometry')
    XC0, YC0, dXC0, dYC0, AngleCS0, AngleSN0, HFacC0, drF0 = read_eccogrid(config_dir)
    Nr0 = HFacC0.shape[0]
       ## Wet cells
    WCmask0 = np.copy(HFacC0)
    WCmask0[WCmask0>0] = 1
    
    ###################################
    ####### Prepare Downscaling #######
    ###################################
    if print_level>=1:
        print('    - Preparing downscaling')
       ## Compute coast mask
    coast = get_coast_indexes(WCmask0)
       ## Compute Delaunay triangulation for interpolation
    if print_level>=1:
        print('        > Calculating Delaunay triangulation')
    Dtri = Delaunay(np.array([XC0.ravel(), YC0.ravel()]).T)
       ## Set the Standard deviation for Gaussian filter 
    if sigmaG == None:
        # Set sigma as the average distance [km] between llc cells in the region
        sigmaG = np.mean([dXC0[((XC0>XC1.min())&(XC0<XC1.max())&(YC0>YC1.min())&(YC0<YC1.max()))],
                          dYC0[((XC0>XC1.min())&(XC0<XC1.max())&(YC0>YC1.min())&(YC0<YC1.max()))]])/1e3
        sigmaG = int(sigmaG)
    else:
        pass # Sigma is already set by the user
    if print_level>=1:
        print('        > Standard deviation for Gaussian filter is %.i:'%sigmaG)
    
    ####################################
    ####### Generate pickup file #######
    ####################################
    #--------- physics condition ---------#
    if print_level>=1:
        print('    - Creating physics initial conditions for the '+model_name+' model from ECCO data')
    pickup_phy = gen_pickup(config_dir, pickup_itr, Nr0, AngleCS0, AngleSN0, WCmask0, drF0, coast,
                            XC1, YC1, Nr1, drF1, HFacC1, Dtri, sigmaG, print_level)
    
    #--------- sea-ice condition ---------#
    if print_level>=1:
        print('    - Creating sea-ice initial conditions for the '+model_name+' model from ECCO data')
    pickup_sic = gen_pickup_seaice(config_dir, pickup_itr, AngleCS0, AngleSN0, WCmask0, coast,
                                   XC1, YC1, HFacC1, Dtri, sigmaG, print_level)
    
    #--------- ggl90 condition ---------#
    if print_level>=1:
        print('    - Creating ggl90 initial conditions for the '+model_name+' model from ECCO data')
    pickup_ggl90 = gen_pickup_1D(config_dir, pickup_itr, Nr0, AngleCS0, AngleSN0, WCmask0, drF0, coast,
                                 XC1, YC1, Nr1, drF1, HFacC1, Dtri, sigmaG, print_level, 'pickup_ggl90')

    if bgc == True:
        #--------- darwin condition ---------#
        if print_level>=1:
            print('    - Creating darwin initial conditions for the '+model_name+' model from ECCO data')
        pickup_darwin = gen_pickup_1D(config_dir, pickup_itr, Nr0, AngleCS0, AngleSN0, WCmask0, drF0, coast,
                                      XC1, YC1, Nr1, drF1, HFacC1, Dtri, sigmaG, print_level, 'pickup_darwin')
        #--------- ptracer condition ---------#
        if print_level>=1:
            print('    - Creating ptracers initial conditions for the '+model_name+' model from ECCO data')
        pickup_ptracers = gen_pickup_ptracers(config_dir, pickup_itr, Nr0, AngleCS0, AngleSN0, WCmask0, drF0, coast,
                                              XC1, YC1, Nr1, drF1, HFacC1, Dtri, sigmaG, print_level, 'pickup_ptracers')
    
    ####################################
    ####### Generate netcdf file #######
    ####################################
    if gennc == True:
        if print_level>=1:
            print('    - Generating the netcdf file')
        gen_nc(XC1, YC1, drF1, Nr1, pickup_phy, pickup_sic, pickup_ggl90, pickup_darwin, pickup_ptracers)

    print("Pickup files successfully generated!")
    return

############################################################
#                          PARSER                          #
############################################################

if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument("-d", "--config_dir", action="store",
                        help="The directory where data files are stored", dest="config_dir",
                        type=str, required=True)
    parser.add_argument("-n", "--reg_nm", action="store",
                        help="Name of the regional cutout", dest="reg_nm",
                        type=str, required=True)
    parser.add_argument("-i", "--itr", action="store",
                        help="iteration of the begining the regional model", dest="itr",
                        type=int, required=True, default=1)
    parser.add_argument("-sg", '--sigma_Gfilt', nargs='?', type=int)
    parser.add_argument("-bgc", "--darwin", action="store_true", default='False',
                        help="generate darwin biogeochemistry pickups")
    parser.add_argument("-v", "--verbose", action="store_true", default='False')
    parser.add_argument("-nc", "--netcdf", action="store_true", default='False',
                        help="generate a verifiaction netcdf file")

    args = parser.parse_args()
    config_dir = args.config_dir
    model_name = args.reg_nm
    pickup_itr = args.itr
    sigmaG = args.sigma_Gfilt
    gennc = args.netcdf
    if args.verbose == True:
    	print_level = 1
    else:
    	print_level = 0
    bgc = args.darwin

    gen_pickup_files(config_dir, model_name, pickup_itr, sigmaG, bgc, print_level, gennc)
























