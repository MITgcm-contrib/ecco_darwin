import os
import argparse
import numpy as np
import xarray as xr
from MITgcmutils import mds, llc

############################################################
#                    READING FUNCTIONS                     #
############################################################

def read_parent(parent_dir, bathy_file):
    grid_dir = os.path.join(parent_dir,'outputs/grid/')
    mitg_dir = os.path.join(parent_dir,'outputs/mitgrid/')

    XC = transp_tiles(mds.rdmds(grid_dir+'XC'))
    YC = transp_tiles(mds.rdmds(grid_dir+'YC'))
    bathy = np.fromfile(os.path.join(grid_dir, bathy_file), '>f4').reshape(XC.shape)
    bathy = transp_tiles(bathy)
    return XC,YC, bathy

def read_ncgrid(config_dir, model_name):
    ds = xr.open_dataset(os.path.join(config_dir,reg_nm+'_ncgrid.nc'))
    XC = ds.XC.values
    YC = ds.YC.values
    Depth = ds.Depth.values
    ds.close()
    return XC, YC, Depth

############################################################
#                  GRID WORK FUNCTIONS                     #
############################################################

def transp_tiles(data, reverse=False):
    nx = data.shape[1]
    tmp = data[7*nx:,::-1]
    if reverse == True:
        transpo = np.concatenate([tmp[-3*nx:],tmp[:3*nx]],axis=1)
        data_out = np.zeros((13*nx,nx))
        data_out[:7*nx] = data[:7*nx]
        data_out[7*nx:][2::3,:] = np.rot90(transpo[2*nx:])
        data_out[7*nx:][1::3,:] = np.rot90(transpo[nx:2*nx])
        data_out[7*nx:][0::3,:] = np.rot90(transpo[:nx])
    else:
        transpo = np.concatenate((tmp[2::3,:].transpose(),tmp[1::3,:].transpose(),tmp[0::3,:].transpose()))
        data_out = np.concatenate((data[:7*nx],np.flipud(transpo[:,:nx]),np.flipud(transpo[:,nx:])))
    return data_out

def get_circle_dist(lon_ref, lat_ref, Lon, Lat):
    earth_radius = 6371000
    lon_ref_radians = np.radians(lon_ref)
    lat_ref_radians = np.radians(lat_ref)
    lons_radians = np.radians(Lon)
    lats_radians = np.radians(Lat)
    lat_diff = lats_radians - lat_ref_radians
    lon_diff = lons_radians - lon_ref_radians
    d = np.sin(lat_diff * 0.5) ** 2 + np.cos(lat_ref_radians) * np.cos(lats_radians) * np.sin(lon_diff * 0.5) ** 2
    h = 2 * earth_radius * np.arcsin(np.sqrt(d))
    return h

############################################################
#                MASK CREATION FUNCTIONS                   #
############################################################

def gen_bnd_mask(XC0, YC0, bathy0, B_coords, resolution):
    ### Initialize
    mask = np.zeros(XC0.shape)
    mask_dict = []
    for n in range(B_coords.shape[0]): 
        ### Get global coordinates closer to regional coordinates
        x = B_coords[n,0]; y = B_coords[n,1]
        dist = get_circle_dist(x, y, XC0, YC0)
        rows, cols = np.where(dist <= resolution*np.sqrt(2))
        ### Store mask information is not land  
        ctr = 1
        for i in range(len(rows)):
            row = rows[i]; col = cols[i]
            if bathy0[row,col]<0:
                if mask[row,col]==0:
                    mask[row,col] = ctr
                    mask_dict.append([row,col])
                    ctr += 1
    mask_dict = np.array(mask_dict)
    ### Collect points information
    npts = len(mask_dict)
    return mask, mask_dict, npts

############################################################
#                        SKELETON                          #
############################################################

def gen_dvmasks(config_dir, reg_nm, bathy_file, boundaries, resolution, print_level):
    
    #### Create the parent folder
    if 'parent' not in os.listdir(config_dir):
        os.mkdir(os.path.join(config_dir,'parent'))
    parent_dir = os.path.join(config_dir,'parent')
    if 'inputs' not in os.listdir(parent_dir):
        os.mkdir(os.path.join(parent_dir,'inputs'))
    save_dir = os.path.join(parent_dir,'inputs')

    #### Read ECCO grid information
    if print_level==1:
         print('- Reading the ECCO global grid information')
    XC0, YC0, bathy0 = read_parent(parent_dir, bathy_file)

    #### Read cut-out grid information
    if print_level==1:
         print(f'- Reading the {reg_nm} model geometry')
    XC1, YC1, bathy1 = read_ncgrid(config_dir, reg_nm)

    #### Get the coordinates of the boundaries 
    S_coords = np.column_stack([XC1[:1, :].ravel(), YC1[:1, :].ravel()])
    N_coords = np.column_stack([XC1[-1:, :].ravel(), YC1[-1:, :].ravel()])
    W_coords = np.column_stack([XC1[:, :1].ravel(), YC1[:, :1].ravel()])
    E_coords = np.column_stack([XC1[:, -1:].ravel(), YC1[:, -1:].ravel()])

    #### Generate masks files
    ## initialize
    masks_nm = []
    masks_grids = []
    masks_dicts = []
    if print_level == 1:
        print('- Creating the boundary masks')
    for i in range(len(boundaries)):
        bnd = boundaries[i]
        ## recover the boundaries needed
        if bnd =='S':
            print('    > Generating SOUTH boundary mask')
            masks_nm.append('south')
        elif bnd =='N':
            print('    > Generating NORTH boundary mask')
            masks_nm.append('north')
        elif bnd =='W':
            print('    > Generating WEST boundary mask')
            masks_nm.append('west')
        elif bnd =='E':
            print('    > Generating EAST boundary mask')
            masks_nm.append('east')
        else:
            raise SystemExit("boundary name "+bnd+" doesn't exist. Please specify a boundary with S,N,W or E")
        ## compute the mask
        tmp = gen_bnd_mask(XC0, YC0, bathy0, locals()[f'{bnd}_coords'], resolution*1e3)
        if print_level == 1:
            print('        The '+masks_nm[i]+' boundary mask has '+str(tmp[2])+' points')
        masks_grids.append(tmp[0])
        masks_dicts.append(tmp[1])

    #### Saving masks files
    if print_level == 1:
        print('- Saving masks binary files')
    for i,nm in enumerate(masks_nm):
        out_fl = os.path.join(save_dir,nm+'_BC_mask.bin')
        mask_save = transp_tiles(masks_grids[i], reverse=True)
        mask_save.ravel(order='C').astype('>f4').tofile(out_fl)     
    return 

############################################################
#                         PARSER                          #
############################################################

if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument("-d", "--config_dir", action="store",
                        help="The directory where parents folder is stored", dest="config_dir",
                        type=str, required=True)
    parser.add_argument("-n", "--reg_nm", action="store",
                        help="Name of the regional cutout", dest="reg_nm",
                        type=str, required=True)
    parser.add_argument("-bfl", "--bathy_fle", action="store",
                        help="The ECCO global model bathymetry file name", dest="bathy_fle",
                        type=str, required=True)
    parser.add_argument("-bnd", "--bnd_nm", action="store",
                        help="boundaries where to generate a mask", dest="bnd_nm",
                        type=str, required=True)
    parser.add_argument("-r", "--reso", action="store",
                        help="Parent grid spacing/ also reprsent the size of the searching area for boundary coordinates", 
                        dest="reso", type=int, required=True)
    parser.add_argument("-v", "--verbose", action="store_true", default='False')

    args = parser.parse_args()
    config_dir = args.config_dir
    reg_nm = args.reg_nm
    bathy_file = args.bathy_fle
    bnd_nm = args.bnd_nm
    reso = args.reso
    if args.verbose == True:
    	print_level = 1
    else:
    	print_level = 0  
    
    gen_dvmasks(config_dir, reg_nm, bathy_file, bnd_nm, reso, print_level)