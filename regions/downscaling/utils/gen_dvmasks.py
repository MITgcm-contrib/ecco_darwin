import os
import argparse
import numpy as np
import xarray as xr

def gen_dict_nc(output_dir,output_name,all_mask_dicts,mask_names_list):
    #### Delete previous file
    if output_name in os.listdir(output_dir):
        os.remove(os.path.join(output_dir,output_name))
    
    #### Generate the xarray
    data = {}; dims = {}
    for i,nm in enumerate(mask_names_list):
        dims.update({nm+'_points': np.arange(len(all_mask_dicts[i]))})
        data.update({nm+'_rows': ([nm+'_points'], np.asarray(all_mask_dicts[i][:,0], dtype=int))})
        data.update({nm+'_cols': ([nm+'_points'], np.asarray(all_mask_dicts[i][:,1], dtype=int))})
    ds = xr.Dataset(data_vars=data, coords=dims)

    #### save into a netcdf file
    ds.to_netcdf(os.path.join(output_dir,output_name), mode='w')
    return

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

def gen_bnd_mask(resolution, XC, YC, Depth, XC_bnd, YC_bnd):
    ### Initialize
    mask_indices = []
    mask_grid = np.zeros(XC.shape)
    counter = 1
    ### Generate the boundary mask
    for i in range(len(XC_bnd)):
        dist = get_circle_dist(XC_bnd[i], YC_bnd[i], XC, YC)
        rows, cols = np.where(dist<resolution*np.sqrt(2))
        for i in range(len(rows)):
            if mask_grid[rows[i],cols[i]]==0 and Depth[rows[i],cols[i]]>0:
                mask_grid[rows[i],cols[i]] = counter
                mask_indices.append([rows[i],cols[i]])
                counter +=1
    mask_indices = np.array(mask_indices)
    return mask_grid, mask_indices

def read_ncgrid(config_dir, model_name):
    ds = xr.open_dataset(os.path.join(config_dir,reg_nm+'_ncgrid.nc'))
    XC = ds.XC.values
    YC = ds.YC.values
    Depth = ds.Depth.values
    ds.close()
    return XC, YC, Depth

def gen_dvmasks(config_dir, reg_nm, boundaries, resolution, print_level):
    
    #### Create the dv folder
    if 'dv' not in os.listdir(os.path.join(config_dir)):
        os.mkdir(os.path.join(config_dir,'dv'))
        
    #### Read grid imformation
    if print_level==1:
         print('- Reading in the model geometry')
    XC, YC, Depth = read_ncgrid(config_dir, reg_nm)
    if print_level==1:
        print('    > The domain shape is ' + str(np.shape(XC)))
        
    #### Compute dv masks
        # Initialize
    mask_names_list = []
    all_mask_grids = []
    all_mask_dicts = []
        # Compute
    if print_level == 1:
        print('- Creating the boundary masks')
    for i in range(len(boundaries)):
        if boundaries[i]=='S':
            print('    > Generating SOUTH boundary mask')
            mask_names_list.append('south')
            XC_bnd = XC[0,:]; YC_bnd = YC[0,:]
        elif boundaries[i]=='N':
            print('    > Generating NORTH boundary mask')
            mask_names_list.append('north')
            XC_bnd = XC[-1,:]; YC_bnd = YC[-1,:]
        elif boundaries[i]=='W':
            print('    > Generating WEST boundary mask')
            mask_names_list.append('west')
            XC_bnd = XC[:,0]; YC_bnd = YC[:,0]
        elif boundaries[i]=='E':
            print('    > Generating EAST boundary mask')
            mask_names_list.append('east')
            XC_bnd = XC[:,-1]; YC_bnd = YC[:,0-1]
        else:
            raise SystemExit("boundary name "+boundaries[i]+" doesn't exist. Please specify a boundary with S,N,W or E")
        mask_grid, mask_indices = gen_bnd_mask(resolution, XC, YC, Depth, XC_bnd, YC_bnd)
        if print_level == 1:
            print('        The '+mask_names_list[i]+' boundary mask has '+str(len(mask_indices))+' points')
        all_mask_grids.append(mask_grid)
        all_mask_dicts.append(mask_indices)

    #### save dv mask into files & dict file
    for i,nm in enumerate(mask_names_list):
        output_file = os.path.join(config_dir,'dv', nm+'_BC_mask.bin')
        mask_grid=np.array(all_mask_grids[i])
        mask_grid.ravel(order='C').astype('>f4').tofile(output_file)
    output_dir = os.path.join(config_dir)
    output_name = reg_nm+'_dv_mask_reference_dict.nc'
    gen_dict_nc(output_dir,output_name,all_mask_dicts,mask_names_list)
    
    return 

if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument("-d", "--config_dir", action="store",
                        help="The directory where .... is stored", dest="config_dir",
                        type=str, required=True)
    parser.add_argument("-n", "--reg_nm", action="store",
                        help="Name of the regional cutout", dest="reg_nm",
                        type=str, required=True)
    parser.add_argument("-b", "--bnd_nm", action="store",
                        help="boundaries where to generate a mask", dest="bnd_nm",
                        type=str, required=True)
    parser.add_argument("-r", "--reso", action="store", help="grid spacing (has to be an integer)", 
                        dest="reso", type=int, required=True)
    parser.add_argument("-v", "--verbose", action="store_true", default='False')

    args = parser.parse_args()
    config_dir = args.config_dir
    reg_nm = args.reg_nm
    bnd_nm = args.bnd_nm
    reso = args.reso
    if args.verbose == True:
        print_level = 1
    else:
        print_level = 0  
    
    gen_dvmasks(config_dir, reg_nm, bnd_nm, reso, print_level)
    print("dv masks generated")