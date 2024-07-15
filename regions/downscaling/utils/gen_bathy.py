import os
import sys
import ast
import argparse
import numpy as np
import netCDF4 as nc4
from scipy.interpolate import griddata

def read_mitgrid(config_dir, model_name, size_dom):
    
    # Get the size of the domain
    n_rows = size_dom[0]
    n_cols = size_dom[1]
    
    # Read mitgrid file
    file_path = os.path.join(config_dir, model_name + '.mitgrid')
    entire_grid = np.fromfile(file_path, dtype='>f8')
    entire_grid = np.reshape(entire_grid, (16, n_cols + 1, n_rows + 1))
    
    # Read coordinates
    Lon_C = entire_grid[0, :, :]; Lon_C = Lon_C[:-1, :-1]
    Lat_C = entire_grid[1, :, :]; Lat_C = Lat_C[:-1, :-1]
    
    return Lon_C, Lat_C

def interp_gebco(bathy_dir, Lon_C, Lat_C):
    
    # Read GEBCO data
    gebco_file = os.path.join(bathy_dir)
    ds = nc4.Dataset(gebco_file)
    lon = ds.variables['lon'][:]
    lat = ds.variables['lat'][:]
    elev = ds.variables['elevation'][:,:]
    ds.close()
    
    # Get the indexes of corners on GEBCO map
    min_lon_index = np.argmin(np.abs(lon - np.min(Lon_C)))
    max_lon_index = np.argmin(np.abs(lon - np.max(Lon_C)))
    min_lat_index = np.argmin(np.abs(lat - np.min(Lat_C)))
    max_lat_index = np.argmin(np.abs(lat - np.max(Lat_C)))
    
    # Isolate the domain on GEBCO map
    lon = lon[min_lon_index:max_lon_index]
    lat = lat[min_lat_index:max_lat_index]
    elev = elev[min_lat_index:max_lat_index, min_lon_index:max_lon_index]
    
    # Interpolate GEBCO on the mitgrid (nearest neigbor)
    lon, lat = np.meshgrid(lon, lat)
    Incoor = np.array([lon.ravel(),lat.ravel()]).T
    Oucoor = np.array([Lon_C.ravel(), Lat_C.ravel()]).T
    bathy = griddata(Incoor, elev.ravel(), Oucoor, method='nearest')
    bathy = np.reshape(bathy,Lon_C.shape)
    bathy[bathy>0] = 0
    
    return bathy

def gen_surf_hFacC(bathy, delR, hFacMin, hFacMinDr):

    # This code is adapted from MITgcm inside ini_masks_etc.F
    # Define grids with same names as those in MITgcm
    RU = -1 * np.cumsum(delR)
    RL = np.concatenate([[0], RU[:-1]])
    drF = delR
    recip_drF = 1 / drF

    #Initialize hFacC grid
    hFacC = np.zeros((np.shape(bathy)[0], np.shape(bathy)[1]))
    # Compute hFacC
    k = 0 #(Surface)
    hFacMnSz = np.max([hFacMin, np.min([hFacMinDr * recip_drF[k], 1])])
    for i in range(np.shape(bathy)[0]):
        for j in range(np.shape(bathy)[1]):
            #      o Non-dimensional distance between grid bound. and domain lower_R bound.
            hFac_loc = (RL[k] - bathy[i, j]) * recip_drF[k]
            #      o Select between, closed, open or partial (0,1,0-1)
            hFac_loc = np.min([np.max([hFac_loc, 0]), 1])
            #      o Impose minimum fraction and/or size (dimensional)
            if hFac_loc <= hFacMnSz * 0.5 or bathy[i, j] >= 0:
                hFacC[i, j] = 0
            else:
                hFacC[i, j] = np.max([hFac_loc, hFacMnSz])

    return hFacC

def gen_connected_mask(start_row, start_col, wet_grid):

    if wet_grid[start_row,start_col]==0:
        raise ValueError(' The start row/col location is  dry')

    rows = np.arange(np.shape(wet_grid)[0])
    cols = np.arange(np.shape(wet_grid)[1])
    Cols,Rows = np.meshgrid(cols,rows)

    mask_grid = 1-np.copy(wet_grid)
    mask_grid[start_row,start_col] = 2
    # in the mask: 0 = unverified; 1 = dry; 2 = wet

    is_remaining = np.logical_and(mask_grid==0,wet_grid==1)
    n_remaining = np.sum(is_remaining)
    continue_iter = True
    for i in range(n_remaining):
        if continue_iter:
            # get the wet rows, cols, and their current mask values
            Wet_Rows = Rows[wet_grid == 1]
            Wet_Cols = Cols[wet_grid == 1]
            Mask_Vals = mask_grid[wet_grid == 1]
            
            # reduce these to the ones that havent been verified yet
            Wet_Rows = Wet_Rows[Mask_Vals == 0]
            Wet_Cols = Wet_Cols[Mask_Vals == 0]
            Mask_Vals = Mask_Vals[Mask_Vals == 0]

            if len(Mask_Vals)>0:

                # for each row/col, see if its connected to one we've verified is connected
                rows_remaining,cols_remaining = np.where(is_remaining)
                for ri in range(n_remaining):
                    row = rows_remaining[ri]
                    col = cols_remaining[ri]
                    
                    # this bit allows for up/down/left/right spreading
                    if row<np.shape(wet_grid)[0]-1:
                        if mask_grid[row+1,col] == 2:
                            mask_grid[row,col] = 2
                    if row > 0:
                        if mask_grid[row - 1, col] == 2:
                            mask_grid[row,col] = 2
                    if col<np.shape(wet_grid)[1]-1:
                        if mask_grid[row,col+1] == 2:
                            mask_grid[row,col] = 2
                    if col > 0:
                        if mask_grid[row, col-1] == 2:
                            mask_grid[row,col] = 2

                is_remaining = np.logical_and(mask_grid == 0, wet_grid == 1)
                n_remaining_now = np.sum(is_remaining)

                if n_remaining_now<n_remaining:
                    n_remaining = n_remaining_now
                else:
                    n_remaining = n_remaining_now
                    continue_iter=False
            else:
                continue_iter = False

    return(mask_grid)

def gen_bathy_file(config_dir, bathy_dir, model_name, size_dom, hFac, Scell_w, Cwet_pts, print_level=1):
    
    sys.path.insert(1, os.path.join(config_dir))
    if print_level>=1:
        print('Creating the bathymetry file')

    # Read domain coordinates
    if print_level>=1:
        print('        - Reading mitgrid file')
    Lon_C, Lat_C = read_mitgrid(config_dir, model_name, size_dom)
    
    # Interpolate with GEBCO
    if print_level >= 1:
        print('        - Interpolating from GEBCO')
    bathy = interp_gebco(bathy_dir, Lon_C, Lat_C)

    # Adjustments the bathymetry to model characteristics (surface cell height)
    hFacMin = hFac[0]; hFacMinDr = hFac[1]
    if print_level >= 1:
        print('        - Lowering shallow points down to a depth of '+str(hFacMinDr)+' m')
    bathy[np.logical_and(bathy>-Scell_w[0],bathy<0)] = -Scell_w[0]

    # Generate the hFacC grid
    if print_level >= 1:
        print('        - Generating a wet grid for the domain')
    hFacC_grid = gen_surf_hFacC(bathy, np.array(Scell_w), hFac[0], hFac[1])
    wet_grid = np.copy(hFacC_grid)
    wet_grid[wet_grid>1] = 1

    if print_level >= 1:
        print('        - Generating a mask for the domain to eliminate unconnected regions')
    Cwet_row = Cwet_pts[0]; Cwet_col = Cwet_pts[1]
    mask_grid = gen_connected_mask(Cwet_row, Cwet_col, wet_grid)
    bathy[mask_grid==0] = 0

    output_file = os.path.join(config_dir,model_name+'_bathymetry.bin')
    bathy.ravel('C').astype('>f4').tofile(output_file)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument("-d", "--config_dir", action="store",
                        help="The directory where the mitgrid is stored and \
                        the bathymetry will be stored", dest="config_dir",
                        type=str, required=True)    
    parser.add_argument("-g", "--gebco_dir", action="store",
                        help="The directory where the gebco netcdf file is stored + the\
                        name of the gebco file", dest="gebco_dir",
                        type=str, required=True)   
    parser.add_argument("-n", "--reg_nm", action="store",
                        help="Name of the regional cutout", dest="reg_nm",
                        type=str, required=True)
    parser.add_argument("-s", "--dom_size", nargs='+', action="store",
                        help="list of the mitgrid shape: \
                        [number of rows, number of colums]", 
                        dest="dom_size", type=int, required=True)
    parser.add_argument("-cs", "--Scell_size", nargs='+', action="store",
                        help="list with the  minimum fraction and dimension size of a \
                        cell: [hFacC_Min, hFacC_MinDir]", 
                        dest="Scell_size", type=float, required=True)
    parser.add_argument("-cw", "--Scell_width", nargs='+', action="store",
                        help="list with the 2 first surface layer height for the \
                        downscaled model: [h1, h2]", 
                        dest="Scell_width", type=float, required=True)
    parser.add_argument("-wp", "--wetC_pt", nargs='+', action="store",
                        help="list with the position of the wet cell situated\
                        in the center of the domain: [row_wet_cell, col_wet_cell]", 
                        dest="wetC_pt", type=int, required=True)
    parser.add_argument("-v", "--verbose", action="store_true", default='False')
    
    args = parser.parse_args()
    config_dir = args.config_dir
    gebco_dir = args.gebco_dir
    reg_nm = args.reg_nm
    dom_size = args.dom_size
    hFac = args.Scell_size
    Scell_w = args.Scell_width
    Cwet_pts = args.wetC_pt
    if args.verbose == True:
    	print_level = 1
    else:
    	print_level = 0

    gen_bathy_file(config_dir, gebco_dir, reg_nm, dom_size, hFac, Scell_w, Cwet_pts, print_level)
    print("Bathymetry file generated")