import os
import glob
import argparse
import numpy as np
import xarray as xr
from MITgcmutils import mds, llc
from scipy.interpolate import interp1d, griddata

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

def gen_bnd_domain(bnd_ls, dv_masks_ls, grid, Nr):
    
    ### Initialize boundary filler information
    var_ls = ['XC','YC','AngleCS','AngleSN','maskC','maskW','maskS']
    boundary_domain_dict = dict.fromkeys(bnd_ls, var_ls)
    for bnd, keys in boundary_domain_dict.items():
        boundary_domain_dict[bnd] = {k: None for k in keys}
    for i in range(len(bnd_ls)):
        npts = np.sum(dv_masks_ls[i]!=0)
        ### make empty arrays to fill in
        for nm in var_ls:
            if nm[:-1]=='mask':
                globals()[nm+'_pts'] = np.zeros((Nr, npts))
            else:
                globals()[nm+'_pts'] = np.zeros((npts, ))
        ### Fill in the interpolater
        dv_rows, dv_cols = np.where(dv_masks_ls[i]!=0)
        if len(dv_rows)>0:
            for row, col in zip(dv_rows, dv_cols):
                for nm in var_ls:
                    if nm[:-1]=='mask':
                        for k in range(Nr):
                            globals()[nm+'_pts'][k,int(dv_masks_ls[i][row,col])-1] = grid["HFac"+nm[-1:]][k,row,col]
                    else:
                        globals()[nm+'_pts'][int(dv_masks_ls[i][row,col])-1] = grid[nm][row,col]
                    boundary_domain_dict[bnd_ls[i]][nm] = globals()[nm+'_pts']
    return boundary_domain_dict

def get_regbnd_info(vnm, bnd, Nr, tstp, grid1):
    
    #---------- Get right mask -----------#
    if vnm == 'UVEL' or vnm == 'UICE':
        mask = grid1['maskW']
    elif vnm == 'VVEL' or vnm == 'VICE':
        mask = grid1['maskS']
    else:
        mask = grid1['maskC']
       
    #-------- Get bnd mask/coord ---------#
    if bnd == 'west':
        XC1 = grid1['XC'][:,:1]
        YC1 = grid1['YC'][:,:1]
        mask1 = mask[:,:,:1]
    elif bnd == 'east':
        XC1 = grid1['XC'][:,-1:]
        YC1 = grid1['YC'][:,-1:]
        mask1 = mask[:,:,-1:]
    elif bnd == 'north':
        XC1 = grid1['XC'][-1:,:]
        YC1 = grid1['YC'][-1:,:]
        mask1 = mask[:,-1:,:]
    elif bnd == 'south':
        XC1 = grid1['XC'][:1,:]
        YC1 = grid1['YC'][:1,:]
        mask1 = mask[:,:1,:]
    else:
        raise ValueError('Boundary '+bnd+' not recognized')

    #----------- surface mask ------------# 
    if vnm in ['AREA','UICE','VICE','HEFF','HSNOW', 'ETAN']:
        mask1 = mask1[:1,:,:]
        obcs_init = np.zeros((tstp, 1, np.size(XC1)))
    else:
        obcs_init = np.zeros((tstp, Nr, np.size(XC1)))
        
    return XC1, YC1, mask1, obcs_init

def gen_obcs(dv_diag, XC0, YC0, mask00, XC1, YC1, mask1):

    #---------- Initialization -----------#
    obcs = np.zeros((mask1.shape))
    coord00 = np.array([XC0, YC0]).T

    #----------- Generaet OBCS -----------#
    for k in range(len(obcs)):
        if np.any(mask1[k] > 0):
            ### initialize interpolation at depth k 
            coord0 = coord00.copy()
            val0 = dv_diag[k].reshape(-1,1)
            mask0 = mask00[k].reshape(-1,1)
            coord0 = coord0[mask0[:,0] != 0, :]
            val0 = val0[mask0[:,0] != 0, :]
            if len(val0)>4:
                val1 = griddata(coord0, val0, (XC1, YC1), method='linear',fill_value=0)
                val1 = val1[:, :, 0]
                if not np.any(val1!=0):
                    val1 = griddata(coord0, val0, (XC1, YC1), method='nearest', fill_value=0)
                    val1 = val1[:, :, 0]
            else:
                val1 = np.zeros_like(XC1).astype(float)
            ### Apply bathymetry masks
            val1[mask1[k, :, :] == 0] = 0
            ### spread the data outward to new wet cells
            val1, n_remaining = Hextrapol(val1, mask1[k, :, :])
            if n_remaining > 0 and k > 0:
                val1 = Vextrapol(obcs, val1, mask1[k, :, :], k, 0)
            # Fill with the nearest neighbor if remaining values to fill
            if n_remaining > 0:
                if len(coord0) > 0:
                    val1_NR = griddata(coord0, val0, (XC1, YC1), method='nearest', fill_value=0)
                    val1_NR = val1_NR[:, :, 0]
                    ids = np.logical_and(val1 == 0, mask1[k, :, :] != 0)
                    val1[ids] = val1_NR[ids]
        else:
            val1 = np.zeros_like(XC1).astype(float)
        obcs[k, :, :] = val1[:, :]

    return obcs

############################################################
#                READING/WRITING FUNCTIONS                 #
############################################################

def read_eccogrid(config_dir, reg_nm):
    grid_dir = os.path.join(config_dir, 'parent/outputs/grid/')
    ####### 1D values #######
    rF = np.asanyarray(mds.rdmds(grid_dir+'RF')[:,0,0], dtype=np.float32)
    drF = np.asanyarray(mds.rdmds(grid_dir+'DRF')[:,0,0], dtype=np.float32)
    ####### 2D values #######
    XC = transp_tiles(mds.rdmds(grid_dir+'XC'))
    YC = transp_tiles(mds.rdmds(grid_dir+'YC'))
    AngleCS = transp_tiles(mds.rdmds(grid_dir+'AngleCS'))
    AngleSN = transp_tiles(mds.rdmds(grid_dir+'AngleSN'))
    ####### 3D values #######
    for nm in ['hFacC','hFacS','hFacW']:
        tmp = mds.rdmds(grid_dir+nm)
        hFac = np.zeros(tmp.shape)
        for i in range(len(hFac)):
            hFac[i] = transp_tiles(tmp[i])
        globals()[nm] = hFac
    return XC, YC, AngleCS, AngleSN, hFacC, hFacS, hFacW, drF, rF

def read_ncgrid(config_dir, reg_nm, var_ls):
    ds = xr.open_dataset(os.path.join(config_dir,reg_nm+'_ncgrid.nc'))
    var_mats = []
    for var in var_ls:
        var_mats.append(ds[var].values)
    ds.close()
    return var_mats

def read_dv_masks(config_dir, boundaries, llc):
    input_dir = os.path.join(config_dir, 'parent/inputs/')
    bnd_ls = []; masks_ls = []
    for nm in boundaries:
        if nm == 'E':
            bnd_ls.append('east')
            tmp = np.fromfile(input_dir+'east_BC_mask.bin', '>f4').reshape((13*llc,llc))
            masks_ls.append(transp_tiles(tmp))
        elif nm == 'W':
            bnd_ls.append('west')
            tmp = np.fromfile(input_dir+'west_BC_mask.bin', '>f4').reshape((13*llc,llc))
            masks_ls.append(transp_tiles(tmp))
        elif nm == 'N':
            bnd_ls.append('north')
            tmp = np.fromfile(input_dir+'north_BC_mask.bin', '>f4').reshape((13*llc,llc))
            masks_ls.append(transp_tiles(tmp))
        elif nm == 'S':
            bnd_ls.append('south')
            tmp = np.fromfile(input_dir+'south_BC_mask.bin', '>f4').reshape((13*llc,llc))
            masks_ls.append(transp_tiles(tmp))
        else:
            raise SystemExit("boundary name "+nm+" doesn't exist. Please specify a boundary with S,N,W or E")
    return bnd_ls, masks_ls

def read_dv_diags(itrs_ls, vnm, bnd, bnd_domain, Nr0, Nr1, grid0, grid1):

    #------------- Initialize ------------#
    dv_diags_dir = os.path.join(config_dir, 'parent/outputs/OBCS/')
    sfx = dv_diags_dir+f"{bnd}_BC_mask_{vnm}"
    itrs = [int(glob.glob(sfx+".*")[i][len(sfx)+1:-4]) for i in range(len(glob.glob(sfx+".*")))]
    if vnm in ['AREA','UICE','VICE','HEFF','HSNOW','ETAN']:
        Nr = 1
    else:
        Nr = Nr0

    #--------- Read dv diagnostic --------#
    for itr in itrs:
        ### Rotate U/V fields
        if 'VEL' in vnm or 'ICE' in vnm:
            comp = vnm[-3:]
            dv_diagU = np.fromfile(sfx[:-4]+"U"+comp+f".{itr:010d}.bin", '>f4')
            dv_diagV = np.fromfile(sfx[:-4]+"V"+comp+f".{itr:010d}.bin", '>f4')
            nm_pt = len(bnd_domain[bnd]['XC'])
            tstp = int(np.size(dv_diagU)/(Nr*nm_pt))
            dv_diagU = np.reshape(dv_diagU, (tstp, Nr, nm_pt))
            dv_diagV = np.reshape(dv_diagV, (tstp, Nr, nm_pt))
            dv_diag = np.zeros((tstp, Nr, nm_pt))
            if vnm in ['UVEL','UICE']:
                for m in range(nm_pt):
                    dv_diag[:,:,m] = bnd_domain[bnd]['AngleCS'][m] * dv_diagU[:,:,m] -\
                                     bnd_domain[bnd]['AngleSN'][m] * dv_diagV[:,:,m]
            elif vnm in ['VVEL','VICE']:
                for m in range(nm_pt):
                    dv_diag[:,:,m] = bnd_domain[bnd]['AngleSN'][m] * dv_diagU[:,:,m] +\
                                     bnd_domain[bnd]['AngleCS'][m] * dv_diagV[:,:,m]
        else:
            dv_diag =  np.fromfile(sfx+f".{itr:010d}.bin", '>f4')
            nm_pt = len(bnd_domain[bnd]['XC'])
            tstp = int(np.size(dv_diag)/(Nr*nm_pt))
            dv_diag = np.reshape(dv_diag, (tstp, Nr, nm_pt)) 
        ### merge timesteps 
        if itr==itrs[0]:
            dv_diag_all = dv_diag
        else:
            dv_diag = np.concatenate([dv_diag_all, dv_diag], axis=0)
        
    #------- Vertical interpolation ------#
    if vnm in ['UVEL','UICE']:
        msk_pts = bnd_domain[bnd]['maskW']
    elif vnm in ['VVEL','VICE']:
        msk_pts = bnd_domain[bnd]['maskS']
    else:
        msk_pts = bnd_domain[bnd]['maskC']

    if vnm not in ['AREA','UICE','VICE','HEFF','HSNOW', 'ETAN']:
        dv_diag_IT, msk_pts_IT = Zinterp(dv_diag, msk_pts, grid0['drF'], grid1['drF'])
    else:
        msk_pts = msk_pts[:1,:]
        dv_diag_IT = dv_diag
        msk_pts_IT = msk_pts

    return dv_diag_IT, msk_pts_IT

############################################################
#               INTERPOLATION FUNCTIONS                    #
############################################################

def Zinterp(dv_diag, msk_pts, delR0, delR1):
    
    #----------- Calculate deph ----------#
    # Parents
    Zbot0 = np.cumsum(delR0)
    Ztop0 = np.concatenate([np.array([0]), Zbot0[:-1]])
    Z0 = (Zbot0 + Ztop0) / 2
    # regions
    Zbot1 = np.cumsum(delR1)
    Ztop1 = np.concatenate([np.array([0]), Zbot1[:-1]])
    Z1 = (Zbot1 + Ztop1) / 2

    #-------- Inteprolate vertical -------#
    ### Diagnostics
    idS = np.where(np.abs(Z1)<np.abs(Z0[0]))[0]
    dv_diag_IT = np.zeros((dv_diag.shape[0], len(Z1), dv_diag.shape[2]))    
    for i in range(dv_diag.shape[2]):
        tmp0 = dv_diag[:,:,i].copy()
        if np.sum(tmp0 != 0) > 1:
            # interpolate on vertical
            tmp0[np.where(tmp0 == 0)] = np.nan
            f = interp1d(Z0, tmp0, axis=1, kind='linear',
                         bounds_error=False, fill_value=np.nan)
            tmp1 = f(Z1)
            # extrapolate surface data
            for j in idS:
                tmp1[:,j] = tmp1[:,idS[-1]+1]
            # handle possible bottom missing values
            if np.size(np.abs(Z0[np.isnan(tmp0[0])])) > 0:
                bottom_depth = np.abs(Z0[np.isnan(tmp0[0])])[0]
                bottom_values = tmp1[:,np.where(~np.isnan(tmp1[0]))[0][-1]]
                idB = np.where(np.logical_and(np.isnan(tmp1[0]), np.abs(Z1) < bottom_depth))[0]
                if len(idB) != 0:
                    tmp1[:,idB[0]] = bottom_values
            # Fill remaining NaNs downward
            tmp1[np.isnan(tmp1)] = 0
            dv_diag_IT[:,:,i] = tmp1
    ### Wet cells mask
    if msk_pts.shape[0] != len(Z1):
        msk_pts_IT = np.zeros((len(Z1), msk_pts.shape[1]))  
        for i in range(msk_pts.shape[1]):
            # interpolate on vertical
            tmp0 = msk_pts[:,i].copy()
            if np.sum(tmp0 != 0) > 1:
                f = interp1d(Z0[tmp0 != 0], tmp0[tmp0 != 0], kind='linear',
                             bounds_error=False, fill_value=np.nan)
                tmp1 = f(Z1)
                # extrapolate surface data
                tmp1[np.abs(Z1) < np.abs(Z0[0])] = tmp1[~np.isnan(tmp1)][0]
                # handle possible bottom missing values
                if np.size(np.abs(Z0[tmp0 == 0])) > 0:
                    bottom_depth = np.abs(Z0[tmp0 == 0])[0]
                    bottom_value = tmp1[~np.isnan(tmp1)][-1]
                msk_pts_IT[:,i] = tmp1
            elif np.sum(tmp0 != 0) == 1:
                msk_pts_IT[0, i] = msk_pts[0, i]
    else:
        msk_pts_IT = msk_pts
    msk_pts_IT[np.isnan(msk_pts_IT)] = 0
    msk_pts_IT = np.round(msk_pts_IT).astype(int)
    
    return dv_diag_IT, msk_pts_IT

############################################################
#               EXTRAPOLATION FUNCTIONS                    #
############################################################

def Hextrapol(var_grid, wet_grid, verbose=False):
    """
    Fill zero values in var_grid using the nearest non-zero neighbors within wet areas.
    Legacy implementation using NumPy.

    Parameters:
        var_grid (ndarray): 2D array with variable data, where zeros represent missing values.
        wet_grid (ndarray): 2D mask where 1 indicates wet cells and 0 indicates dry.

    Returns:
        tuple: Updated var_grid and number of remaining unfilled wet cells.
    """
    rows = np.arange(np.shape(var_grid)[0])
    cols = np.arange(np.shape(var_grid)[1])
    Cols, Rows = np.meshgrid(cols, rows)

    is_remaining = np.logical_and(var_grid == 0, wet_grid == 1)
    n_remaining = np.sum(is_remaining)
    continue_iter = True

    for _ in range(n_remaining):
        if verbose:
            if n_remaining>0:
                print(f"         - Remaining cells to fill: {n_remaining}")
        if continue_iter:
            Wet_Rows = Rows[wet_grid == 1]
            Wet_Cols = Cols[wet_grid == 1]
            Wet_Vals = var_grid[wet_grid == 1]
            Wet_Rows = Wet_Rows[Wet_Vals != 0]
            Wet_Cols = Wet_Cols[Wet_Vals != 0]
            Wet_Vals = Wet_Vals[Wet_Vals != 0]

            if len(Wet_Vals) > 0:
                rows_remaining, cols_remaining = np.where(is_remaining)
                for ri in range(n_remaining):
                    row = rows_remaining[ri]
                    col = cols_remaining[ri]
                    row_col_dist = ((Wet_Rows.astype(float) - row) ** 2 + (Wet_Cols.astype(float) - col) ** 2) ** 0.5
                    closest_index = np.argmin(row_col_dist)
                    if row_col_dist[closest_index] < np.sqrt(2):
                        var_grid[row, col] = Wet_Vals[closest_index]

                is_remaining = np.logical_and(var_grid == 0, wet_grid == 1)
                n_remaining_now = np.sum(is_remaining)
                if n_remaining_now < n_remaining:
                    n_remaining = n_remaining_now
                else:
                    continue_iter = False
            else:
                continue_iter = False

    return var_grid, n_remaining

def Vextrapol(full_grid, level_grid, wet_grid, level, mean_vertical_difference):
    """
    Spreads values from a given depth level vertically downward into wet grid cells where no data exists.

    Parameters:
    - full_grid (ndarray): The full variable grid from which to source values.
    - level_grid (ndarray): The grid to be filled at the current level.
    - wet_grid (ndarray): Mask indicating wet points (1 = wet, 0 = dry).
    - level (int): Current level (must be > 0).
    - mean_vertical_difference (float): Constant to add while copying downward.

    Returns:
    - level_grid (ndarray): The updated grid with missing values filled vertically.
    """

    if level == 0:
        bad_row, bad_col = np.where(np.logical_and(level_grid == 0, wet_grid == 1))
        raise ValueError(
            'Cannot spread vertically in the surface layer e.g. at row=' + str(bad_row[0]) + ', col=' + str(bad_col[0]))

    # Identify locations that are wet but uninitialized in level_grid
    is_remaining = np.logical_and(level_grid == 0, wet_grid == 1)
    rows_remaining, cols_remaining = np.where(is_remaining)

    for ri in range(len(rows_remaining)):
        row = rows_remaining[ri]
        col = cols_remaining[ri]
        # Fill the level grid by copying from one level above and adding a vertical offset
        level_grid[row, col] = full_grid[level - 1, row, col] + mean_vertical_difference

    return level_grid

############################################################
#                     MAIN FUNCTION                        #
############################################################

def gen_obcs_files(config_dir, reg_nm, boundaries, itrs, bgc, print_level):

    grd_ls = ['XC', 'YC', 'AngleCS', 'AngleSN', 'HFacC', 'HFacS', 'HFacW', 'drF']
    ###############################################
    #######   Read parent/region grid info  #######
    ###############################################
    #--------------- Parent --------------#
    if print_level>=1:
        print('> Reading in the parent model tile geometry')
    tmp = read_eccogrid(config_dir, reg_nm)
    grid0 = dict(zip(grd_ls+['rF'], tmp))
    Nr0 = tmp[grd_ls.index("HFacC")].shape[0]
    llc = tmp[grd_ls.index("XC")].shape[1]
    del tmp
    #--------------- region --------------#
    if print_level>=1:
        print('> Reading in the regional model tile geometry')
    tmp = read_ncgrid(config_dir, reg_nm, grd_ls)
    tmp[grd_ls.index("HFacS")] = tmp[grd_ls.index("HFacS")][:,:-1,:]
    tmp[grd_ls.index("HFacW")] = tmp[grd_ls.index("HFacW")][:,:,:-1]
    grid1 = dict(zip(grd_ls, tmp))
    Nr1 = tmp[grd_ls.index("HFacC")].shape[0]
    del tmp

    #################################################
    #######         Read/Generate masks       #######
    #################################################
    bnd_ls, dv_masks_ls = read_dv_masks(config_dir, boundaries, llc)
    for nm in ['C','S','W']:
        tmp = np.copy(grid1["HFac"+nm]); tmp[tmp>0] = 1
        grid1.update({"mask"+nm: tmp}); del tmp

    #################################################
    #######        Generate OBCS files        #######
    #################################################
    out_dir = os.path.join(config_dir,'forcings/OBCS/')
    #-------- Get boundary domain --------#
    if print_level>=1:
        print('> Getting boundary domain')
    bnd_domain = gen_bnd_domain(bnd_ls, dv_masks_ls, grid0, Nr0) 
    #--------- physics condition ---------#
    if print_level>=1:
        print('> Creating physics OBCS for the '+reg_nm+' model')
    vnms = ['THETA','SALT','UVEL','VVEL', 'ETAN']
    for vnm in vnms:
        if print_level>=1:
            print(f'    - Processing {vnm} file')
        for bnd in bnd_ls:
            output_file = vnm+'_'+bnd+'.bin'
            ### Read diagnotic_vec output
            dv_diag, msk_pts = read_dv_diags(itrs, vnm, bnd, bnd_domain, Nr0, Nr1, grid0, grid1)
            tstp = len(dv_diag)
            ### Get regional boundaary info
            XC1, YC1, mask1, obcs = get_regbnd_info(vnm, bnd, Nr1, tstp, grid1)
            ### Horizontal intrpolation
            for t in range(tstp):
                obcs_tmp = gen_obcs(dv_diag[t], bnd_domain[bnd]['XC'], bnd_domain[bnd]['YC'], msk_pts, XC1, YC1, mask1)
                if bnd in ['east','west']:
                    obcs_tmp = obcs_tmp[:,:,0]
                elif bnd in ['north','south']:
                    obcs_tmp = obcs_tmp[:,0,:]
                obcs[t] = obcs_tmp
            if print_level>=1:
                print(f'       * saving the {bnd} OBCS file')
            obcs.ravel('C').astype('>f4').tofile( os.path.join(out_dir,output_file))
     
    #--------- sea-ice condition ---------#
    if print_level>=1:
        print('> Creating sea-ice OBCS for the '+reg_nm+' model')
    vnms = ['AREA','HEFF','HSNOW','UICE','VICE']
    for vnm in vnms:
        if print_level>=1:
            print(f'    - Processing {vnm} file')
        for bnd in bnd_ls:
            output_file = vnm+'_'+bnd+'.bin'
            ### Read diagnotic_vec output
            dv_diag, msk_pts = read_dv_diags(itrs, vnm, bnd, bnd_domain, Nr0, Nr1, grid0, grid1)
            tstp = len(dv_diag)
            ### Get regional boundaary info
            XC1, YC1, mask1, obcs = get_regbnd_info(vnm, bnd, Nr1, tstp, grid1)
            ### Horizontal intrpolation
            for t in range(tstp):
                obcs_tmp = gen_obcs(dv_diag[t], bnd_domain[bnd]['XC'], bnd_domain[bnd]['YC'], msk_pts, XC1, YC1, mask1)
                if bnd in ['east','west']:
                    obcs_tmp = obcs_tmp[:,:,0]
                elif bnd in ['north','south']:
                    obcs_tmp = obcs_tmp[:,0,:]
                obcs[t] = obcs_tmp
            if print_level>=1:
                print(f'       * saving the {bnd} OBCS file')
            obcs.ravel('C').astype('>f4').tofile( os.path.join(out_dir,output_file))
     
    #--------- ptracer condition ---------#
    if bgc == True:
        if print_level>=1:
            print('> Creating biogeochemical OBCS for the '+reg_nm+' model')
        sfx = os.path.join(config_dir, 'parent/outputs/OBCS/')+f"{bnd_ls[0]}_BC_mask_PTRACE*"
        ptr_ls = glob.glob(sfx)
        ptr_ls = np.unique([int(ptr_ls[i][len(sfx)-1:-15]) for i in range(len(ptr_ls))])
        vnms = [f'PTRACE{i:02d}' for i in ptr_ls]
        for vnm in vnms:
            if print_level>=1:
                print(f'    - Processing {vnm} file')
            for bnd in bnd_ls:
                output_file = vnm+'_'+bnd+'.bin'
                ### Read diagnotic_vec output
                dv_diag, msk_pts = read_dv_diags(itrs, vnm, bnd, bnd_domain, Nr0, Nr1, grid0, grid1)
                tstp = len(dv_diag)
                ### Get regional boundaary info
                XC1, YC1, mask1, obcs = get_regbnd_info(vnm, bnd, Nr1, tstp, grid1)
                ### Horizontal intrpolation
                for t in range(tstp):
                    obcs_tmp = gen_obcs(dv_diag[t], bnd_domain[bnd]['XC'], bnd_domain[bnd]['YC'], msk_pts, XC1, YC1, mask1)
                    if bnd in ['east','west']:
                        obcs_tmp = obcs_tmp[:,:,0]
                    elif bnd in ['north','south']:
                        obcs_tmp = obcs_tmp[:,0,:]
                    obcs[t] = obcs_tmp
                if print_level>=1:
                    print(f'       * saving the {bnd} OBCS file')
                obcs.ravel('C').astype('>f4').tofile( os.path.join(out_dir,output_file))
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
    parser.add_argument("-bnd", "--bnd_nm", action="store",
                        help="boundaries where to generate a mask", dest="bnd_nm",
                        type=str, required=True)
    parser.add_argument("-i", "--itr", nargs='+', action="store",
                        help="iteration of the begining the regional model", dest="itr",
                        type=int, required=True)
    parser.add_argument("-bgc", "--darwin", action="store_true", default='False',
                        help="generate darwin biogeochemistry pickups")
    parser.add_argument("-v", "--verbose", action="store_true", default='False')
    

    args = parser.parse_args()
    config_dir = args.config_dir
    reg_nm = args.reg_nm
    bds = args.bnd_nm
    itrs = args.itr
    bgc = args.darwin
    if args.verbose == True:
        print_level = 1
    else:
        print_level = 0

    gen_obcs_files(config_dir, reg_nm, bds, itrs, bgc, print_level)
























