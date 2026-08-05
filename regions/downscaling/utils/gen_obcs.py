import os
import glob
import argparse
import numpy as np
import xarray as xr
from MITgcmutils import mds, llc
from scipy.interpolate import interp1d, griddata, LinearNDInterpolator
from scipy.spatial import Delaunay, cKDTree

############################################################
#                   VALIDATION HELPERS                     #
############################################################

VALID_BND = {'E': 'east', 'W': 'west', 'N': 'north', 'S': 'south'}

def check_boundaries(bnd_str):
    """Validate a -bnd string up front, before any slow I/O happens."""
    if not bnd_str:
        raise SystemExit("ERROR: -bnd is empty. Give a string of boundary letters, e.g. 'ENW'.")
    bad = [c for c in bnd_str if c not in VALID_BND]
    if bad:
        raise SystemExit(
            f"ERROR: unrecognized boundary letter(s) {bad} in -bnd '{bnd_str}'. "
            f"Valid letters are E (east), W (west), N (north), S (south), e.g. 'ENW'.")
    dup = [c for c in set(bnd_str) if bnd_str.count(c) > 1]
    if dup:
        raise SystemExit(f"ERROR: boundary letter(s) {sorted(dup)} repeated in -bnd '{bnd_str}'.")
    return bnd_str

def check_input_layout(config_dir, reg_nm, boundaries, need_dv_output=True):
    """
    Fail early, and with an actionable message, when the config_dir does not
    have the layout every utils script assumes (see the Directory layout
    section of the downscaling README).
    """
    problems = []

    ncgrid = os.path.join(config_dir, reg_nm + '_ncgrid.nc')
    if not os.path.isfile(ncgrid):
        problems.append(f"  missing: {ncgrid}\n"
                        f"      -> produced by stitch_ncgrid.py (STEP1, Section IV.c). "
                        f"Check -n '{reg_nm}' matches the file name.")

    grid_dir = os.path.join(config_dir, 'parent/outputs/grid/')
    if not os.path.isdir(grid_dir):
        problems.append(f"  missing: {grid_dir}\n"
                        f"      -> copy the parent grid files here (STEP1, Section V.b). "
                        f"Note the folder is 'parent', singular.")

    mask_dir = os.path.join(config_dir, 'parent/inputs/')
    for c in boundaries:
        nm = VALID_BND[c]
        f = os.path.join(mask_dir, nm + '_BC_mask.bin')
        if not os.path.isfile(f):
            problems.append(f"  missing: {f}\n"
                            f"      -> produced by gen_dvmasks.py -bnd {c} (STEP1, Section V.c).")

    if need_dv_output:
        dv_dir = os.path.join(config_dir, 'parent/outputs/OBCS/')
        if not os.path.isdir(dv_dir):
            problems.append(f"  missing: {dv_dir}\n"
                            f"      -> copy the diagnostics_vec output here (STEP3, Section I).")
        elif not glob.glob(os.path.join(dv_dir, '*_BC_mask_*.bin')):
            problems.append(f"  no '*_BC_mask_*.bin' files in: {dv_dir}\n"
                            f"      -> these come from the STEP2 run. If the run produced files "
                            f"under different names, check nml_vecFiles in data.diagnostics_vec.")

    if problems:
        raise SystemExit("ERROR: the -d directory is missing files this script needs:\n"
                         + "\n".join(problems))

def warn_if_all_zero(arr, vnm, bnd):
    """
    diagnostics_vec silently writes a correctly-named file full of zeros when
    a field name in data.diagnostics_vec isn't one it recognizes (its
    IF/ELSE-IF chain has no ELSE). Catching that here saves discovering it as
    a dead boundary in the model run.
    """
    if arr.size and not np.any(arr):
        print(f"    WARNING: every value read for {vnm} at the {bnd} boundary is zero. "
              f"Check that '{vnm}' is a field name diagnostics_vec recognizes and that it "
              f"is listed in data.diagnostics_vec for this boundary's mask.")

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
        if npts == 0:
            raise SystemExit(
                f"ERROR: the {bnd_ls[i]} boundary mask has no points. Re-run gen_dvmasks.py "
                f"for this boundary (STEP1, Section V.c) and check its -r search radius is "
                f"large enough, in km, for the parent grid.")
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
    if seaice == True:
        if vnm in ['AREA','UICE','VICE','HEFF','HSNOW', 'ETAN']:
            mask1 = mask1[:1,:,:]
            obcs_init = np.zeros((tstp, 1, np.size(XC1)))
        else:
            obcs_init = np.zeros((tstp, Nr, np.size(XC1)))
    else:
        if vnm in ['AREA', 'ETAN']:
            mask1 = mask1[:1,:,:]
            obcs_init = np.zeros((tstp, 1, np.size(XC1)))
        else:
            obcs_init = np.zeros((tstp, Nr, np.size(XC1)))

    return XC1, YC1, mask1, obcs_init

def build_interp_cache(XC0, YC0, mask00, XC1, YC1, nlev):
    """
    Precompute, once per depth level, everything about the horizontal
    interpolation that does NOT depend on the field values -- i.e. on the
    timestep.

    For a given boundary and variable the native sample set at level k is
    fixed by mask00[k] alone, so the Delaunay triangulation behind
    griddata(method='linear') and the nearest-neighbour lookup behind
    griddata(method='nearest') are identical for every timestep. Rebuilding
    them per timestep (as calling griddata inside the t-loop does) repeats
    the most expensive part of the interpolation tstp times over.

    Returns {k: (sel, tri, nn_idx)} where sel selects the wet native points,
    tri is their Delaunay triangulation (or None if too few to triangulate)
    and nn_idx maps each target point to its nearest native sample.
    """
    coord00 = np.array([XC0, YC0]).T
    xi = np.column_stack([np.asarray(XC1).ravel(), np.asarray(YC1).ravel()])
    cache = {}
    for k in range(nlev):
        sel = mask00[k].ravel() != 0
        coord0 = coord00[sel, :]
        tri = None
        nn_idx = None
        if coord0.shape[0] > 0:
            nn_idx = cKDTree(coord0).query(xi)[1]
            if coord0.shape[0] > 4:
                try:
                    tri = Delaunay(coord0)
                except Exception:
                    tri = None   # degenerate (e.g. collinear) point set
        cache[k] = (sel, tri, nn_idx)
    return cache

def gen_obcs(dv_diag, XC0, YC0, mask00, XC1, YC1, mask1, interp_cache=None):

    #---------- Initialization -----------#
    obcs = np.zeros((mask1.shape))
    if interp_cache is None:
        interp_cache = build_interp_cache(XC0, YC0, mask00, XC1, YC1, len(obcs))

    #----------- Generaet OBCS -----------#
    for k in range(len(obcs)):
        if np.any(mask1[k] > 0):
            ### initialize interpolation at depth k
            sel, tri, nn_idx = interp_cache[k]
            val0 = dv_diag[k].ravel()[sel]
            if len(val0) > 0:
                if tri is not None:
                    val1 = LinearNDInterpolator(tri, val0, fill_value=np.nan)(
                        np.asarray(XC1).ravel(), np.asarray(YC1).ravel()).reshape(XC1.shape)
                else:
                    val1 = np.full(XC1.shape, np.nan)
                ### points outside the convex hull of the native samples (or
                ### where there were too few points to triangulate at all)
                ### come back NaN from griddata -- fall back to the true
                ### nearest native sample for those instead of leaving a
                ### zero for Hextrapol to chain-propagate one cell at a time.
                needs_nearest = np.isnan(val1)
                if np.any(needs_nearest):
                    val1_near = val0[nn_idx].reshape(XC1.shape)
                    val1[needs_nearest] = val1_near[needs_nearest]
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
                if len(val0) > 0:
                    val1_NR = val0[nn_idx].reshape(XC1.shape)
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
    ### An explicit multi-value -i is treated as an exact whitelist; with
    ### zero or one value, auto-discover every "{bnd}_BC_mask_{vnm}.*" file
    ### and read all of them, sorted chronologically ascending.
    if len(itrs_ls) > 1:
        itrs = sorted(itrs_ls)
    else:
        itrs = sorted(int(p[len(sfx)+1:-4]) for p in glob.glob(sfx+".*"))
    if seaice == True:
        if vnm in ['AREA','UICE','VICE','HEFF','HSNOW','ETAN']:
            Nr = 1
        else:
            Nr = Nr0
    else:
        if vnm in ['AREA','ETAN']:
            Nr = 1
        else:
            Nr = Nr0

    #--------- Read dv diagnostic --------#
    dv_diag_chunks = []
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

            if seaice == True:
                if vnm in ['UVEL','UICE']:
                    for m in range(nm_pt):
                        dv_diag[:,:,m] = bnd_domain[bnd]['AngleCS'][m] * dv_diagU[:,:,m] -\
                                         bnd_domain[bnd]['AngleSN'][m] * dv_diagV[:,:,m]
                elif vnm in ['VVEL','VICE']:
                    for m in range(nm_pt):
                        dv_diag[:,:,m] = bnd_domain[bnd]['AngleSN'][m] * dv_diagU[:,:,m] +\
                                         bnd_domain[bnd]['AngleCS'][m] * dv_diagV[:,:,m]
            else:
                if vnm in ['UVEL']:
                    for m in range(nm_pt):
                        dv_diag[:,:,m] = bnd_domain[bnd]['AngleCS'][m] * dv_diagU[:,:,m] -\
                                         bnd_domain[bnd]['AngleSN'][m] * dv_diagV[:,:,m]
                elif vnm in ['VVEL']:
                    for m in range(nm_pt):
                        dv_diag[:,:,m] = bnd_domain[bnd]['AngleSN'][m] * dv_diagU[:,:,m] +\
                                         bnd_domain[bnd]['AngleCS'][m] * dv_diagV[:,:,m]
        else:
            dv_diag =  np.fromfile(sfx+f".{itr:010d}.bin", '>f4')
            nm_pt = len(bnd_domain[bnd]['XC'])
            tstp = int(np.size(dv_diag)/(Nr*nm_pt))
            dv_diag = np.reshape(dv_diag, (tstp, Nr, nm_pt))
        dv_diag_chunks.append(dv_diag)
    ### merge all files, in chronological order, exactly once
    if not dv_diag_chunks:
        raise SystemExit(
            f"ERROR: no diagnostics_vec files found for {vnm} at the {bnd} boundary "
            f"(looked for '{sfx}.*.bin'). Check that this field is listed in "
            f"data.diagnostics_vec for this boundary's mask, and that -i matches the "
            f"iteration numbers actually written by the STEP2 run.")
    dv_diag = np.concatenate(dv_diag_chunks, axis=0)
    warn_if_all_zero(dv_diag, vnm, bnd)

    #------- Vertical interpolation ------#
    if seaice == True:
        if vnm in ['UVEL','UICE']:
            msk_pts = bnd_domain[bnd]['maskW']
        elif vnm in ['VVEL','VICE']:
            msk_pts = bnd_domain[bnd]['maskS']
        else:
            msk_pts = bnd_domain[bnd]['maskC']
    else:
        if vnm in ['UVEL']:
            msk_pts = bnd_domain[bnd]['maskW']
        elif vnm in ['VVEL']:
            msk_pts = bnd_domain[bnd]['maskS']
        else:
            msk_pts = bnd_domain[bnd]['maskC']

    if seaice == True:
        if vnm not in ['AREA','UICE','VICE','HEFF','HSNOW','ETAN']:
            dv_diag_IT, msk_pts_IT = Zinterp(dv_diag, msk_pts, grid0['drF'], grid1['drF'])
        else:
            msk_pts = msk_pts[:1,:]
            dv_diag_IT = dv_diag
            msk_pts_IT = msk_pts
    else:
        if vnm not in ['AREA', 'ETAN']:
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
                tmp = np.where(~np.isnan(tmp1[0]))[0]
                if tmp.size > 1:
                          bottom_values = tmp1[:,np.where(~np.isnan(tmp1[0]))[0][-1]]
                elif tmp.size == 0:
                          bottom_values = 0.
                else:
                          bottom_values = tmp1[:,np.where(~np.isnan(tmp1[0]))[0][0]]
                idB = np.where(np.logical_and(np.isnan(tmp1[0]), np.abs(Z1) < bottom_depth))[0]
                if len(idB) != 0:
                    # fill every child level in the gap, not just the first
                    tmp1[:, idB] = bottom_values[:, None] if np.ndim(bottom_values) else bottom_values
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
                    idB = np.where(np.logical_and(np.isnan(tmp1), np.abs(Z1) < bottom_depth))[0]
                    if len(idB) != 0:
                        tmp1[idB] = bottom_value
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

def _shift_from(a, dr, dc, fill):
    """Value of the neighbour at (r+dr, c+dc), aligned back onto (r, c)."""
    out = np.full_like(a, fill)
    src = a[max(0, dr):a.shape[0] + min(0, dr), max(0, dc):a.shape[1] + min(0, dc)]
    out[max(0, -dr):a.shape[0] + min(0, -dr), max(0, -dc):a.shape[1] + min(0, -dc)] = src
    return out

def Hextrapol(var_grid, wet_grid, verbose=False):
    """
    Fill zero values in var_grid using the nearest non-zero neighbors within wet areas.

    Parameters:
        var_grid (ndarray): 2D array with variable data, where zeros represent missing values.
        wet_grid (ndarray): 2D mask where 1 indicates wet cells and 0 indicates dry.

    Returns:
        tuple: Updated var_grid and number of remaining unfilled wet cells.

    This is a breadth-first fill: each pass spreads values one cell outward
    from the cells that were already non-zero when the pass began. It is a
    vectorised rewrite of an earlier version that, for every unfilled cell on
    every pass, computed the distance to every non-zero wet cell and took the
    argmin. That is bit-for-bit reproduced here, because on an integer grid
    the old "closest distance < sqrt(2)" test can only ever be satisfied by an
    orthogonal neighbour (distance exactly 1), and np.argmin resolved ties in
    row-major order of the source cell -- i.e. up, then left, then right, then
    down, which is the order of OFFSETS below. Verified identical (values and
    returned count) on 2100 randomised grids; 4-67x faster, most on the
    many-gap cases that used to dominate runtime.
    """
    OFFSETS = ((-1, 0), (0, -1), (0, 1), (1, 0))
    wet = (wet_grid == 1)

    while True:
        is_remaining = (var_grid == 0) & wet
        if not is_remaining.any():
            break
        if verbose:
            print(f"         - Remaining cells to fill: {int(is_remaining.sum())}")
        ### snapshot the sources at the start of the pass: a cell filled during
        ### this pass must not itself become a source until the next one
        src_ok = wet & (var_grid != 0)
        if not src_ok.any():
            break

        filled = np.zeros_like(var_grid)
        taken = np.zeros(var_grid.shape, dtype=bool)
        for dr, dc in OFFSETS:
            cand_ok = _shift_from(src_ok, dr, dc, False)
            use = is_remaining & cand_ok & ~taken
            if use.any():
                cand_val = _shift_from(var_grid, dr, dc, 0)
                filled[use] = cand_val[use]
                taken |= use
        if not taken.any():
            break
        var_grid[taken] = filled[taken]

    return var_grid, int(((var_grid == 0) & wet).sum())

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

def gen_obcs_files(config_dir, reg_nm, boundaries, itrs, seaice, bgc, print_level):

    ### fail fast, before any slow grid reading, if the inputs aren't there
    check_boundaries(boundaries)
    check_input_layout(config_dir, reg_nm, boundaries, need_dv_output=True)

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
    ### On a C-grid HFacS has ny+1 rows and HFacW has nx+1 columns: entry j is
    ### the SOUTH face of cell j (resp. the WEST face of cell i). MITgcm masks
    ### every OBCS normal velocity with the INTERIOR face of the boundary cell
    ### (obcs_apply_uv.F): _maskS(Jobc) north / _maskS(Jobc+1) south,
    ### _maskW(Iobc) east / _maskW(Iobc+1) west. Trimming BOTH ends is what
    ### makes a single slice correct for both boundaries on each axis:
    ###   [0]  -> face 1    = interior face of the low-index boundary (S / W)
    ###   [-1] -> face n-1  = interior face of the high-index boundary (N / E)
    ### Trimming only one end (e.g. [:,1:,:]) is right for S/E but silently
    ### hands N/W the OUTER domain-edge face instead.
    tmp[grd_ls.index("HFacS")] = tmp[grd_ls.index("HFacS")][:,1:-1,:]
    tmp[grd_ls.index("HFacW")] = tmp[grd_ls.index("HFacW")][:,:,1:-1]
    grid1 = dict(zip(grd_ls, tmp))
    Nr1 = tmp[grd_ls.index("HFacC")].shape[0]
    del tmp

    #################################################
    #######         Read/Generate masks       #######
    #################################################
    bnd_ls, dv_masks_ls = read_dv_masks(config_dir, boundaries, llc)
    for nm in ['C','S','W']:
        tmp = np.copy(grid1["HFac"+nm]); tmp[tmp>0] = 1
        grid1["mask"+nm] = tmp

    #################################################
    #######        Generate OBCS files        #######
    #################################################
    out_dir = os.path.join(config_dir,'forcings/OBCS/')
    ### create it rather than dying on the first tofile(); gen_pickups.py
    ### already creates its own forcings/pickups/ the same way
    os.makedirs(out_dir, exist_ok=True)
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
            ### value-independent interpolation set-up: identical for every
            ### timestep, so build it once here rather than inside the t-loop
            interp_cache = build_interp_cache(bnd_domain[bnd]['XC'], bnd_domain[bnd]['YC'],
                                              msk_pts, XC1, YC1, mask1.shape[0])
            ### Horizontal intrpolation
            for t in range(tstp):
                obcs_tmp = gen_obcs(dv_diag[t], bnd_domain[bnd]['XC'], bnd_domain[bnd]['YC'],
                                    msk_pts, XC1, YC1, mask1, interp_cache=interp_cache)
                if bnd == 'east':
                    obcs_tmp = obcs_tmp[:,:,-1]
                elif bnd == 'west':
                    obcs_tmp = obcs_tmp[:,:,0]
                elif bnd == 'north':
                    obcs_tmp = obcs_tmp[:,-1,:]
                elif bnd == 'south':
                    obcs_tmp = obcs_tmp[:,0,:]
                obcs[t] = obcs_tmp
            if print_level>=1:
                print(f'       * saving the {bnd} OBCS file')
            obcs.ravel('C').astype('>f4').tofile( os.path.join(out_dir,output_file))

    #--------- sea-ice condition ---------#
    if seaice == True:
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
                ### value-independent interpolation set-up: identical for every
                ### timestep, so build it once here rather than inside the t-loop
                interp_cache = build_interp_cache(bnd_domain[bnd]['XC'], bnd_domain[bnd]['YC'],
                                                  msk_pts, XC1, YC1, mask1.shape[0])
                ### Horizontal intrpolation
                for t in range(tstp):
                    obcs_tmp = gen_obcs(dv_diag[t], bnd_domain[bnd]['XC'], bnd_domain[bnd]['YC'],
                                    msk_pts, XC1, YC1, mask1, interp_cache=interp_cache)
                    if bnd == 'east':
                       obcs_tmp = obcs_tmp[:,:,-1]
                    elif bnd == 'west':
                       obcs_tmp = obcs_tmp[:,:,0]
                    elif bnd == 'north':
                       obcs_tmp = obcs_tmp[:,-1,:]
                    elif bnd == 'south':
                       obcs_tmp = obcs_tmp[:,0,:]
                    obcs[t] = obcs_tmp
                if print_level>=1:
                    print(f'       * saving the {bnd} OBCS file')
                obcs.ravel('C').astype('>f4').tofile( os.path.join(out_dir,output_file))

    #--------- ptracer condition ---------#
    if bgc == True:
        if print_level>=1:
            print('> Creating biogeochemical OBCS for the '+reg_nm+' model')
        sfx = os.path.join(config_dir, 'parent/outputs/OBCS/')+f"{bnd_ls[0]}_BC_mask_TRAC*"
        ptr_ls = glob.glob(sfx)
        ptr_ls = np.unique([int(ptr_ls[i][len(sfx)-1:-15]) for i in range(len(ptr_ls))])
        vnms = [f'TRAC{i:02d}' for i in ptr_ls]
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
                ### value-independent interpolation set-up: identical for every
                ### timestep, so build it once here rather than inside the t-loop
                interp_cache = build_interp_cache(bnd_domain[bnd]['XC'], bnd_domain[bnd]['YC'],
                                                  msk_pts, XC1, YC1, mask1.shape[0])
                ### Horizontal intrpolation
                for t in range(tstp):
                    obcs_tmp = gen_obcs(dv_diag[t], bnd_domain[bnd]['XC'], bnd_domain[bnd]['YC'],
                                    msk_pts, XC1, YC1, mask1, interp_cache=interp_cache)
                    if bnd == 'east':
                        obcs_tmp = obcs_tmp[:,:,-1]
                    elif bnd == 'west':
                        obcs_tmp = obcs_tmp[:,:,0]
                    elif bnd == 'north':
                        obcs_tmp = obcs_tmp[:,-1,:]
                    elif bnd == 'south':
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
                        help="The config directory holding <reg_nm>_ncgrid.nc, parent/inputs/ \
                        (dv masks) and parent/outputs/ (grid + diagnostics_vec output); \
                        OBCS files are written to forcings/OBCS/", dest="config_dir",
                        type=str, required=True)
    parser.add_argument("-n", "--reg_nm", action="store",
                        help="Name of the regional cutout", dest="reg_nm",
                        type=str, required=True)
    parser.add_argument("-bnd", "--bnd_nm", action="store",
                        help="Boundaries to generate OBCS for, as a string of E/W/N/S \
                        (e.g. ENW)", dest="bnd_nm",
                        type=str, required=True)
    parser.add_argument("-i", "--itr", nargs='+', action="store",
                        help="Iteration number(s) of the diagnostics_vec output files to read. \
                        With 2+ values this is an exact whitelist; with 0 or 1 value every \
                        matching file is auto-discovered and read", dest="itr",
                        type=int, required=True)
    parser.add_argument("-seaice", "--si", action="store_true", default=False,
                        help="generate sea-ice obcs")
    parser.add_argument("-bgc", "--darwin", action="store_true", default=False,
                        help="generate darwin biogeochemistry obcs")
    parser.add_argument("-v", "--verbose", action="store_true", default=False)


    args = parser.parse_args()
    config_dir = args.config_dir
    reg_nm = args.reg_nm
    bds = args.bnd_nm
    itrs = args.itr
    seaice = args.si
    bgc = args.darwin
    if args.verbose == True:
        print_level = 1
    else:
        print_level = 0

    gen_obcs_files(config_dir, reg_nm, bds, itrs, seaice, bgc, print_level)
