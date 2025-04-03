import os
import ast
import argparse
import numpy as np
import xarray as xr

def gen_nc(Hvnm, Hgrd_st, Zvnm,Zgrd_st):
    ### Create coordinates dictionary
    dims = dict(X=np.arange(Hgrd_st[0].shape[1]), Y=np.arange(Hgrd_st[0].shape[0]),
                Xp1=np.arange(Hgrd_st[7].shape[1]), Yp1=np.arange(Hgrd_st[9].shape[0]),
                z=np.arange(Zgrd_st[4].shape[0]))
    ### Create dictionary
        # Horizontal variables
    data = {}
    for i,nm in enumerate(Hvnm):
        if i >=7 and i<9:
            data.update({nm: (["Y","Xp1"], Hgrd_st[i])})
        elif i>=9 and i<11:
            data.update({nm: (["Yp1","X"], Hgrd_st[i])})
        else:
            data.update({nm: (["Y","X"], Hgrd_st[i])})   
        # Verical variables
    data.update({Zvnm[0]: (["Y","X"], Zgrd_st[0])})
    data.update({Zvnm[1]: (["Z","Y","X"], Zgrd_st[1])})
    data.update({Zvnm[2]: (["Z","Yp1","X"], Zgrd_st[2])})
    data.update({Zvnm[3]: (["Z","Y","Xp1"], Zgrd_st[3])})
    data.update({Zvnm[4]: (["Z"], Zgrd_st[4])})
    ### Generate the xarra
    ds = xr.Dataset(data_vars=data, coords=dims)
    return ds

def stitch_Hgrd(config_dir, model_name, sNx, sNy, tiles, vnm, grd_st):
    mnc_dir = os.path.join(config_dir, 'mncs/')
    N = tiles.shape[0]*tiles.shape[1]
    ### Stictch the tiles
    for r in range(tiles.shape[0]):
        for c in range(tiles.shape[1]):
            tile_number = int(tiles[r][c])
            for n in range(N):
                if 'grid.t'+'{:03d}'.format(tile_number)+'.nc' in os.listdir(os.path.join(mnc_dir,'mnc_'+'{:04d}'.format(n+1))):
                    ds = xr.open_dataset(os.path.join(mnc_dir, 'mnc_' + '{:04d}'.format(n + 1),
                                                      'grid.t'+'{:03d}'.format(tile_number) + '.nc'))
                    for i,nm in enumerate(vnm):
                        if i>=2 and i<4:
                            grd_st[i][r*sNy:(r+1)*sNy, c*sNx:(c+1)*sNx] = ds[nm].data
                            grd_st[i][r*sNy:(r+1)*sNy, c*sNx:(c+1)*sNx][:,-1:] = grd_st[i][r*sNy:(r+1)*sNy, c*sNx:(c+1)*sNx][:,-2:-1]
                            grd_st[i][r*sNy:(r+1)*sNy, c*sNx:(c+1)*sNx][-1:] = grd_st[i][r*sNy:(r+1)*sNy, c*sNx:(c+1)*sNx][-2:-1]
                        elif i>=5 and i<7:
                            grd_st[i][r*sNy:(r+1)*sNy, c*sNx:(c+1)*sNx] = ds[nm][:-1,:-1].data
                        elif i>=7 and i<9:
                            grd_st[i][r*sNy:(r+1)*sNy, c*sNx:(c+1)*sNx+1] = ds[nm].data
                        elif i>=9 and i<11:
                            grd_st[i][r*sNy:(r+1)*sNy+1, c*sNx:(c+1)*sNx] = ds[nm].data
                        else:
                            grd_st[i][r*sNy:(r+1)*sNy, c*sNx:(c+1)*sNx] = ds[nm].data
    ### Correcting articfacts
    grd_st[7][:,:1] = grd_st[7][:,1:2]; grd_st[9][:1] = grd_st[9][1:2]
    return grd_st

def stitch_Zgrd(config_dir, model_name, sNx, sNy, tiles, vnm, grd_st):
    mnc_dir = os.path.join(config_dir, 'mncs/')
    N = tiles.shape[0]*tiles.shape[1]
    ### Stictch the tiles
    for r in range(tiles.shape[0]):
        for c in range(tiles.shape[1]):
            tile_number = int(tiles[r][c])
            for n in range(N):
                if 'grid.t'+'{:03d}'.format(tile_number)+'.nc' in os.listdir(os.path.join(mnc_dir,'mnc_'+'{:04d}'.format(n+1))):
                    ds = xr.open_dataset(os.path.join(mnc_dir, 'mnc_' + '{:04d}'.format(n + 1),
                                                      'grid.t'+'{:03d}'.format(tile_number) + '.nc'))
                    grd_st[0][r*sNy:(r+1)*sNy, c*sNx:(c+1)*sNx] = ds[vnm[0]].data
                    grd_st[1][:, r*sNy:(r+1)*sNy, c*sNx:(c+1)*sNx] = ds[vnm[1]].data
                    grd_st[2][:, r*sNy:(r+1)*sNy+1, c*sNx:(c+1)*sNx] = ds[vnm[2]].data
                    grd_st[3][:, r*sNy:(r+1)*sNy, c*sNx:(c+1)*sNx+1] = ds[vnm[3]].data
    grd_st[4] = ds[vnm[4]].data
    return grd_st

def gen_grd(config_dir, model_name, Nr, sNx, sNy, Htiles, Ztiles):    
    ### Prepare Horizontal grid
    Hgrd_st = []
    Hvnm = ['XC','YC','AngleCS','AngleSN','rA',
           'XG','YG','dxC','dyG','dyC','dxG']
    Htiles = np.array(Htiles)
    for i,nm in enumerate(Hvnm):
        if i>=7 and i<9:
            Hgrd_st.append(np.zeros((sNy*Htiles.shape[0],sNx*Htiles.shape[1]+1)))
        elif i>=9 and i<11:   
            Hgrd_st.append(np.zeros((sNy*Htiles.shape[0]+1,sNx*Htiles.shape[1])))
        else:
            Hgrd_st.append(np.zeros((sNy*Htiles.shape[0],sNx*Htiles.shape[1])))
    ### Prepare Vertical grid varibale
    Zgrd_st = []
    Zvnm = ['Depth','HFacC','HFacS','HFacW','drF']
    Ztiles = np.array(Ztiles)
    Zgrd_st.append(np.zeros((sNy*Ztiles.shape[0],sNx*Ztiles.shape[1])))
    Zgrd_st.append(np.zeros((Nr,sNy*Ztiles.shape[0],sNx*Ztiles.shape[1])))
    Zgrd_st.append(np.zeros((Nr,sNy*Ztiles.shape[0]+1,sNx*Ztiles.shape[1])))
    Zgrd_st.append(np.zeros((Nr,sNy*Ztiles.shape[0],sNx*Ztiles.shape[1]+1)))
    Zgrd_st.append(np.zeros(Nr))
    #### Stitch the grid files
    Hgrd_st = stitch_Hgrd(config_dir, model_name, sNx, sNy, Htiles, Hvnm, Hgrd_st)
    Zgrd_st = stitch_Zgrd(config_dir, model_name, sNx, sNy, Ztiles, Zvnm, Zgrd_st)
    ### Generate a nc file
    ds = gen_nc(Hvnm, Hgrd_st, Zvnm, Zgrd_st)
    return ds

def stitching(config_dir, rg_nm, z_num, size_dom, size_proc):
    print('Stitching the nc grid files')
    Nr = z_num
    sNx = size_proc[0]; sNy = size_proc[1]
    nPx = int(size_dom[0]/sNx); nPy = int(size_dom[1]/sNy)
    Htiles = np.cumsum(np.ones(nPx*nPy)).reshape(nPx,nPy).T
    Ztiles = np.cumsum(np.ones(nPx*nPy)).reshape(nPy,nPx)
    ### Read the grid information and generate the netcdf file
    ds = gen_grd(config_dir, rg_nm, Nr, sNx, sNy, Htiles, Ztiles)
    ### Save the netcdf file
    ds.to_netcdf(os.path.join(config_dir, rg_nm+'_ncgrid.nc'), mode='w')
    return


if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument("-d", "--config_dir", action="store",
                        help="The directory where .... is stored", dest="config_dir",
                        type=str, required=True)
    parser.add_argument("-n", "--reg_nm", action="store",
                        help="Name of the regional cutout", dest="reg_nm",
                        type=str, required=True)
    parser.add_argument("-z", "--z_num", action="store",
                        help="Number of vertical levels: Nr", 
                        dest="z_num", type=int, required=True)
    parser.add_argument("-s", "--size_dom", nargs='+', action="store",
                        help="list of the mitgrid shape: \
                        [number of rows, number of colums]", 
                        dest="size_dom", type=int, required=True)
    parser.add_argument("-p", "--size_proc", nargs='+', action="store",
                        help="list of the processor size: \
                        [sNx, SNy] from Size.h file", 
                        dest="size_proc", type=int, required=True)

    args = parser.parse_args()
    config_dir = args.config_dir
    reg_nm = args.reg_nm
    z_num = args.z_num
    size_dom = args.size_dom
    size_proc = args.size_proc

    stitching(config_dir, reg_nm, z_num, size_dom, size_proc)
    print('nc grid files generated')