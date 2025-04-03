import os
import numpy as np
import netCDF4 as nc4
import matplotlib.pyplot as plt
import argparse
import sys
import ast


def create_mitgrid(config_dir, print_level):

    L1_model_name = 'L1_mac_delta'

    sys.path.insert(1, os.path.join(config_dir, 'utils', 'init_file_creation'))

    x0 = 206.8565
    y0 = 68.4727
    dx = 0.0287
    dy = 0.0091
    Nx = 1224
    Ny = 744

    xg = np.arange(x0,x0+Nx*dx+dx,dx)
    yg = np.arange(y0,y0+Ny*dy+dy,dy)
    XG, YG = np.meshgrid(xg,yg)

    xc = np.arange(x0+dx/2, x0 + Nx * dx, dx)
    yc = np.arange(y0+dy/2, y0 + Ny * dy, dy)
    XC, YC = np.meshgrid(xc, yc)

    # be sure to flip these upside down because simplegrid is weirdly indexed
    XC = np.flipud(XC)
    YC = np.flipud(YC)
    XG = np.flipud(XG)
    YG = np.flipud(YG)

    # ds = nc4.Dataset(os.path.join(config_dir,'nc_grids',L1_model_name+'_grid_old.nc'))
    # XG_test = ds.variables['XG'][:, :]
    # YG_test = ds.variables['YG'][:, :]
    # XC_test = ds.variables['XC'][:, :]
    # YC_test = ds.variables['YC'][:, :]
    # ds.close()
    #
    # plt.subplot(2,2,1)
    # plt.imshow(XG-XG_test,origin='lower')
    # plt.title('Max Difference: '+str(np.max(np.abs(XG-XG_test))))
    # plt.subplot(2, 2, 2)
    # plt.imshow(YG - YG_test, origin='lower')
    # plt.title('Max Difference: ' + str(np.max(np.abs(YG - YG_test))))
    # plt.subplot(2, 2, 3)
    # plt.imshow(XC - XC_test, origin='lower')
    # plt.title('Max Difference: ' + str(np.max(np.abs(XC - XC_test))))
    # plt.subplot(2, 2, 4)
    # plt.imshow(YC - YC_test, origin='lower')
    # plt.title('Max Difference: ' + str(np.max(np.abs(YC - YC_test))))
    # plt.show()

    import create_mitgrid_from_XY as cm
    cm.create_mitgrid_file(config_dir,L1_model_name,XC,YC,XG,YG,print_level)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument("-d", "--config_dir", action="store",
                        help="The directory where the L1, L1, and L1 configurations are stored.", dest="config_dir",
                        type=str, required=True)

    args = parser.parse_args()
    config_dir = args.config_dir

    create_mitgrid(config_dir, print_level=5)


