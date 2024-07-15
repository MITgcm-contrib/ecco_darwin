import os
import sys
import argparse
import numpy as np
import simplegrid as sg
from pyproj import Transformer


def gen_XY(min_x, max_x, min_y, max_y, reso_x, reso_y, espg, print_level):

    # Handle the drop from 180 to -180
    if min_x > max_x:
        if espg!=4326: 
            raise ValueError("Can't handle the [-180;180] drop in a different ESPG than WGS84")
        if print_level >= 1:
            print('    - Handling the discontinuity in the longitude')
        init_max_x = max_x
        max_x += 360
    else:
        init_max_x = 0
    
    # Generate XG & YG matrices
    xg = np.arange(min_x, max_x+reso_x, reso_x)
    yg = np.arange(min_y, max_y+reso_y, reso_y)
    XG, YG = np.meshgrid(xg,yg)

    # Generate XC & YC matrices
    xc = np.arange(min_x+(reso_x/2), max_x+(3*reso_x/2), reso_x)
    if len(xg)-len(xc) == 0:
        xc = xc[:-1]
    yc = np.arange(min_y+(reso_y/2), max_y+(3*reso_y/2), reso_y)
    if len(yg)-len(yc) == 0:
        yc = yc[:-1]
    XC, YC = np.meshgrid(xc,yc)

    # correct the drop from 180 to -180 if exists
    if init_max_x != 0:
        XG[XG>180] -= 360
        XC[XC>180] -= 360
    
    return(XG, YG, XC, YC)

def gen_mat(domain_coords):
    mitgrid_matrices = dict()
    mitgrid_matrices['XG'] = domain_coords[2]
    mitgrid_matrices['YG'] = domain_coords[3]
    mitgrid_matrices['XC'] = domain_coords[0]
    mitgrid_matrices['YC'] = domain_coords[1]
    return(mitgrid_matrices)

def gen_file(output_file,mitgrid_matrices,XG,YG,factor):
    mg_new, n_rows, n_cols = sg.regrid.regrid(mitgrid_matrices=mitgrid_matrices,
                                              lon_subscale=factor, lat_subscale=factor,
                                              lon1=XG[0, 0], lat1=YG[0, 0], 
                                              lon2=XG[-1, -1], lat2=YG[-1, -1],
                                              verbose=False, outfile = output_file)
    
    print('    - Output shape: ('+str(n_rows)+', '+str(n_cols)+')')

    sg.gridio.write_mitgridfile(output_file, mg_new, n_rows, n_cols)

def gen_mitgrid_file(config_dir,reg_nm,XC,YC,XG,YG,print_level):

    domain_coords = [XC, YC, XG, YG]
    mitgrid_matrices = gen_mat(domain_coords)

    if print_level >=1:
        print('    - Generating mitgrid')
    output_file = os.path.join(config_dir,reg_nm+'.mitgrid')\

    gen_file(output_file, mitgrid_matrices, XG, YG, factor=1)

def gen_mitgrid(config_dir, reg_nm, corners, reso, espg, print_level):

    sys.path.insert(1, os.path.join(config_dir))
    if print_level >= 1:
        print('Creating the mitgrid')
        print('    - Generating the grid')
    
    # Create coordinate matrices
    XG, YG, XC, YC = gen_XY(corners[0], corners[1], corners[2], corners[3],
                            reso[0], reso[1], espg, print_level)

    # Reproject the grid to lon, lat coordinates
    if espg != 4326:
        if print_level >= 1:
            print('    - Reprojecting the grid to lat/lon')
        transformer = Transformer.from_crs('EPSG:' + str(espg), 'EPSG:' + str(4326))
        
        Lat_G, Lon_G = transformer.transform(XG.ravel(), YG.ravel())
        Lat_G = np.reshape(Lat_G, np.shape(XG))
        Lon_G = np.reshape(Lon_G, np.shape(XG))
        
        Lat_C, Lon_C = transformer.transform(XC.ravel(), YC.ravel())
        Lat_C = np.reshape(Lat_C,np.shape(XC))
        Lon_C = np.reshape(Lon_C, np.shape(XC))
    else:
        Lat_G, Lon_G = YG, XG
        Lat_C, Lon_C = YC, XC
    
    # flip coordinates upside down for simplegrid
    Lon_G = np.flipud(Lon_G); Lat_G = np.flipud(Lat_G)
    Lon_C = np.flipud(Lon_C); Lat_C = np.flipud(Lat_C)

    gen_mitgrid_file(config_dir,reg_nm,Lon_C,Lat_C,Lon_G, Lat_G,print_level)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument("-d", "--config_dir", action="store",
                        help="The directory where to store the mitgrid", dest="config_dir",
                        type=str, required=True)
    parser.add_argument("-n", "--reg_nm", action="store",
                        help="Name of the regional cutout", dest="reg_nm",
                        type=str, required=True)
    parser.add_argument("-c", "--corners", nargs='+', action="store",
                        help="list of the coordinates of the corners: \
                        left longitude right longitude left latitude right latitude]", 
                        dest="corners", type=float, required=True)
    parser.add_argument("-r", "--reso", nargs='+', action="store",
                        help="list of the longitude and latitude increment: \
                        longitude increment latitude increment]", 
                        dest="reso", type=float, required=True)
    parser.add_argument("-e", "--espg", action="store",
                        help="ESPG of the input coordinates (default is WGS84; ESPG:4326)", 
                        dest="espg", type=int, required=False, default=4326)
    parser.add_argument("-v", "--verbose", action="store_true", default='False')
    
    args = parser.parse_args()
    config_dir = args.config_dir
    reg_nm = args.reg_nm
    corners = args.corners
    reso = args.reso
    espg = args.espg
    if args.verbose == True:
    	print_level = 1
    else:
    	print_level = 0

    gen_mitgrid(config_dir, reg_nm, corners, reso, espg, print_level)
    print("mitgrid file generated")