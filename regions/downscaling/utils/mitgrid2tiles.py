import os
import argparse
import numpy as np

def split_mitgrid(config_dir, model_name, size_dom, size_proc):

    print('Splitting the mitgrid file')
    # Read parameters
    n_rows = size_dom[0]; rows = size_proc[0]
    n_cols = size_dom[1]; cols = size_proc[1]

    # read in the grid subset to interpolate onto
    mitgrid_file = os.path.join(config_dir,model_name+'.mitgrid')
    entire_grid = np.fromfile(mitgrid_file, dtype='>f8')
    entire_grid = np.reshape(entire_grid, (16, n_cols + 1, n_rows + 1))

    #Create a tiles folder
    if 'tiles' not in os.listdir(config_dir):
        os.mkdir(os.path.join(config_dir,'tiles'))
        
    #Generate the tiles tiled
    nr = int(n_rows/rows)
    nc = int(n_cols/cols)
    counter = 1
    for ri in range(nr):
        for ci in range(nc):
            tile_subset = entire_grid[:,cols*ci:cols*(ci+1)+1
                                       ,rows*ri:rows*(ri+1)+1]
            output_file = os.path.join(config_dir,'tiles','tile'+'{:03d}'.format(counter)+'.mitgrid')
            tile_subset.ravel('C').astype('>f8').tofile(output_file)
            counter+=1

if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument("-d", "--config_dir", action="store",
                        help="The directory where to store the mitgrid", dest="config_dir",
                        type=str, required=True)
    parser.add_argument("-n", "--reg_nm", action="store",
                        help="Name of the regional cutout", dest="reg_nm",
                        type=str, required=True)
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
    size_dom = args.size_dom
    size_proc = args.size_proc

    split_mitgrid(config_dir, reg_nm, size_dom, size_proc)
    print("mitgrid files generated")