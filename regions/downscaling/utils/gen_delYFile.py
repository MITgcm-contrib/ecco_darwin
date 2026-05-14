import os
import argparse
import numpy as np
import xarray as xr
from MITgcmutils import mds, llc

############################################################
#                    READING FUNCTIONS                     #
############################################################

def gen_delYFile(config_dir, model_name, n_rows, n_cols):
    mitgrid_file = os.path.join(config_dir,model_name+'.mitgrid')
    entire_grid = np.fromfile(mitgrid_file, dtype='>f8')
    entire_grid = np.reshape(entire_grid, (16, n_cols + 1, n_rows + 1))
    YG = entire_grid[6, :, :]
    delY = np.diff(YG[:, 0]).reshape(-1, 1)
    print(delY)
    print(delY.shape)
    output_file = os.path.join(config_dir,'delYFile')
    delY.astype('>f8').tofile(output_file)

############################################################
#                         PARSER                          #
############################################################

if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument("-d", "--config_dir", action="store",
                        help="The directory where parents folder is stored", dest="config_dir",
                        type=str, required=True)
    parser.add_argument("-n", "--model_name", action="store",
                        help="Name of the regional cutout", dest="model_name",
                        type=str, required=True)
    parser.add_argument("-x", "--n_rows", action="store",
                        help="grid size in x direction", dest="n_rows",
                        type=int, required=True)
    parser.add_argument("-y", "--n_cols", action="store",
                        help="grid size in y direction", dest="n_cols",
                        type=int, required=True)
    
    args = parser.parse_args()
    config_dir = args.config_dir
    model_name = args.model_name
    n_rows = args.n_rows
    n_cols = args.n_cols
    
    gen_delYFile(config_dir, model_name, n_rows, n_cols)
