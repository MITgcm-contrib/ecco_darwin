import os
import argparse
import numpy as np
import xarray as xr
from MITgcmutils import mds, llc

############################################################
#                    READING FUNCTIONS                     #
############################################################

def gen_delYFile(config_dir, model_name):
    ds = xr.open_dataset(os.path.join(config_dir,model_name+'_ncgrid.nc'))
    YG = ds.YG.values
    ds.close()
    delY = np.diff(YG[:, 0]).reshape(-1, 1)
    print(delY)
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

    args = parser.parse_args()
    config_dir = args.config_dir
    model_name = args.model_name
    
    gen_delYFile(config_dir, model_name)
