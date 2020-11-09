#!/usr/bin/python
###############################################################################
#  Function to convert global llc270 binary files into Mackenzie llc270
#  binary files zoomed on the region
#
# Clément BERTIN, Université de La Rochelle
# clement.bertin1@univ-lr.fr
# April 17, 2020
###############################################################################

import numpy as np

###############################################################################
#Read Global llc 270 data file for a variable and return data centered on the
#Mackenzie setup region
# > Input :
#   - fname = Path to the llc binary file
#   - nt = number of iterations in binary file
#
# > Output :
#   - fld_Mac = Varible matrix wirh (nt, ny, nx) shape
###############################################################################
def RnConv(fname, nt):
    ########################################
    # Fix selection constants and vectors :
    ########################################
    nX=270; nY=nX*13; nZ=50 #SIze to read global data
    I3 = np.arange(nX*6,(nX*7)+1).astype('int') # Idexes to select face 3
    f=8; # tile number
    I2 = np.arange(0,(nX*3),3).astype('int')+7*nX+f-8 # Idexes to select face 4 / tile 8
    ########################################
    # Read binary file :
    ########################################
    with open(fname, 'rb') as fid:
        tmp = np.fromfile(fid, '>f4')
    fld = tmp.reshape((nt,nY,nX)) # Create a data array
    fld = np.transpose(fld, (2,1,0)) # Reshape the array
    ########################################
    # Select data on Mackenzie area  :
    ########################################
    st7=fld[:,I3,:]
    st8=fld[:,I2,:]
    # create Region of 46x68 fields
    fld_Mac = np.concatenate((st7[231:270,202:270,:],st8[0:7,202:270,:]),axis=0).T
    return(fld_Mac)
