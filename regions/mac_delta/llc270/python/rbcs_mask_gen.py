#!/usr/bin/env python
# coding: utf-8

###############################################################################
# Make RBCS mask
#
# Clément BERTIN, Université de La Rochelle
# clement.bertin1@univ-lr.fr
# June 2nd, 2019
###############################################################################
import numpy as np
from struct import *

####################
# > Parameters
####################
nx=46
ny=68
mask=np.zeros([nx,ny])
relaxMaskAFile="Mac_rbcs_mask.bin"
#from data.obcs
seaiceSpongeThickness = 8
Arelaxobcsinner = 691200.0
Arelaxobcsbound = 86400.0
#from data.rbcs
tauRelaxA=691200.

####################
# > Compute
####################
dt = (Arelaxobcsinner-Arelaxobcsbound)/(seaiceSpongeThickness-1)
Teff = Arelaxobcsbound + np.arange(0,seaiceSpongeThickness)*dt
mm = tauRelaxA/Teff

####################
# > Create mask
####################
# West Sponge
for i in range(seaiceSpongeThickness):
    mask[i,:]=mm[i]
# South Sponge
for i in range(seaiceSpongeThickness):
    mask[i:,i]=mm[i]
# Modify mask dtype
mask = np.asfortranarray(mask.T, dtype=np.float32)

print(mask)

####################
# > Save binary file
####################
newFile = open(relaxMaskAFile, "wb")
bitearray = pack(">%sf" % len(mask.flatten()),*mask.flatten())
newFile.write(bitearray);
