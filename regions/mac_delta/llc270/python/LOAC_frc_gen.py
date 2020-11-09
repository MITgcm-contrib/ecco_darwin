#!/usr/bin/env python
# coding: utf-8

###############################################################################
# Generate Mac270 forcing files from llc270 forcing file
#
# Clément BERTIN, Université de La Rochelle
# clement.bertin1@univ-lr.fr
# Octobre 8th, 2020
###############################################################################

# Basic Tools
import os
import sys
import time
import numpy as np
from struct import *
from calendar import monthrange
from progressbar import ProgressBar
from llc2Mack import RnConv

############################################################
# Functions
############################################################
def leap_year(year):
    if monthrange(year, 2)[1] == 29:
        day_nb = 366
    else:
        day_nb = 365
    return(day_nb)

def bin_save(fnm, array):
    newFile = open(fnm, "wb")
    bitearray = pack(">%sf" % len(array.flatten()),*array.flatten())
    newFile.write(bitearray);

############################################################
# Read llc270 files (by year) and generate Mack270
# corresponding file
############################################################

# Set path to data
inpath = input("Set the path to llc270 Yearly files\n")
outpath = inpath+"Mack270/"
fnm = input("Set generic name of file (w/o year)\n")
fnm_out = input("Set generic name of output file (w/o year)\n")

# Set date window
st_dt = int(input("Set starting date\n"))
end_dt = int(input("Set ending date\n"))

# Ceate a folder to save Mack270 files
try:
    os.mkdir(outpath)
except OSError:
    print ("Creation of the directory %s failed" % outpath)
else:
    print ("Successfully created the directory %s" % outpath)

# Read file names in llc270 directory
# (Make sure there is only llc270 files to reshape in the folder)
yer_ls = np.arange(st_dt, end_dt+1)
itr_nb = len(yer_ls); itr = 0

# Read and reshape llc270 file
with ProgressBar(max_value=itr_nb) as bar:
    for i in range(itr_nb):
        # Leap year verification
        ndy = leap_year(yer_ls[i])
        # Read and zoom llc270 files
        Mac_data = RnConv(inpath+fnm+str(yer_ls[i]), ndy)
        # Save JRA55 Mack270 files
        bin_save(outpath+fnm_out+str(yer_ls[i]), Mac_data)
        # Update progressbar
        time.sleep(0.1)
        itr+=1
        bar.update(itr)
