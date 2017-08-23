"""
Created on Thurs Jul 14 2016

@author: rachel
"""

import os, errno
from netCDF4 import Dataset
import netCDF4
import numpy as np
import datetime as dt
import Ngl
import xray
import constants as const
import matplotlib.pyplot as plt
from scipy.interpolate import griddata

startyr = ["0001", "00010"]
endyr = ["0035","0049"]

Expertype = "CAM4POP_f45g37"
ExperimentCTL = ""

Dir = "/home/disk/rachel/CESM_outfiles/CAM4POP_CTL/atm/hist/"
filenameD = Dir + "/AnnAvg_CAM4POP_CTL.cam2.h0.0001-0011.nc"

print filenameD
filein= xray.open_dataset(filenameD,decode_times=False)
TLON = filein['TLONG'].values
TLAT = filein['TLAT'].values

nlats = filein.dims['nlat']
nlons = filein.dims['nlon']
print nlats,nlons

SurfArea = filein['TAREA'] / (100.0 * 100.0) 	# convert from cm2 to m2

MaskArea = filein['HBLT'].values

OceanMask = np.where(np.isnan(MaskArea),0,1)

totarea = np.nansum(SurfArea)
#check:
if np.nansum(totarea) < 5.05e+14 or np.nansum(totarea) > 5.06e+14:
	exit("Total Earth Surface Area not 5.05e+14m2")

OceanArea = np.squeeze(OceanMask * SurfArea.values)

#check:
if np.nansum(OceanArea)/totarea < 0.7 or np.nansum(OceanArea)/totarea > 0.72:
	exit("ocean fraction is not 0.71!")
#For NA:
# TLAT >= 25.0 and TLAT <= 55
# TLON >= 270 and TLON <=360
NAarea = 0
for ilat in range(0,nlats):
	for ilon in range(0,nlons):
		if TLAT[ilat,ilon] >=25 and TLAT[ilat,ilon] <=55:
			if TLON[ilat,ilon] >= 270 and TLON[ilat,ilon] <=360:
				NAarea += OceanArea[ilat,ilon]

print "NA area, m2: ", NAarea
print "NA total PW, for 10W/m2 = ", NAarea * 10.0 / 1.0e15


# For SH:
# TLAT >= -75 and TLAT <=0.0
# TLON >= 0 and TLON <=360

SHarea = 0
for ilat in range(0,nlats):
        for ilon in range(0,nlons):
                if TLAT[ilat,ilon] >=-75 and TLAT[ilat,ilon] <=0:
                        if TLON[ilat,ilon] >= 0 and TLON[ilat,ilon] <=360:
                                SHarea += OceanArea[ilat,ilon]
print "SH area, m2: ", SHarea
print "SH total PW, for 10W/m2 = ", SHarea * 10.0 / 1.0e15

# For SO:
# TLAT >=-75 and TLAT <=45
# TLON >= 0 and TLON <=360
SOarea = 0
for ilat in range(0,nlats):
        for ilon in range(0,nlons):
                if TLAT[ilat,ilon] >=-75 and TLAT[ilat,ilon] <=-45:
                        if TLON[ilat,ilon] >= 0 and TLON[ilat,ilon] <=360:
                                SOarea += OceanArea[ilat,ilon]
print "SO area, m2: ", SOarea
print "SO total PW, for 10W/m2 = ", SOarea * 10.0 / 1.0e15



