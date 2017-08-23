# -*- coding: utf-8 -*- 
"""
Created on Mon July 11 2016

@author: rachel rachel.2hite@cantab.net
"""


import os, errno
from netCDF4 import Dataset
import netCDF4
import numpy as np
import datetime as dt
import pandas
import xray
import Ngl
from scipy import stats
import math

pi = 3.14159
rE = 6.371E6
Experiment = "CAM4POP_f45g37"

Dir = "/home/disk/eos4/rachel/CESM_outfiles/HYAK/" + Experiment + "/atm/hist/"
filenameD = Dir + "/" + Experiment + ".cam2.h0.0001-01.nc"

print filenameD

FileIn = xray.open_dataset(filenameD)

lats = FileIn['lat']
lons = FileIn['lon']

lats_s = FileIn['slat']

weights = FileIn['w_stag']

#convert to radians
latsr = lats * pi / 180.0
latssr = lats_s * pi / 180.0
lonsr = lons * pi / 180.0

nlats = len(lats)
nlons = len(lons)

area = np.zeros([nlats,nlons],np.float)
lonvalue = np.zeros(nlons,np.float)
lonweight = np.zeros(nlons,np.float)

# Almost all grids will have equal longitude spacing, but just in case:
for ilon in range(0,nlons-1):
        lonvalue[ilon] = abs(lonsr[ilon+1] - lonsr[ilon])
        lonweight[ilon] = abs(lonsr[ilon+1] - lonsr[ilon])/2.0 * pi

lonvalue[nlons-1] = abs(lonsr[nlons-1] - lonsr[nlons-2])
lonweight[nlons-1] = abs(lonsr[nlons-1] - lonsr[nlons-2])/2.0*pi

"""
for ilat in range(0,nlats):
	latweight = weights[ilat]
	for ilon in range(0,nlons):
		area[ilat,ilon] = abs(pi * rE * rE * latweight * lonweight[ilon])
""" 

for ilat in range(1,nlats-1):
        print ilat
        # Based on: area above a latitude lat = 2piR^2(1 - sin(lat)
        # Thus area between two latitudes: 2piR^2(sin(lat1) - sin(lat2))
        # Break into 2pi and multiply by difference between lons: R^2(sin(lat1)-sin(lat2)) * (lon1 - lon2)
        latvalue = abs(np.sin(latssr[ilat]) - np.sin(latssr[ilat-1]))

        for ilon in range(0,nlons):
                area[ilat,ilon] = abs(rE * rE * latvalue * lonvalue[ilon])

# ilat == 0 and nlats-1
latvalue = abs(np.sin((latssr[0] - latsr[0])))
for ilon in range(0,nlons):
	area[0,ilon] = abs(rE * rE * latvalue * lonvalue[ilon])
	area[nlats-1,ilon] = abs(rE * rE * latvalue * lonvalue[ilon])


ncfile = Dataset(Dir + "/SurfaceArea.nc", 'w')
ncfile.createDimension('lon', nlons)
ncfile.createDimension('lat', nlats)


SurfA = ncfile.createVariable('SurfaceArea','f4',('lat','lon'),fill_value=-9999)
setattr(SurfA,'Extra Info','Gridded surface area in m2')
SurfA[...] = area[...]

print np.sum(SurfA)

ncfile.close()

Ngl.end()


