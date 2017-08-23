# -*- coding: utf-8 -*-
"""
Created on Mon Jul 11 2016

The basis for many of these functions have been taken from Brian Rose's climlab github page:
https://github.com/brian-rose/climlab

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

from rhwhitepackages.readwrite import XrayOpen
from rhwhitepackages.plots import initcontourplot,setplotrange,initCScontourplot
import argparse
import copy

parser = argparse.ArgumentParser(description="generic ocean plots")
parser.add_argument('--expN',metavar='exp',type=str,nargs='?',default="CAM4POP_f19g16C_noTopo",help='experiment WITHOUT orography')
parser.add_argument('--res',metavar='res',type=str,nargs='?',default="f19g16",help='resolution')


args = parser.parse_args()
print "here's what I have as arguments: "
print args

preciplat = 10
ExpN = args.expN
res = args.res

if ExpN == "CAM4POP_f19g16C_noTopo":
        startyr = 211
        endyr = 240
elif ExpN == "CAM4POP_f19g16C_noR":
        startyr = 411
        endyr = 440

ExperimentM = "CAM4POP_B1850_NCAR"

FigDir = '/home/disk/eos4/rachel/Figures/AMOC/'

Dir = "/home/disk/rachel/CESM_outfiles/"
filenameM = Dir + "/" + ExperimentM + "/atm/hist/ClimAnn_b40.1850.track1.2deg.003.cam.h0.500-529.nc"
filenameN = Dir + "/" + ExpN + "/atm/hist/ClimAnn_" + ExpN + ".cam.h0." + str(startyr) + "-" + str(endyr) + ".nc"

fileinM= XrayOpen(filenameM,False)
fileinN =XrayOpen(filenameN,False)
lons = fileinM['lon'].values
lats = fileinM['lat'].values
lats2 = fileinN['lat'].values

gwA = fileinM['gw']
nlats = len(lats)

filenameOM = Dir + "/" + ExperimentM + "/ocn/hist/ClimAnn_b40.1850.track1.2deg.003.pop.h.500-529.nc"
filenameON = Dir + "/" + ExpN + "/ocn/hist/ClimAnn_" + ExpN + ".pop.h." + str(startyr) + "-" + str(endyr) + ".nc"

fileinOM = XrayOpen(filenameOM,False)
fileinON =XrayOpen(filenameON,False)
latsO = fileinOM['lat_aux_grid'].values
gwO = fileinOM['TAREA']
nlatsO = len(latsO)
print nlatsO

if nlatsO != nlats:
	print "ocean lats not the same as atmos lats"

def inferred_heat_transport( energy_in, lat_deg ):
    '''Returns the inferred heat transport (in PW) by integrating the net energy imbalance from pole to pole.'''
    from scipy import integrate
    lat_rad = np.deg2rad( lat_deg )
    return ( 1E-15 * 2 * np.math.pi * const.a**2 * 
            integrate.cumtrapz( np.cos(lat_rad)*energy_in,
            x=lat_rad, initial=0. ) )

def CESM_heat_transport(ncdata):
    lat = ncdata.variables['lat'][:]
    # TOA radiation
    OLR = np.mean(ncdata.variables['FLNT'][:], axis=2)
    ASR = np.mean(ncdata.variables['FSNT'][:], axis=2)
    Rtoa = ASR - OLR  # net downwelling radiation
    #  surface fluxes  (all positive UP)
    LHF = np.mean(ncdata.variables['LHFLX'][:], axis=2)  # latent heat flux (evaporation)
    SHF = np.mean(ncdata.variables['SHFLX'][:], axis=2) # sensible heat flux
    LWsfc = np.mean(ncdata.variables['FLNS'][:], axis=2)  # net longwave radiation at surface
    SWsfc = -np.mean(ncdata.variables['FSNS'][:], axis=2)  # net shortwave radiation at surface
    #  energy flux due to snowfall
    SnowFlux =  np.mean(ncdata.variables['PRECSC'][:]+
                        ncdata.variables['PRECSL'][:], axis=2)*const.rho_w*const.Lhfus
    #  hydrological cycle
    Evap = np.mean(ncdata.variables['QFLX'][:], axis=2)  # kg/m2/s or mm/s
    Precip = np.mean(ncdata.variables['PRECC'][:]+
                    ncdata.variables['PRECL'][:], axis=2)*const.rho_w  # kg/m2/s or mm/s
    EminusP = Evap - Precip  # kg/m2/s or mm/s
    PrecipWatts = np.mean(ncdata.variables['PRECC'][:]+
                    ncdata.variables['PRECL'][:], axis=2)*const.rho_w*const.Lhvap
    SurfaceRadiation = LWsfc + SWsfc  # net upward radiation from surface
    SurfaceHeatFlux = SurfaceRadiation + LHF + SHF + SnowFlux  # net upward surface heat flux
    Fatmin = Rtoa + SurfaceHeatFlux  # net heat flux in to atmosphere
    # heat transport terms
    HTann = {}
    HTann['total'] = np.squeeze(inferred_heat_transport(Rtoa, lats))
    HTann['atm'] = np.squeeze(inferred_heat_transport(Fatmin, lats))
    HTann['ocean'] = np.squeeze(inferred_heat_transport(-SurfaceHeatFlux, lats))
    HTann['latent'] = np.squeeze(inferred_heat_transport(EminusP*const.Lhvap, lats)) # atm. latent heat transport from moisture imbal.
    HTann['dse'] = HTann['atm'] - HTann['latent']  # dry static energy transport as residual
    HTann['precip'] = np.squeeze(PrecipWatts)

    return HTann

HT_control = CESM_heat_transport(fileinM)
HT_ocean_ctl = fileinOM['N_HEAT'].values[0][0][0]
HT_exper = CESM_heat_transport(fileinN)
HT_ocean_exper = fileinON['N_HEAT'].values[0][0][0]
runs = [HT_control,HT_exper]
N = len(runs) + 1

# Get equatorial values
CTL_CrE_atm = griddata(lats,HT_control['atm'],[0])
EXP_CrE_atm = griddata(lats,HT_exper['atm'],[0])

CTL_CrE_ocn = griddata(latsO,HT_ocean_ctl,[0])
EXP_CrE_ocn = griddata(latsO,HT_ocean_exper,[0])

# Get precip asymmetry values
precipNH = 0
precipSH = 0
precipNHE = 0
precipSHE = 0
for ilat in range(0,nlats):
        if lats[ilat] < preciplat and lats[ilat] > 0:
                precipNH +=HT_control['precip'][ilat]
		precipNHE += HT_exper['precip'][ilat]
        elif lats[ilat] > -preciplat and lats[ilat] < 0:
                precipSH +=HT_control['precip'][ilat]
		precipSHE += HT_exper['precip'][ilat]

CTL_Pcp_NHSH = precipNH.values - precipSH.values
EXP_Pcp_NHSH = precipNHE.values - precipSHE.values

ticks = [-90, -60, -30, 0, 30, 60, 90]
titles = [ExperimentM,ExpN,"CTL-" + ExpN] 

#Plot
for n, HT in enumerate(runs):
    fig = plt.figure(figsize=(10,4))

    ax = fig.add_subplot(1, 1, 1)
    ax.plot(lats, griddata(latsO,HT_ocean_ctl,lats) + HT['atm'], 'k-', label='atm+ocean', linewidth=2)
    ax.plot(lats, HT['atm'], 'r-', label='atm', linewidth=2)
    ax.plot(lats, HT['dse'], 'r--', label='dry')
    ax.plot(lats, HT['latent'], 'r:', label='latent')
    ax.plot(lats, HT['precip']*np.amax(HT['atm'])/np.amax(HT['precip']), 'g:', label='normalized precip')
    ax.plot(lats, HT['ocean'], 'b-', label='ocean', linewidth=2)
    ax.plot(lats, griddata(latsO,HT_ocean_ctl,lats), 'b--', label='oceanDIAG', linewidth=2)

    plt.ylabel('Northward heat transport, PW')
    plt.title(titles[n] + str(startyr) + 'to' + str(endyr) + 'yr avg')
    ax.set_xlim(-90,90)
    ax.set_xticks(ticks)
    ax.legend(loc='upper left')
    ax.grid()
#Plot differences
    figtitle = titles[n] + "NorthwardHeatFlux_Diag"
    fig.tight_layout(pad=2)        # Put more white space around plots
    plt.savefig(FigDir + figtitle + ".eps",format='eps',dpi=1200., facecolor='w', edgecolor='w')


fig = plt.figure(figsize=(10,4))

ax = fig.add_subplot(1, 1 , 1)
ax.plot(lats, griddata(latsO,HT_ocean_ctl-HT_ocean_exper,lats) + HT_control['atm']-HT_exper['atm'], 'k-', label='atm+ocn', linewidth=2)
ax.plot(lats, HT_control['atm']-HT_exper['atm'], 'r-', label='atm', linewidth=2)
ax.plot(lats, HT_control['dse']-HT_exper['dse'], 'r--', label='dry')
ax.plot(lats, HT_control['latent']-HT_exper['latent'], 'r:', label='latent')
ax.plot(lats, (HT_control['precip']-HT_exper['precip']) * np.amax(HT_control['atm']-HT_exper['atm'])/np.amax((HT_control['precip']-HT_exper['precip'])), 'g:', label='normalized precip')
ax.plot(lats, griddata(latsO,HT_ocean_ctl-HT_ocean_exper,lats), 'b-', label='ocean', linewidth=2)


plt.ylabel('Northward heat transport, PW')
plt.title(titles[2] + 'clim avg')
ax.set_xlim(-90,90)
ax.set_xticks(ticks)
ax.legend(loc='upper left')
ax.grid()


figtitle = "CTL-" + ExpN + "_NorthwardHeatFlux" 
fig.tight_layout(pad=2)        # Put more white space around plots
plt.savefig(FigDir + figtitle + ".eps",format='eps',dpi=1200., facecolor='w', edgecolor='w')

NetTOPWm2 = np.squeeze(fileinM["FSNT"] - fileinM["FLNT"]).sum(dim='lon').values
NetSURFWm2 = np.squeeze(-fileinM["FSNS"] + fileinM["FLNS"] + fileinM["LHFLX"] + fileinM["SHFLX"]).sum(dim='lon').values

NetATMOSWm2 = NetTOPWm2 + NetSURFWm2

NetTOPWm2N = np.squeeze(fileinN["FSNT"] - fileinN["FLNT"]).sum(dim='lon').values
NetSURFWm2N = np.squeeze(-fileinN["FSNS"] + fileinN["FLNS"] + fileinN["LHFLX"] + fileinN["SHFLX"]).sum(dim='lon').values
NetATMOSWm2N = NetTOPWm2N + NetSURFWm2N

print np.nansum(gwA)
# Broadcast lat array
# Create SurfArea mask for NH and SH
LatsNH = np.where(lats <= 0,0,lats)
SurfAreaNH = np.where(lats < 0,0,gwA)
print np.nansum(SurfAreaNH)
LatsSH = np.where(lats >= 0,0,lats)
SurfAreaSH = np.where(lats > 0,0,gwA)
print np.nansum(SurfAreaSH)

SurfArea = gwA

#print NetTOPWm2N * SurfAreaSH

print "area-weighted net TOM: ", np.nanmean(NetTOPWm2 * SurfArea), ' vs ', np.nanmean(NetTOPWm2N * SurfArea)
print "area-weighted net out of surface: ", np.nanmean(NetSURFWm2 * SurfArea), ' vs ',np.nanmean(NetSURFWm2N * SurfArea)


print "NH area-weighted net TOM: ", np.nanmean(NetTOPWm2 * SurfAreaNH), ' vs ', np.nanmean(NetTOPWm2N * SurfAreaNH), " W/m2"
print "SH area-weighted net TOM: ", np.nanmean(NetTOPWm2 * SurfAreaSH), ' vs ', np.nanmean(NetTOPWm2N * SurfAreaSH), " W/m2"

print "NH area-weighted net out of surface: ", np.nanmean(NetSURFWm2 * SurfAreaNH), ' vs ', np.nanmean(NetSURFWm2N * SurfAreaNH)
print "SH area-weighted net out of surface: ", np.nanmean(NetSURFWm2 * SurfAreaSH), ' vs ', np.nanmean(NetSURFWm2N * SurfAreaSH)


print "area-weighted INTO atmosphere:"
print np.nanmean(NetATMOSWm2 * SurfArea)

NetATMOS_NH = np.nanmean(NetATMOSWm2 * SurfAreaNH)
NetATMOS_SH = np.nanmean(NetATMOSWm2 * SurfAreaSH)

NetATMOS_NHN = np.nanmean(NetATMOSWm2N * SurfAreaNH)
NetATMOS_SHN = np.nanmean(NetATMOSWm2N * SurfAreaSH)


print "Net Atmosphere in NH: ", NetATMOS_NH, ' vs ', NetATMOS_NHN
print "Net Atmosphere in SH: ", NetATMOS_SH, ' vs ', NetATMOS_SHN


print "STATISTICS ON ANNUAL AVERAGES"
print " mean atmosphere: ", CTL_CrE_atm, EXP_CrE_atm, 'diff: ', CTL_CrE_atm-EXP_CrE_atm
print "mean ocean: ", CTL_CrE_ocn, EXP_CrE_ocn, 'diff: ', CTL_CrE_ocn - EXP_CrE_ocn
print "mean total: ", CTL_CrE_atm + CTL_CrE_ocn, EXP_CrE_atm + EXP_CrE_ocn, 'diff: ', CTL_CrE_atm + CTL_CrE_ocn - (EXP_CrE_atm + EXP_CrE_ocn)
print "mean precip: ", CTL_Pcp_NHSH, EXP_Pcp_NHSH, 'diff: ', CTL_Pcp_NHSH - EXP_Pcp_NHSH






