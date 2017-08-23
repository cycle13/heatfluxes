# -*- coding: utf-8 -*-
"""
Created on Mon Jul 11 2016

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

test = 0
startyr = [5,6]
endyr = [105,95]
spinup = 60
Expertype = "CAM4POP_f45g37"
ExperimentCTL = ""
#Experiments = np.array(["_10All_75S-0S_0"])
#Experiments = np.array(["_10All_75S-45S_0"])
Experiments = np.array(["_10At25_55N_1","_10At25_55N_2","_10At25_55N_3","_10At25_55N_4"])
Long = np.array([0,0,0,0])
RunAvgs = 10

if test == 1:
	nexps = 1
else:
	nexps = len(Experiments)


def XrayOpen(filenamein,decodetimes=True):

        try:
                if decodetimes:
                        filein= xray.open_dataset(filenamein)
                else:
                        filein=xray.open_dataset(filenamein,decode_times=False)
        except RuntimeError:
                print filenamein
                exit("couldn't find file")
        return filein



FigDir = '/home/disk/eos4/rachel/Figures/HeatFluxes/'

Dir = "/home/disk/eos4/rachel/CESM_outfiles/HYAK/FluxEXP/"
filenameCTL = Dir + "/" + Expertype + ExperimentCTL + "/atm/ClimAvg_" + Expertype + ExperimentCTL + ".cam2.h0." + '{:04d}'.format(startyr[0]) + "-" + '{:04d}'.format(endyr[0]) + ".nc"
print filenameCTL
filein= XrayOpen(filenameCTL)
lons = filein['lon'].values
lats = filein['lat'].values
nlats = len(lats)

filenameOCTL = Dir + "/" + Expertype + ExperimentCTL + "/ocn/ClimAvg_" + Expertype + ExperimentCTL + ".pop.h." + '{:04d}'.format(startyr[0]) + "-" + '{:04d}'.format(endyr[0]) + ".nc"
print filenameOCTL

fileinO= XrayOpen(filenameOCTL,False)
latsO = fileinO['lat_aux_grid'].values

nlatsO = len(latsO)

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
    PrecipE = np.mean(ncdata.variables['PRECC'][:]+
                    ncdata.variables['PRECL'][:], axis=2)*const.rho_w  # kg/m2/s or mm/s
    EminusP = Evap - PrecipE  # kg/m2/s or mm/s
    Precip = np.mean(ncdata.variables['PRECT'][:],axis=2) * 1000 * 60*60*24  # from m/s to mm/day

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
    HTann['precip'] = np.squeeze(Precip)

    return HTann


# Open Atmos surface area
SurfAreaFile = Dataset(Dir + "/" + Expertype + "_SurfaceArea.nc")
SurfArea = SurfAreaFile["SurfaceArea"]

# Calculate values for Control:
HT_control = CESM_heat_transport(filein)
HT_ocean_ctl = fileinO['N_HEAT'].values[0][0][0]

# Get equatorial values
CTL_CrE_atm = griddata(lats,HT_control['atm'],[0])
print CTL_CrE_atm
CTL_CrE_ocn = griddata(latsO,HT_ocean_ctl,[0])
print CTL_CrE_ocn

# Get precip asymmetry values
precipNH = 0
precipSH = 0
for ilat in range(0,nlats):
	if lats[ilat] < 30 and lats[ilat] > 0:
		precipNH +=HT_control['precip'][ilat]
	elif lats[ilat] > -30 and lats[ilat] < 0:
                precipSH +=HT_control['precip'][ilat]

CTL_Pcp_NHSH = precipNH - precipSH

# Loop through all ensemble members
if np.any(Long == 1):
        nyears = endyr[2]-startyr[2]
else:
        nyears = endyr[1]-startyr[1]

EXP_CrE_atm = np.empty([nexps,nyears]) * np.nan
EXP_CrE_ocn = np.empty([nexps,nyears]) * np.nan
EXP_Pcp_NHSH = np.empty([nexps,nyears]) * np.nan
years = np.zeros(nyears)

for iexp in range(0,nexps):
        if Long[iexp] == 1:
                start = startyr[2]
                end =endyr[2]
        else:
                start = startyr[1]
                end = endyr[1]
	yearcount = 0
	for iyear in range(start,end):
		if RunAvgs > 0:
			filenameE = Dir + "/" + Expertype + Experiments[iexp] + "/atm/RunAvg" + str(RunAvgs) + "_" + Expertype + Experiments[iexp] + ".cam2.h0." + '{:04d}'.format(iyear) + ".nc"
			filenameEO = Dir + "/" + Expertype + Experiments[iexp] + "/ocn/RunAvg" + str(RunAvgs) + "_" + Expertype + Experiments[iexp] + ".pop.h." + '{:04d}'.format(iyear) + ".nc"
		else:
                        filenameE = Dir + "/" + Expertype + Experiments[iexp] + "/atm/AnnAvg_" + Expertype + Experiments[iexp] + ".cam2.h0." + '{:04d}'.format(iyear) + ".nc"
                        filenameEO = Dir + "/" + Expertype + Experiments[iexp] + "/ocn/" + Expertype + Experiments[iexp] + ".pop.h." + '{:04d}'.format(iyear) + ".nc"

                fileinE =XrayOpen(filenameE)
                fileinEO =XrayOpen(filenameEO,False)

		HT_exper = CESM_heat_transport(fileinE)
		HT_ocean_exper = fileinEO['N_HEAT'].values[0][0][0]

		EXP_CrE_atm[iexp,yearcount] = griddata(lats,HT_exper['atm'],[0]) - CTL_CrE_atm
                EXP_CrE_ocn[iexp,yearcount] = griddata(latsO,HT_ocean_exper,[0]) - CTL_CrE_ocn
		years[yearcount] = iyear	

		precipNH = 0
		precipSH = 0
		for ilat in range(0,nlats):
			if lats[ilat] < 30 and lats[ilat] > 0:
				precipNH +=HT_exper['precip'][ilat]
			elif lats[ilat] > -30 and lats[ilat] < 0:
				precipSH +=HT_exper['precip'][ilat]

		EXP_Pcp_NHSH[iexp,yearcount] = precipNH - precipSH - CTL_Pcp_NHSH
		
		yearcount += 1
		
ticks = [-90, -60, -30, 0, 30, 60, 90]
#titles = [ExperimentCTL,Experiment,"CTL-" + Experiment] 


# Plot

fig = plt.figure(figsize=(10,4))

ax = fig.add_subplot(1,1,1)
ax2 = ax.twinx()


for iexp in range(0,nexps):
	if iexp == 0:
		ax.plot(years, np.nanmean(EXP_CrE_atm[:,:],axis=0), 'r-', label='atm', linewidth=2)
		ax.plot(years, np.nanmean(EXP_CrE_ocn[:,:],axis=0), 'b-', label='ocn', linewidth=2)
		ax2.plot(years, np.nanmean(EXP_Pcp_NHSH[:,:],axis=0),'g-',label='0-30 NH-SH precip',linewidth=2)
        
	ax.plot(years, EXP_CrE_atm[iexp,:], 'r-', label='_nolegend_', linewidth=0.5)
        ax.plot(years, EXP_CrE_ocn[iexp,:], 'b-', label='_nolegend_', linewidth=0.5)
        ax2.plot(years, EXP_Pcp_NHSH[iexp,:],'g-',label='_nolegend',linewidth=0.5)


ax.set_ylabel(r'$\Delta$ c-eq heat transport, PW')
ax2.set_ylabel(r'$\Delta$ 0-30 NH-SH precip, mm/day')
ax.set_xlabel('Years')
if Experiments[0] == "_10At25_55N_1":
	plt.title("North Atlantic 25-55N heating")
else:
	plt.title(Expertype + Experiments[0] + " - " + Expertype + " CTL, ensemble")
#ax.axhline(0, linestyle='-', color='k')
#maxy1 = int(max(np.amax(abs(EXP_CrE_atm)),np.amax(abs(EXP_CrE_ocn)))*1.25)
#maxy2 = int(np.amax(abs(EXP_Pcp_NHSH))
#ax.set_ylim(-maxy1,maxy1)

ax.legend(loc="upper left")
ax2.legend(loc="lower right")
ax.grid()
figtitle = Expertype + "CTL-" + Experiments[0] + "Ens_CrossEqHeatFlux_Diag"
fig.tight_layout(pad=2)        # Put more white space around plots
plt.savefig(FigDir + figtitle + ".eps",format='eps',dpi=1200., facecolor='w', edgecolor='w')

exit()
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
    plt.title(Expertype + titles[n] + startyr[n] + 'to' + endyr[n] + 'yr avg')
    ax.set_xlim(-90,90)
    ax.set_xticks(ticks)
    ax.legend(loc='upper left')
    ax.grid()
#Plot differences
    figtitle = Expertype + titles[n] + "NorthwardHeatFlux_Diag"
    fig.tight_layout(pad=2)        # Put more white space around plots
    plt.savefig(FigDir + figtitle + ".eps",format='eps',dpi=1200., facecolor='w', edgecolor='w')


fig = plt.figure(figsize=(10,4))

ax = fig.add_subplot(1, 1 , 1)
ax.plot(lats, griddata(latsO,HT_ocean_exper-HT_ocean_ctl,lats) + HT_exper['atm']-HT_control['atm'], 'k-', label='atm+ocn', linewidth=2)
ax.plot(lats, HT_exper['atm']-HT_control['atm'], 'r-', label='atm', linewidth=2)
ax.plot(lats, HT_exper['dse']-HT_control['dse'], 'r--', label='dry')
ax.plot(lats, HT_exper['latent']-HT_control['latent'], 'r:', label='latent')
ax.plot(lats, (HT_exper['precip']-HT_control['precip']) * np.amax(HT_exper['atm']-HT_control['atm'])/np.amax((HT_exper['precip']-HT_control['precip'])), 'g:', label='normalized precip')
ax.plot(lats, griddata(latsO,HT_ocean_exper-HT_ocean_ctl,lats), 'b-', label='ocean', linewidth=2)


plt.ylabel('Northward heat transport, PW')
plt.title(Expertype + titles[2] + 'clim avg')
ax.set_xlim(-90,90)
ax.set_xticks(ticks)
ax.legend(loc='upper left')
ax.grid()


figtitle = Expertype + Experiment + "-CTL_NorthwardHeatFlux" 
fig.tight_layout(pad=2)        # Put more white space around plots
plt.savefig(FigDir + figtitle + ".eps",format='eps',dpi=1200., facecolor='w', edgecolor='w')


exit()
FSDT = filein["FSDTOA"]
NetTOPWm2 = filein["FSNT"] - filein["FLNT"]
NetSURFWm2 = -filein["FSNS"] + filein["FLNS"] + filein["LHFLX"] + filein["SHFLX"]

NetATMOSWm2 = NetTOPWm2 + NetSURFWm2

SurfAreaFile = Dataset(Dir + "/" + Expertype + "_SurfaceArea.nc")
SurfArea = SurfAreaFile["SurfaceArea"]
print np.nansum(SurfArea)
print SurfArea[:,0]
# Broadcast lat array
latsarray = np.resize(lats,SurfArea.shape)
# Create SurfArea mask for NH and SH
LatsNH = np.where(latsarray <= 0,0,latsarray)
print np.nansum(LatsNH)
SurfAreaNH = np.where(latsarray < 0,0,SurfArea)
#SurfAreaNH = np.where(latsarray == 0,SurfArea[:,:]/2.0,SurfAreaNH)

LatsSH = np.where(latsarray >= 0,0,latsarray)
print np.nansum(LatsSH)
SurfAreaSH = np.where(latsarray > 0,0,SurfArea)
#SurfAreaSH = np.where(latsarray == 0,SurfArea[:,:]/2.0,SurfAreaSH)


print "area-weighted downwelling at TOA:", np.nansum(FSDT * SurfArea) / np.nansum(SurfArea)
print "area-weighted net TOM: ", np.nansum(NetTOPWm2 * SurfArea)/np.nansum(SurfArea)
print "area-weighted net out of surface: ", np.nansum(NetSURFWm2 * SurfArea)/np.nansum(SurfArea)


print "NH area-weighted net TOM: ", np.nansum(NetTOPWm2 * SurfAreaNH)/np.nansum(SurfAreaNH)
print "SH area-weighted net TOM: ", np.nansum(NetTOPWm2 * SurfAreaSH)/np.nansum(SurfAreaSH)

print "NH area-weighted net out of surface: ", np.nansum(NetSURFWm2 * SurfAreaNH)/np.nansum(SurfAreaNH)
print "SH area-weighted net out of surface: ", np.nansum(NetSURFWm2 * SurfAreaSH)/np.nansum(SurfAreaSH)


print "area-weighted INTO atmosphere:"
print np.nansum(NetATMOSWm2 * SurfArea)/np.nansum(SurfArea)

NetATMOS_NH = np.nansum(NetATMOSWm2 * SurfAreaNH)/np.nansum(SurfAreaNH)

NetATMOS_SH = np.nansum(NetATMOSWm2 * SurfAreaSH)/np.nansum(SurfAreaSH)

print np.nansum(SurfAreaSH)
print np.nansum(SurfAreaNH)

print "Net Atmosphere in NH: ", NetATMOS_NH
print "Net Atmosphere in SH: ", NetATMOS_SH

#print np.nansum(NetATMOSW)

#Now calculate net for NH and SH separately



filein.close()
SurfAreaFile.close()




