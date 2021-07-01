# -*- coding: utf-8 -*-
"""
Created on Thu Jun 16 13:24:45 2016
This version frozen on Thu Jul 1 5;43 2021

@author: Ent00002
@author: pkeys
"""

#%% Import libraries
import os
os.chdir('/Users/patrickkeys/rams_gdrive/work/research/models/WAM2layersPython-master/WAM2layersPython3/')

import sys

import numpy as np
from netCDF4 import Dataset
import scipy.io as sio
import calendar
from getconstants import getconstants
#from extrafunctions import extrafunctions
from timeit import default_timer as timer
import datetime
import copy 

import matplotlib.pyplot as plt

import cartopy as ct
from mpl_toolkits.axes_grid1 import make_axes_locatable

import matplotlib as mpl

import warnings
warnings.filterwarnings("ignore", category=RuntimeWarning)

#%%
def plotMap(axes, data, lats, lons, cmap='coolwarm', vmin=None, vmax=None):
    image = plt.pcolormesh(lons,lats,data,transform=ct.crs.PlateCarree(),vmin=vmin,vmax=vmax,cmap=cmap, shading='flat')
    axes.coastlines(color='black', linewidth=1.2) 
    divider = make_axes_locatable(axes)
    ax_cb = divider.new_horizontal(size="2%", pad=0.1, axes_class=plt.Axes)  
    plt.gcf().add_axes(ax_cb)    
    cb = plt.colorbar(image,cax=ax_cb)    
    plt.sca(axes)    
    return cb,image
## example of how to plot stuff
#plt.figure(figsize=(10,6))
#ax = plt.subplot(1,1,1,projection=ct.crs.PlateCarree())
#cb,image = plotMap(ax,sp[0,:,:]/100.,latitude,longitude)
#cb.ax.set_ylabel('hPa')
#plt.savefig('figures/testFigure', dpi=400, bbox_inches='tight')
#plt.show()

#%% var names

### deepcore
input_folder2d = '/Volumes/jammer/merra_2/data_files/merra2_2d_files/'
input_folder2d_biasP = '/Volumes/onenight/merra_2/data_files/merra2_2d_files/'
input_folderVertStates = '/Volumes/catfish/merra_2/data_files/merra2_2d_files/'
input_folderVertFluxes = '/Volumes/onenight/merra_2/data_files/merra2_2d_files/'
input_folder2d_LandEP = '/Volumes/catfish/merra_2/data_files/merra2_2d_files/'
# input_folder3d = '/Volumes/catfish/merra_2/data_files/merra2_3d_files/'
input_folder3d = '/Volumes/catfish/merra_2/data_files/merra2_3d_files/'
invariant_data  = '/Users/patrickkeys/rams_gdrive/work/research/models/WAM2layersPython-master/WAM2layersPython3/data/constants/MERRA2_101.const_2d_asm_Nx.00000000.nc4'
# output_fluxes_states_folder = '/Volumes/coffey/merra_2/fluxes_and_states_PRECTOTCORR_mod/'
output_fluxes_states_folder = '/Volumes/coffey/merra_2/fluxes_and_states_2020_75/'


#%% variables names

spVarName = 'PS'
humVarName = 'QV'
uVarName = 'U'
vVarName = 'V'
evapVarName = 'EVAP'
evapLandVarName = 'EVLAND'
# precipVarName = 'PRECTOT'
precipVarName = 'PRECTOTCORR'
tclwVarName = 'TQL'
tciwVarName = 'TQI'
tcwvVarName = 'TQV'
elwVarName = 'UFLXQL' #eastward_flux_of_atmospheric_liquid_water
ewvVarName = 'UFLXQV' #eastward_flux_of_atmospheric_water_vapor
eiVarName = 'UFLXQI' #eastward_flux_of_atmospheric_ice
nlwVarName = 'VFLXQL' #northward_flux_of_atmospheric_liquid_water
nwvVarName = 'VFLXQV' #northward_flux_of_atmospheric_water_vapor
niVarName = 'VFLXQI' #northward_flux_of_atmospheric_ice


#%% pressure levels

# MERRA2 documentation: http://wiki.seas.harvard.edu/geos-chem/index.php/GEOS-Chem_vertical_grids#72-layer_vertical_grid 
numPresLevels= 72
# Reported here as TOP OF ATMOS to the SURFACE
grabLevs = np.array([26,34,38,44,49,53,56,59,62,65,68,70,71])
edgeLevs = np.array([0,31,37,42,47,52,55,58,61,64,67,70,71,72])
boundary = 8
#%%

timeincPS = np.arange(0,24,6) #assuming 24 timesteps per day z(one every hour)
timeincU = np.arange(0,4,1) #assuming 4 timesteps per day (one every 6 hours)
timeincV = np.arange(0,4,1) #assuming 4 timesteps per day (one every 6 hours)
timeincQ = np.arange(0,4,1) #assuming 4 timesteps per day (one every 6 hours)
timeincEP = np.arange(0,24,1) #assuming 24 timesteps per day (one every hour)
timeincTCW = np.arange(0,24,6) #assuming 24 timesteps per day (one every hour)
timeincCorrFlux = np.arange(0,24,6) #assuming 24 timesteps per day (one every hour)

timeindexTomorrow = 0

#%%BEGIN OF INPUT (FILL THIS IN)
years = np.arange(1980,1981) #fill in the years
start_iDOY = 0#start day of year index for initial year in "years" vector
end_iDOY = 366#end day of year ; PWK NOTE = 366 is a full run for 0 indexed I think? 

divt = 48 # division of the timestep, 24 means a calculation timestep of 6/24 = 0.25 hours (numerical stability purposes)
count_time = 4 # number of indices to get data from (for six hourly data this means everytime one day), note that this code assumes the count_time for precipitation to be double this number

# Manage the extent of your dataset (FILL THIS IN)
# Define the latitude and longitude cell numbers to consider and corresponding lakes that should be considered part of the land
latnrs = np.arange(20,341)#(0,361)
lonnrs = np.arange(0,576)
isglobal = 1 # fill in 1 for global computations (i.e. Earth round), fill in 0 for a local domain with boundaries

lake_mask = np.zeros((1,2))+10

#END OF INPUT

#%% function to determine what the NASA hundred number is
def get_NASA_wildNumStr(year):
    if(year>=1980 and year<1992):
        NASA_wildNumStr = '100'
    elif(year>=1990 and year<2001):
        NASA_wildNumStr = '200'
    elif(year>=2000 and year<2011):
        NASA_wildNumStr = '300'
    elif(year>=2010 and year<2021):
        NASA_wildNumStr = '400'
    else:
        sys.exit("NASA wildcard number string could not be determined. Check your input year.")
    
    return NASA_wildNumStr     
        
#%% Datapaths (FILL THIS IN)

# other scripts use exactly this sequence, do not change it unless you change it also in the scripts
def data_path(iDOY,yearnumber):   
        
    aDate = datetime.date(yearnumber,1,1) + datetime.timedelta(int(iDOY))
    dayMonth = f"{aDate:%Y%m%d}"[4:]

    aDateTmrw = datetime.date(yearnumber,1,1) + datetime.timedelta(int(iDOY)+1)
    dayMonthTmrw = f"{aDateTmrw:%Y%m%d}"[4:]
    
    if(dayMonth[0:2]=='12' and dayMonthTmrw[0:2]=='01'):
        yearnumberTmrw = yearnumber + 1
    else:
        yearnumberTmrw = yearnumber

    NASA_wildNumStr = get_NASA_wildNumStr(yearnumber)
    NASA_wildNumStrTmrw = get_NASA_wildNumStr(yearnumberTmrw)

    varName = 'ps' # Datapath 0, 1
    sp_dataPath = os.path.join(input_folder2d + varName + '/', 'MERRA2_' + NASA_wildNumStr + '.inst1_2d_asm_Nx.' + str(yearnumber) + dayMonth + '.nc4.nc4') #surface pressure
    sp_dataPath_tmrw = os.path.join(input_folder2d + varName + '/', 'MERRA2_' + NASA_wildNumStrTmrw + '.inst1_2d_asm_Nx.' + str(yearnumberTmrw) + dayMonthTmrw + '.nc4.nc4') #surface pressure Tomorrow

    varName = 'uvq' # Datapath 2, 3
    humwind_dataPath = os.path.join(input_folder3d + varName + '/', 'MERRA2_' + NASA_wildNumStr + '.inst6_3d_ana_Nv.' + str(yearnumber) + dayMonth + '.nc4.nc') # u v q
    humwind_dataPath_tmrw = os.path.join(input_folder3d + varName + '/', 'MERRA2_' + NASA_wildNumStrTmrw + '.inst6_3d_ana_Nv.'  + str(yearnumberTmrw) + dayMonthTmrw + '.nc4.nc') # u v q Tomorrow

    varName = 'vert_integ_states' # Datapath 4, 5
    tcw_dataPath = os.path.join(input_folderVertStates + varName + '/', 'MERRA2_' + NASA_wildNumStr + '.inst1_2d_int_Nx.' + str(yearnumber) + dayMonth + '.nc4.nc4') #integrated water column variables
    tcw_dataPath_tmrw = os.path.join(input_folderVertStates + varName + '/', 'MERRA2_' + NASA_wildNumStrTmrw + '.inst1_2d_int_Nx.' + str(yearnumberTmrw) + dayMonthTmrw + '.nc4.nc4') #integrated water column variables Tomorrow

    varName = 'ep' # Datapath 6
    ep_dataPath = os.path.join(input_folder2d + varName + '/', 'MERRA2_' + NASA_wildNumStr + '.tavg1_2d_flx_Nx.' + str(yearnumber) + dayMonth + '.nc4.nc4') #evap and precip
    
    varName = 'vert_integ_fluxes' # Datapath 7, 8
    corrflux_dataPath = os.path.join(input_folderVertFluxes + varName + '/', 'MERRA2_' + NASA_wildNumStr + '.tavg1_2d_int_Nx.' + str(yearnumber) + dayMonth + '.nc4.nc4') #integrated water column variables
    corrflux_dataPath_tmrw = os.path.join(input_folderVertFluxes + varName + '/', 'MERRA2_' + NASA_wildNumStrTmrw + '.tavg1_2d_int_Nx.' + str(yearnumberTmrw) + dayMonthTmrw + '.nc4.nc4') #integrated water column variables Tomorrow
    
    varName = 'biasP' # Datapath 9
    biasP_dataPath = os.path.join(input_folder2d_biasP + varName + '/', 'MERRA2_' + NASA_wildNumStr + '.tavg1_2d_flx_Nx.' + str(yearnumber) + dayMonth + '.nc4.nc4') #evap and precip

    varName = 'land_e_p' # Datapath 10
    landEP_dataPath = os.path.join(input_folder2d_LandEP + varName + '/', 'MERRA2_' + NASA_wildNumStr + '.tavg1_2d_lnd_Nx.' + str(yearnumber) + dayMonth + '.nc4.nc4') #evap and precip

    # Datapath 11
    save_path = os.path.join(output_fluxes_states_folder, str(yearnumber) + '-' + str(iDOY) + 'fluxes_storages.mat')
      
    return sp_dataPath, sp_dataPath_tmrw, humwind_dataPath, humwind_dataPath_tmrw, tcw_dataPath, tcw_dataPath_tmrw, ep_dataPath, corrflux_dataPath, corrflux_dataPath_tmrw, biasP_dataPath, landEP_dataPath, save_path


#%% Code (no need to look at this for running)

def load_sp(a,yearnumber,spNow):
    
    print('loading sp')    

    if(len(spNow)<1):
        print('loading today spNow')
        # surface pressure (state at 00.00, 06.00, 12.00, 18.00 and 00.00 the next day)
        ncid = Dataset(datapath[0], mode = 'r')
        spNow = ncid.variables[spVarName][timeincPS,latnrs,lonnrs]
        #swap MERRA grid to ERA grid
        spNow = spNow[:,ilat_SwapGridtoERA,:]
        spNow = spNow[:,:,ilon_SwapGridtoERA]   
        # done swapping grid
        ncid.close()
    
    ncid = Dataset(datapath[1], mode = 'r')
    sp_tmrw = ncid.variables[spVarName][timeincPS,latnrs,lonnrs]
    #swap MERRA grid to ERA grid
    sp_tmrw = sp_tmrw[:,ilat_SwapGridtoERA,:]
    sp_tmrw = sp_tmrw[:,:,ilon_SwapGridtoERA]   
    # done swapping grid
    sp = np.insert(spNow,[4],(sp_tmrw[timeindexTomorrow,:,:]), axis = 0) #Pa
    ncid.close()
    
    return sp, sp_tmrw 

def load_q(a,yearnumber,qNow):
    
    print('loading q')
    
    if(len(qNow)<1):
        print('loading today qNow')        
        # specific humidity (state at 00.00 06.00 12.00 18.00)
        ncid = Dataset(datapath[2], mode = 'r')
        qNow = ncid.variables[humVarName][timeincQ,grabLevs,latnrs,lonnrs]
        #swap MERRA grid to ERA grid
        qNow = qNow[:,:,ilat_SwapGridtoERA,:]
        qNow = qNow[:,:,:,ilon_SwapGridtoERA]   
        # done swapping grid
        ncid.close()

    ncid = Dataset(datapath[3], mode = 'r')
    q_tmrw = ncid.variables[humVarName][timeincQ,grabLevs,latnrs,lonnrs]
    #swap MERRA grid to ERA grid
    q_tmrw = q_tmrw[:,:,ilat_SwapGridtoERA,:]
    q_tmrw = q_tmrw[:,:,:,ilon_SwapGridtoERA]   
    # done swapping grid    
    q = np.insert(qNow,[4],(q_tmrw[timeindexTomorrow,:,:,:]), axis = 0) #kg/kg
    ncid.close()
    
    return q, q_tmrw  

def load_uv(a,yearnumber, uNow, vNow):
    # u stands for wind in zonal direction = west-east
    # v stands for wind in meridional direction = south-north 

    print('loading uv')

    if(len(uNow)<1 or len(vNow)<1):
        print('loading today uNow, vNow')        
        ncid = Dataset(datapath[2], mode = 'r')
        uNow = ncid.variables[uVarName][timeincU,grabLevs,latnrs,lonnrs]
        vNow = ncid.variables[vVarName][timeincV,grabLevs,latnrs,lonnrs]
        #swap MERRA grid to ERA grid
        uNow = uNow[:,:,ilat_SwapGridtoERA,:]
        uNow = uNow[:,:,:,ilon_SwapGridtoERA] 
        vNow = vNow[:,:,ilat_SwapGridtoERA,:]
        vNow = vNow[:,:,:,ilon_SwapGridtoERA]        
        # done swapping grid           
        ncid.close()

    # get tomorrow's data
    
    # read the u-wind data    
    ncid = Dataset(datapath[3], mode = 'r')
    U_tmrw = ncid.variables[uVarName][timeincU,grabLevs,latnrs,lonnrs]
    #swap MERRA grid to ERA grid
    U_tmrw = U_tmrw[:,:,ilat_SwapGridtoERA,:]
    U_tmrw = U_tmrw[:,:,:,ilon_SwapGridtoERA]        
    # done swapping grid      
    U = np.insert(uNow,[4],(U_tmrw[timeindexTomorrow,:,:,:]), axis = 0)
    
    
    # read the v-wind data
    V_tmrw = ncid.variables[vVarName][timeincV,grabLevs,latnrs,lonnrs]
    #swap MERRA grid to ERA grid
    V_tmrw = V_tmrw[:,:,ilat_SwapGridtoERA,:]
    V_tmrw = V_tmrw[:,:,:,ilon_SwapGridtoERA]        
    # done swapping grid    
    V = np.insert(vNow,[4],(V_tmrw[timeindexTomorrow,:,:,:]), axis = 0)
    ncid.close()

    return U,V,U_tmrw,V_tmrw

def load_tcw(a,yearnumber,tcwNow):
    #tclwVarName = 'TQL'
    #tciwVarName = 'TQI'
    #tcwvVarName = 'TQV'
    
    print('loading tcw')    

    if(len(tcwNow)<1):
        print('loading today tcwNow')
        ncid = Dataset(datapath[4], mode = 'r')
        tcwNow = ncid.variables[tclwVarName][timeincTCW,latnrs,lonnrs] \
            + ncid.variables[tciwVarName][timeincTCW,latnrs,lonnrs] \
            + ncid.variables[tcwvVarName][timeincTCW,latnrs,lonnrs]
        #swap MERRA grid to ERA grid
        tcwNow = tcwNow[:,ilat_SwapGridtoERA,:]
        tcwNow = tcwNow[:,:,ilon_SwapGridtoERA] 
        # done swapping grid             
        ncid.close()
    
    ncid = Dataset(datapath[5], mode = 'r')
    tcw_tmrw = ncid.variables[tclwVarName][timeincTCW,latnrs,lonnrs] \
            + ncid.variables[tciwVarName][timeincTCW,latnrs,lonnrs] \
            + ncid.variables[tcwvVarName][timeincTCW,latnrs,lonnrs]
    #swap MERRA grid to ERA grid
    tcw_tmrw = tcw_tmrw[:,ilat_SwapGridtoERA,:]
    tcw_tmrw = tcw_tmrw[:,:,ilon_SwapGridtoERA] 
    # done swapping grid            
    tcw = np.insert(tcwNow,[4],(tcw_tmrw[timeindexTomorrow,:,:]), axis = 0) #Pa
    ncid.close()
    
    # print('shape of tcw')
    # print(np.shape(tcw))
    # print('shape of tcw_tmrw')
    # print(np.shape(tcw_tmrw))    
    # exit()
    
    return tcw, tcw_tmrw 

def load_corrflux(a,yearnumber,ecorrfluxNow,ncorrfluxNow):
    #elwVarName = 'UFLXQL' #eastward_flux_of_atmospheric_liquid_water
    #ewvVarName = 'UFLXQV' #eastward_flux_of_atmospheric_water_vapor
    #eiVarName = 'UFLXQI' #eastward_flux_of_atmospheric_ice
    #nlwVarName = 'VFLXQL' #northward_flux_of_atmospheric_liquid_water
    #nwvVarName = 'VFLXQV' #northward_flux_of_atmospheric_water_vapor
    #niVarName = 'VFLXQI' #northward_flux_of_atmospheric_ice

    print('loading corrective fluxes')    

    if(len(ecorrfluxNow)<1 or len(ncorrfluxNow)<1):
        print('loading today corrfluxNow')
        ncid = Dataset(datapath[7], mode = 'r')
        ecorrfluxNow = \
            + ncid.variables[elwVarName][timeincCorrFlux,latnrs,lonnrs] \
            + ncid.variables[ewvVarName][timeincCorrFlux,latnrs,lonnrs] \
            + ncid.variables[eiVarName][timeincCorrFlux,latnrs,lonnrs]  
        ncorrfluxNow = \
            + ncid.variables[nlwVarName][timeincCorrFlux,latnrs,lonnrs] \
            + ncid.variables[nwvVarName][timeincCorrFlux,latnrs,lonnrs] \
            + ncid.variables[niVarName][timeincCorrFlux,latnrs,lonnrs]
        
        #swap MERRA grid to ERA grid
        ecorrfluxNow = ecorrfluxNow[:,ilat_SwapGridtoERA,:]
        ecorrfluxNow = ecorrfluxNow[:,:,ilon_SwapGridtoERA] 
        ncorrfluxNow = ncorrfluxNow[:,ilat_SwapGridtoERA,:]
        ncorrfluxNow = ncorrfluxNow[:,:,ilon_SwapGridtoERA] 
        # done swapping grid                 
        ncid.close()
    
    ncid = Dataset(datapath[8], mode = 'r')
    ecorrflux_tmrw = \
            + ncid.variables[elwVarName][timeincCorrFlux,latnrs,lonnrs] \
            + ncid.variables[ewvVarName][timeincCorrFlux,latnrs,lonnrs] \
            + ncid.variables[eiVarName][timeincCorrFlux,latnrs,lonnrs]  
    ncorrflux_tmrw = \
            + ncid.variables[nlwVarName][timeincCorrFlux,latnrs,lonnrs] \
            + ncid.variables[nwvVarName][timeincCorrFlux,latnrs,lonnrs] \
            + ncid.variables[niVarName][timeincCorrFlux,latnrs,lonnrs]
    #swap MERRA grid to ERA grid
    ecorrflux_tmrw = ecorrflux_tmrw[:,ilat_SwapGridtoERA,:]
    ecorrflux_tmrw = ecorrflux_tmrw[:,:,ilon_SwapGridtoERA] 
    ncorrflux_tmrw = ncorrflux_tmrw[:,ilat_SwapGridtoERA,:]
    ncorrflux_tmrw = ncorrflux_tmrw[:,:,ilon_SwapGridtoERA] 
    # done swapping grid  

    ecorrflux = np.insert(ecorrfluxNow,[4],(ecorrflux_tmrw[timeindexTomorrow,:,:]), axis = 0) 
    ncorrflux = np.insert(ncorrfluxNow,[4],(ncorrflux_tmrw[timeindexTomorrow,:,:]), axis = 0)     
    ncid.close()
    
#    print('shape of ecorrflux')
#    print(np.shape(ecorrflux))
#    print('shape of ecorrflux_tmrw')
#    print(np.shape(ecorrflux_tmrw))    
    
    #units of kg*m-1*s-1
    return ecorrflux, ncorrflux, ecorrflux_tmrw, ncorrflux_tmrw

    
#%%
# The model level as downloaded from ERA-Interim are hard-defined within this code. So it has to be changed if other model levels are downloaded
def getW(q,sp,tcwRaw,latnrs,lonnrs,count_time,density_water,latitude,longitude,g,A_gridcell,boundary):
    
    # Reported here as TOP OF ATMOS to the SURFACE
    A = np.array([1., 2., 3.27, 4.758501, 6.600001, 8.934502, 11.9703, 15.9495, 21.1349, 27.85261, 36.50411, 47.58061, 61.67791, 79.51341, 101.944, 130.051, 165.079, 208.497, 262.0211, 327.6431, 407.6571, 504.6801, 621.6801, 761.9842, 929.2942, 1127.69, 1364.34, 1645.71, 1979.16, 2373.041, 2836.781, 3381.001, 4017.541, 4764.391, 5638.791, 6660.341, 7851.231, 9236.572, 10866.3, 12783.7, 15039.3, 17693., 20119.2, 21686.5, 22436.3, 22389.8, 21877.6, 21215., 20325.9, 19309.7, 18161.9, 16960.9, 15626., 14291., 12869.59, 11895.86, 10918.17, 9936.521, 8909.992, 7883.422, 7062.198, 6436.264, 5805.321, 5169.611, 4533.901, 3898.201, 3257.081, 2609.201, 1961.311, 1313.48, 659.3752, 4.804826, 0.])
    B = np.array([0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 8.175413e-09, 0.006960025, 0.02801004, 0.06372006, 0.1136021, 0.1562241, 0.2003501, 0.2467411, 0.2944031, 0.3433811, 0.3928911, 0.4437402, 0.4945902, 0.5463042, 0.5810415, 0.6158184, 0.6506349, 0.6858999, 0.721166, 0.7493782, 0.7706375, 0.7919469, 0.8133039, 0.8346609, 0.856018, 0.877429, 0.898908, 0.920387, 0.941865, 0.963406, 0.984952, 1.])

    p_swap = np.zeros((len(edgeLevs), count_time+1, len(latitude), len(longitude))) #Pa
    for m in range(len(edgeLevs)):
        p_swap[m] =  A[edgeLevs[m]] + B[edgeLevs[m]] * sp[:,:,:]
        
    # make cwv_Calc vector
    q_swap = np.swapaxes(q,0,1)
    cwCalc_swap = np.zeros((np.shape(q_swap))) #kg/m2    
    for n in range(len(grabLevs)):
        cwCalc_swap[n] = (np.squeeze(q_swap[n])*np.squeeze(p_swap[n+1] - p_swap[n])) / g # column water vapor = specific humidity * pressure levels length / g [kg/m2]
        # plotDiagnostics(np.squeeze(q_swap[n,1,:,:]),'q_swap'+str(n))
    cwv_Calc = np.swapaxes(cwCalc_swap,0,1)

    # make tcwv_Calc vector
    tcwv_Calc = np.squeeze(np.sum(cwv_Calc,1)) #total column water vapor, cwv_Calc is summed over the vertical [kg/m2]

#    # make cw vector
    cwCorrected_swap = np.zeros((np.shape(cwCalc_swap)))
    for n in range(len(grabLevs)):
        cwCorrected_swap[n] = (tcwRaw / tcwv_Calc) * np.squeeze(cwCalc_swap[n])
    cwCorrected = np.swapaxes(cwCorrected_swap,0,1)
    
    # just a test, return this variable when running tests
    # tcwCorrected = np.squeeze(np.sum(cwCorrected,1)) #total column water, cw is summed over the vertical [kg/m2]    
    testvar = tcwRaw - tcwv_Calc # should be around zero for most cells

    # put A_gridcell on a 3D grid
    A_gridcell2D = np.tile(A_gridcell,[1,len(longitude)])
    A_gridcell_1_2D = np.reshape(A_gridcell2D, [1,len(latitude),len(longitude)])
    A_gridcell_plus3D = np.tile(A_gridcell_1_2D,[count_time+1,1,1])
    
    # water volumes        
    vapor_top = np.squeeze(np.sum(cwv_Calc[:,0:boundary,:,:],1))
    vapor_down = np.squeeze(np.sum(cwv_Calc[:,boundary:,:,:],1))
    vapor = vapor_top + vapor_down
    W_top = tcwRaw * (vapor_top / vapor) * A_gridcell_plus3D / density_water #m3
    W_down = tcwRaw * (vapor_down / vapor) * A_gridcell_plus3D / density_water #m3

    return cwCorrected, W_top, W_down, p_swap,testvar 

#%% Code
def getFa(latnrs,lonnrs,boundary,cw,U,V,ecorrflux,ncorrflux,count_time,yearnumber):
    
    #eastward and northward fluxes
    Fa_E_p = U * cw
    Fa_N_p = V * cw

    # uncorrected down and top fluxes
    Fa_E_down_uncorr = np.squeeze(np.sum(Fa_E_p[:,boundary:,:,:],1)) #kg*m-1*s-1
    Fa_N_down_uncorr = np.squeeze(np.sum(Fa_N_p[:,boundary:,:,:],1)) #kg*m-1*s-1
    Fa_E_top_uncorr = np.squeeze(np.sum(Fa_E_p[:,0:boundary,:,:],1)) #kg*m-1*s-1
    Fa_N_top_uncorr = np.squeeze(np.sum(Fa_N_p[:,0:boundary,:,:],1)) #kg*m-1*s-1
    
    # # ##############################################################
    # Using EAB's method which is faster
    corr_E = np.zeros([count_time+1,len(latitude),len(longitude)])   
    corr_N = np.zeros([count_time+1,len(latitude),len(longitude)])
    corr_E = np.minimum(2,np.maximum(0,ecorrflux/(Fa_E_down_uncorr + Fa_E_top_uncorr)))
    corr_N = np.minimum(2,np.maximum(0,ncorrflux/(Fa_N_down_uncorr + Fa_N_top_uncorr)))
    Fa_E_down = corr_E * Fa_E_down_uncorr #kg*m-1*s-1
    Fa_N_down = corr_N * Fa_N_down_uncorr #kg*m-1*s-1
    Fa_E_top = corr_E * Fa_E_top_uncorr #kg*m-1*s-1
    Fa_N_top = corr_N * Fa_N_top_uncorr #kg*m-1*s-1

    # make the fluxes during the timestep
    Fa_E_down = 0.5*(Fa_E_down[0:-1,:,:]+Fa_E_down[1:,:,:]);
    Fa_N_down = 0.5*(Fa_N_down[0:-1,:,:]+Fa_N_down[1:,:,:]);
    Fa_E_top = 0.5*(Fa_E_top[0:-1,:,:]+Fa_E_top[1:,:,:]);
    Fa_N_top = 0.5*(Fa_N_top[0:-1,:,:]+Fa_N_top[1:,:,:]);
    
    return Fa_E_top,Fa_N_top,Fa_E_down,Fa_N_down,corr_E,corr_N


#%% Code
def getEP(latnrs,lonnrs,yearnumber,count_time,latitude,longitude,A_gridcell,lsm):
        
    #(accumulated after the forecast at 00.00 and 12.00 by steps of 3 hours in time
    # EVAPORATION FROM TURBULENCE
    ncid = Dataset(datapath[6], mode = 'r')
    turbEvaporation24 = ncid.variables[evapVarName][timeincEP,latnrs,lonnrs] #m
    ncid.close()
    
    # LAND EVAP ONLY
    ncid = Dataset(datapath[10], mode = 'r')
    landEvaporation24 = ncid.variables[evapLandVarName][timeincEP,latnrs,lonnrs] #m
    ncid.close()

    #swap MERRA grid to ERA grid
    turbEvaporation24 = turbEvaporation24[:,ilat_SwapGridtoERA,:]
    turbEvaporation24 = turbEvaporation24[:,:,ilon_SwapGridtoERA] 
    turbEvaporation_orig = turbEvaporation24
    landEvaporation24 = landEvaporation24[:,ilat_SwapGridtoERA,:]
    landEvaporation24 = landEvaporation24[:,:,ilon_SwapGridtoERA]
    # done swapping grid          

    
    # RAW PRECIPITATION DATA 
    # ncid = Dataset(datapath[6], mode = 'r')
    # precipitation24 = ncid.variables[precipVarName][timeincEP,latnrs,lonnrs] #m, SEP 1, 2020: Commenting this out in lieu of PRECTOTCORR
    # ncid.close()
    # BIAS CORRECTED PRECIPITATION DATA
    ncid = Dataset(datapath[9], mode = 'r')
    precipitation24 = ncid.variables[precipVarName][timeincEP,latnrs,lonnrs] #m
    ncid.close()
    
    #swap MERRA grid to ERA grid
    precipitation24 = precipitation24[:,ilat_SwapGridtoERA,:]
    precipitation24 = precipitation24[:,:,ilon_SwapGridtoERA] 
    # done swapping grid     
    
    
    landEvaporation24[landEvaporation24 > 10000 ] = 0    
    landEvaporation24 = landEvaporation24 * lsm
    
    # create ocean mask
    ocean_mask = np.abs(lsm-1)

    turbEvaporation24 = turbEvaporation24 * ocean_mask
    
    evaporation24 = turbEvaporation24 + landEvaporation24
    
    ##***multiplying evaporation by -1 for MERRA2***
    evaporation24 = -1.0*evaporation24

#    # *******************************    
    # data comes in kg/m^2/s average rate over a 1 hour period
    # we want the average rate every 3 hours, which is still in kg/m^2/s
#    # *******************************    
        
    evaporation = np.zeros((8,np.shape(evaporation24)[1],np.shape(evaporation24)[2]))
    precipitation = np.zeros((8,np.shape(precipitation24)[1],np.shape(precipitation24)[2]))
    
    for x in np.arange(0,8):
        evaporation[x,:,:] = np.mean(evaporation24[x*3:x*3+2+1,:,:],axis=0)
        precipitation[x,:,:] = np.mean(precipitation24[x*3:x*3+2+1,:,:],axis=0)
    
    #delete and transfer negative values, change sign convention to all positive
    precipitation = np.reshape(np.maximum(np.reshape(precipitation, (np.size(precipitation))) + np.maximum(np.reshape(evaporation, (np.size(evaporation))),0.0),0.0),
                        (np.int(count_time*2),len(latitude),len(longitude))) 
    evaporation = np.reshape(np.abs(np.minimum(np.reshape(evaporation, (np.size(evaporation))),0.0)),(np.int(count_time*2),len(latitude),len(longitude)))   
            
    #convert from kg m-2 s-1 TO meters per day, for which the code was originally written (m d-1)
    evaporation_converted = evaporation*10800./1000. # 10800 = number of seconds in three hours; 1000 is the number of kg in a cubic meter of water
    precipitation_converted = precipitation*10800./1000.
    
    #calculate volumes
    A_gridcell2D = np.tile(A_gridcell,[1,len(longitude)])
    A_gridcell_1_2D = np.reshape(A_gridcell2D, [1,len(latitude),len(longitude)])
    A_gridcell_max3D = np.tile(A_gridcell_1_2D,[count_time*2,1,1])

    E = evaporation_converted * A_gridcell_max3D
    P = precipitation_converted * A_gridcell_max3D
    
    return E, P

#%% Code
def getrefined(Fa_E_top,Fa_N_top,Fa_E_down,Fa_N_down,W_top,W_down,E,P,divt,count_time,latitude,longitude):
    
    #for 3 hourly information
    divt2 = divt/2.
    oddvector2 = np.zeros((1,np.int(count_time*2*divt2)))
    partvector2 = np.zeros((1,np.int(count_time*2*divt2)))
    da = np.arange(1,divt2)

    
    for o in np.arange(0,np.int(count_time*2*divt2),24): #NOTE PWK: Changed on Sept 15, 2020. This was hard-coded as 12, which is half what the former divt of 24 was. I changed divt today, to 72, to see whether I could create a more temporally resolved model. so I manually changed this to 36. 
        for i in range(len(da)):
            oddvector2[0,o+i]    = (divt2-da[i])/divt2
            partvector2[0,o+i+1] = da[i]/divt2

    E_small = np.nan*np.zeros((np.int(count_time*2*divt2),len(latitude),len(longitude)))
    for t in range(1,np.int(count_time*2*divt2)+1):
        E_small[t-1] = (1./divt2) * E[np.int(t/divt2+oddvector2[0,t-1]-1)]
    E = E_small            
    
    P_small = np.nan*np.zeros((np.int(count_time*2*divt2),len(latitude),len(longitude)))
    for t in range(1,np.int(count_time*2*divt2)+1):
        P_small[t-1] = (1./divt2) * P[np.int(t/divt2+oddvector2[0,t-1]-1)]
    P = P_small
    
    # for 6 hourly info
    oddvector = np.zeros((1,np.int(count_time*divt)))
    partvector = np.zeros((1,np.int(count_time*divt)))
    da = np.arange(1,divt) 
    divt = np.float(divt)
    for o in np.arange(0,np.int(count_time*divt),np.int(divt)):
        for i in range(len(da)):
            oddvector[0,o+i]    = (divt-da[i])/divt
            partvector[0,o+i+1] = da[i]/divt    
        
    W_top_small = np.nan*np.zeros((np.int(count_time*divt+1),len(latitude),len(longitude)))
    for t in range(1,np.int(count_time*divt)+1):
        W_top_small[t-1] = W_top[np.int(t/divt+oddvector[0,t-1]-1)] + partvector[0,t-1] * (W_top[np.int(t/divt+oddvector[0,t-1])] - W_top[np.int(t/divt+oddvector[0,t-1]-1)])
        W_top_small[-1] = W_top[-1]
    W_top = W_top_small

    W_down_small = np.nan*np.zeros((np.int(count_time*divt+1),len(latitude),len(longitude)))
    for t in range(1,np.int(count_time*divt)+1):
        W_down_small[t-1] = W_down[np.int(t/divt+oddvector[0,t-1]-1)] + partvector[0,t-1] * (W_down[np.int(t/divt+oddvector[0,t-1])] - W_down[np.int(t/divt+oddvector[0,t-1]-1)])
        W_down_small[-1] = W_down[-1]
    W_down = W_down_small

    Fa_E_down_small = np.nan*np.zeros((np.int(count_time*divt),len(latitude),len(longitude)))
    Fa_N_down_small = np.nan*np.zeros((np.int(count_time*divt),len(latitude),len(longitude)))
    Fa_E_top_small = np.nan*np.zeros((np.int(count_time*divt),len(latitude),len(longitude)))
    Fa_N_top_small = np.nan*np.zeros((np.int(count_time*divt),len(latitude),len(longitude)))
    
    for t in range(1,np.int(count_time*divt)+1):
        Fa_E_down_small[t-1] = Fa_E_down[np.int(t/divt+oddvector[0,t-1]-1)]
        Fa_N_down_small[t-1] = Fa_N_down[np.int(t/divt+oddvector[0,t-1]-1)]
        Fa_E_top_small[t-1] = Fa_E_top[np.int(t/divt+oddvector[0,t-1]-1)]
        Fa_N_top_small[t-1] = Fa_N_top[np.int(t/divt+oddvector[0,t-1]-1)]

    Fa_E_down = Fa_E_down_small
    Fa_N_down = Fa_N_down_small
    Fa_E_top = Fa_E_top_small
    Fa_N_top = Fa_N_top_small
    
    return Fa_E_top,Fa_N_top,Fa_E_down,Fa_N_down,E,P,W_top,W_down


#%% Code
def get_stablefluxes(W_top,W_down,Fa_E_top_1,Fa_E_down_1,Fa_N_top_1,Fa_N_down_1,
                                   timestep,divt,L_EW_gridcell,density_water,L_N_gridcell,L_S_gridcell,latitude):
    
    #redefine according to units
    Fa_E_top_kgpmps = Fa_E_top_1
    Fa_E_down_kgpmps = Fa_E_down_1
    Fa_N_top_kgpmps = Fa_N_top_1
    Fa_N_down_kgpmps = Fa_N_down_1
    
    #convert to m3
    Fa_E_top = Fa_E_top_kgpmps * timestep/np.float(divt) * L_EW_gridcell / density_water # [kg*m^-1*s^-1*s*m*kg^-1*m^3]=[m3]
    Fa_E_down = Fa_E_down_kgpmps * timestep/np.float(divt) * L_EW_gridcell / density_water # [s*m*kg*m^-1*s^-1*kg^-1*m^3]=[m3]

    Fa_N_top_swap = np.zeros((len(latitude),np.int(count_time*np.float(divt)),len(longitude)))
    Fa_N_down_swap = np.zeros((len(latitude),np.int(count_time*np.float(divt)),len(longitude)))
    Fa_N_top_kgpmps_swap = np.swapaxes(Fa_N_top_kgpmps,0,1)
    Fa_N_down_kgpmps_swap = np.swapaxes(Fa_N_down_kgpmps,0,1)
    for c in range(len(latitude)):
        Fa_N_top_swap[c] = Fa_N_top_kgpmps_swap[c] * timestep/np.float(divt) * 0.5 *(L_N_gridcell[c]+L_S_gridcell[c]) / density_water # [s*m*kg*m^-1*s^-1*kg^-1*m^3]=[m3]
        Fa_N_down_swap[c] = Fa_N_down_kgpmps_swap[c] * timestep/np.float(divt) * 0.5*(L_N_gridcell[c]+L_S_gridcell[c]) / density_water # [s*m*kg*m^-1*s^-1*kg^-1*m^3]=[m3]

    Fa_N_top = np.swapaxes(Fa_N_top_swap,0,1) 
    Fa_N_down = np.swapaxes(Fa_N_down_swap,0,1) 
    
    #find out where the negative fluxes are
    Fa_E_top_posneg = np.ones(np.shape(Fa_E_top))
    Fa_E_top_posneg[Fa_E_top < 0] = -1
    Fa_N_top_posneg = np.ones(np.shape(Fa_E_top))
    Fa_N_top_posneg[Fa_N_top < 0] = -1
    Fa_E_down_posneg = np.ones(np.shape(Fa_E_top))
    Fa_E_down_posneg[Fa_E_down < 0] = -1
    Fa_N_down_posneg = np.ones(np.shape(Fa_E_top))
    Fa_N_down_posneg[Fa_N_down < 0] = -1
    
    #make everything absolute   
    Fa_E_top_abs = np.abs(Fa_E_top)
    Fa_E_down_abs = np.abs(Fa_E_down)
    Fa_N_top_abs = np.abs(Fa_N_top)
    Fa_N_down_abs = np.abs(Fa_N_down)
    
    # stabilize the outfluxes / influxes 
    stab = 1./2.  # during the reduced timestep the water cannot move further than 1/x * the gridcell, 
                    #in other words at least x * the reduced timestep is needed to cross a gridcell
    Fa_E_top_stable = np.reshape(np.minimum(np.reshape(Fa_E_top_abs, (np.size(Fa_E_top_abs))), (np.reshape(Fa_E_top_abs, (np.size(Fa_E_top_abs)))  / 
                                    (np.reshape(Fa_E_top_abs, (np.size(Fa_E_top_abs)))  + np.reshape(Fa_N_top_abs, (np.size(Fa_N_top_abs))))) * stab 
                                            * np.reshape(W_top[:-1,:,:], (np.size(W_top[:-1,:,:])))),(np.int(count_time*np.float(divt)),len(latitude),len(longitude)))
    Fa_N_top_stable = np.reshape(np.minimum(np.reshape(Fa_N_top_abs, (np.size(Fa_N_top_abs))), (np.reshape(Fa_N_top_abs, (np.size(Fa_N_top_abs)))  / 
                                    (np.reshape(Fa_E_top_abs, (np.size(Fa_E_top_abs)))  + np.reshape(Fa_N_top_abs, (np.size(Fa_N_top_abs))))) * stab 
                                            * np.reshape(W_top[:-1,:,:], (np.size(W_top[:-1,:,:])))),(np.int(count_time*np.float(divt)),len(latitude),len(longitude)))
    Fa_E_down_stable = np.reshape(np.minimum(np.reshape(Fa_E_down_abs, (np.size(Fa_E_down_abs))), (np.reshape(Fa_E_down_abs, (np.size(Fa_E_down_abs)))  / 
                                    (np.reshape(Fa_E_down_abs, (np.size(Fa_E_down_abs)))  + np.reshape(Fa_N_down_abs, (np.size(Fa_N_down_abs))))) * stab 
                                             * np.reshape(W_down[:-1,:,:], (np.size(W_down[:-1,:,:])))),(np.int(count_time*np.float(divt)),len(latitude),len(longitude)))
    Fa_N_down_stable = np.reshape(np.minimum(np.reshape(Fa_N_down_abs, (np.size(Fa_N_down_abs))), (np.reshape(Fa_N_down_abs, (np.size(Fa_N_down_abs)))  / 
                                    (np.reshape(Fa_E_down_abs, (np.size(Fa_E_down_abs)))  + np.reshape(Fa_N_down_abs, (np.size(Fa_N_down_abs))))) * stab 
                                             * np.reshape(W_down[:-1,:,:], (np.size(W_down[:-1,:,:])))),(np.int(count_time*np.float(divt)),len(latitude),len(longitude)))
    
    #get rid of the nan values
    Fa_E_top_stable[np.isnan(Fa_E_top_stable)] = 0
    Fa_N_top_stable[np.isnan(Fa_N_top_stable)] = 0
    Fa_E_down_stable[np.isnan(Fa_E_down_stable)] = 0
    Fa_N_down_stable[np.isnan(Fa_N_down_stable)] = 0
    
    #redefine
    Fa_E_top = Fa_E_top_stable * Fa_E_top_posneg
    Fa_N_top = Fa_N_top_stable * Fa_N_top_posneg
    Fa_E_down = Fa_E_down_stable * Fa_E_down_posneg
    Fa_N_down = Fa_N_down_stable * Fa_N_down_posneg
    
    
    return Fa_E_top,Fa_E_down,Fa_N_top,Fa_N_down

#%% Code
def getFa_Vert(Fa_E_top,Fa_E_down,Fa_N_top,Fa_N_down,E,P,W_top,W_down,divt,count_time,latitude,longitude,boundary):
    
    #total moisture in the column
    W = W_top + W_down
    
    #define the horizontal fluxes over the boundaries
    # fluxes over the eastern boundary
    Fa_E_top_boundary = np.zeros(np.shape(Fa_E_top))
    Fa_E_top_boundary[:,:,:-1] = 0.5 * (Fa_E_top[:,:,:-1] + Fa_E_top[:,:,1:])
    if isglobal == 1:
        Fa_E_top_boundary[:,:,-1] = 0.5 * (Fa_E_top[:,:,-1] + Fa_E_top[:,:,0])
    Fa_E_down_boundary = np.zeros(np.shape(Fa_E_down))
    Fa_E_down_boundary[:,:,:-1] = 0.5 * (Fa_E_down[:,:,:-1] + Fa_E_down[:,:,1:])
    if isglobal == 1:
        Fa_E_down_boundary[:,:,-1] = 0.5 * (Fa_E_down[:,:,-1] + Fa_E_down[:,:,0])

    # find out where the positive and negative fluxes are
    Fa_E_top_pos = np.ones(np.shape(Fa_E_top))
    Fa_E_down_pos = np.ones(np.shape(Fa_E_down))
    Fa_E_top_pos[Fa_E_top_boundary < 0] = 0
    Fa_E_down_pos[Fa_E_down_boundary < 0] = 0
    Fa_E_top_neg = Fa_E_top_pos - 1
    Fa_E_down_neg = Fa_E_down_pos - 1

    # separate directions west-east (all positive numbers)
    Fa_E_top_WE = Fa_E_top_boundary * Fa_E_top_pos
    Fa_E_top_EW = Fa_E_top_boundary * Fa_E_top_neg
    Fa_E_down_WE = Fa_E_down_boundary * Fa_E_down_pos
    Fa_E_down_EW = Fa_E_down_boundary * Fa_E_down_neg

    # fluxes over the western boundary
    Fa_W_top_WE = np.nan*np.zeros(np.shape(P))
    Fa_W_top_WE[:,:,1:] = Fa_E_top_WE[:,:,:-1]
    Fa_W_top_WE[:,:,0] = Fa_E_top_WE[:,:,-1]
    Fa_W_top_EW = np.nan*np.zeros(np.shape(P))
    Fa_W_top_EW[:,:,1:] = Fa_E_top_EW[:,:,:-1]
    Fa_W_top_EW[:,:,0] = Fa_E_top_EW[:,:,-1]
    Fa_W_down_WE = np.nan*np.zeros(np.shape(P))
    Fa_W_down_WE[:,:,1:] = Fa_E_down_WE[:,:,:-1]
    Fa_W_down_WE[:,:,0] = Fa_E_down_WE[:,:,-1]
    Fa_W_down_EW = np.nan*np.zeros(np.shape(P))
    Fa_W_down_EW[:,:,1:] = Fa_E_down_EW[:,:,:-1]
    Fa_W_down_EW[:,:,0] = Fa_E_down_EW[:,:,-1]    

    # fluxes over the northern boundary
    Fa_N_top_boundary = np.nan*np.zeros(np.shape(Fa_N_top))
    Fa_N_top_boundary[:,1:,:] = 0.5 * ( Fa_N_top[:,:-1,:] + Fa_N_top[:,1:,:] )
    Fa_N_down_boundary = np.nan*np.zeros(np.shape(Fa_N_down));
    Fa_N_down_boundary[:,1:,:] = 0.5 * ( Fa_N_down[:,:-1,:] + Fa_N_down[:,1:,:] )

    # find out where the positive and negative fluxes are
    Fa_N_top_pos = np.ones(np.shape(Fa_N_top))
    Fa_N_down_pos = np.ones(np.shape(Fa_N_down))
    Fa_N_top_pos[Fa_N_top_boundary < 0] = 0
    Fa_N_down_pos[Fa_N_down_boundary < 0] = 0
    Fa_N_top_neg = Fa_N_top_pos - 1
    Fa_N_down_neg = Fa_N_down_pos - 1

    # separate directions south-north (all positive numbers)
    Fa_N_top_SN = Fa_N_top_boundary * Fa_N_top_pos
    Fa_N_top_NS = Fa_N_top_boundary * Fa_N_top_neg
    Fa_N_down_SN = Fa_N_down_boundary * Fa_N_down_pos
    Fa_N_down_NS = Fa_N_down_boundary * Fa_N_down_neg

    # fluxes over the southern boundary
    Fa_S_top_SN = np.nan*np.zeros(np.shape(P))
    Fa_S_top_SN[:,:-1,:] = Fa_N_top_SN[:,1:,:]
    Fa_S_top_NS = np.nan*np.zeros(np.shape(P))
    Fa_S_top_NS[:,:-1,:] = Fa_N_top_NS[:,1:,:]
    Fa_S_down_SN = np.nan*np.zeros(np.shape(P))
    Fa_S_down_SN[:,:-1,:] = Fa_N_down_SN[:,1:,:]
    Fa_S_down_NS = np.nan*np.zeros(np.shape(P))
    Fa_S_down_NS[:,:-1,:] = Fa_N_down_NS[:,1:,:]
    
    # check the water balance
    Sa_after_Fa_down = np.zeros([1,len(latitude),len(longitude)])
    Sa_after_Fa_top = np.zeros([1,len(latitude),len(longitude)])
    Sa_after_all_down = np.zeros([1,len(latitude),len(longitude)])
    Sa_after_all_top = np.zeros([1,len(latitude),len(longitude)])
    residual_down = np.zeros(np.shape(P)) # residual factor [m3]
    residual_top = np.zeros(np.shape(P)) # residual factor [m3]

    for t in range(np.int(count_time*divt)):
        # down: calculate with moisture fluxes:
        Sa_after_Fa_down[0,1:-1,:] = (W_down[t,1:-1,:] - Fa_E_down_WE[t,1:-1,:] + Fa_E_down_EW[t,1:-1,:] + Fa_W_down_WE[t,1:-1,:] - Fa_W_down_EW[t,1:-1,:] - Fa_N_down_SN[t,1:-1,:] 
                                     + Fa_N_down_NS[t,1:-1,:] + Fa_S_down_SN[t,1:-1,:] - Fa_S_down_NS[t,1:-1,:])

        # top: calculate with moisture fluxes:
        Sa_after_Fa_top[0,1:-1,:] = (W_top[t,1:-1,:]- Fa_E_top_WE[t,1:-1,:] + Fa_E_top_EW[t,1:-1,:] + Fa_W_top_WE[t,1:-1,:] - Fa_W_top_EW[t,1:-1,:] - Fa_N_top_SN[t,1:-1,:] 
                                     + Fa_N_top_NS[t,1:-1,:] + Fa_S_top_SN[t,1:-1,:]- Fa_S_top_NS[t,1:-1,:])
    
        # down: substract precipitation and add evaporation
        Sa_after_all_down[0,1:-1,:] = Sa_after_Fa_down[0,1:-1,:] - P[t,1:-1,:] * (W_down[t,1:-1,:] / W[t,1:-1,:]) + E[t,1:-1,:]
    
        # top: substract precipitation
        Sa_after_all_top[0,1:-1,:] = Sa_after_Fa_top[0,1:-1,:] - P[t,1:-1,:] * (W_top[t,1:-1,:] / W[t,1:-1,:])

        # down: calculate the residual
        residual_down[t,1:-1,:] = W_down[t+1,1:-1,:] - Sa_after_all_down[0,1:-1,:]
    
        # top: calculate the residual
        residual_top[t,1:-1,:] = W_top[t+1,1:-1,:] - Sa_after_all_top[0,1:-1,:]


    # compute the resulting vertical moisture flux
    Fa_Vert_raw = W_down[1:,:,:] / W[1:,:,:] * (residual_down + residual_top) - residual_down # the vertical velocity so that the new residual_down/W_down =  residual_top/W_top (positive downward)

    # find out where the negative vertical flux is
    Fa_Vert_posneg = np.ones(np.shape(Fa_Vert_raw))
    Fa_Vert_posneg[Fa_Vert_raw < 0] = -1

    # make the vertical flux absolute
    Fa_Vert_abs = np.abs(Fa_Vert_raw)

    # stabilize the outfluxes / influxes
    stab = 1./4. #during the reduced timestep the vertical flux can maximally empty/fill 1/x of the top or down storage
    
    Fa_Vert_stable = np.reshape(np.minimum(np.reshape(Fa_Vert_abs, (np.size(Fa_Vert_abs))), np.minimum(stab*np.reshape(W_top[1:,:,:], (np.size(W_top[1:,:,:]))), 
                                        stab*np.reshape(W_down[1:,:,:], (np.size(W_down[1:,:,:]))))),(np.int(count_time*np.float(divt)),len(latitude),len(longitude)))
                
    # redefine the vertical flux
    Fa_Vert = Fa_Vert_stable * Fa_Vert_posneg;
    
    return Fa_Vert_raw, Fa_Vert

# #### End of code

#%% Runtime & Results

start1 = timer()

# obtain the constants
latitude,longitude,ilat_SwapGridtoERA,ilon_SwapGridtoERA,lsm,g,density_water,timestep,A_gridcell,L_N_gridcell,L_S_gridcell,L_EW_gridcell,gridcellLat,gridcellLon = getconstants(latnrs,lonnrs,lake_mask,invariant_data, lsmName='FRLAND')

def plotDiagnostics(data, figName, lats=latitude, lons=longitude,vmin=None, vmax=None):
    plt.figure(figsize=(10,6))
    ax = plt.subplot(1,1,1,projection=ct.crs.PlateCarree())
    cb,image = plotMap(ax,data,lats,lons,vmin=vmin,vmax=vmax)
    plt.title(figName)
    plt.savefig('figures/' + figName, dpi=400, bbox_inches='tight')
    plt.show()
    
# make area map
area_grid = np.tile(A_gridcell, (1,len(longitude)))

# delcare U,V,q as empty
U_tmrw, V_tmrw, q_tmrw, sp_tmrw, tcw_tmrw, ecorrflux_tmrw, ncorrflux_tmrw = [], [], [], [], [], [], []

# loop through the years
for yearnumber in years:
    
    daysinYear = (datetime.date(yearnumber,12,31) - datetime.date(yearnumber,1,1)).days+1
    print(daysinYear)
    for iDOY in np.arange(start_iDOY,end_iDOY):  #np.arange(0,daysinYear):
        start = timer() # start timing for each day

        datapath = data_path(iDOY,yearnumber) # global variable
                
        #2 get q, ps and integrate specific humidity to get the (total) column water (vapor)
        sp, sp_tmrw = load_sp(iDOY,yearnumber,sp_tmrw)
        q, q_tmrw = load_q(iDOY,yearnumber,q_tmrw)
        tcw, tcw_tmrw = load_tcw(iDOY,yearnumber,tcw_tmrw)        
        cw,W_top,W_down, presTest, testvar = getW(q,sp,tcw,latnrs,lonnrs,count_time,density_water,latitude,longitude,g,A_gridcell,boundary)
        del q, sp
        plotDiagnostics(testvar[0,:,:],'testvar')
        plotDiagnostics(W_top[0,:,:]/W_down[0,:,:],'Wtop div by Wdown')


        #3 get winds and calculate horizontal moisture fluxes with corrections
        U, V, U_tmrw, V_tmrw = load_uv(iDOY, yearnumber, U_tmrw, V_tmrw)      
        ecorrflux, ncorrflux, ecorrflux_tmrw, ncorrflux_tmrw = load_corrflux(iDOY, yearnumber, ecorrflux_tmrw, ncorrflux_tmrw)
        Fa_E_top,Fa_N_top,Fa_E_down,Fa_N_down,corr_E,corr_N = getFa(latnrs,lonnrs,boundary,cw,U,V,ecorrflux, ncorrflux, count_time,yearnumber)
        plotDiagnostics(corr_E[0,:,:],'corr E')
        plotDiagnostics(corr_N[0,:,:],'corr N')
        del U,V,ecorrflux,ncorrflux,corr_E,corr_N
        
        #4 evaporation and precipitation
        E,P = getEP(latnrs,lonnrs,yearnumber,count_time,latitude,longitude,A_gridcell,lsm)
        
        #5 put data on a smaller time step
        Fa_E_top_1,Fa_N_top_1,Fa_E_down_1,Fa_N_down_1,E,P,W_top,W_down = getrefined(Fa_E_top,Fa_N_top,Fa_E_down,Fa_N_down,W_top,W_down,E,P,divt,count_time,latitude,longitude)

        #6 stabilize horizontal fluxes and get everything in (m3 per smaller timestep)
        Fa_E_top,Fa_E_down,Fa_N_top,Fa_N_down = get_stablefluxes(W_top,W_down,Fa_E_top_1,Fa_E_down_1,Fa_N_top_1,Fa_N_down_1,
                               timestep,divt,L_EW_gridcell,density_water,L_N_gridcell,L_S_gridcell,latitude)
        
        #7 determine the vertical moisture flux
        Fa_Vert_raw,Fa_Vert = getFa_Vert(Fa_E_top,Fa_E_down,Fa_N_top,Fa_N_down,E,P,W_top,W_down,divt,count_time,latitude,longitude,boundary)

        #8 save the data        
        sio.savemat(datapath[11], {'Fa_E_top':Fa_E_top, 'Fa_N_top':Fa_N_top, 'Fa_E_down':Fa_E_down,'Fa_N_down':Fa_N_down, 'Fa_Vert':Fa_Vert, 'E':E, 'P':P,  
                                                                                'W_top':W_top, 'W_down':W_down}, do_compression=True)
        del E, P
        # alternative, but slower and more spacious
        # np.savez_compressed(datapath[23],Fa_E_top,Fa_N_top,Fa_E_down,Fa_N_down,Fa_Vert,E,P,W_top,W_down)
 
        # output the timer results
        end = timer()
        print('================================================')
        print('Runtime fluxes_and_storages for day ' + str(iDOY+1) + ' in year ' + str(yearnumber) + ' is',(end - start),' seconds.')
        print('================================================')

    start_iDOY = 0 # restart initial day of year for next year if not zero so that next year begins Jan. 1        

end1 = timer()
print('')
print('')
print('================================================')
print('================================================')
print('The total runtime is',(end1-start1),' seconds.')
print('================================================')
print('================================================')
#%%







