# -*- coding: utf-8 -*-
"""
Created on Tue Jun 21 12:38:45 2016

@author: Ent00002
"""
import numpy as np
from netCDF4 import Dataset

landfrac = 0.5

def getconstants(latnrs,lonnrs,lake_mask,invariant_data,latName='lat',lonName='lon', lsmName='lsm'): # def getconstants in Python is the same as function in MATLAB. 
    
    # load the latitude and longitude from the invariants file
    latitude = Dataset(invariant_data, mode = 'r').variables[latName][latnrs] # [degrees north]    
    longitude = Dataset(invariant_data, mode = 'r').variables[lonName][lonnrs] # [degrees east]
    
    # flip lats for MERRA2
    ilat_SwapGridtoERA = np.argsort(-1*latitude)
    latitude[:] = latitude[ilat_SwapGridtoERA]
    
    # cut longitudes to run from 0 to 360
    longitude[longitude<-0.0001] = longitude[longitude<-0.0001]+360
    ilon_SwapGridtoERA = np.argsort(longitude)
    longitude[:] = longitude[ilon_SwapGridtoERA]
    
    # Create land-sea-mask (in this model lakes are considered part of the land)
    lsm = np.squeeze(Dataset(invariant_data, mode = 'r').variables[lsmName][0,latnrs,lonnrs]) # 0 = sea, 1 = land
    #swap MERRA grid to ERA grid
    lsm = lsm[ilat_SwapGridtoERA, :]
    lsm = lsm[:, ilon_SwapGridtoERA]    
    # done swapping grid      
    
    lsm[lsm<landfrac] = 0
    lsm[lsm>0] = 1
    lsm = np.array(lsm.astype('int'))
    
    # Constants 
    g = 9.80665 # [m/s2] from ERA-interim archive
    density_water = 1000 # [kg/m3]
    dg = 111089.56 # [m] length of 1 degree latitude
    timestep = 6*3600 # [s] timestep in the ERA-interim archive (watch out! P & E have 3 hour timestep)
    Erad = 6.371e6 # [m] Earth radius
    
    # Semiconstants
    gridcellLon = np.abs(longitude[1] - longitude[0]) # [degrees] grid cell size
    gridcellLat = np.abs(latitude[1] - latitude[0]) # [degrees] grid cell size
    
    # new area size calculation:
    lat_n_bound = np.minimum(90.0 , latitude + 0.5*gridcellLat)
    lat_s_bound = np.maximum(-90.0 , latitude - 0.5*gridcellLat)
    
    A_gridcell = np.zeros([len(latitude),1])
    A_gridcell[:,0] = (np.pi/180.0)*Erad**2 * abs( np.sin(lat_s_bound*np.pi/180.0) - np.sin(lat_n_bound*np.pi/180.0) ) * gridcellLon
  
    L_N_gridcell = gridcellLat * np.cos((latitude + gridcellLat / 2.0) * np.pi / 180.0) * dg # [m] length northern boundary of a cell
    L_S_gridcell = gridcellLat * np.cos((latitude - gridcellLat / 2.0) * np.pi / 180.0) * dg # [m] length southern boundary of a cell
    L_EW_gridcell = gridcellLon * dg # [m] length eastern/western boundary of a cell 
    
    return latitude , longitude, ilat_SwapGridtoERA, ilon_SwapGridtoERA, lsm , g , density_water , timestep , A_gridcell , L_N_gridcell , L_S_gridcell , L_EW_gridcell , gridcellLat, gridcellLon
