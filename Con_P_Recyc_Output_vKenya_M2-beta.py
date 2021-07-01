# -*- coding: utf-8 -*-
"""
Created on Thu Jun 16 13:24:45 2016

@author: Ent00002
"""
# #%% PAUSING CODE 
# import time
# print('pausing code for 37332 seconds')
# time.sleep(37332)
 
#%% Import libraries
import os
os.chdir('/Users/patrickkeys/rams_gdrive/work/research/models/WAM2layersPython-master/WAM2layersPython3/')

import numpy as np
import sys
from netCDF4 import Dataset
import scipy.io as sio
import calendar
import datetime
from getconstants import getconstants
from timeit import default_timer as timer
import copy
import time
import matplotlib.pyplot as plt

import cartopy as ct
from mpl_toolkits.axes_grid1 import make_axes_locatable

import matplotlib as mpl

import warnings
warnings.filterwarnings("ignore", category=RuntimeWarning)


#%% BEGIN OF INPUT1 (FILL THIS IN)
years = np.arange(1980,2020) #fill in the years
yearpart = np.arange(0,364) # for a full (leap)year fill in np.arange(0,366)

# Manage the extent of your dataset (FILL THIS IN)
# Define the latitude and longitude cell numbers to consider and corresponding lakes that should be considered part of the land
latnrs = np.arange(20,341)#(0,361)
lonnrs = np.arange(0,576)# Full MERRA2 range is (0,576)
isglobal = 1 # fill in 1 for global computations (i.e. Earth round), fill in 0 for a local domain with boundaries

lake_mask = np.zeros((1,2))

# obtain the constants
invariant_data = '/Users/patrickkeys/rams_gdrive/work/research/models/WAM2layersPython-master/WAM2layersPython3/data/constants/MERRA2_101.const_2d_asm_Nx.00000000.nc4'
latitude,longitude,ilat_SwapGridtoERA,ilon_SwapGridtoERA,lsm,g,density_water,timestep,A_gridcell,L_N_gridcell,L_S_gridcell,L_EW_gridcell,gridcellLat,gridcellLon = getconstants(latnrs,lonnrs,lake_mask,invariant_data, lsmName='FRLAND')

# BEGIN OF INPUT 2 (FILL THIS IN)
# Region = lsm
kenya_region_mat = sio.loadmat('/Users/patrickkeys/rams_gdrive/work/research/models/WAM2layersPython-master/WAM2layersPython3/kenya_grid.mat')['kenya_grid']
Region = kenya_region_mat
Region = Region[latnrs,:]
Region = Region[:,lonnrs]
Region = Region[ilat_SwapGridtoERA,:]
Region = Region[:,ilon_SwapGridtoERA]

Kvf = 3 # vertical dispersion factor (advection only is 0, dispersion the same size of the advective flux is 1, for stability don't make this more than 3)
daily = 0 # 1 for writing out daily data, 0 for only monthly data
timetracking = 0 # 0 for not tracking time and 1 for tracking time
veryfirstrun = 1 # type '1' if no run has been done before from which can be continued, otherwise type '0'

# # #interdata_folder = r'C:\Users\bec\Desktop\WAM2\interdata' #must be an existing folder, existence is not checked
output_dir = '/Volumes/onenight/merra_2/output/forward_tracking/2020_72/'

# # #interdata_folder = r'C:\Users\bec\Desktop\WAM2\interdata' #must be an existing folder, existence is not checked
# # nostromo
interdata_folder = '/Volumes/lindsey/merra_2/forward_tracking/Sa_track_files/' #must be an existing folder, existence is not checked
sub_interdata_folder = os.path.join(interdata_folder, '2020_72/')        

def get_FS_dir_wildNumStr(year):
    if(year>=1980 and year<1990):
        fluxes_and_states_dir = '/Volumes/coffey/merra_2/fluxes_and_states_2020_70/' #must be an existing folder, existence is not checked
    elif(year>=1990 and year<=2020):
        fluxes_and_states_dir = '/Volumes/kane/merra_2/fluxes_and_states_2020_70/' #must be an existing folder, existence is not checked
    else:
        sys.exit("Fluxes and States wildcard number string could not be determined. Check your input year.")
    
    return fluxes_and_states_dir

output_dir = '/Volumes/onenight/merra_2/output/forward_tracking/2020_72/'
#END OF INPUT

#%% Datapaths (FILL THIS IN)

def data_path(y,a):
    load_Sa_track = os.path.join(sub_interdata_folder, str(y) + '-' + str(a) + 'Sa_track.mat')
    
    load_Sa_time = os.path.join(sub_interdata_folder, str(y) + '-' + str(a) + 'Sa_time.mat')
    
    load_fluxes_and_storages = os.path.join(fluxes_and_states_dir, str(y) + '-' + str(a) + 'fluxes_storages.mat')
    
    save_path = os.path.join(output_dir, 'P_track_continental_full' + str(years[0]) + '-' + str(years[-1]) + '-timetracking' + str(timetracking) + '.mat')
    
    save_path_daily = os.path.join(output_dir, 'P_track_continental_daily_full' + str(y) + '-timetracking' + str(timetracking) + '.mat')
    
    return load_Sa_track,load_Sa_time,load_fluxes_and_storages,save_path,save_path_daily

#%% Runtime & Results

start1 = timer()
startyear = years[0]

E_per_year_per_month = np.zeros((len(years),12,len(latitude),len(longitude)))
P_track_per_year_per_month = np.zeros((len(years),12,len(latitude),len(longitude)))
P_per_year_per_month = np.zeros((len(years),12,len(latitude),len(longitude)))
Sa_track_down_per_year_per_month = np.zeros((len(years),12,len(latitude),len(longitude)))
Sa_track_top_per_year_per_month = np.zeros((len(years),12,len(latitude),len(longitude)))
W_down_per_year_per_month = np.zeros((len(years),12,len(latitude),len(longitude)))
W_top_per_year_per_month = np.zeros((len(years),12,len(latitude),len(longitude)))
north_loss_per_year_per_month = np.zeros((len(years),12,1,len(longitude)))
south_loss_per_year_per_month = np.zeros((len(years),12,1,len(longitude)))
down_to_top_per_year_per_month = np.zeros((len(years),12,len(latitude),len(longitude)))
top_to_down_per_year_per_month = np.zeros((len(years),12,len(latitude),len(longitude)))
water_lost_per_year_per_month = np.zeros((len(years),12,len(latitude),len(longitude)))
# shadow_storage_per_year_per_month = np.zeros((len(years),12,len(latitude),len(longitude)))

if timetracking == 1:
    Sa_time_down_per_year_per_month = np.zeros((len(years),12,len(latitude),len(longitude)))
    Sa_time_top_per_year_per_month = np.zeros((len(years),12,len(latitude),len(longitude)))
    P_time_per_year_per_month = np.zeros((len(years),12,len(latitude),len(longitude)))
       
for i in range(len(years)):
    y = years[i]
    ly = int(calendar.isleap(y))
    final_time = 364+ly
    
    fluxes_and_states_dir = get_FS_dir_wildNumStr(y)    
    
    E_per_day = np.zeros((365+ly,len(latitude),len(longitude)))
    P_track_per_day = np.zeros((365+ly,len(latitude),len(longitude)))
    P_per_day = np.zeros((365+ly,len(latitude),len(longitude)))
    Sa_track_down_per_day = np.zeros((365+ly,len(latitude),len(longitude)))
    Sa_track_top_per_day = np.zeros((365+ly,len(latitude),len(longitude)))
    W_down_per_day = np.zeros((365+ly,len(latitude),len(longitude)))
    W_top_per_day = np.zeros((365+ly,len(latitude),len(longitude)))
    north_loss_per_day = np.zeros((365+ly,1,len(longitude)))
    south_loss_per_day = np.zeros((365+ly,1,len(longitude)))
    down_to_top_per_day = np.zeros((365+ly,len(latitude),len(longitude)))
    top_to_down_per_day = np.zeros((365+ly,len(latitude),len(longitude)))
    water_lost_per_day = np.zeros((365+ly,len(latitude),len(longitude)))
    if timetracking == 1:
        Sa_time_down_per_day = np.zeros((365+ly,len(latitude),len(longitude)))
        Sa_time_top_per_day = np.zeros((365+ly,len(latitude),len(longitude)))
        P_time_per_day = np.zeros((365+ly,len(latitude),len(longitude)))
        
    for j in range(len(yearpart)):
        start = timer()
        a = yearpart[j]
        datapath = data_path(y,a)
        if a > final_time: # a = 365 (366th index) and not a leapyear\
            pass
        else:
            # load tracked data
            loading_ST = sio.loadmat(datapath[0],verify_compressed_data_integrity=False)
            Sa_track_top = loading_ST['Sa_track_top']
            Sa_track_down = loading_ST['Sa_track_down']
            north_loss = loading_ST['north_loss']
            south_loss = loading_ST['south_loss']
            down_to_top = loading_ST['down_to_top']
            top_to_down = loading_ST['top_to_down']
            water_lost = loading_ST['water_lost']
            Sa_track = Sa_track_top + Sa_track_down
            if timetracking == 1:
                loading_STT = sio.loadmat(datapath[1],verify_compressed_data_integrity=False)
                Sa_time_top = loading_STT['Sa_time_top']
                Sa_time_down = loading_STT['Sa_time_down']
            
            # load the total moisture data
            loading_FS = sio.loadmat(datapath[2],verify_compressed_data_integrity=False)
            Fa_E_top = loading_FS['Fa_E_top']
            Fa_N_top = loading_FS['Fa_N_top']
            Fa_E_down = loading_FS['Fa_E_down']
            Fa_N_down = loading_FS['Fa_N_down']
            Fa_Vert = loading_FS['Fa_Vert']
            E = loading_FS['E']
            P = loading_FS['P']
            W_top = loading_FS['W_top']
            W_down = loading_FS['W_down']

            W = W_top + W_down
            
            # compute tracked precipitation
            P_track = P[:,:,:] * (Sa_track[:-1,:,:] / W[:-1,:,:])
            
            # save per day
            E_per_day[a,:,:] = np.sum(E, axis =0)
            P_track_per_day[a,:,:] = np.sum(P_track, axis =0)
            P_per_day[a,:,:] = np.sum(P, axis =0)
            Sa_track_down_per_day[a,:,:] = np.mean(Sa_track_down[:-1,:,:], axis =0)
            Sa_track_top_per_day[a,:,:] = np.mean(Sa_track_top[:-1,:,:], axis =0)
            W_down_per_day[a,:,:] = np.mean(W_down[:-1,:,:], axis =0)
            W_top_per_day[a,:,:] = np.mean(W_top[:-1,:,:], axis =0)

            north_loss_per_day[a,:,:] = np.sum(north_loss, axis =0)
            south_loss_per_day[a,:,:] = np.sum(south_loss, axis =0)
            down_to_top_per_day[a,:,:] = np.sum(down_to_top, axis =0)
            top_to_down_per_day[a,:,:] = np.sum(top_to_down, axis =0)
            water_lost_per_day[a,:,:] = np.sum(water_lost, axis =0)

            if timetracking == 1:
                # compute tracked precipitation time
                P_track_down = P[:,:,:] * (Sa_track_down[:-1,:,:] / W[:-1,:,:])
                P_track_top = P[:,:,:] * (Sa_track_top[:-1,:,:] / W[:-1,:,:])
                P_time_down = 0.5 * ( Sa_time_down[:-1,:,:] + Sa_time_down[1:,:,:] ) # seconds
                P_time_top = 0.5 * ( Sa_time_top[:-1,:,:] + Sa_time_top[1:,:,:] ) # seconds
                
                # save per day
                Sa_time_down_per_day[a,:,:] = (np.mean(Sa_time_down[:-1,:,:] * Sa_track_down[:-1,:,:], axis = 0) 
                    / Sa_track_down_per_day[a,:,:]) # seconds
                Sa_time_top_per_day[a,:,:] = (np.mean(Sa_time_top[:-1,:,:] * Sa_track_top[:-1,:,:], axis = 0) 
                    / Sa_track_top_per_day[a,:,:]) # seconds
                P_time_per_day[a,:,:] = (np.sum((P_time_down * P_track_down + P_time_top * P_track_top), axis = 0) 
                    / P_track_per_day[a,:,:]) # seconds
                    
                # remove nans
                where_are_NaNs = np.isnan(P_time_per_day)
                P_time_per_day[where_are_NaNs] = 0
                            
        end = timer()
        print('Runtime output for day ' + str(a+1) + ' in year ' + str(y) + ' is',(end - start),' seconds.')
        
    if daily == 1:
        if timetracking == 0: # create dummy values
            Sa_time_down_per_day = 0
            Sa_time_top_per_day = 0
            P_time_per_day = 0
            
        sio.savemat(datapath[4],
                    {'E_per_day':E_per_day,'P_track_per_day':P_track_per_day,'P_per_day':P_per_day,
                     'Sa_track_down_per_day':Sa_track_down_per_day,'Sa_track_top_per_day':Sa_track_top_per_day, 
                     'Sa_time_down_per_day':Sa_time_down_per_day,'Sa_time_top_per_day':Sa_time_top_per_day, 
                     'W_down_per_day':W_down_per_day,'W_top_per_day':W_top_per_day,
                     'P_time_per_day':P_time_per_day},do_compression=True)   
                        
    # values per month
    for m in range(12):
        first_day = int(datetime.date(y,m+1,1).strftime("%j"))
        last_day = int(datetime.date(y,m+1,calendar.monthrange(y,m+1)[1]).strftime("%j"))
        days = np.arange(first_day,last_day+1)-1 # -1 because Python is zero-based
        
        E_per_year_per_month[y-startyear,m,:,:] = (np.squeeze(np.sum(E_per_day[days,:,:], axis = 0)))
        P_track_per_year_per_month[y-startyear,m,:,:] = (np.squeeze(np.sum(P_track_per_day[days,:,:], axis = 0)))
        P_per_year_per_month[y-startyear,m,:,:] = (np.squeeze(np.sum(P_per_day[days,:,:], axis = 0)))
        Sa_track_down_per_year_per_month[y-startyear,m,:,:] = (np.squeeze(np.mean(Sa_track_down_per_day[days,:,:], axis = 0)))
        Sa_track_top_per_year_per_month[y-startyear,m,:,:] = (np.squeeze(np.mean(Sa_track_top_per_day[days,:,:], axis = 0)))
        W_down_per_year_per_month[y-startyear,m,:,:] = (np.squeeze(np.mean(W_down_per_day[days,:,:], axis = 0)))
        W_top_per_year_per_month[y-startyear,m,:,:] = (np.squeeze(np.mean(W_top_per_day[days,:,:], axis = 0)))
        north_loss_per_year_per_month[y-startyear,m,:,:] = (np.squeeze(np.sum(north_loss_per_day[days,:,:], axis = 0)))
        south_loss_per_year_per_month[y-startyear,m,:,:] = (np.squeeze(np.sum(south_loss_per_day[days,:,:], axis = 0)))
        down_to_top_per_year_per_month[y-startyear,m,:,:] = (np.squeeze(np.sum(down_to_top_per_day[days,:,:], axis = 0)))
        top_to_down_per_year_per_month[y-startyear,m,:,:] = (np.squeeze(np.sum(top_to_down_per_day[days,:,:], axis = 0)))
        water_lost_per_year_per_month[y-startyear,m,:,:] = (np.squeeze(np.sum(water_lost_per_day[days,:,:], axis = 0)))
        
        if timetracking == 1:
            Sa_time_down_per_year_per_month[y-startyear,m,:,:] = ( np.squeeze( np.mean( Sa_time_down_per_day[days,:,:]
                * Sa_track_down_per_day[days,:,:], axis = 0)) 
                / np.squeeze(Sa_track_down_per_year_per_month[y-startyear,m,:,:]) )
            Sa_time_top_per_year_per_month[y-startyear,m,:,:] = ( np.squeeze( np.mean( Sa_time_top_per_day[days,:,:] 
                * Sa_track_top_per_day[days,:,:],axis = 0)) 
                / np.squeeze(Sa_track_top_per_year_per_month[y-startyear,m,:,:]) )
            P_time_per_year_per_month[y-startyear,m,:,:] = ( np.squeeze( np.sum( P_time_per_day[days,:,:] 
                * P_track_per_day[days,:,:], axis = 0)) 
                / np.squeeze(P_track_per_year_per_month[y-startyear,m,:,:]) )
        elif timetracking == 0: # dummy values
            Sa_time_down_per_year_per_month = 0
            Sa_time_top_per_year_per_month = 0
            P_time_per_year_per_month = 0
        
if timetracking == 1:
    where_are_NaNs = np.isnan(P_time_per_year_per_month)
    P_time_per_year_per_month[where_are_NaNs] = 0

# save monthly data
sio.savemat(datapath[3],
            {'E_per_year_per_month':E_per_year_per_month,'P_track_per_year_per_month':P_track_per_year_per_month,'P_per_year_per_month':P_per_year_per_month,
             'Sa_track_down_per_year_per_month':Sa_track_down_per_year_per_month,'Sa_track_top_per_year_per_month':Sa_track_top_per_year_per_month, 
             'Sa_time_down_per_year_per_month':Sa_time_down_per_year_per_month,'Sa_time_top_per_year_per_month':Sa_time_top_per_year_per_month, 
             'P_time_per_year_per_month':P_time_per_year_per_month,
             'W_down_per_year_per_month':W_down_per_year_per_month,'W_top_per_year_per_month':W_top_per_year_per_month,
             'north_loss_per_year_per_month':north_loss_per_year_per_month,'south_loss_per_year_per_month':south_loss_per_year_per_month,
             'down_to_top_per_year_per_month':down_to_top_per_year_per_month,'top_to_down_per_year_per_month':top_to_down_per_year_per_month,
             'water_lost_per_year_per_month':water_lost_per_year_per_month},do_compression=True)  
            
end1 = timer()
print('The total runtime of Con_P_Recyc_Output is',(end1-start1),' seconds.')
