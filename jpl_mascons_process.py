# -*- coding: utf-8 -*-
"""
Script to read JPL mascon NetCDF file, fill temporal data gaps with interpolation and generate output file of required format

Input arguments
mascon_file: JPL mascon file with full path
tws_op_file: File name with full path of the NetCDF file to be generated
    
Outputs
Single NetCDF file which contains time interpolated TWS values in appropriate format from JPL mascons

This script calls the function create_nc_file which is also provided in the directory.
"""
# Importing modules
import netCDF4 as nc
import matplotlib.pyplot as plt
import numpy as np
from scipy import interpolate
from datetime import datetime
from create_nc_file import create_nc_file

# Specifying input information
mascon_file = "D:\IISC_COURSES\GEODETIC_SIGNAL_PROCESSING\PROJECT\TWS\DATA\GRCTellus.JPL.200204_202207.GLO.RL06.1M.MSCNv03CRI.nc"
tws_op_file = "D:\IISC_COURSES\GEODETIC_SIGNAL_PROCESSING\PROJECT\TWS\OUTPUTS\JPL_CRI_TWS_0P5_DEG.nc"

# Reading the mascon file and getting the required variables
ds = nc.Dataset(mascon_file)
tws_data_file = ds['lwe_thickness'][:]
time_var_file = ds['time']
time_data_file = time_var_file[:]
uncertainty_data_file = ds['uncertainty'][:]
lat_file = ds['lat'][:]
lat_bounds_file = ds['lat_bounds'][:]
lon_file = ds['lon'][:]
lon_bounds_file = ds['lon_bounds'][:]
time_bounds_var_file = ds['time_bounds']
time_bounds_file = time_bounds_var_file[:]

# Getting some information about the sizes
n_times = np.size(time_data_file)
n_latitudes = np.size(lat_file)
n_longitudes = np.size(lon_file)

# Converting to mm units from cm
tws_data_mm = tws_data_file * 1e1
uncertainty_data_mm = uncertainty_data_file * 1e1

# Spatial flipping of data to get latitudes from 90 to -90 and longitude from -180 to 180
tws_flipped = np.zeros(np.shape(tws_data_mm))
uncertainty_flipped = np.zeros(np.shape(uncertainty_data_mm))
half_n_lon = int(n_longitudes/2)
for i in range(0,n_times):
    tws_cur = tws_data_mm[i]
    tws_latflip = np.flipud(tws_cur)
    tws_flipped[i,:,0:half_n_lon] = tws_latflip[:,half_n_lon:]
    tws_flipped[i,:,half_n_lon:] = tws_latflip[:,0:half_n_lon]
    
    uncertainty_cur = uncertainty_data_mm[i]
    uncertainty_latflip = np.flipud(uncertainty_cur)
    uncertainty_flipped[i,:,0:half_n_lon] = uncertainty_latflip[:,half_n_lon:]
    uncertainty_flipped[i,:,half_n_lon:] = uncertainty_latflip[:,0:half_n_lon]

# Computing latitude values of the flipped array
lat_flipped = np.flipud(lat_file)
lat_bounds_flipped = np.flipud(lat_bounds_file)

# Computing the longitude values of the flipped array
lon_flipped = np.zeros(np.shape(lon_file))
lon_flipped[0:half_n_lon] = lon_file[half_n_lon:]
lon_flipped[half_n_lon:] = lon_file[0:half_n_lon]
lon_flipped[lon_flipped>180] = lon_flipped[lon_flipped>180] - 360

lon_bounds_flipped = np.zeros(np.shape(lon_bounds_file))
lon_bounds_flipped[0:half_n_lon,:] = lon_bounds_file[half_n_lon:,:]
lon_bounds_flipped[half_n_lon:,:] = lon_bounds_file[0:half_n_lon,:]
lon_bounds_flipped[lon_bounds_flipped>180] = lon_bounds_flipped[lon_bounds_flipped>180] - 360
lon_bounds_flipped[0,0] = -180

# Getting the missing times from data
missing_times = ds.__dict__['months_missing']
missing_times = missing_times.split(';')
##### There is a typo in the months_missing attribute. Last element is given as 2018-08-2018-09 but it should have been 2018-08;2018-09. Fixing that here
missing_times = missing_times[0:np.size(missing_times)-1] + ['2018-08'] + ['2018-09']
missing_times = np.array(missing_times)
n_missing_times = np.size(missing_times)

# Generating the required time object in proper format for the missing times
missing_time_value = np.zeros(n_missing_times)
missing_time_bounds = np.zeros((n_missing_times,2))
for i in range(0,n_missing_times):
    year_cur = int(missing_times[i][0:4])
    month_cur = int(missing_times[i][5:7])
    
    missing_time_value[i] = nc.date2num(datetime(year_cur,month_cur,15,0),time_var_file.units)
    missing_time_bounds[i,0] = nc.date2num(datetime(year_cur,month_cur,1,0),time_var_file.units)
    missing_time_bounds[i,1] = missing_time_value[i] + 10
    
# Generating the array of continuous times
time_continuous = np.sort(np.concatenate((time_data_file,missing_time_value)))
time_bounds_continuous_unsorted = np.concatenate((time_bounds_file,missing_time_bounds),axis=0)
time_bounds_continuous = time_bounds_continuous_unsorted[time_bounds_continuous_unsorted[:,0].argsort(),:]
n_times_continuous = np.size(time_continuous)

# Getting the indices where data is missing and data is available
index_data_missing = np.where(time_continuous == missing_time_value[0])[0]
for i in range(1,n_missing_times):
    ind = np.where(time_continuous == missing_time_value[i])[0]
    index_data_missing = np.concatenate((index_data_missing,ind))
index_data_avail = np.delete(np.arange(0,n_times_continuous),index_data_missing)

# Interpolating the TWS for missing times
tws_interp = np.zeros((n_times_continuous,n_latitudes,n_longitudes))
for i in range(0,n_latitudes):
    for j in range(0,n_longitudes):
        tck = interpolate.splrep(index_data_avail, tws_flipped[:,i,j], s=0)
        y_interp = interpolate.splev(index_data_missing, tck, der=0)
        tws_interp[index_data_missing,i,j] = y_interp
tws_interp[index_data_avail] = tws_flipped

# Set uncertainities at data missing times as the uncertainities in the previous month where data is available
uncertainty_interp = np.zeros((n_times_continuous,n_latitudes,n_longitudes))
uncertainty_interp[index_data_avail] = uncertainty_flipped
for i in range(0,n_missing_times):
    uncertainty_interp[index_data_missing[i]] = uncertainty_interp[index_data_missing[i]-1]

# Inflating uncertainities from the native 3 degree resolution to 0.5 degree resolution
uncertainty_inf = uncertainty_interp * 35/6

# Creating the UTC string and year-month array
time_continuous_utc = nc.num2date(time_continuous,time_var_file.units)
time_cont_utc_str = [str(time_continuous_utc[0].year) + '-' + str(time_continuous_utc[0].month).zfill(2) + '-' + str(time_continuous_utc[0].day).zfill(2) + 'T' + str(time_continuous_utc[0].hour).zfill(2) + ':' + str(time_continuous_utc[0].minute).zfill(2) + ':' + str(time_continuous_utc[0].second).zfill(2) + '.000']
year_month = [time_cont_utc_str[0][0:7]]
for i in range(1,n_times_continuous):
    time_utc_cur = [str(time_continuous_utc[i].year) + '-' + str(time_continuous_utc[i].month).zfill(2) + '-' + str(time_continuous_utc[i].day).zfill(2) + 'T' + str(time_continuous_utc[i].hour).zfill(2) + ':' + str(time_continuous_utc[i].minute).zfill(2) + ':' + str(time_continuous_utc[i].second).zfill(2) + '.000']
    time_cont_utc_str = np.concatenate((time_cont_utc_str,time_utc_cur))
    year_month = np.concatenate((year_month,[time_cont_utc_str[i][0:7]]))

# Writing to netCDF file
flag = create_nc_file(tws_op_file,tws_interp,uncertainty_inf,lat_flipped,lon_flipped,time_continuous,lat_bounds_flipped,lon_bounds_flipped,time_bounds_continuous,time_var_file.units,time_cont_utc_str,year_month)
