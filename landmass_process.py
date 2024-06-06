# -*- coding: utf-8 -*-
"""
Script to read CSR/GFZ landmass data product and generate NetCDF files of TWS in the required format

Input arguments
landmass_folder: Path to directory containing CSR/GFZ landmass files for the entire duration of interest
scalefactor_file: Filename with full path of the scale factor file
tws_op_file: File name with full path of the NetCDF file to be generated

Outputs
Single NetCDF file which contains time interpolated TWS values in appropriate format from the monthly landmass files

This script calls the function create_nc_file which is also provided in the directory.
"""
# Importing modules
import netCDF4 as nc
import matplotlib.pyplot as plt
import numpy as np
from scipy import interpolate
from datetime import datetime
from create_nc_file import create_nc_file
import glob

# Specifying some input information
landmass_folder = "D:\IISC_COURSES\GEODETIC_SIGNAL_PROCESSING\PROJECT\TWS\DATA\CSR_LANDMASS"
scalefactor_file = "D:\IISC_COURSES\GEODETIC_SIGNAL_PROCESSING\PROJECT\TWS\DATA\CLM4.SCALE_FACTOR.DS.G300KM.RL05.DSTvSCS1409.nc"
tws_op_file = "D:\IISC_COURSES\GEODETIC_SIGNAL_PROCESSING\PROJECT\TWS\OUTPUTS\CSR_TWS_1_DEG.nc"

# Getting the list of land mass files
flist = glob.glob(landmass_folder+'\*.nc')
flist = np.sort(flist)
n_files = len(flist)

# Getting the required input information from the landmass files
ds = nc.Dataset(flist[0])
lat_file = ds['lat'][:]
lat_bounds_file = ds['lat_bounds'][:]
lon_file = ds['lon'][:]
lon_bounds_file = ds['lon_bounds'][:]
time_units = ds['time'].units
ds.close()

n_latitudes = np.size(lat_file)
n_longitudes = np.size(lon_file)
n_times = n_files

tws_data_file = np.zeros((n_files,n_latitudes,n_longitudes))
uncertainty_data_file = np.zeros((n_files,n_latitudes,n_longitudes))
time_data_file = np.zeros(n_files)
time_bounds_file = np.zeros((n_files,2))
for i in range(0,n_files):
    ds = nc.Dataset(flist[i])
    
    # Getting the TWS data and uncertainty
    tws_data_file[i] = ds['lwe_thickness'][:]
    uncertainty_data_file[i] = ds['uncertainty'][:]
    
    # Getting the time related information
    time_data_file[i] = ds['time'][:]
    time_bounds_file[i] = ds['time_bounds'][:]
    
# Converting to mm units from m
tws_data_mm = tws_data_file * 1e3
uncertainty_data_mm = uncertainty_data_file * 1e3

# Applying scale factors to the data
scale_ds = nc.Dataset(scalefactor_file)
scale_factor = scale_ds['SCALE_FACTOR'][:]
tws_data_mm = np.multiply(tws_data_mm,scale_factor)

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

# Identifying the times where data is missing by comparing with the array of continuous times in the period of interest
time_utc_data = nc.num2date(time_data_file,time_units)
month_year_data = [str(time_utc_data[0].year) + '-' + str(time_utc_data[0].month).zfill(2)]
for i in range(1,n_times):
    month_year_cur = [str(time_utc_data[i].year) + '-' + str(time_utc_data[i].month).zfill(2)]
    month_year_data = np.concatenate((month_year_data,month_year_cur))
############ Fixing the month-year identification issues for 2011-11 and 2015-05
month_year_data[110] = '2011-11'
month_year_data[144] = '2015-05'

month_year_full_match = []
start_year = int(month_year_data[0][0:4])
end_year = int(month_year_data[-1][0:4])
years_interest = np.arange(start_year,end_year+1)
months_interest = np.arange(1,13)
for i in range(0,np.size(years_interest)):
    for j in range(0,np.size(months_interest)):
        month_year_full_match = month_year_full_match + [str(years_interest[i]) + '-' + str(months_interest[j]).zfill(2)]

start_month = int(month_year_data[0][5:7])
index_start = np.where(months_interest == start_month)[0][0]
end_month = int(month_year_data[-1][5:7])
index_end = 12 - np.where(months_interest == end_month)[0][0]
month_year_full = month_year_full_match[index_start:np.size(month_year_full_match)-index_end+1]
month_year_full = np.array(month_year_full)

missing_times = np.sort(np.array(list(set(month_year_full) - set(month_year_data))))
n_missing_times = np.size(missing_times)

# Generating the required time object in proper format for the missing times
missing_time_value = np.zeros(n_missing_times)
missing_time_bounds = np.zeros((n_missing_times,2))
for i in range(0,n_missing_times):
    year_cur = int(missing_times[i][0:4])
    month_cur = int(missing_times[i][5:7])
    
    missing_time_value[i] = nc.date2num(datetime(year_cur,month_cur,15,0),time_units)
    missing_time_bounds[i,0] = nc.date2num(datetime(year_cur,month_cur,1,0),time_units)
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

# Creating the UTC string and year-month array
time_continuous_utc = nc.num2date(time_continuous,time_units)
time_cont_utc_str = [str(time_continuous_utc[0].year) + '-' + str(time_continuous_utc[0].month).zfill(2) + '-' + str(time_continuous_utc[0].day).zfill(2) + 'T' + str(time_continuous_utc[0].hour).zfill(2) + ':' + str(time_continuous_utc[0].minute).zfill(2) + ':' + str(time_continuous_utc[0].second).zfill(2) + '.000']
year_month = [time_cont_utc_str[0][0:7]]
for i in range(1,n_times_continuous):
    time_utc_cur = [str(time_continuous_utc[i].year) + '-' + str(time_continuous_utc[i].month).zfill(2) + '-' + str(time_continuous_utc[i].day).zfill(2) + 'T' + str(time_continuous_utc[i].hour).zfill(2) + ':' + str(time_continuous_utc[i].minute).zfill(2) + ':' + str(time_continuous_utc[i].second).zfill(2) + '.000']
    time_cont_utc_str = np.concatenate((time_cont_utc_str,time_utc_cur))
    year_month = np.concatenate((year_month,[time_cont_utc_str[i][0:7]]))

# Writing to netCDF file
flag = create_nc_file(tws_op_file,tws_interp,uncertainty_interp,lat_flipped,lon_flipped,time_continuous,lat_bounds_flipped,lon_bounds_flipped,time_bounds_continuous,time_units,time_cont_utc_str,year_month)
