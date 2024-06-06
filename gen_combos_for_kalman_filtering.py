# -*- coding: utf-8 -*-
"""
Script to generate the ET estimates at each time for Kalman Smoother

Input arguements
twsc_path: Path to the directory where TWSC files are stored. 
twsc_files: Array of TWSC files. These are generated after applying central difference scheme. Each element of this array corresponds to one combination.
precip_path: Path to the directory where Precipitation files are stored
precip_files: Array of precipitation files that are generated after applying moving window averaging. Each element of this array corresponds to one combination.
runoff_path: Path to the directory where Runoff files are stored.
runoff_files: Array of runoff files that are generated after applying moving window averaging. Each element of this array corresponds to one combination.
time_start_interest: datetime.datetime object of start time of analysis
time_end_interest: datetime.datetime object of end time of analysis
op_file: Filename with full path of the output NetCDF file to be generated that contains the different ET estimates at each time step.

Outputs
Single NetCDF file which contains ET estimates at each time step
"""
##### Importing modules
import netCDF4 as nc
import matplotlib.pyplot as plt
import numpy as np
import glob
import datetime as dt

##### Specifying inputs
# Providing information about datasets
twsc_path = 'D:/IISC_COURSES/GEODETIC_SIGNAL_PROCESSING/PROJECT/TIME_SERIES/TWSC/'
twsc_files = ['TWSC_CSR.nc','TWSC_GFZ.nc','TWSC_JPL.nc']

precip_path = 'D:/IISC_COURSES/GEODETIC_SIGNAL_PROCESSING/PROJECT/TIME_SERIES/PRECIPITATION/'
precip_files = ['PRECIP_COBRA.nc','PRECIP_ERA5.nc','PRECIP_GPCC.nc','PRECIP_IMERG.nc']

runoff_path = 'D:/IISC_COURSES/GEODETIC_SIGNAL_PROCESSING/PROJECT/TIME_SERIES/RUNOFF/'
runoff_files = ['RUNOFF_CLSM.nc','RUNOFF_GRUN.nc','RUNOFF_GRUN_CRUST.nc','RUNOFF_GRUN_era5.nc','RUNOFF_GRUN_merra2.nc','RUNOFF_GRUN_wfde5.nc','RUNOFF_NOAH.nc','RUNOFF_VIC.nc'] # GRUN considered

# Providing information about time of interest
time_start_interest = dt.datetime(2002,12,31,23,59,59)
time_end_interest = dt.datetime(2007,12,31,23,59,59)
op_file = 'D:/IISC_COURSES/GEODETIC_SIGNAL_PROCESSING/PROJECT/KALMAN_FILTERING/TRIAL4/ET_MEASUREMENT_COMBOS_JAN2003_DEC2007.nc'

##### Reading the information from the time series files and subsetting to time of interest
### TWSC files
n_twsc_files = np.size(twsc_files)

# Reading a sample TWSC file to get the information of dimensions
twsc_file_cur = glob.glob(twsc_path + twsc_files[0])[0]
ds_twsc = nc.Dataset(twsc_file_cur)
twsc_lat = ds_twsc['lat'][:]
twsc_lon = ds_twsc['lon'][:]
twsc_time = ds_twsc['time'][:]
twsc_time_units = ds_twsc['time'].units
twsc_mm = ds_twsc['twsc_mm'][:]
twsc_unc_mm = ds_twsc['uncertainty_mm'][:]

n_lat = np.size(twsc_lat)
n_lon = np.size(twsc_lon)

# Subsetting the data to time of interest
time_start_interest_num = nc.date2num(time_start_interest,twsc_time_units)
time_end_interest_num = nc.date2num(time_end_interest,twsc_time_units)

index_interest = np.where((twsc_time > time_start_interest_num) & (twsc_time < time_end_interest_num))[0]
n_months_interest = np.size(index_interest)

twsc_mm_interest = np.zeros((n_twsc_files,n_months_interest,n_lat,n_lon))
twsc_unc_mm_interest = np.zeros((n_twsc_files,n_months_interest,n_lat,n_lon))

twsc_mm_interest[0] = twsc_mm[index_interest]
twsc_unc_mm_interest[0] = twsc_unc_mm[index_interest]

# Looping to get information from the rest of the files
for i in range(1,n_twsc_files):
    twsc_file_cur = glob.glob(twsc_path + twsc_files[i])[0]
    ds_twsc = nc.Dataset(twsc_file_cur)
    twsc_time = ds_twsc['time'][:]
    twsc_time_units = ds_twsc['time'].units
    twsc_mm = ds_twsc['twsc_mm']
    twsc_unc_mm = ds_twsc['uncertainty_mm'][:]
    
    time_start_interest_num = nc.date2num(time_start_interest,twsc_time_units)
    time_end_interest_num = nc.date2num(time_end_interest,twsc_time_units)
    index_interest = np.where((twsc_time > time_start_interest_num) & (twsc_time < time_end_interest_num))[0]
    
    twsc_mm_interest[i] = twsc_mm[index_interest]
    twsc_unc_mm_interest[i] = twsc_unc_mm[index_interest]

### PRECIPITATION FILES
n_precip_files = np.size(precip_files)

# Reading a sample precipitation file to get the information of dimensions
precip_file_cur = glob.glob(precip_path + precip_files[0])[0]
ds_precip = nc.Dataset(precip_file_cur)
precip_lat = ds_precip['lat'][:]
precip_lon = ds_precip['lon'][:]
precip_time = ds_precip['time'][:]
precip_time_units = ds_precip['time'].units
precip_mm = ds_precip['precip_avg_mm'][:]
precip_unc_mm = ds_precip['uncertainty_mm'][:]

n_lat = np.size(precip_lat)
n_lon = np.size(precip_lon)

# Subsetting the data to time of interest
time_start_interest_num = nc.date2num(time_start_interest,precip_time_units)
time_end_interest_num = nc.date2num(time_end_interest,precip_time_units)

index_interest = np.where((precip_time > time_start_interest_num) & (precip_time < time_end_interest_num))[0]
n_months_interest = np.size(index_interest)

precip_mm_interest = np.zeros((n_precip_files,n_months_interest,n_lat,n_lon))
precip_unc_mm_interest = np.zeros((n_precip_files,n_months_interest,n_lat,n_lon))

precip_mm_interest[0] = precip_mm[index_interest]
precip_unc_mm_interest[0] = precip_unc_mm[index_interest]

# Looping to get information from the rest of the files
for i in range(1,n_precip_files):
    precip_file_cur = glob.glob(precip_path + precip_files[i])[0]
    ds_precip = nc.Dataset(precip_file_cur)
    precip_time = ds_precip['time'][:]
    precip_time_units = ds_precip['time'].units
    precip_mm = ds_precip['precip_avg_mm']
    precip_unc_mm = ds_precip['uncertainty_mm'][:]
    
    time_start_interest_num = nc.date2num(time_start_interest,precip_time_units)
    time_end_interest_num = nc.date2num(time_end_interest,precip_time_units)
    index_interest = np.where((precip_time > time_start_interest_num) & (precip_time < time_end_interest_num))[0]
    
    precip_mm_interest[i] = precip_mm[index_interest]
    precip_unc_mm_interest[i] = precip_unc_mm[index_interest]

### RUNOFF FILES
n_runoff_files = np.size(runoff_files)

# Reading a sample runoff file to get the information of dimensions
runoff_file_cur = glob.glob(runoff_path + runoff_files[0])[0]
ds_runoff = nc.Dataset(runoff_file_cur)
runoff_lat = ds_runoff['lat'][:]
runoff_lon = ds_runoff['lon'][:]
runoff_time = ds_runoff['time'][:]
runoff_time_units = ds_runoff['time'].units
runoff_mm = ds_runoff['runoff_avg_mm'][:]
runoff_unc_mm = ds_runoff['uncertainty_mm'][:]

n_lat = np.size(runoff_lat)
n_lon = np.size(runoff_lon)

# Subsetting the data to time of interest
time_start_interest_num = nc.date2num(time_start_interest,runoff_time_units)
time_end_interest_num = nc.date2num(time_end_interest,runoff_time_units)

index_interest = np.where((runoff_time > time_start_interest_num) & (runoff_time < time_end_interest_num))[0]
n_months_interest = np.size(index_interest)

runoff_mm_interest = np.zeros((n_runoff_files,n_months_interest,n_lat,n_lon))
runoff_unc_mm_interest = np.zeros((n_runoff_files,n_months_interest,n_lat,n_lon))

runoff_mm_interest[0] = runoff_mm[index_interest]
runoff_unc_mm_interest[0] = runoff_unc_mm[index_interest]

# Looping to get information from the rest of the files
for i in range(1,n_runoff_files):
    runoff_file_cur = glob.glob(runoff_path + runoff_files[i])[0]
    ds_runoff = nc.Dataset(runoff_file_cur)
    runoff_time = ds_runoff['time'][:]
    runoff_time_units = ds_runoff['time'].units
    runoff_mm = ds_runoff['runoff_avg_mm']
    runoff_unc_mm = ds_runoff['uncertainty_mm'][:]
        
    time_start_interest_num = nc.date2num(time_start_interest,runoff_time_units)
    time_end_interest_num = nc.date2num(time_end_interest,runoff_time_units)
    index_interest = np.where((runoff_time > time_start_interest_num) & (runoff_time < time_end_interest_num))[0]
    
    runoff_mm_interest[i] = runoff_mm[index_interest]
    runoff_unc_mm_interest[i] = runoff_unc_mm[index_interest]

#### Making different combinations for ET measurements
n_et_measurement = n_twsc_files*n_precip_files*n_runoff_files
combination_name = []
et_measurement = np.zeros((n_et_measurement,n_months_interest,n_lat,n_lon),dtype=np.float32)
et_unc = np.zeros(np.shape(et_measurement),dtype=np.float32)
flag = 0
for i in range(0,n_runoff_files):
    for j in range(0,n_precip_files):
        for k in range(0,n_twsc_files):
            combo_cur = runoff_files[i].split('.')[0] + '_' + precip_files[j].split('.')[0] + '_' + twsc_files[k].split('.')[0]
            combination_name = combination_name + [combo_cur]
            et_measurement[flag] = precip_mm_interest[j] - runoff_mm_interest[i] - twsc_mm_interest[k]
            et_unc[flag] = np.sqrt((precip_unc_mm_interest[j])**2 + (runoff_unc_mm_interest[i])**2 + (twsc_unc_mm_interest[k])**2)
            flag = flag + 1
combination_name = np.array(combination_name)

# Randomizing the order of the combinations
index_combos = np.arange(0,n_et_measurement,1)
np.random.shuffle(index_combos)
combination_name = combination_name[index_combos]
et_measurement = et_measurement[index_combos]
et_unc = et_unc[index_combos]

#### Saving the combinations in the netCDF file
# Creating NetCDF file
ncfile = nc.Dataset(op_file,mode='w')

# Defining the dimensions in the NetCDF file
lat_dim = ncfile.createDimension('lat', n_lat)
lon_dim = ncfile.createDimension('lon', n_lon)
time_dim = ncfile.createDimension('time', n_months_interest)
combo_dim = ncfile.createDimension('combo',n_et_measurement)

# Defining the variables to be stored in the netCDF file along with their long names and appropriate units
lat = ncfile.createVariable('lat', np.float32, ('lat',))
lat.units = 'degrees_north'
lat.long_name = 'latitude'
lon = ncfile.createVariable('lon', np.float32, ('lon',))
lon.units = 'degrees_east'
lon.long_name = 'longitude'
time = ncfile.createVariable('time', np.float64, ('time',))
time.units = runoff_time_units
time.long_name = 'time'
combo = ncfile.createVariable('combo', np.int16, ('combo',))
et = ncfile.createVariable('et_mm',np.float64,('combo','time','lat','lon'))
et.units = 'mm'
et.long_name = 'Measurement of ET values'
unc = ncfile.createVariable('uncertainty_mm',np.float64,('combo','time','lat','lon'))
unc.units = 'mm'
unc.long_name = 'Uncertainities in Measured ET values'
combo_name = ncfile.createVariable('combo_name', str, ('combo',))
combo_name.long_name = 'Name of the combinations'

# Assigning values to the variables
lat[:] = runoff_lat
lon[:] = runoff_lon
time[:] = runoff_time[index_interest]
combo[:] = np.arange(0,n_et_measurement)
et[:] = et_measurement
unc[:] = et_unc
combo_name[:] = np.array(combination_name, dtype='object')

# Closing the netCDF file
ncfile.close()