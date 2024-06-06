# -*- coding: utf-8 -*-
"""
Function that takes the interpolated mascon fields/landmass fields and uncertainties, and creates an output nc file with proper dimensions and variables

Function usage
flag = create_nc_file(ncfile_name,ewh_data,uncertainty,delta,n_values_lat,n_values_lon,utc_time,year_month)

Input arguments
ncfile: Name with full path of the NC file that needs to be created (dtype=string)
tws_data: 3D cube of TWS (in mm) obtained after interpolation (dtype=float, size=n_timebins X n_latitudes X n_longitudes)
uncertainty: 3D cube of uncertainities (in mm) appropriately interpolated in data gaps (dtype=float,size=same as tws_data)
lat_values: Latitude values in degrees north (dtype=float, size=n_latitudes)
lon_values: Longitude values in degrees east (dtype=float, size=n_longitudes)
time_values: Time values in appropriate units (dtype=float, size=n_timebins)
lat_bounds: Latitude boundary values in degrees north (dtype=float, size=n_latitudes X 2)
lon_bounds: Longitude boundary values in degrees east (dtype=float, size=n_longitudes X 2)
time_bounds: Time boundary values in appropriate units (dtype=float, size=n_timebins X 2)
time_units: Units in which time axis is defined (dtype=string)
utc_time: UTC time in YYYY-MM-DDTHH:mm:SS.zzz. When data is available, these are the values from mascons or landmass files (typically mid point of the month), and in datagaps they have been set to the mid point of the month (dtype=string, size=n_timebins)
year_month: String of year and month in YYYY-MM format (dtype=string,size=n_timebins)

Output arguments
flag: Just to return an output from the function
"""
def create_nc_file(ncfile_name,tws_data,uncertainty,lat_values,lon_values,time_values,lat_bounds_values,lon_bounds_values,time_bounds_values,time_units,utc_time,year_month):
    
    # Importing modules
    import numpy as np
    import netCDF4 as nc
    
    # Extracting size related info
    n_values_lat = np.size(lat_values)
    n_values_lon = np.size(lon_values)
    n_values_time = np.size(time_values)
        
    # Defining UTC time and year month strings in object datatype to write in NetCDF file
    utc_time_str = np.array(utc_time, dtype='object')
    year_month_str = np.array(year_month,dtype='object')

    # Creating NetCDF file
    ncfile = nc.Dataset(ncfile_name,mode='w')
    
    # Defining the dimensions in the NetCDF file
    lat_dim = ncfile.createDimension('lat', n_values_lat)
    lon_dim = ncfile.createDimension('lon', n_values_lon)
    time_dim = ncfile.createDimension('time', n_values_time)
    bounds_dim = ncfile.createDimension('bounds', 2)

    # Defining the variables to be stored in the NetCDF file along with their long names and appropriate units
    lat = ncfile.createVariable('lat', np.float32, ('lat',))
    lat.units = 'degrees_north'
    lat.long_name = 'latitude'
    lon = ncfile.createVariable('lon', np.float32, ('lon',))
    lon.units = 'degrees_east'
    lon.long_name = 'longitude'
    time = ncfile.createVariable('time', np.float64, ('time',))
    time.units = time_units
    time.long_name = 'time'
    
    lat_bounds = ncfile.createVariable('lat_bounds', np.float32, ('lat','bounds'))
    lat_bounds.units = 'degrees_north'
    lat_bounds.long_name = 'latitude boundaries'
    lon_bounds = ncfile.createVariable('lon_bounds', np.float32, ('lon','bounds'))
    lon_bounds.units = 'degrees_east'
    lon_bounds.long_name = 'longitude boundaries'
    time_bounds = ncfile.createVariable('time_bounds', np.float64, ('time','bounds'))
    time_bounds.units = time_units
    time_bounds.long_name = 'time boundaries'

    tws = ncfile.createVariable('tws_mm',np.float64,('time','lat','lon'))
    tws.units = 'mm'
    tws.long_name = 'Total Water Storage'
    unc = ncfile.createVariable('uncertainty_mm',np.float64,('time','lat','lon'))
    unc.units = 'mm'
    unc.long_name = 'Uncertainities in TWS'

    utc = ncfile.createVariable('UTC_time', str, ('time',))
    utc.long_name = 'UTC time'
    ym = ncfile.createVariable('yearmonth', str, ('time',))
    ym.long_name = 'Year-Month'
    
    lat[:] = lat_values
    lon[:] = lon_values
    time[:] = time_values
    
    lat_bounds[:] = lat_bounds_values
    lon_bounds[:] = lon_bounds_values
    time_bounds[:] = time_bounds_values
    
    tws[:,:,:] = tws_data
    unc[:,:,:] = uncertainty
      
    utc[:] = utc_time_str
    ym[:] = year_month_str

    ncfile.close()
    flag = 1
    return flag

