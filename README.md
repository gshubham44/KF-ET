# KF-ET
This package contains codes used for development of Kalman Filter based Global Gridded Evapotranspiration product using water budget approach (under revision in ERL)
The package contains following python scripts:
1.	create_nc_file.py: Function that takes the interpolated mascon fields/landmass fields and uncertainties and creates an output .nc file with proper dimensions and variables.
2.	gen_combos_for_kalman_filtering.py: Script to generate the ET estimates at each time for Kalman Filtering.
3.	jpl_mascons_process.py: Script to read JPL mascon NetCDF file, fill temporal data gaps with interpolation and generate output file of required format.
4.	landmass_process.py Script to read CSR/GFZ landmass data product and generate NetCDF files of TWS in the required format.
5.	run_kalman_filtering.py: Script to run Kalman Filter to generate the best ET estimate at each time step.

#General comments

1. Random order of combinations.

2. No uncertainty inflation for negative ET values

##### *withGRUN

Time period: Jan-2003 to Nov-2014

Datasets used:

TWS : JPL mascons, CSR and GFZ landmass products (3)

Precipitation : COBRA, ERA5, GPCC, IMERG (4)

Runoff : CLSM, GRUN, GRUN_CRUST, GRUN_era5, GRUN_merra2, GRUN_wfde5, NOAH, VIC (For CLSM, VIC and NOAH uncertainties are taken as 100% of the value) (8)

ET Initial : ERA5

Total no of combinations = 96 (8x4x3)

##### *withoutGRUN

Time period: Dec-2014 to Dec-20116

Datasets used:

TWS : JPL mascons, CSR and GFZ landmass products (3)

Precipitation : COBRA, ERA5, GPCC, IMERG (4)

Runoff : CLSM, GRUN_CRUST, GRUN_era5, GRUN_merra2, GRUN_wfde5, NOAH, VIC (For CLSM, VIC and NOAH uncertainties are taken as 100% of the value) (8)

ET Initial : ERA5

Total no of combinations = 84 (7x4x3)

Code Developer:
Netra S. Pillai
[Centre for Earth Science, IISc Bangalore, India &
Indian Space Research Organization (ISRO), Bangalore, India]
