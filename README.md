# KF-ET
A novel and meaningful ET using Kalman Filter and Water Budget

General comments

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
