# -*- coding: utf-8 -*-
"""
Created on Tue Oct 18 15:09:26 2022

@author: veenstra
"""

import os
import subprocess
import pandas as pd
from pathlib import Path
import requests


def download_ERA5(varkey,
                  date_min, date_max,
                  longitude_min, longitude_max, latitude_min, latitude_max, 
                  dir_output='.', overwrite=False):
    
    #TODO: describe something about the .cdsapirc file
    import cdsapi # Import cdsapi and create a Client instance # https://cds.climate.copernicus.eu/api-how-to #TODO: move to top of script? (then make dependency of dfm_tools)
    c = cdsapi.Client()
    
    #dictionary with ERA5 variables
    variables_dict = {'ssr':'surface_net_solar_radiation',
                      'sst':'sea_surface_temperature',
                      'strd':'surface_thermal_radiation_downwards',
                      'slhf':'surface_latent_heat_flux',
                      'sshf':'surface_sensible_heat_flux',
                      'str':'surface_net_thermal_radiation',
                      'chnk':'charnock',
                      'd2m':'2m_dewpoint_temperature',
                      't2m':'2m_temperature',
                      'tcc':'total_cloud_cover',
                      'mslp':'mean_sea_level_pressure',
                      'u10':'10m_u_component_of_wind',
                      'u10n':'10m_u_component_of_neutral_wind',
                      'v10':'10m_v_component_of_wind',
                      'v10n':'10m_v_component_of_neutral_wind',
                      'mer':'mean_evaporation_rate',
                      'mtpr':'mean_total_precipitation_rate',
                      }
    if varkey not in variables_dict.keys():
        raise Exception(f'"{varkey}" not available, choose from: {list(variables_dict.keys())}')
    
    period_range = pd.period_range(date_min,date_max,freq='M')
    print(f'retrieving data from {period_range[0]} to {period_range[-1]} (freq={period_range.freq})')
    
    for date in period_range:
        name_output = f'era5_{varkey}_{date.strftime("%Y-%m")}.nc'
        file_out = Path(dir_output,name_output)
        if file_out.is_file() and not overwrite:
            print(f'"{name_output}" found and overwrite=False, continuing.')
            continue
        print (f'retrieving ERA5 data for variable "{varkey}" and month {date.strftime("%Y-%m")} (YYYY-MM)')

        
        request_dict = {'product_type':'reanalysis',
                        'variable':variables_dict[varkey],
                        'year': date.strftime('%Y'),
                        'month':date.strftime('%m'),
                        #'month':[f'{x:02d}' for x in range(1,12+1)], #all months, but instead retrieving per month
                        'day':[f'{x:02d}' for x in range(1,31+1)], #all days
                        'time':[f'{x:02d}:00' for x in range(0,23+1)], #all times/hours
                        'area':[latitude_max,longitude_min,latitude_min,longitude_max], # north, west, south, east. default: global - option not available through the Climate Data Store (CDS) web interface
                        #'grid': [1.0, 1.0], # latitude/longitude grid: east-west (longitude) and north-south resolution (latitude). default: 0.25 x 0.25 - option not available through the Climate Data Store (CDS) web interface
                        'format':'netcdf'}
        
        c.retrieve(name='reanalysis-era5-single-levels', request=request_dict, target=file_out)
    return


def download_CMEMS(varkey,
                   date_min, date_max,
                   longitude_min, longitude_max, latitude_min, latitude_max, 
                   dir_output='.', overwrite=False,
                   source_combination=None, #str from source_dict.keys(). Or provide motu_url/service/product arguments
                   motu_url=None, service=None, product=None, #or provide source_combination argument
                   credentials=None, #credentials=['username','password'], or create "%USER%/motucredentials.txt" with username on line 1 and password on line 2. Register at: https://resources.marine.copernicus.eu/registration-form'
                   timeout=30): #in seconds #TODO: set timeout back to 300?
    
    """
    How to get motu_url/service/product:
        - find your service+product on https://resources.marine.copernicus.eu/products
        - go to the data access tab, e.g.: https://resources.marine.copernicus.eu/product-detail/GLOBAL_MULTIYEAR_PHY_001_030/DATA-ACCESS
        - click on the API button in the SUBS column, e.g.: https://my.cmems-du.eu/motu-web/Motu?action=describeproduct&service=GLOBAL_MULTIYEAR_PHY_001_030-TDS&product=cmems_mod_glo_phy_my_0.083_P1D-m
        - from this url you see the motu_url (https://my.cmems-du.eu), the service (GLOBAL_MULTIYEAR_PHY_001_030-TDS) and the product (cmems_mod_glo_phy_my_0.083_P1D-m)
        - some example combinations are available in source_dict
    Some examples can be found in source_dict
    """
    
    import motuclient #used in motu_commands, so has to be importable. conda install -c conda-forge motuclient #TODO: move to top of script (then make dependency of dfm_tools)
    #TODO: consider opendap instead?
    
    #get credentials
    if credentials is None:
        file_credentials = f'{os.path.expanduser("~")}/motucredentials.txt'
        if not os.path.exists(file_credentials):
            raise Exception(f'credentials argument not supplied and file_credentials not available ({file_credentials})')
        with open(file_credentials) as fc:
            username = fc.readline().strip()
            password = fc.readline().strip()
    else:
        username,password = credentials
       
    source_dict =   {'multiyear_physchem':{'motu_url':'https://my.cmems-du.eu', # multiyear reanalysis data (01-01-1993 12:00 till 31-05-2020 12:00)
                                           'service': 'GLOBAL_MULTIYEAR_PHY_001_030-TDS',
                                           'product': 'cmems_mod_glo_phy_my_0.083_P1D-m'},
                     'multiyear_bio':     {'motu_url':'https://my.cmems-du.eu',
                                           'service': 'GLOBAL_MULTIYEAR_BGC_001_029-TDS',
                                           'product': 'cmems_mod_glo_bgc_my_0.25_P1D-m'},
                     'forecast_physchem': {'motu_url':'https://nrt.cmems-du.eu', # operational forecast data (01-01-2019 12:00 till now + several days)
                                           'service': 'GLOBAL_ANALYSIS_FORECAST_PHY_001_024-TDS',
                                           'product': 'global-analysis-forecast-phy-001-024'},
                     'forecast_bio':      {'motu_url':'https://nrt.cmems-du.eu',
                                           'service': 'GLOBAL_ANALYSIS_FORECAST_BIO_001_028-TDS',
                                           'product': 'global-analysis-forecast-bio-001-028-daily'},
                     }
    
    if source_combination is not None:
        if (motu_url is not None) or (service is not None) or (product is not None):
            raise Exception('motu_url, service and/or product arguments provided while source_combination is also provided.')
        if not source_combination in source_dict.keys():
            raise Exception(f'provided source_combination={source_combination} argument is not valid, options are {list(source_dict.keys())}. Alternatively provide motu_url/service/product arguments.')
        motu_url = source_dict[source_combination]['motu_url']
        service = source_dict[source_combination]['service']
        product = source_dict[source_combination]['product']
    
    #test if supplied motu_url is valid
    requests.get(motu_url)
    
    period_range = pd.period_range(date_min,date_max,freq='D')
    print(f'retrieving data from {period_range[0]} to {period_range[-1]} (freq={period_range.freq})')
    
    for date in period_range:
        date_str = date.strftime('%Y-%m-%d')
        name_output = f'cmems_{varkey}_{date_str}.nc'
        file_out = Path(dir_output,name_output)
        if file_out.is_file() and not overwrite:
            print(f'"{name_output}" found and overwrite=False, continuing.')
            continue
        print(f'retrieving variable {varkey} for {date_str}: {name_output}')
        #TODO: alternatively use opendap, like half way this webpage: https://help.marine.copernicus.eu/en/articles/5182598-how-to-consume-the-opendap-api-and-cas-sso-using-python
        motu_command = ' '.join(['motuclient', '--motu', f'{motu_url}/motu-web/Motu', '--service-id', service, '--product-id', product,
                                 '--longitude-min', str(longitude_min), '--longitude-max', str(longitude_max),
                                 '--latitude-min', str(latitude_min), '--latitude-max', str(latitude_max),
                                 '--date-min', date_str, '--date-max', date_str, #+' 12:00:00',
                                 '--depth-min', '0', '--depth-max', '2e31',
                                 '--variable', str(varkey),
                                 '--out-dir', dir_output, '--out-name', name_output,
                                 '--user', username, '--pwd', password])
        
        try: #this try/except is necessary to capture/raise the motuclient logging in the 'finally' of the loop. Otherwise user only gets unclear error message: "CalledProcessError: Command [] returned non-zero exit status 1." 
            out = subprocess.run(f'python -m {motu_command}', capture_output=True, check=True, universal_newlines=True, timeout=timeout)
        except Exception as e: # capture error and save as out for the finally loop
            out = e 
            raise Exception(e)
        finally: #to capture also ERRORS/WARNINGS in logging that are not raised by motuclient
            if ('ERROR' in out.stdout) or ('WARNING' in out.stdout): #catch all other errors, and the relevant information in TimeoutExpired and CalledProcessError
                if 'variable not found' in out.stdout:
                    raise Exception(f'ERROR/WARNING found in motuclient logging:\nOUT: {out.stdout}\nERR: {out.stderr}\n\nNetCDF variable {varkey} not found: check available variables at: {motu_url}/motu-web/Motu?action=describeproduct&service={service}&product={product}')
                else:
                    raise Exception(f'ERROR/WARNING found in motuclient logging:\nOUT: {out.stdout}\nERR: {out.stderr}')
            #else:
            #    print(f'OUT: {out.stdout}\nERR: {out.stderr}')

    print('done')
    return