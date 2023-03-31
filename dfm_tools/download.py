# -*- coding: utf-8 -*-
"""
Created on Tue Oct 18 15:09:26 2022

@author: veenstra
"""

import os
import pandas as pd
from pathlib import Path
import xarray as xr
from pydap.client import open_url
from pydap.cas.get_cookies import setup_session
import warnings


def download_ERA5(varkey,
                  longitude_min, longitude_max, latitude_min, latitude_max, 
                  date_min, date_max,
                  dir_output='.', overwrite=False):
    
    #TODO: describe something about the .cdsapirc file
    #TODO: make this function cdsapi generic, instead of ERA5 hardcoded (make flexible for product_type/name/name_output) (variables_dict is not used actively anymore, so this is possible)
    
    import cdsapi # Import cdsapi and create a Client instance # https://cds.climate.copernicus.eu/api-how-to #TODO: move to top of script? (then make dependency of dfm_tools)
    c = cdsapi.Client()
    
    #dictionary with ERA5 variables #this is not actively used
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
                      'msl':'mean_sea_level_pressure',
                      'u10':'10m_u_component_of_wind',
                      'u10n':'10m_u_component_of_neutral_wind',
                      'v10':'10m_v_component_of_wind',
                      'v10n':'10m_v_component_of_neutral_wind',
                      'mer':'mean_evaporation_rate',
                      'mtpr':'mean_total_precipitation_rate',
                      }
    if varkey not in variables_dict.keys(): #TODO: how to get list of available vars? mean_sea_level_pressure and msl both return a dataset with msl varkey, but standard_name air_pressure_at_mean_sea_level returns an error
        raise KeyError(f'"{varkey}" not available, choose from: {list(variables_dict.keys())}')
    
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
                        'year':date.strftime('%Y'),
                        'month':date.strftime('%m'),
                        #'month':[f'{x:02d}' for x in range(1,12+1)], #all months, but instead retrieving per month
                        'day':[f'{x:02d}' for x in range(1,31+1)], #all days
                        'time':[f'{x:02d}:00' for x in range(0,23+1)], #all times/hours
                        'area':[latitude_max,longitude_min,latitude_min,longitude_max], # north, west, south, east. default: global - option not available through the Climate Data Store (CDS) web interface (for cmip data)
                        #'grid': [1.0, 1.0], # latitude/longitude grid: east-west (longitude) and north-south resolution (latitude). default: 0.25 x 0.25 - option not available through the Climate Data Store (CDS) web interface
                        'format':'netcdf'}
        
        c.retrieve(name='reanalysis-era5-single-levels', request=request_dict, target=file_out)
    return


def open_OPeNDAP_xr(dataset_url,
                    credentials=None): #for CMEMS. credentials=['username','password'], or create "%USERPROFILE%/CMEMS_credentials.txt" with username on line 1 and password on line 2. Register at: https://resources.marine.copernicus.eu/registration-form'
    
    """
    How to get the opendap dataset_url (CMEMS example):
        - https://data.marine.copernicus.eu/products
        - go to the data access tab of a product, e.g.: https://data.marine.copernicus.eu/product/GLOBAL_MULTIYEAR_PHY_001_030/services
        - click the opendap link of the dataset of your choice
        - copy the dataset_url from the adress bar (excl .html), e.g.: https://my.cmems-du.eu/thredds/dodsC/cmems_mod_glo_phy_my_0.083_P1D-m
        
        Example multiyear phys/chem:
            https://data.marine.copernicus.eu/product/GLOBAL_MULTIYEAR_PHY_001_030/services
            https://my.cmems-du.eu/thredds/dodsC/cmems_mod_glo_phy_my_0.083_P1D-m
            ['bottomT','mlotst','siconc','sithick','so','thetao','uo','usi','vo','vsi','zos']
        Example multiyear bio:
            https://data.marine.copernicus.eu/product/GLOBAL_MULTIYEAR_BGC_001_029/services
            https://my.cmems-du.eu/thredds/dodsC/cmems_mod_glo_bgc_my_0.25_P1D-m
            ['chl','no3','nppv','o2','po4','si']
        Example forecast phys/chem:
            https://data.marine.copernicus.eu/product/GLOBAL_ANALYSIS_FORECAST_PHY_001_024/services
            https://nrt.cmems-du.eu/thredds/dodsC/global-analysis-forecast-phy-001-024
            ['bottomT','mlotst','siconc','sithick','so','thetao','uo','usi','vo','vsi','zos']
        Example forecast bio:
            https://data.marine.copernicus.eu/product/GLOBAL_ANALYSIS_FORECAST_BIO_001_028/services
            https://nrt.cmems-du.eu/thredds/dodsC/global-analysis-forecast-bio-001-028-daily
            ['chl','fe','no3','nppv','o2','ph','phyc','po4','si','spco2']
    
    How to get the opendap dataset_url (HYCOM example):
        - https://www.hycom.org/dataserver
        - Select a product and search for THREDDS, e.g.: https://www.hycom.org/dataserver/gofs-3pt1/analysis
        - find an opendap dataset_url, it depends per product/run where to find it.
        Some examples:
            https://tds.hycom.org/thredds/dodsC/GLBu0.08/expt_19.1/2010
            https://tds.hycom.org/thredds/dodsC/GLBy0.08/expt_93.0
    """
    
    def copernicusmarine_datastore(dataset_url, username, password):
        #https://help.marine.copernicus.eu/en/articles/5182598-how-to-consume-the-opendap-api-and-cas-sso-using-python
        cas_url = 'https://cmems-cas.cls.fr/cas/login'
        session = setup_session(cas_url, username, password)
        cookies_dict = session.cookies.get_dict()
        if not 'CASTGC' in cookies_dict.keys():
            raise Exception('CASTGC key missing from session cookies_dict, probably authentication failure')
        session.cookies.set("CASTGC", cookies_dict['CASTGC'])
        #TODO: add check for wrong dataset_id (now always "AttributeError: You cannot set the charset when no content-type is defined")
        DAP_dataset = open_url(dataset_url, session=session)#, user_charset='utf-8') # TODO: user_charset needs PyDAP >= v3.3.0 see https://github.com/pydap/pydap/pull/223/commits 
        data_store = xr.backends.PydapDataStore(DAP_dataset)
        return data_store
    
    if isinstance(dataset_url,list):
        dataset_url_one = dataset_url[0]
    else:
        dataset_url_one = dataset_url
    
    if 'cmems-du.eu' in dataset_url_one:
        if isinstance(dataset_url,list):
            raise Exception('no list supported by opendap method used for cmems')
        
        #parce credentials to username/password #TODO: now CMEMS specific, make more generic
        if credentials is None:
            file_credentials = f'{os.path.expanduser("~")}/CMEMS_credentials.txt'
            if not os.path.exists(file_credentials):
                raise Exception(f'credentials argument not supplied and file_credentials not available ({file_credentials})')
            with open(file_credentials) as fc:
                username = fc.readline().strip()
                password = fc.readline().strip()
        else:
            username,password = credentials

        print(f'opening pydap connection to opendap dataset: {dataset_url}.html')
        data_store = copernicusmarine_datastore(dataset_url=dataset_url, username=username, password=password)
        print('xarray opening opendap dataset')
        data_xr = xr.open_dataset(data_store)
    elif 'hycom.org' in dataset_url_one:
        if isinstance(dataset_url,list):
            print(f'xarray opening opendap dataset like: {dataset_url[0]}.html ({len(dataset_url)} urls/years)')
            data_xr = xr.open_mfdataset(dataset_url,decode_times=False) #TODO: for some reason decode_times does not work: "ValueError: unable to decode time units 'hours since analysis' with 'the default calendar'."
        else:
            print(f'xarray opening opendap dataset: {dataset_url}.html')
            data_xr = xr.open_dataset(dataset_url,decode_times=False) #TODO: for some reason decode_times does not work: "ValueError: unable to decode time units 'hours since analysis' with 'the default calendar'."
        import cftime
        data_xr['time'] = cftime.num2date(data_xr.time,units=data_xr.time.units,calendar=data_xr.time.calendar)
        data_xr = data_xr.rename({'lon':'longitude','lat':'latitude'})
    else:
        print(f'unspecified dataset_url, might fail: {dataset_url}')
        if isinstance(dataset_url,list):
            print(f'xarray opening opendap dataset like: {dataset_url[0]}.html ({len(dataset_url)} urls/years)')
            data_xr = xr.open_mfdataset(dataset_url)
        else:
            print(f'xarray opening opendap dataset: {dataset_url}.html')
            data_xr = xr.open_dataset(dataset_url)
        if 'lon' in data_xr.dims:
            data_xr = data_xr.rename({'lon':'longitude','lat':'latitude'})
        
    return data_xr


def download_OPeNDAP(dataset_url,
                     varkey,
                     longitude_min, longitude_max, latitude_min, latitude_max, 
                     date_min, date_max, freq='D',
                     dir_output='.', file_prefix='', overwrite=False,
                     credentials=None):
    
    data_xr = open_OPeNDAP_xr(dataset_url=dataset_url, credentials=credentials)
    
    print(f'xarray subsetting data (variable \'{varkey}\' and lon/lat extents)')
    if varkey not in data_xr.data_vars:
        raise Exception(f'{varkey} not found in dataset, available are: {list(data_xr.data_vars)}')
    data_xr_var = data_xr[[varkey]]
    data_xr_var = data_xr_var.sel(longitude=slice(longitude_min,longitude_max), #TODO: add depth selection?
                                  latitude=slice(latitude_min,latitude_max))
    data_xr_times = data_xr_var.time.to_series()
    print(f'available time range in dataset from {data_xr_times.index[0]} to {data_xr_times.index[-1]}')
    period_range = pd.period_range(date_min,date_max,freq=freq)
    
    #check if date_min/date_max are available in dataset
    if not (data_xr_times.index[0] <= period_range[0].to_timestamp() <= data_xr_times.index[-1]):
        raise Exception(f'date_min ({period_range[0]}) is outside available time range in dataset: {data_xr_times.index[0]} to {data_xr_times.index[-1]}')
    if not (data_xr_times.index[0] <= period_range[-1].to_timestamp() <= data_xr_times.index[-1]):
        raise Exception(f'date_max ({period_range[-1]}) is outside available time range in dataset: {data_xr_times.index[0]} to {data_xr_times.index[-1]}')
    
    for date in period_range:
        date_str = str(date)
        name_output = f'{file_prefix}{varkey}_{date_str}.nc'
        file_out = Path(dir_output,name_output)
        if file_out.is_file() and not overwrite:
            print(f'"{name_output}" found and overwrite=False, continuing.')
            continue
        
        print(f'xarray subsetting data per {period_range.freq}: {date_str}')
        data_xr_var_seltime = data_xr_var.sel(time=slice(date_str,date_str)) #+' 12:00:00', 
        
        print(f'xarray writing netcdf file: {name_output}')
        data_xr_var_seltime.to_netcdf(os.path.join(dir_output,name_output)) #TODO: add chunks={'time':1} or only possible with opening?
        data_xr_var_seltime.close()
    
    print('done')
    return
