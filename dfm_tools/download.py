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
from dfm_tools.errors import OutOfRangeError
import cdsapi
import cftime
import getpass


def download_ERA5(varkey,
                  longitude_min, longitude_max, latitude_min, latitude_max, 
                  date_min, date_max,
                  dir_output='.', overwrite=False):
    """
    empty docstring
    """
    
    #TODO: describe something about the .cdsapirc file
    #TODO: make this function cdsapi generic, instead of ERA5 hardcoded (make flexible for product_type/name/name_output) (variables_dict is not used actively anymore, so this is possible)
    
    # create $HOME/.cdsapirc if it does not exist
    cds_credentials()
    
    c = cdsapi.Client() # import cdsapi and create a Client instance # https://cds.climate.copernicus.eu/api-how-to
    
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
                      'p140209':'air_density_over_the_oceans', # TODO: paramID might be replaced with shortname rhoao: https://jira.ecmwf.int/plugins/servlet/desk/portal/4/SD-82050
                      }
    if varkey not in variables_dict.keys(): #TODO: how to get list of available vars? mean_sea_level_pressure and msl both return a dataset with msl varkey, but standard_name air_pressure_at_mean_sea_level returns an error
        raise KeyError(f'"{varkey}" not available, choose from: {list(variables_dict.keys())}')
    
    period_range = pd.period_range(date_min,date_max,freq='M')
    print(f'retrieving data from {period_range[0]} to {period_range[-1]} (freq={period_range.freq})')
    
    #make sure the data fully covers the desired spatial extent. Download 1 additional grid cell (resolution is 1/4 degrees) in both directions
    longitude_min -= 1/4
    longitude_max += 1/4
    latitude_min  -= 1/4
    latitude_max  += 1/4
    
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


def download_CMEMS(varkey,
                   longitude_min, longitude_max, latitude_min, latitude_max, 
                   date_min, date_max, freq='D',
                   dir_output='.', file_prefix='', overwrite=False):
    """
    empty docstring
    """

    date_min, date_max = round_timestamp_to_outer_noon(date_min,date_max)
    
    global product #set product as global variable, so it only has to be retreived once per download run (otherwise once per variable)
    if 'product' not in globals():
        print('retrieving time range of CMEMS reanalysis and forecast products') #assuming here that physchem and bio reanalyisus/multiyear datasets have the same enddate, this seems safe
        reanalysis_tstart, reanalysis_tstop = get_OPeNDAP_xr_ds_timerange(dataset_url='https://my.cmems-du.eu/thredds/dodsC/cmems_mod_glo_phy_my_0.083_P1D-m')
        forecast_tstart, forecast_tstop = get_OPeNDAP_xr_ds_timerange(dataset_url='https://nrt.cmems-du.eu/thredds/dodsC/cmems_mod_glo_phy_anfc_0.083deg_P1D-m')
        if (date_min >= reanalysis_tstart) & (date_max <= reanalysis_tstop):
            product = 'reanalysis'
            print(f"The CMEMS '{product}' product will be used.")
        elif (date_min >= forecast_tstart) & (date_max <= forecast_tstop):
            product = 'analysisforecast'
            print(f"The CMEMS '{product}' product will be used.")
        else:
            raise ValueError(f'Requested timerange ({date_min} to {date_max}) is not fully within timerange of reanalysis product ({reanalysis_tstart} to {reanalysis_tstop}) or forecast product ({forecast_tstart} to {forecast_tstop}).')
    
    Path(dir_output).mkdir(parents=True, exist_ok=True)
    if varkey in ['bottomT','tob','mlotst','siconc','sithick','so','thetao','uo','vo','usi','vsi','zos']: #for physchem
        buffer = 2/12 # resolution is 1/12 degrees
        if product == 'analysisforecast': #forecast: https://data.marine.copernicus.eu/product/GLOBAL_ANALYSISFORECAST_PHY_001_024/description
            if varkey in ['uo','vo']: #anfc datset is splitted over multiple urls, construct the correct one here.
                varkey_name = 'phy-cur'
            elif varkey in ['so','thetao']:
                varkey_name = 'phy-'+varkey
            else:
                varkey_name = 'phy'
            dataset_url = f'https://nrt.cmems-du.eu/thredds/dodsC/cmems_mod_glo_{varkey_name}_anfc_0.083deg_P1D-m'#.html' #TODO: also PT6H-i timeresolution available but not for all variables and not for reanalysis
        else: #reanalysis: https://data.marine.copernicus.eu/product/GLOBAL_MULTIYEAR_PHY_001_030/description
            dataset_url = 'https://my.cmems-du.eu/thredds/dodsC/cmems_mod_glo_phy_my_0.083_P1D-m'
    else: #for bio
        buffer = 2/4 # resolution is 1/4 degrees
        if product == 'analysisforecast': #forecast: https://data.marine.copernicus.eu/product/GLOBAL_ANALYSIS_FORECAST_BIO_001_028/description
            dataset_url = 'https://nrt.cmems-du.eu/thredds/dodsC/global-analysis-forecast-bio-001-028-daily' #contains ['chl','fe','no3','nppv','o2','ph','phyc','po4','si','spco2']
        else: #https://data.marine.copernicus.eu/product/GLOBAL_MULTIYEAR_BGC_001_029/description
            dataset_url = 'https://my.cmems-du.eu/thredds/dodsC/cmems_mod_glo_bgc_my_0.25_P1D-m' #contains ['chl','no3','nppv','o2','po4','si']
    
    #make sure the data fully covers the desired spatial extent. Download 2 additional grid cells (resolution is 1/12 degrees, but a bit more/less in alternating cells) in each direction
    longitude_min -= buffer
    longitude_max += buffer
    latitude_min  -= buffer
    latitude_max  += buffer
    
    download_OPeNDAP(dataset_url=dataset_url,
                     varkey=varkey,
                     longitude_min=longitude_min, longitude_max=longitude_max, latitude_min=latitude_min, latitude_max=latitude_max,
                     date_min=date_min, date_max=date_max,
                     dir_output=dir_output, file_prefix=file_prefix, overwrite=overwrite)


def cds_credentials():
    """
    create $HOME/.cdsapirc via getpass if necessary
    """
    #TODO: put this in a PR at https://github.com/ecmwf/cdsapi (https://github.com/ecmwf/cdsapi/blob/master/cdsapi/api.py#L303)
    file_credentials = f'{os.path.expanduser("~")}/.cdsapirc'
    if os.path.exists(file_credentials):
        print('found CDS apikey')
    else:
        print("Downloading CDS/ERA5 data requires a CDS API key, copy the key from https://cds.climate.copernicus.eu/api-how-to (first register and sign in) ")
        apikey = getpass.getpass("\nEnter your CDS API Key: ")
        if ":" not in apikey:
            raise ValueError('A ":"-symbol is missing in the provided key (format should be "UID:apikey"). Please check the instructions above again.')
        with open(file_credentials,'w') as fc:
            fc.write('url: https://cds.climate.copernicus.eu/api/v2\n')
            fc.write(f'key: {apikey}')


def copernicusmarine_credentials():
    """
    get CMEMS username/password from file or via getpass
    """
    file_credentials = f'{os.path.expanduser("~")}/CMEMS_credentials.txt'
    if os.path.exists(file_credentials):
        print('found CMEMS credentials')
        with open(file_credentials) as fc:
            username = fc.readline().strip()
            password = fc.readline().strip()
    else: #query username and password with getpass
        print("Downloading CMEMS data requires a Copernicus Marine username and password, sign up for free at: https://data.marine.copernicus.eu/register.")
        username = getpass.getpass("Enter your Copernicus Marine username: ")
        password = getpass.getpass("Enter your Copernicus Marine password: ")
        userpass_save = input("Do you want to save your credentials as plain text? [y/n]: ")
        if userpass_save == 'y':
            with open(file_credentials,'w') as fc:
                fc.write(f'{username}\n{password}\n')
    return username, password


def copernicusmarine_datastore(dataset_url):
    """
    Setting up a copernicus marine PydapDataStore with authentication via username and password.
    """
    if 'session' not in globals():
        username, password = copernicusmarine_credentials()
        #setting up a session and making the variable global so we do not have to repeat it
        #https://help.marine.copernicus.eu/en/articles/5182598-how-to-consume-the-opendap-api-and-cas-sso-using-python
        cas_url = 'https://cmems-cas.cls.fr/cas/login'
        global session
        session = setup_session(cas_url, username, password)
    
    cookies_dict = session.cookies.get_dict()
    if 'CASTGC' not in cookies_dict.keys():
        raise KeyError('CASTGC key missing from session cookies_dict, probably authentication failure')
    session.cookies.set("CASTGC", cookies_dict['CASTGC'])
    #TODO: add check for wrong dataset_id (now always "AttributeError: You cannot set the charset when no content-type is defined")
    dap_dataset = open_url(dataset_url, session=session, user_charset='utf-8')
    data_store = xr.backends.PydapDataStore(dap_dataset)
    return data_store


def open_OPeNDAP_xr(dataset_url):
    """
    How to get the opendap dataset_url (CMEMS example):
        - https://data.marine.copernicus.eu/products
        - go to the data access tab of a product, e.g.: https://data.marine.copernicus.eu/product/GLOBAL_MULTIYEAR_PHY_001_030/services
        - click the opendap link of the dataset of your choice
        - copy the dataset_url from the adress bar (excl .html), e.g.: https://my.cmems-du.eu/thredds/dodsC/cmems_mod_glo_phy_my_0.083_P1D-m
    
    How to get the opendap dataset_url (HYCOM example):
        - https://www.hycom.org/dataserver
        - Select a product and search for THREDDS, e.g.: https://www.hycom.org/dataserver/gofs-3pt1/analysis
        - find an opendap dataset_url, it depends per product/run where to find it.
        Some examples:
            https://tds.hycom.org/thredds/dodsC/GLBu0.08/expt_19.1/2010
            https://tds.hycom.org/thredds/dodsC/GLBy0.08/expt_93.0

    """
        
    if isinstance(dataset_url,list):
        dataset_url_one = dataset_url[0]
    else:
        dataset_url_one = dataset_url
    
    if 'cmems-du.eu' in dataset_url_one:
        if isinstance(dataset_url,list):
            raise TypeError('list not supported by opendap method used for cmems')
        
        print(f'opening pydap connection to opendap dataset and opening with xarray: {dataset_url}.html')
        data_store = copernicusmarine_datastore(dataset_url=dataset_url)
        data_xr = xr.open_dataset(data_store)
    elif 'hycom.org' in dataset_url_one:
        if isinstance(dataset_url,list):
            print(f'xarray opening opendap dataset like: {dataset_url[0]}.html ({len(dataset_url)} urls/years)')
            data_xr = xr.open_mfdataset(dataset_url,decode_times=False) #TODO: for some reason decode_times does not work: "ValueError: unable to decode time units 'hours since analysis' with 'the default calendar'."
        else:
            print(f'xarray opening opendap dataset: {dataset_url}.html')
            data_xr = xr.open_dataset(dataset_url,decode_times=False) #TODO: for some reason decode_times does not work: "ValueError: unable to decode time units 'hours since analysis' with 'the default calendar
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


def download_OpenDAP_gettimes(dataset_url):
    ds = open_OPeNDAP_xr(dataset_url=dataset_url)
    ds_time = ds.time.to_series()
    return ds_time


def download_OPeNDAP(dataset_url,
                     varkey,
                     longitude_min, longitude_max, latitude_min, latitude_max, 
                     date_min, date_max, freq='D',
                     dir_output='.', file_prefix='', overwrite=False):
    """
    

    Parameters
    ----------
    dataset_url : TYPE
        DESCRIPTION.
    varkey : TYPE
        DESCRIPTION.
    longitude_min : TYPE
        DESCRIPTION.
    longitude_max : TYPE
        DESCRIPTION.
    latitude_min : TYPE
        DESCRIPTION.
    latitude_max : TYPE
        DESCRIPTION.
    date_min : TYPE
        DESCRIPTION.
    date_max : TYPE
        DESCRIPTION.
    freq : TYPE, optional
        DESCRIPTION. The default is 'D'.
    dir_output : TYPE, optional
        DESCRIPTION. The default is '.'.
    file_prefix : TYPE, optional
        DESCRIPTION. The default is ''.
    overwrite : TYPE, optional
        DESCRIPTION. The default is False.
    
    Raises
    ------
    KeyError
        DESCRIPTION.
    OutOfRangeError
        DESCRIPTION.

    Returns
    -------
    None.

    """
    
    data_xr = open_OPeNDAP_xr(dataset_url=dataset_url)
    
    print(f'xarray subsetting data (variable \'{varkey}\' and lon/lat extents)')
    if varkey not in data_xr.data_vars:
        raise KeyError(f'"{varkey}" not found in dataset, available are: {list(data_xr.data_vars)}')
    data_xr_var = data_xr[[varkey]]
    data_xr_var = data_xr_var.sel(longitude=slice(longitude_min,longitude_max), #TODO: add depth selection?
                                  latitude=slice(latitude_min,latitude_max))
    data_xr_times = data_xr_var.time.to_series()
    print(f'available time range in dataset from {data_xr_times.index[0]} to {data_xr_times.index[-1]}')
    period_range = pd.period_range(date_min,date_max,freq=freq)
    
    #check if date_min/date_max are available in dataset
    if not (data_xr_times.index[0] <= period_range[0].to_timestamp() <= data_xr_times.index[-1]):
        raise OutOfRangeError(f'date_min ({period_range[0]}) is outside available time range in dataset: {data_xr_times.index[0]} to {data_xr_times.index[-1]}')
    if not (data_xr_times.index[0] <= period_range[-1].to_timestamp() <= data_xr_times.index[-1]):
        raise OutOfRangeError(f'date_max ({period_range[-1]}) is outside available time range in dataset: {data_xr_times.index[0]} to {data_xr_times.index[-1]}')
    
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


def get_OPeNDAP_xr_ds_timerange(dataset_url):
    ds = open_OPeNDAP_xr(dataset_url=dataset_url)
    ds_times = ds.time.to_series()
    ds_tstart, ds_tstop = ds_times.iloc[0], ds_times.iloc[-1]
    return ds_tstart, ds_tstop


def round_timestamp_to_outer_noon(date_min, date_max):
    """
    Since the CMEMS dataset only contains noon-values, we often need to round to the previous or next noon timestep to download enough data.
    
    """
    
    td_12h = pd.Timedelta(hours=12)
    date_min = (pd.Timestamp(date_min) + td_12h).floor('1d') - td_12h
    date_max = (pd.Timestamp(date_max) - td_12h).ceil('1d') + td_12h
    return date_min, date_max

