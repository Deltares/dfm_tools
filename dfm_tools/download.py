import os
import pandas as pd
from pathlib import Path
import xarray as xr
from dfm_tools.errors import OutOfRangeError
import cdsapi
import copernicusmarine
import cftime
import getpass
import shutil
import subprocess
import sys

__all__ = [
    "download_ERA5",
    "download_CMEMS",
    "download_OPeNDAP",
]


def download_ERA5(varkey,
                  longitude_min, longitude_max, latitude_min, latitude_max, 
                  date_min, date_max,
                  dir_output='.', overwrite=False):
    """
    empty docstring
    """
    
    #TODO: describe something about the .cdsapirc file
    #TODO: make this function cdsapi generic, instead of ERA5 hardcoded (make flexible for product_type/name/name_output) (variables_dict is not used actively anymore, so this is possible)
    
    # create ~/.cdsapirc if it does not exist
    cds_credentials()
    
    # import cdsapi and create a Client instance # https://cds.climate.copernicus.eu/api-how-to
    c = cds_client_withargs()

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
    buffer = 2/4
    longitude_min -= buffer
    longitude_max += buffer
    latitude_min  -= buffer
    latitude_max  += buffer
    
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


def cds_get_file():
    file_cds_credentials = os.environ.get("CDSAPI_RC", os.path.expanduser("~/.cdsapirc"))
    return file_cds_credentials


def cds_credentials():
    """
    get cdsapikey from environment variables or file or query via getpass if necessary
    """
    #TODO: put this in a PR at https://github.com/ecmwf/cdsapi (https://github.com/ecmwf/cdsapi/blob/master/cdsapi/api.py#L303)
    cds_url = os.environ.get("CDSAPI_URL", "https://cds.climate.copernicus.eu/api/v2")
    # set default/provided CDSAPI_URL back to environ for platforms like EDITO 
    # that depend on environ (only CDSAPI_KEY has to be set in that case)
    os.environ["CDSAPI_URL"] = cds_url
    cds_uid_apikey = os.environ.get("CDSAPI_KEY")
    
    # read credentials from file if it exists. This has higher precedence over env vars
    file_cds_credentials = cds_get_file()
    if os.path.isfile(file_cds_credentials):
        config = cdsapi.api.read_config(file_cds_credentials)
        cds_url = config["url"]
        cds_uid_apikey = config["key"]
    
    try:
        # checks whether CDS apikey is in environment variable or ~/.cdsapirc file and if it is in correct format
        c = cds_client_withargs()
        # checks whether credentials (uid and apikey) are correct
        c.retrieve(name="dummy", request={})
    except Exception as e:
        if "Missing/incomplete configuration file" in str(e):
            # to catch "Exception: Missing/incomplete configuration file"
            # query uid and apikey if not present
            print("Downloading CDS/ERA5 data requires a CDS API key, copy your UID and API-key from https://cds.climate.copernicus.eu/user (first register, login and accept the terms). ")
            cds_uid = getpass.getpass("\nEnter your CDS UID (six digits): ")
            cds_apikey = getpass.getpass("\nEnter your CDS API-key (string with dashes): ")
            cds_uid_apikey = f"{cds_uid}:{cds_apikey}"
            os.environ["CDSAPI_URL"] = cds_url
            os.environ["CDSAPI_KEY"] = cds_uid_apikey
            with open(file_cds_credentials,'w') as fc:
                fc.write(f'url: {cds_url}\n')
                fc.write(f'key: {cds_uid_apikey}')
            cds_credentials()
        elif "not the correct format" in str(e):
            # to catch "AssertionError: The cdsapi key provided is not the correct format, please ensure it conforms to: <UID>:<APIKEY>."
            cds_remove_credentials()
            raise Exception(f"{e}. The CDS apikey environment variables were deleted. Try again.")
        elif "Authorization Required" in str(e):
            cds_remove_credentials()
            raise Exception("Authorization failed. The CDS apikey environment variables were deleted. Try again.")
        elif "Resource dummy not found" in str(e):
            # catching incorrect name, but authentication was successful
            print('found CDS credentials and authorization successful')
        else:
            raise e


def cds_remove_credentials():
    """
    remove CDS url and uid:apikey environment variables and ~/.cdsapirc file
    environment variables defined in https://github.com/ecmwf/cdsapi/blob/main/cdsapi/api.py
    """
    
    keys_toremove = ["CDSAPI_URL",
                     "CDSAPI_KEY"]
    for key in keys_toremove:
        if key in os.environ.keys():
            os.environ.pop(key)
    
    file_cds_credentials = cds_get_file()
    if os.path.isfile(file_cds_credentials):
        os.remove(file_cds_credentials)


def cds_client_withargs():
    """
    Initialize a csdapi client with url and key from os environment variables
    These are the default values for the url/key arguments, but somehow this 
    is necessary to use the key/url from environment variables.
    """
    c = cdsapi.Client(url=os.environ.get("CDSAPI_URL"),
                      key=os.environ.get("CDSAPI_KEY"))
    return c


def download_CMEMS(varkey,
                   longitude_min, longitude_max, latitude_min, latitude_max, 
                   date_min, date_max, freq='D',
                   dataset_id=None, buffer=None,
                   dir_output='.', file_prefix='', overwrite=False):
    """
    https://help.marine.copernicus.eu/en/articles/8283072-copernicus-marine-toolbox-api-subset
    """
    copernicusmarine_remove_manual_credentials_file()
    copernicusmarine_credentials()
    
    # We floor/ceil the input timestamps to make sure we
    # subset enough data in case of data with daily timesteps.
    date_min = pd.Timestamp(date_min).floor('1d')
    date_max = pd.Timestamp(date_max).ceil('1d')

    if dataset_id is None:
        dataset_id = copernicusmarine_get_dataset_id(varkey, date_min, date_max)
    if buffer is None:
        buffer = copernicusmarine_get_buffer(dataset_id)
    # date_range with same start as stoptime is a bit tricky so we limit freqs: https://github.com/Deltares/dfm_tools/issues/720
    if freq not in ["D","M"]:
        raise ValueError(f"freq should be 'D' or 'M', not {freq}")
    
    # make sure the data fully covers more than the desired spatial extent
    longitude_min -= buffer
    longitude_max += buffer
    latitude_min  -= buffer
    latitude_max  += buffer

    dataset = copernicusmarine.open_dataset(
         dataset_id = dataset_id,
         variables = [varkey],
         minimum_longitude = longitude_min,
         maximum_longitude = longitude_max,
         minimum_latitude = latitude_min,
         maximum_latitude = latitude_max,
         start_datetime = date_min,
         end_datetime = date_max,
    )
    
    Path(dir_output).mkdir(parents=True, exist_ok=True)
    
    if freq is None:
        date_str = f"{date_min.strftime('%Y%m%d')}_{date_max.strftime('%Y%m%d')}"
        name_output = f'{file_prefix}{varkey}_{date_str}.nc'
        output_filename = Path(dir_output,name_output)
        if output_filename.is_file() and not overwrite:
            print(f'"{name_output}" found and overwrite=False, returning.')
            return
        print(f'xarray writing netcdf file: {name_output}')
        dataset.to_netcdf(output_filename)
    else:
        period_range = pd.period_range(date_min,date_max,freq=freq)
        for date in period_range:
            date_str = str(date)
            name_output = f'{file_prefix}{varkey}_{date_str}.nc'
            output_filename = Path(dir_output,name_output)
            if output_filename.is_file() and not overwrite:
                print(f'"{name_output}" found and overwrite=False, continuing.')
                continue
            dataset_perperiod = dataset.sel(time=slice(date_str, date_str))
            print(f'xarray writing netcdf file: {name_output}')
            dataset_perperiod.to_netcdf(output_filename)


def copernicusmarine_get_product(date_min, date_max):
    # time extents as global variables, so they only has to be retreived once per download run (otherwise once per variable)
    global reanalysis_tstart, reanalysis_tstop
    global forecast_tstart, forecast_tstop
    if 'reanalysis_tstart' not in globals():
        print('retrieving time range of CMEMS reanalysis and forecast products') #assuming here that physchem and bio reanalyisus/multiyear datasets have the same enddate, this seems safe
        reanalysis_tstart, reanalysis_tstop = copernicusmarine_dataset_timerange(dataset_id="cmems_mod_glo_phy_my_0.083deg_P1D-m")
        forecast_tstart, forecast_tstop = copernicusmarine_dataset_timerange(dataset_id="cmems_mod_glo_phy_anfc_0.083deg_P1D-m")
    
    if (date_min >= reanalysis_tstart) & (date_max <= reanalysis_tstop):
        product = 'reanalysis'
    elif (date_min >= forecast_tstart) & (date_max <= forecast_tstop):
        product = 'analysisforecast'
    else:
        raise ValueError(f'Requested timerange ({date_min} to {date_max}) is not fully within timerange of reanalysis product ({reanalysis_tstart} to {reanalysis_tstop}) or forecast product ({forecast_tstart} to {forecast_tstop}).')
    print(f"The CMEMS '{product}' product will be used.")
    return product


def copernicusmarine_get_dataset_id(varkey, date_min, date_max):
    #TODO: maybe get dataset_id from 'copernicusmarine describe --include-datasets --contains <search_token>'
    
    product = copernicusmarine_get_product(date_min, date_max)
    
    if varkey in ['bottomT','tob','mlotst','siconc','sithick','so','thetao','uo','vo','usi','vsi','zos']: #for physchem
        # resolution is 1/12 degrees in lat/lon dimension, but a bit more/less in alternating cells
        if product == 'analysisforecast': #forecast: https://data.marine.copernicus.eu/product/GLOBAL_ANALYSISFORECAST_PHY_001_024/description
            if varkey in ['uo','vo']: #anfc datset is splitted over multiple urls
                dataset_id = 'cmems_mod_glo_phy-cur_anfc_0.083deg_P1D-m'
            elif varkey in ['so']:
                dataset_id = 'cmems_mod_glo_phy-so_anfc_0.083deg_P1D-m'
            elif varkey in ['thetao']:
                dataset_id = 'cmems_mod_glo_phy-thetao_anfc_0.083deg_P1D-m'
            else:
                dataset_id = 'cmems_mod_glo_phy_anfc_0.083deg_P1D-m'
        else: #reanalysis: https://data.marine.copernicus.eu/product/GLOBAL_MULTIYEAR_PHY_001_030/description
            dataset_id = 'cmems_mod_glo_phy_my_0.083deg_P1D-m'
    else: #for bio (resolution is 1/4 degrees)
        if product == 'analysisforecast': #forecast: https://data.marine.copernicus.eu/product/GLOBAL_ANALYSISFORECAST_BGC_001_028/description
            dataset_id = 'global-analysis-forecast-bio-001-028-daily'
        else: #https://data.marine.copernicus.eu/product/GLOBAL_MULTIYEAR_BGC_001_029/description
            dataset_id = 'cmems_mod_glo_bgc_my_0.25_P1D-m'
    return dataset_id


def copernicusmarine_get_buffer(dataset_id):
    ds = copernicusmarine.open_dataset(dataset_id=dataset_id)
    try:
        resolution = ds.latitude.attrs["step"]
        buffer = 2 * resolution
    except (AttributeError, KeyError, TypeError):
        print("failed to automatically derive a buffer from the dataset, using buffer=0")
        buffer = 0
    return buffer


def copernicusmarine_remove_manual_credentials_file():
    cmems_file_old = os.path.expanduser("~/CMEMS_credentials.txt")
    if os.path.isfile(cmems_file_old):
        os.remove(cmems_file_old)


def copernicusmarine_credentials():
    """
    To login at copernicusmarine if this was not the case yet.
    If the credentials file is present, the function returns None.
    If the credentials file is not present, get_and_check_username_password first
    checks env vars and if not present it prompts the user for credentials.
    Feeding the returned credentials to copernicusmarine.login() generates the credentials file.
    If the file is available, it gets the credentials from the file.
    Either way, the credentials are returned for use in e.g. ftp login
    """
    from copernicusmarine.core_functions.credentials_utils import (
        DEFAULT_CLIENT_CREDENTIALS_FILEPATH,
        InvalidUsernameOrPassword,
        get_and_check_username_password,
    )
    login_kwargs = dict(username=None, password=None, credentials_file=DEFAULT_CLIENT_CREDENTIALS_FILEPATH, no_metadata_cache=False)
    if not DEFAULT_CLIENT_CREDENTIALS_FILEPATH.is_file():
        print("Downloading CMEMS data requires a Copernicus Marine username and password, sign up for free at: https://data.marine.copernicus.eu/register.")
        username, password = get_and_check_username_password(**login_kwargs)
        success = copernicusmarine.login(username, password)
        if not success:
            raise InvalidUsernameOrPassword("invalid credentials")
    else:
        username, password = get_and_check_username_password(**login_kwargs)
    return username, password


def copernicusmarine_reset(remove_folder=False, overwrite_cache=True, update_package=False):
    if remove_folder:
        dir_copernicusmarine = os.path.expanduser("~/.copernicusmarine")
        print(f"reset copernicusmarine: removing {dir_copernicusmarine}, you will have to login again.")
        shutil.rmtree(dir_copernicusmarine, ignore_errors=True)
    if overwrite_cache:
        print("reset copernicusmarine: overwriting copernicusmarine metadata cache (takes some time)")
        subprocess.run("copernicusmarine describe --overwrite-metadata-cache")
    if update_package:
        print("reset copernicusmarine: updating copernicusmarine")
        subprocess.check_call(f"{sys.executable} -m pip install copernicusmarine -U")


def copernicusmarine_dataset_timerange(dataset_id):
    ds = copernicusmarine.open_dataset(dataset_id=dataset_id)
    ds_tstart = pd.Timestamp(ds.time.isel(time=0).values)
    ds_tstop = pd.Timestamp(ds.time.isel(time=-1).values)
    return ds_tstart, ds_tstop


def open_OPeNDAP_xr(dataset_url):
    """
    
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
    
    if 'hycom.org' in dataset_url_one:
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
