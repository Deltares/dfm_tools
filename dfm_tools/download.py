import os
import pandas as pd
from pathlib import Path
import xarray as xr
from dfm_tools.errors import OutOfRangeError
import cdsapi
import copernicusmarine
from copernicusmarine.core_functions.credentials_utils import InvalidUsernameOrPassword
import cftime
import getpass
import shutil
import subprocess
import sys
from requests import HTTPError

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
                      'p140209':'air_density_over_the_oceans', # TODO: paramID might be replaced with shortname rhoao: https://jira.ecmwf.int/plugins/servlet/desk/portal/4/SD-82050
                      }
    if varkey not in variables_dict.keys(): #TODO: how to get list of available vars? mean_sea_level_pressure and msl both return a dataset with msl varkey, but standard_name air_pressure_at_mean_sea_level returns an error
        raise KeyError(f'"{varkey}" not available, choose from: {", ".join(variables_dict.keys())}')
    
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


def cds_credentials():
    """
    get cdsapikey from environment variables or file or query via getpass if necessary
    """
    # TODO: put this in a PR at https://github.com/ecmwf/cdsapi
    
    if "CDSAPI_KEY" in os.environ.keys():
        # in case of envvars, also set CDSAPI_URL envvar so it does not have to be supplied
        cds_url = os.environ.get("CDSAPI_URL", "https://cds.climate.copernicus.eu/api")
        os.environ["CDSAPI_URL"] = cds_url
    
    try:
        # gets url/key from env vars or ~/.cdsapirc file
        c = cdsapi.Client()
        cds_url = c.url
        cds_apikey = c.key
    except Exception as e:
        if "Missing/incomplete configuration file" in str(e):
            # query apikey if not present in file or envvars
            print("Downloading CDS/ERA5 data requires a ECMWF API-key, copy "
                  "your API-key from https://cds.climate.copernicus.eu/profile "
                  "(first register, login and accept the terms). "
                  "More info in https://forum.ecmwf.int/t/3743.")
            cds_apikey = getpass.getpass("\nEnter your ECMWF API-key (string with dashes): ")
            cds_set_credentials(cds_url, cds_apikey)
        else:
            raise e
    
    # remove cdsapirc file or env vars if the url/apikey are according to old format
    old_urls = [
        "https://cds.climate.copernicus.eu/api/v2",
        "https://cds-beta.climate.copernicus.eu/api",
        ]
    if cds_url in old_urls:
        # to avoid "HTTPError: 401 Client Error: Unauthorized for url"
        cds_remove_credentials_raise(reason='Old CDS URL found')
    if ":" in cds_apikey:
        # to avoid "AttributeError: 'Client' object has no attribute 'client'"
        cds_remove_credentials_raise(reason='Old CDS API-key found (with :)')
    
    # check if the authentication works
    try:
        c = cdsapi.Client()
        c.client.check_authentication()
        print('found ECMWF API-key and authorization successful')
    except HTTPError as e:
        cds_remove_credentials_raise(reason=str(e))


def cds_get_file():
    file_cds_credentials = os.environ.get("CDSAPI_RC", os.path.expanduser("~/.cdsapirc"))
    return file_cds_credentials

    
def cds_set_credentials_rcfile(cds_url, cds_apikey):
    file_cds_credentials = cds_get_file()
    with open(file_cds_credentials,'w') as fc:
        fc.write(f'url: {cds_url}\n')
        fc.write(f'key: {cds_apikey}')
    
    
def cds_set_credentials(cds_url, cds_apikey):
    # set env vars
    os.environ["CDSAPI_URL"] = cds_url
    os.environ["CDSAPI_KEY"] = cds_apikey
    
    # set ~/.cdsapirc file
    cds_set_credentials_rcfile(cds_url, cds_apikey)


def cds_remove_credentials_raise(reason=''):
    """
    remove CDS url and uid:apikey environment variables and ~/.cdsapirc file
    environment variables defined in https://github.com/ecmwf/cdsapi/blob/main/cdsapi/api.py
    """
    
    keys_toremove = ["CDSAPI_URL",
                     "CDSAPI_KEY",
                     "CDSAPI_RC"]
    for key in keys_toremove:
        if key in os.environ.keys():
            os.environ.pop(key)

    file_cds_credentials = cds_get_file()
    if os.path.isfile(file_cds_credentials):
        os.remove(file_cds_credentials)
    
    raise ValueError(f"{reason}. The CDS/ECMWF apikey file (~/.cdsapirc) was deleted. "
                     "Please try again and provide your apikey when prompted.")


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


def copernicusmarine_get_product(date_min, date_max, vartype):
    assert vartype in ['phy','bio']
    
    # the time extents between phy and bio can be different so we have to retrieve them both
    
    # time extents as global variables, so they only has to be retreived once per download run (otherwise once per variable)
    global phy_reanalysis_tstart, phy_reanalysis_tstop
    global phy_reanalysis_int_tstart, phy_reanalysis_int_tstop
    global phy_forecast_tstart, phy_forecast_tstop
    global bio_reanalysis_tstart, bio_reanalysis_tstop
    global bio_reanalysis_int_tstart, bio_reanalysis_int_tstop
    global bio_forecast_tstart, bio_forecast_tstop
    
    # retrieve times
    if vartype=='phy' and 'phy_reanalysis_tstart' not in globals():
        print('retrieving time range of CMEMS reanalysis, reanalysis-interim and forecast products (phy)') #assuming here that physchem and bio reanalyisus/multiyear datasets have the same enddate, this seems safe
        phy_reanalysis_tstart, phy_reanalysis_tstop = copernicusmarine_dataset_timerange(dataset_id="cmems_mod_glo_phy_my_0.083deg_P1D-m")
        phy_reanalysis_int_tstart, phy_reanalysis_int_tstop = copernicusmarine_dataset_timerange(dataset_id="cmems_mod_glo_phy_myint_0.083deg_P1D-m")
        phy_forecast_tstart, phy_forecast_tstop = copernicusmarine_dataset_timerange(dataset_id="cmems_mod_glo_phy_anfc_0.083deg_P1D-m")
    if vartype=='bio' and 'bio_reanalysis_tstart' not in globals():
        print('retrieving time range of CMEMS reanalysis, reanalysis-interim and forecast products (bio)') #assuming here that physchem and bio reanalyisus/multiyear datasets have the same enddate, this seems safe
        bio_reanalysis_tstart, bio_reanalysis_tstop = copernicusmarine_dataset_timerange(dataset_id="cmems_mod_glo_bgc_my_0.25deg_P1D-m")
        bio_reanalysis_int_tstart, bio_reanalysis_int_tstop = copernicusmarine_dataset_timerange(dataset_id="cmems_mod_glo_bgc_myint_0.25deg_P1D-m")
        bio_forecast_tstart, bio_forecast_tstop = copernicusmarine_dataset_timerange(dataset_id="cmems_mod_glo_bgc-pft_anfc_0.25deg_P1D-m")
    
    # set current start/stop times dependent on whether we request phy/bio
    if vartype=='phy':
        reanalysis_tstart, reanalysis_tstop = phy_reanalysis_tstart, phy_reanalysis_tstop
        reanalysis_int_tstart, reanalysis_int_tstop = phy_reanalysis_int_tstart, phy_reanalysis_int_tstop
        forecast_tstart, forecast_tstop = phy_forecast_tstart, phy_forecast_tstop
    if vartype=='bio' :
        reanalysis_tstart, reanalysis_tstop = bio_reanalysis_tstart, bio_reanalysis_tstop
        reanalysis_int_tstart, reanalysis_int_tstop = bio_reanalysis_int_tstart, bio_reanalysis_int_tstop
        forecast_tstart, forecast_tstop = bio_forecast_tstart, bio_forecast_tstop
        
    if (date_min >= reanalysis_tstart) & (date_max <= reanalysis_tstop):
        product = 'reanalysis'
    elif (date_min >= reanalysis_int_tstart) & (date_max <= reanalysis_int_tstop):
        product = 'reanalysis-interim'
    elif (date_min >= forecast_tstart) & (date_max <= forecast_tstop):
        product = 'analysisforecast'
    else:
        raise ValueError(f'The requested timerange ({date_min} to {date_max}) is not fully within the timerange of one of the following datasets:\n'
                         f'- reanalysis: {reanalysis_tstart} to {reanalysis_tstop}\n'
                         f'- reanalysis-interim: {reanalysis_int_tstart} to {reanalysis_int_tstop}\n'
                         f'- analysisforecast: {forecast_tstart} to {forecast_tstop}\n'
                         'Please adjust the requested timerange and try again.')
    print(f"The CMEMS '{product}' product will be used.")
    return product


def copernicusmarine_get_dataset_id(varkey, date_min, date_max):
    #TODO: maybe get dataset_id from 'copernicusmarine describe --include-datasets --contains <search_token>'
    
    if varkey in ['bottomT','tob','mlotst','siconc','sithick','so','thetao','uo','vo','usi','vsi','zos']: #for physchem
        vartype = 'phy'
    elif varkey in ['nppv','o2','talk','dissic','ph','spco2','no3','po4','si','fe','chl','phyc']: # for bio
        vartype = 'bio'
    else:
        raise KeyError(f"unknown varkey for cmems: {varkey}")
    
    product = copernicusmarine_get_product(date_min, date_max, vartype)
    
    if vartype == 'phy': #for physchem
        # resolution is 1/12 degrees in lat/lon dimension, but a bit more/less in alternating cells
        if product == 'analysisforecast': # forecast: https://data.marine.copernicus.eu/product/GLOBAL_ANALYSISFORECAST_PHY_001_024/description
            if varkey in ['uo','vo']: # anfc datset is splitted over multiple urls
                dataset_id = 'cmems_mod_glo_phy-cur_anfc_0.083deg_P1D-m'
            elif varkey in ['so']:
                dataset_id = 'cmems_mod_glo_phy-so_anfc_0.083deg_P1D-m'
            elif varkey in ['thetao']:
                dataset_id = 'cmems_mod_glo_phy-thetao_anfc_0.083deg_P1D-m'
            else:
                dataset_id = 'cmems_mod_glo_phy_anfc_0.083deg_P1D-m'
        elif product == 'reanalysis-interim': # reanalysis-interim: https://data.marine.copernicus.eu/product/GLOBAL_MULTIYEAR_BGC_001_029/description
            dataset_id = 'cmems_mod_glo_phy_myint_0.083deg_P1D-m'
        else: # reanalysis: https://data.marine.copernicus.eu/product/GLOBAL_MULTIYEAR_PHY_001_030/description
            dataset_id = 'cmems_mod_glo_phy_my_0.083deg_P1D-m'
    elif vartype == 'bio': # for bio
        # resolution is 1/4 degrees
        if product == 'analysisforecast': # forecast: https://data.marine.copernicus.eu/product/GLOBAL_ANALYSISFORECAST_BGC_001_028/description
            if varkey in ['nppv','o2']:
                dataset_id = 'cmems_mod_glo_bgc-bio_anfc_0.25deg_P1D-m'
            elif varkey in ['talk','dissic','ph']:
                dataset_id = 'cmems_mod_glo_bgc-car_anfc_0.25deg_P1D-m'
            elif varkey in ['spco2']:
                dataset_id = 'cmems_mod_glo_bgc-co2_anfc_0.25deg_P1D-m'
            elif varkey in ['no3','po4','si','fe']:
                dataset_id = 'cmems_mod_glo_bgc-nut_anfc_0.25deg_P1D-m'
            elif varkey in ['chl','phyc']:
                dataset_id = 'cmems_mod_glo_bgc-pft_anfc_0.25deg_P1D-m'
        elif product == 'reanalysis-interim': # reanalysis-interim: https://data.marine.copernicus.eu/product/GLOBAL_MULTIYEAR_BGC_001_029/description
            dataset_id = 'cmems_mod_glo_bgc_myint_0.25deg_P1D-m'
        else: # reanalysis: https://data.marine.copernicus.eu/product/GLOBAL_MULTIYEAR_BGC_001_029/description
            dataset_id = 'cmems_mod_glo_bgc_my_0.25deg_P1D-m'
    else:
        raise ValueError(f"unknown vartype for cmems: {vartype}")
    return dataset_id


def copernicusmarine_get_buffer(dataset_id):
    ds = copernicusmarine.open_dataset(dataset_id=dataset_id)
    try:
        resolution = ds.latitude.attrs["step"]
        buffer = 2 * resolution
    except (AttributeError, KeyError, TypeError):
        print("failed to automatically derive a buffer from the dataset, using buffer=0.5")
        buffer = 0.5
    return buffer


def copernicusmarine_remove_manual_credentials_file():
    cmems_file_old = os.path.expanduser("~/CMEMS_credentials.txt")
    if os.path.isfile(cmems_file_old):
        os.remove(cmems_file_old)


def copernicusmarine_credentials():
    """
    Login at copernicusmarine if user not logged in yet.
    Works via prompt, environment variables or credentials file.
    """
    print("Downloading CMEMS data requires a Copernicus Marine username and password, "
          "sign up for free at: https://data.marine.copernicus.eu/register.")
    success = copernicusmarine.login(skip_if_user_logged_in=True)
    if not success:
        raise InvalidUsernameOrPassword("Invalid credentials, please try again")


def copernicusmarine_reset(update_package=False, remove_folder=False, overwrite_cache=True):
    if update_package:
        print("reset copernicusmarine: updating copernicusmarine")
        subprocess.check_call(f"{sys.executable} -m pip install copernicusmarine -U")
    if remove_folder:
        dir_copernicusmarine = os.path.expanduser("~/.copernicusmarine")
        print(f"reset copernicusmarine: removing {dir_copernicusmarine}, you will have to login again.")
        shutil.rmtree(dir_copernicusmarine, ignore_errors=True)
    if overwrite_cache:
        print("reset copernicusmarine: overwriting copernicusmarine metadata cache (takes some time)")
        subprocess.run("copernicusmarine describe --overwrite-metadata-cache")


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
