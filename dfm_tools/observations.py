import re
import pandas as pd
import numpy as np
from urllib.request import urlopen
import unicodedata
import geopandas as gpd
from shapely import Point
import os
import xarray as xr
from ftplib import FTP
from dfm_tools.download import copernicusmarine_credentials
from erddapy import ERDDAP
import requests
from zipfile import ZipFile
from io import BytesIO
import functools
import tempfile
import ddlpy

__all__ = ["ssh_catalog_subset",
           "ssh_catalog_toxynfile",
           "ssh_retrieve_data",
           ]


def _make_hydrotools_consistent(ds):
    """
    to make resulting netcdf file consistent with hydro_tools matlab post-processing
    """
    
    # assert presence of time variable/dim (case sensitive)
    assert "time" in ds.variables
    assert "time" in ds.dims
    
    # assert presence and units of waterlevel variable
    #TODO: add/check for standard_name attr?
    assert "waterlevel" in ds.data_vars
    assert hasattr(ds["waterlevel"], "units")
    assert ds["waterlevel"].attrs["units"] == "m"
    
    ds["station_name"] = xr.DataArray(ds.attrs["station_name"]).astype("S64")
    ds["station_id"] = xr.DataArray(ds.attrs["station_id"]).astype("S64")
    x_attrs = dict(standard_name = 'longitude',
                   units = 'degrees_east')
    y_attrs = dict(standard_name = 'latitude',
                   units = 'degrees_north')
    ds["station_x_coordinate"] = xr.DataArray(ds.attrs["longitude"]).assign_attrs(x_attrs)
    ds["station_y_coordinate"] = xr.DataArray(ds.attrs["latitude"]).assign_attrs(y_attrs)


def _remove_accents(input_str):
    nfkd_form = unicodedata.normalize('NFKD', input_str)
    only_ascii = nfkd_form.encode('ASCII', 'ignore').decode('ASCII')
    return only_ascii


def _check_ssc_groups_valid(groups):
    list_validgroups = ['psmsl','ioc','ptwc','gloss','uhslc']
    if not isinstance(groups, list):
        groups = [groups]
    for groupname in groups:
        if groupname not in list_validgroups:
            raise ValueError(f"groupname should be one of {list_validgroups}, '{groupname}' is not valid")


def ssc_sscid_from_otherid(group_id, groupname):
    """
    sscid_fromcatalog = get_sscid_fromIOCcatalog(group_id=347, groupname='uhslc')
    fname_fromcatalog = IOC_catalog_pd.loc[sscid_fromcatalog,'station_name_fname']
    """
    group_id = str(group_id)
    ssc_catalog_pd = ssc_ssh_read_catalog()
    _check_ssc_groups_valid(groupname)
        
    bool_strinseries = ssc_catalog_pd[groupname].apply(lambda x: group_id in x)
    if bool_strinseries.sum() < 1:
        raise ValueError('sscid not found for id %s in group %s'%(group_id, groupname))
    if bool_strinseries.sum() > 1:
        raise ValueError('More than 1 sscid found for id %s in group %s:\n%s'%(group_id, groupname, ssc_catalog_pd.loc[bool_strinseries,['name','country', 'geo:lat', 'geo:lon',groupname]]))
    
    sscid = ssc_catalog_pd.loc[bool_strinseries].index[0]
    return sscid


def ssc_ssh_subset_groups(groups, ssc_catalog_gpd=None):
    
    if ssc_catalog_gpd is None:
        ssc_catalog_gpd = ssc_ssh_read_catalog()
    
    if not isinstance(groups, list):
        groups = [groups]
    _check_ssc_groups_valid(groups)
    bool_ingroup = ssc_catalog_gpd[groups].apply(lambda x: x.str.len()).sum(axis=1)!=0
    ssc_catalog_gpd = ssc_catalog_gpd[bool_ingroup]

    return ssc_catalog_gpd


def ssc_ssh_read_catalog():
    """
    The SSC catalog contains e.g. UHSLC and GLOSS ids, this can be used
    for unique station naming across several observation datasets

    info: https://www.ioc-sealevelmonitoring.org/ssc.php
    station list: https://www.ioc-sealevelmonitoring.org/ssc/
    """
    
    url_json = 'https://www.ioc-sealevelmonitoring.org/ssc/service.php?format=json'
    ssc_catalog_pd = pd.read_json(url_json)
    
    #TODO: country column has 2-digit codes instead of 3
    
    #convert all cells with ids to list of strings or NaN
    for colname in ['psmsl','ioc','ptwc','gloss','uhslc','sonel_gps','sonel_tg']:
        ssc_catalog_pd[colname] = ssc_catalog_pd[colname].apply(lambda x: x if isinstance(x,list) else [] if x is np.nan else [x])
    
    ssc_catalog_pd['station_name'] = ssc_catalog_pd['name']
    ssc_catalog_pd['station_id'] = ssc_catalog_pd['ssc_id']
    
    #generate station_name_fname (remove all non numeric/letter characters and replace by -, strip '-' from begin/end, set ssc_id as prefix)
    ssc_catalog_pd['station_name_fname'] = ssc_catalog_pd['name'].str.replace('ø','o') # necessary since otherwise these letters are dropped with remove_accents()
    ssc_catalog_pd['station_name_fname'] = ssc_catalog_pd['station_name_fname'].apply(lambda x: _remove_accents(x)) #remove accents from letters
    bool_somethingdropped = (ssc_catalog_pd['station_name_fname'].str.len() != ssc_catalog_pd['name'].str.len())
    if bool_somethingdropped.any():
        raise Exception('lengths mismatch, characters were dropped:\n%s'%(ssc_catalog_pd.loc[bool_somethingdropped,['name','station_name_fname']]))
    ssc_catalog_pd['station_name_fname'] = ssc_catalog_pd['station_name_fname'].apply(lambda x: re.sub("[^0-9a-zA-Z]+", "-", x)) #replace comma and brackets with dash
    ssc_catalog_pd['station_name_fname'] = ssc_catalog_pd['station_name_fname'].str.strip('-') #remove first/last dash from name if present
    col = ssc_catalog_pd.pop('station_name_fname')
    ssc_catalog_pd.insert(3, col.name, col)
    ssc_catalog_pd['station_name_unique'] = ssc_catalog_pd['ssc_id']+'_'+ssc_catalog_pd['station_name_fname']
    
    #set ssc_id as index
    ssc_catalog_pd = ssc_catalog_pd.set_index('ssc_id',drop=False)
    
    #convert datetimestrings to datetime values, check for new ones
    # last_retrieve = dt.date(2023,2,22)
    # SSC_catalog_pd['dcterms:modified'] = pd.to_datetime(SSC_catalog_pd['dcterms:modified'])#.dt.to_pydatetime()
    # bool_newstations = SSC_catalog_pd['dcterms:modified'] > pd.Timestamp(last_retrieve,tzinfo=pytz.UTC)#,tzinfo=dt.datetime.tzname('GMT'))
    # num_newstations = bool_newstations.sum()
    # if num_newstations > 0:
    #     raise Exception('Caution, %i new stations since last retrieve on %s:\n%s'%(num_newstations, last_retrieve, SSC_catalog_pd.loc[bool_newstations.values,['name','dcterms:modified']]))
    
    # remove 64 DART stations (Deep-ocean Assessment and Reporting of Tsunamis)
    # SSC_catalog_pd = SSC_catalog_pd.loc[~SSC_catalog_pd['name'].str.startswith('DART ')]
    
    # generate geom and geodataframe
    geom = [Point(x["geo:lon"], x["geo:lat"]) for irow, x in ssc_catalog_pd.iterrows()]
    ssc_catalog_gpd = gpd.GeoDataFrame(data=ssc_catalog_pd, geometry=geom, crs='EPSG:4326')
    ssc_catalog_gpd = ssc_catalog_gpd.drop(["geo:lon","geo:lat"], axis=1)
    
    # compare coordinates of station metadata with coordinates of IOC/UHSLC linked stations
    if 0: 
        for station_ssc_id, row in ssc_catalog_gpd.iterrows():
            idx = ssc_catalog_gpd.index.tolist().index(station_ssc_id)
            print(f'station {idx+1} of {len(ssc_catalog_gpd)}: {station_ssc_id}')
            ssc_catalog_pd_stat_ioc_uhslc = ssc_catalog_gpd.loc[station_ssc_id,['ioc','uhslc']]
    
            url_station = f'https://www.ioc-sealevelmonitoring.org/ssc/stationdetails.php?id={station_ssc_id}'
            url_response = urlopen(url_station)
            url_response_read = url_response.read()
            
            station_meta_lon = ssc_catalog_gpd.loc[station_ssc_id].geometry.x
            station_meta_lat = ssc_catalog_gpd.loc[station_ssc_id].geometry.y
    
            if (ssc_catalog_pd_stat_ioc_uhslc.str.len()==0).all(): #skip station if no IOC/UHSLC id present (after retrieval of precise coordinate)
                continue
            # fix html, fetch last matched table and set row with codes/location/lat/lon/sensors as column names
            # TODO: report missing <tr> to VLIZ?
            url_response_read_fixed = url_response_read.replace(b' colspan="100%"',b'').replace(b'<td><a href',b'<tr><td><a href')
            table2 = pd.read_html(url_response_read_fixed, match='Linked codes', header=0)
            tab_linked = table2[-1]
            tab_linked.columns = tab_linked.iloc[0]

            #loop over IOC/UHSLC linked stations
            bool_tocheck = tab_linked["Codes"].str.contains('UHSLC') | tab_linked["Codes"].str.contains('IOC')
            if bool_tocheck.sum()==0:
                continue
            tab_linked_tocheck = tab_linked.loc[bool_tocheck]
            station_check_dict = {}
            station_check_dist_all = []
            for _, row in tab_linked_tocheck.iterrows():
                station_check_lat = float(row['Latitude'])
                station_check_lon = float(row['Longitude'])
                station_check_dist = np.sqrt((station_meta_lat-station_check_lat)**2+(station_meta_lon-station_check_lon)**2)
                station_check_dist_all.append(station_check_dist)
                station_check_dict[row['Codes']] = [station_check_dist,station_check_lat,station_check_lon]
            ssc_catalog_gpd.loc[station_ssc_id,'dist_dict'] = [station_check_dict]
            ssc_catalog_gpd.loc[station_ssc_id,'dist_min'] = np.min(station_check_dist_all)
            ssc_catalog_gpd.loc[station_ssc_id,'dist_max'] = np.max(station_check_dist_all)

    return ssc_catalog_gpd


def get_cmems_params(source):
    params_cmems_my = {"host": "my.cmems-du.eu",
                       "cwd": "Core/INSITU_GLO_PHY_SSH_DISCRETE_MY_013_053/cmems_obs-ins_glo_phy-ssh_my_na_PT1H"}
    params_cmems_nrt = {"host": "nrt.cmems-du.eu",
                        "cwd": "Core/INSITU_GLO_PHYBGCWAV_DISCRETE_MYNRT_013_030/cmems_obs-ins_glo_phybgcwav_mynrt_na_irr"}
    params_dict = {"cmems": params_cmems_my,
                   "cmems-nrt": params_cmems_nrt}
    params = params_dict[source]
    return params


def cmems_my_ssh_read_catalog():
    cmems_catalog_gpd = cmems_ssh_read_catalog(source="cmems")
    return cmems_catalog_gpd
    

def cmems_nrt_ssh_read_catalog():
    cmems_catalog_gpd = cmems_ssh_read_catalog(source="cmems-nrt")
    return cmems_catalog_gpd


def cmems_ssh_read_catalog(source):
    cmems_params = get_cmems_params(source)
    
    # setup ftp connection
    ftp = FTP(host=cmems_params["host"])
    username, password = copernicusmarine_credentials()
    ftp.login(user=username, passwd=password)
    ftp.cwd(cmems_params["cwd"])
    
    # read index
    fname = 'index_history.txt'
    with open(fname, 'wb') as fp:
        ftp.retrbinary(f'RETR {fname}', fp.write)
    with open(fname, 'r') as f:
        for line in f:
            if line.startswith('#'):
                header = line
            else:
                break #stop when there are no more #
    colnames = header.strip('#').strip().split(',')
    index_history_pd = pd.read_csv(fname,comment='#',names=colnames)
    os.remove(fname) # remove the local file again
    
    # filter only history tidegauges (TG) containing SLEV variable, relevant for nrt dataset
    # TODO: why are there non-SLEV files in TG folder? Cleanup possible?
    bool_tidegauge = index_history_pd["file_name"].str.contains("/history/TG/")
    bool_slev = index_history_pd["parameters"].str.contains("SLEV")
    index_history_pd = index_history_pd.loc[bool_tidegauge & bool_slev]
    
    # drop andratx station, lat/lon vary over time in nrt dataset
    # TODO: remove this exception when the CMEMS nrt dataset is cleaned up
    bool_moving = index_history_pd["file_name"].str.contains("MO_TS_TG_ANDRATX")
    index_history_pd = index_history_pd.loc[~bool_moving]
    
    # generate geom and geodataframe
    assert (index_history_pd["geospatial_lon_min"] == index_history_pd["geospatial_lon_max"]).all()
    assert (index_history_pd["geospatial_lat_min"] == index_history_pd["geospatial_lat_max"]).all()
    geom = [Point(x["geospatial_lon_min"], x["geospatial_lat_min"]) for irow, x in index_history_pd.iterrows()]
    index_history_gpd = gpd.GeoDataFrame(data=index_history_pd, geometry=geom, crs='EPSG:4326')
    drop_list = ["geospatial_lon_min", "geospatial_lon_max",
                 "geospatial_lat_min", "geospatial_lat_max"]
    index_history_gpd = index_history_gpd.drop(drop_list, axis=1)
    
    # add dummy country column for the test to pass
    # TODO: hopefully country metadata is available via cmems API in the future
    index_history_gpd["country"] = ""
    stat_ids = index_history_gpd["file_name"].apply(lambda x: os.path.basename(x).replace(".nc",""))
    stat_names = stat_ids.str.split("_").str[3] # corresponds to the cmems platform_code
    index_history_gpd["station_name"] = stat_names
    index_history_gpd["station_id"] = stat_names
    index_history_gpd["station_name_unique"] = stat_ids
    
    # rename columns
    rename_dict = {'time_coverage_start':'time_min',
                   'time_coverage_end':'time_max'}
    index_history_gpd = index_history_gpd.rename(rename_dict, axis=1)
    index_history_gpd["time_min"] = pd.to_datetime(index_history_gpd["time_min"])
    index_history_gpd["time_max"] = pd.to_datetime(index_history_gpd["time_max"])
    
    index_history_gpd["time_ndays"] = (index_history_gpd['time_max'] - index_history_gpd['time_min']).dt.total_seconds()/3600/24
    return index_history_gpd


def uhslc_ssh_read_catalog(source):
    # TODO: country is "New Zealand" and country_code is 554. We would like country/country_code=NZL
    # TODO: maybe use min of rqds and max of fast for time subsetting
    # TODO: maybe enable merging of datasets?
    uhslc_gpd = gpd.read_file("https://uhslc.soest.hawaii.edu/data/meta.geojson")
    
    for drop_col in ["rq_basin", "rq_versions"]:
        if drop_col in uhslc_gpd.columns:
            uhslc_gpd = uhslc_gpd.drop(drop_col, axis=1)
    
    uhslc_gpd = uhslc_gpd.set_index('uhslc_id', drop=False)
    
    # shift from 0to360 to -180to180
    from shapely import Point
    geom_shift = [Point(((pnt.x + 180)%360 - 180), pnt.y) for pnt in uhslc_gpd.geometry]
    uhslc_gpd.geometry = geom_shift
    
    timespan_dict = {"uhslc-fast":"fd_span", "uhslc-rqds":"rq_span"}
    timespan_var = timespan_dict[source]
    time_min = uhslc_gpd[timespan_var].apply(lambda x: x["oldest"])
    time_max = uhslc_gpd[timespan_var].apply(lambda x: x["latest"])
    uhslc_gpd = uhslc_gpd.loc[~time_min.isnull()].copy()
    uhslc_gpd["time_min"] = pd.to_datetime(time_min)
    uhslc_gpd["time_max"] = pd.to_datetime(time_max)
    
    uhslc_gpd["station_name"] = uhslc_gpd['name']
    uhslc_gpd["station_id"] = uhslc_gpd['uhslc_id']
    stat_names = source + "-" + uhslc_gpd['uhslc_id'].apply(lambda x: f"{x:03d}")
    uhslc_gpd["station_name_unique"] = stat_names
    
    uhslc_gpd["time_ndays"] = (uhslc_gpd['time_max'] - uhslc_gpd['time_min']).dt.total_seconds()/3600/24
    return uhslc_gpd


def uhslc_rqds_ssh_read_catalog():
    uhslc_gpd = uhslc_ssh_read_catalog(source="uhslc-rqds")
    return uhslc_gpd


def uhslc_fast_ssh_read_catalog():
    uhslc_gpd = uhslc_ssh_read_catalog(source="uhslc-fast")
    return uhslc_gpd


def ioc_ssh_read_catalog(drop_uhslc=True, drop_dart=True, drop_nonutc=True):
    """
    Generates a list of all active IOC stations (showall=a).
    showall=all is all known and returns more stations, but also returns nonexistent stations.
    The stations that are already in UHSLC are dropped,
    as well as DART stations and non-UTC stations.
    """
    #TODO: "Code" contains more station codes than "code", what is the difference?
    #TODO: "Location" contains full name, but contains spaces etcetera, retrieve from SSC instead?
    
    url_json = 'https://www.ioc-sealevelmonitoring.org/service.php?query=stationlist&showall=a'
    resp = requests.get(url_json)
    if resp.status_code==404: #continue to next station if not found
        raise Exception(f'url 404: {resp.text}')    
    resp_json = resp.json()
    ioc_catalog_pd = pd.DataFrame.from_dict(resp_json)
    
    #set ssc_id as index
    ioc_catalog_pd = ioc_catalog_pd.set_index('Code',drop=False)
    
    #derive start/stop times indications from metadata
    ioc_catalog_pd["time_min"] = pd.to_datetime(ioc_catalog_pd["date_created"])
    ioc_catalog_pd["time_max"] = pd.to_datetime(ioc_catalog_pd["lasttime"])
    
    # generate geom and geodataframe and remove the old columns
    geom = [Point(x["lon"], x["lat"]) for irow, x in ioc_catalog_pd.iterrows()]
    ioc_catalog_gpd = gpd.GeoDataFrame(data=ioc_catalog_pd, geometry=geom, crs='EPSG:4326')
    drop_list = ["lon","lat"]
    ioc_catalog_gpd = ioc_catalog_gpd.drop(drop_list, axis=1)
    
    ioc_catalog_gpd["station_name"] = ioc_catalog_gpd['Code']
    ioc_catalog_gpd["station_id"] = ioc_catalog_gpd['Code']
    stat_names = "ioc-" + ioc_catalog_gpd['Code'] + "-" + ioc_catalog_gpd['code'].astype(str)
    ioc_catalog_gpd["station_name_unique"] = stat_names

    if drop_uhslc:
        # filter non-UHSLC stations
        ssc_catalog_pd = ssh_catalog_subset(source='ssc')
        ioc_catalog_gpd['UHSLC'] = False
        ioc_catalog_gpd['inSSClist'] = False
        for ioc_code in ioc_catalog_gpd['Code']:
            iocstat_uhslcid = ssc_catalog_pd.loc[ssc_catalog_pd['ioc'].apply(lambda x: ioc_code in x),'uhslc']
            if len(iocstat_uhslcid)==0: #ioc_code not in the SSC list
                continue
            else:
                ioc_catalog_gpd.loc[ioc_code,'inSSClist'] = True
                if len(iocstat_uhslcid.iloc[0])==0: #station contains no UHSLC code
                    continue
                else:
                    ioc_catalog_gpd.loc[ioc_code,'UHSLC'] = True
        ioc_catalog_gpd = ioc_catalog_gpd.loc[~ioc_catalog_gpd['UHSLC']]
    
    if drop_dart:
        bool_dart = ioc_catalog_gpd["Location"].str.startswith("DART ")
        ioc_catalog_gpd = ioc_catalog_gpd.loc[~bool_dart]
    
    if drop_nonutc:
        # filter out all non-UTC stations
        ioc_catalog_gpd = ioc_catalog_gpd.loc[ioc_catalog_gpd['UTCOffset']==0]
    
    ioc_catalog_gpd["time_ndays"] = (ioc_catalog_gpd['time_max'] - ioc_catalog_gpd['time_min']).dt.total_seconds()/3600/24
    return ioc_catalog_gpd


def psmsl_gnssir_ssh_read_catalog():
    # https://psmsl.org/data/gnssir/metadatainfo.php
    # https://psmsl.org/data/gnssir/useful_files.php
    # url = "https://psmsl.org/data/gnssir/data/maplayers/good_sites.json"
    # TODO: use only good_sites instead of all (request field in json)
    # TODO: request time extents in json instead of with psmsl_gnssir_ssh_read_catalog_gettimes()
    url = "https://psmsl.org/data/gnssir/data/sites.json"
    station_list_pd = pd.read_json(url).T
    
    rename_dict = {"CountryCode":"country"}
    station_list_pd = station_list_pd.rename(rename_dict, axis=1)
    
    # generate geom and geodataframe and remove the old columns
    geom = [Point(x["Longitude"], x["Latitude"]) for irow, x in station_list_pd.iterrows()]
    station_list_gpd = gpd.GeoDataFrame(data=station_list_pd, geometry=geom, crs='EPSG:4326')
    drop_list = ["Longitude","Latitude"]
    station_list_gpd = station_list_gpd.drop(drop_list, axis=1)
    
    station_list_gpd['station_name'] = station_list_gpd['Code']
    station_list_gpd["station_id"] = station_list_gpd.index
    stat_names = "psmsl-gnssir-" + station_list_gpd['station_id'].astype(str) + "-" + station_list_gpd['station_name']
    station_list_gpd["station_name_unique"] = stat_names
    return station_list_gpd


def psmsl_gnssir_ssh_read_catalog_gettimes(station_list_gpd):
    # the catalog json does not contain time ranges so we derive it from daily csv files
    station_list_gpd["time_min"] = pd.NaT
    station_list_gpd["time_max"] = pd.NaT
    print(f"retrieving psmsl-gnssir time extents for {len(station_list_gpd)} stations:", end=" ")
    for station_id in station_list_gpd.index:
        irow = station_list_gpd.index.tolist().index(station_id)
        print(irow+1, end=" ")
        url = f"https://psmsl.org/data/gnssir/data/daily/{station_id}_daily.csv"
        data_daily = pd.read_csv(url)
        time_min = pd.Timestamp(data_daily["time"].iloc[0])
        time_max = pd.Timestamp(data_daily["time"].iloc[-1])
        station_list_gpd.loc[station_id, "time_min"] = time_min
        station_list_gpd.loc[station_id, "time_max"] = time_max
        station_list_gpd.loc[station_id, "time_ndays"] = (time_max - time_min).total_seconds()/3600/24
    print()

    return station_list_gpd


def gesla3_ssh_read_catalog(file_gesla3_meta=None, only_coastal=True):
    if file_gesla3_meta is None:
        file_gesla3_meta = r"p:\1230882-emodnet_hrsm\data\GESLA3\GESLA3_ALL 2.csv"
    
    if not os.path.isfile(file_gesla3_meta):
        raise FileNotFoundError(f"The 'file_gesla3_meta' file '{file_gesla3_meta}' was not found. "
                                "You can download it from https://gesla787883612.wordpress.com/downloads and provide the path")
    
    station_list_pd = pd.read_csv(file_gesla3_meta)
    station_list_pd.columns = [c.replace(" ", "_").replace("(", "").replace(")", "").replace("/", "_").lower() for c in station_list_pd.columns]
    station_list_pd = station_list_pd.set_index('file_name', drop=False)
    station_list_pd["start_date_time"] = pd.to_datetime(station_list_pd["start_date_time"])
    station_list_pd["end_date_time"] = pd.to_datetime(station_list_pd["end_date_time"])
    
    # drop non-coastal
    if only_coastal:
        station_list_pd = station_list_pd.loc[station_list_pd['gauge_type']=='Coastal']
        
    # generate geom and geodataframe and remove the old columns
    geom = [Point(x["longitude"], x["latitude"]) for irow, x in station_list_pd.iterrows()]
    station_list_gpd = gpd.GeoDataFrame(data=station_list_pd, geometry=geom, crs='EPSG:4326')
    drop_list = ["longitude","latitude"]
    station_list_gpd = station_list_gpd.drop(drop_list, axis=1)
    
    stat_names = station_list_gpd["file_name"]
    station_list_gpd["station_name_unique"] = stat_names
    station_list_gpd["station_id"] = station_list_gpd["file_name"]
    station_list_gpd["station_name"] = station_list_gpd["site_name"]
    
    # rename columns
    rename_dict = {'start_date_time':'time_min',
                   'end_date_time':'time_max'}
    station_list_gpd = station_list_gpd.rename(rename_dict, axis=1)
    station_list_gpd["time_ndays"] = (station_list_gpd['time_max'] - station_list_gpd['time_min']).dt.total_seconds()/3600/24
    return station_list_gpd


def rwsddl_ssh_meta_dict():
    # combination for measured waterlevels
    meta_dict = {'Grootheid.Code':'WATHTE', 'Groepering.Code':'NVT'}
    return meta_dict


def rwsddl_ssh_get_time_max(locations):
    locations = locations.copy() # copy to prevent SettingWithCopyWarning
    
    print(f"getting time_max for {len(locations)} locations: ", end="")
    dtstart = pd.Timestamp.now()
    list_time_latest = []
    for stat_code, location in locations.iterrows():
        try:
            meas_latest = ddlpy.measurements_latest(location)
            time_latest = meas_latest.index.max().tz_convert(None)
        except ddlpy.ddlpy.NoDataException:
            time_latest = pd.NaT
        list_time_latest.append(time_latest)
    locations["time_max"] = list_time_latest
    print(f'{(pd.Timestamp.now()-dtstart).total_seconds():.2f} sec')
    return locations


def rwsddl_ssh_read_catalog(meta_dict=None):
    """
    convert LocatieLijst to geopandas dataframe
    """
    if meta_dict is None:
        meta_dict = rwsddl_ssh_meta_dict()
    
    locations = ddlpy.locations()
    selected = locations.copy()
    for key in meta_dict.keys():
        value = meta_dict[key]
        try:
            bool_sel = selected[key].isin([value])
            selected = selected.loc[bool_sel]
        except KeyError:
            cols_with_code = [x for x in selected.columns if ".Code" in x]
            raise KeyError(f"ddlpy.locations() cannot subset for '{key}', available are {cols_with_code}")            
    
    # add "Code" index as column and reset the index
    selected = selected.reset_index()
    
    xcoords = selected["X"]
    ycoords = selected["Y"]
    epsg_all = selected["Coordinatenstelsel"]
    epsg_uniq = epsg_all.unique()
    if len(epsg_uniq)>1:
        raise ValueError(f"multiple EPSG codes in one LocatieLijst not supported: {epsg_uniq.tolist()}")
    epsg = epsg_uniq[0]
    geom_points = [Point(x,y) for x,y in zip(xcoords,ycoords)]
    ddl_slev_gdf = gpd.GeoDataFrame(selected, geometry=geom_points, crs=epsg)
    
    # convert coordinates to wgs84
    ddl_slev_gdf = ddl_slev_gdf.to_crs(4326)
    
    ddl_slev_gdf["station_name"] = ddl_slev_gdf["Naam"]
    ddl_slev_gdf["station_name_unique"] = ddl_slev_gdf["Code"]
    ddl_slev_gdf["station_id"] = ddl_slev_gdf["Code"]
    ddl_slev_gdf["country"] = "NLD"
    
    return ddl_slev_gdf


@functools.lru_cache
def cmems_ftp_login(host, dir_data):
    """
    cache ftp connection so it does not matter that we call it for each station
    """
    # setup ftp connection
    ftp = FTP(host=host)
    username, password = copernicusmarine_credentials()
    ftp.login(user=username, passwd=password)
    ftp.cwd(dir_data)
    return ftp


def cmems_ssh_retrieve_data(row, dir_output, time_min=None, time_max=None):
    """
    Retrieve data from FTP
    Can only retrieve entire files, subsetting is done during reconstruction
    
    """
    
    # get source from gdf, uniqueness is checked in ssh_retrieve_data
    source = row['source']
    cmems_params = get_cmems_params(source)

    # get cached ftp connection
    host = cmems_params["host"]
    dir_data = os.path.dirname(row["file_name"]).split(host)[1]
    ftp = cmems_ftp_login(host, dir_data)
    
    fname = os.path.basename(row["file_name"])
    tempdir = tempfile.gettempdir()
    fname_out_raw = os.path.join(tempdir, "dfmtools_cmems_ssh_retrieve_data_temporary_file.nc")
    with open(fname_out_raw, 'wb') as fp:
        ftp.retrbinary(f'RETR {fname}', fp.write)
    
    # reconstruct this dataset (including time subsetting) and write again
    ds = xr.open_dataset(fname_out_raw)
    
    # reduce DEPTH dimension if present (nrt dataset)
    if "DEPTH" in ds.dims:
        ds = ds.max(dim="DEPTH", keep_attrs=True)

    ds = ds.rename(TIME="time")
    
    # keep only data with "good_data" quality flag
    ds["SLEV"] = ds.SLEV.where(ds.SLEV_QC==1)
    ds = ds.rename_vars(SLEV="waterlevel")
    
    # drop all unnecessary vars
    ds = ds.drop_vars(["SLEV_QC","TIME_QC","DEPH_QC","LATITUDE","LONGITUDE","STATION","POSITION"], errors="ignore")
    
    # check order of time variable
    if not ds.time.to_pandas().index.is_monotonic_increasing:
        # TODO: happens in some MO_TS_TG_RMN-* stations in NRT dataset, asked to fix
        # Genova, Imperia, LaSpezia, Livorno, Ravenna, Venice
        stat_name = row.loc["station_name_unique"]
        print(f"[{stat_name} NOT MONOTONIC] ", end="")
        del ds
        os.remove(fname_out_raw)
        return
        # ds = ds.sortby("time")
    
    # slice on time extent
    ds = ds.sel(time=slice(time_min, time_max))
    if len(ds.time) == 0:
        print("[NODATA] ", end="")
        del ds
        os.remove(fname_out_raw)
        return
    
    return ds


def uhslc_ssh_retrieve_data(row, dir_output, time_min=None, time_max=None):
    
    # TODO: maybe merge rqds/fast datasets automatically
    # TODO: 5 meter offset?
    
    # docs from https://ioos.github.io/erddapy/ and https://ioos.github.io/erddapy/02-extras-output.html#
    
    # setup server connection, this takes no time so does not have to be cached
    server = "https://uhslc.soest.hawaii.edu/erddap"
    e = ERDDAP(server=server, protocol="tabledap", response="nc") #opendap is way slower than nc/csv/html
    
    dataset_id_dict = {"uhslc-fast":"global_hourly_fast",
                       "uhslc-rqds":"global_hourly_rqds"}
    
    uhslc_id = row.name
    source = row["source"]
    dataset_id = dataset_id_dict[source]
    e.dataset_id = dataset_id
    
    # set erddap constraints
    e.constraints = {}
    e.constraints["uhslc_id="] = uhslc_id
    if time_min is not None:
        e.constraints["time>="] = pd.Timestamp(time_min)
    if time_max is not None:
        e.constraints["time<="] = pd.Timestamp(time_max)
    
    from httpx import HTTPError
    try:
        ds = e.to_xarray()
    except HTTPError:
        # no data so early return
        return
    
    # change longitude from 0to360 to -180to180
    ds["longitude"] = (ds.longitude + 180)%360 - 180
    ds = ds.sortby("longitude")
    
    # change units to meters
    with xr.set_options(keep_attrs=True):
        assert ds["sea_level"].attrs["units"] == "millimeters"
        ds["sea_level"] = ds["sea_level"]/1000
        ds["sea_level"] = ds["sea_level"].assign_attrs(units="m")
    ds = ds.rename_vars(sea_level="waterlevel")
    
    # set time index
    ds = ds.set_index(obs="time").rename(obs="time")
    ds['time'] = ds.time.dt.round('s') #round to seconds
    return ds


@functools.lru_cache
def gesla3_cache_zipfile(file_gesla3_data=None):
    if file_gesla3_data is None:
        file_gesla3_data = r"p:\1230882-emodnet_hrsm\data\GESLA3\GESLA3.0_ALL.zip"
    
    if not os.path.isfile(file_gesla3_data):
        raise FileNotFoundError(f"The 'file_gesla3_data' file '{file_gesla3_data}' was not found. "
                                "You can download it from https://gesla787883612.wordpress.com/downloads and provide the path")
    
    gesla3_zip = ZipFile(file_gesla3_data)
    return gesla3_zip


def gesla3_ssh_retrieve_data(row, dir_output, time_min=None, time_max=None,
                             file_gesla3_data=None):
    
    # get cached gesla3 zipfile instance
    gesla3_zip = gesla3_cache_zipfile(file_gesla3_data=file_gesla3_data)
    
    file_gesla = row.name
    with gesla3_zip.open(file_gesla, "r") as f:
        data = pd.read_csv(f, comment='#', delim_whitespace = True,
                           names=["date", "time", "sea_level", "qc_flag", "use_flag"],
                           parse_dates=[[0, 1]], index_col=0)
    
    # clean up time duplicates
    data.index.name = 'time'
    bool_duplicate = data.index.duplicated()
    if bool_duplicate.sum() > 0:
        data = data.loc[~bool_duplicate]
        print(f"[{bool_duplicate.sum()} duplicate timestamps removed] ", end="")
    
    # convert to xarray and add metadata
    meta_sel = row.copy()
    geometry = meta_sel.pop('geometry')
    meta_sel["time_min"] = str(meta_sel["time_min"])
    meta_sel["time_max"] = str(meta_sel["time_max"])
    meta_sel["longitude"] = geometry.x
    meta_sel["latitude"] = geometry.y
    ds = data.to_xarray().assign_attrs(meta_sel.to_dict())
    ds['site_name'] = xr.DataArray([meta_sel.site_name], dims=('stations'))
    ds['latitude'] = xr.DataArray([meta_sel.latitude], dims=('stations'))
    ds['longitude'] = xr.DataArray([meta_sel.longitude], dims=('stations'))
    ds['sea_level'] = ds['sea_level'].assign_attrs(units="m")
    ds = ds.rename_vars(sea_level="waterlevel")
    
    # subset time
    ds = ds.sel(time=slice(time_min, time_max))
    if len(ds.time) == 0:
        return
    
    # filter bad quality data
    ds = ds.where(ds.qc_flag==1)
    return ds


def ioc_ssh_retrieve_data(row, dir_output, time_min, time_max, subset_hourly=False):
    
    # https://www.ioc-sealevelmonitoring.org/service.php?query=help
    
    station_code = row.name
    period_range = pd.period_range(time_min, time_max, freq="1M")
    
    if time_min is None or time_max is None:
        raise ValueError("cannot supply None for 'time_min' or 'time_max' to 'ioc_ssh_retrieve_data()'")
    
    results_list = []
    for date in period_range:
        print('.', end="")
        year = date.year
        month = date.month
        starttime = date.to_timestamp()
        if month==12:
            endtime = pd.Timestamp(year+1,1,1)
        else:
            endtime = pd.Timestamp(year,month+1,1)
        url_json = (f"https://www.ioc-sealevelmonitoring.org/service.php?query=data&"
                    f"code={station_code}&format=json&"
                    f"timestart={starttime.isoformat()}&timestop={endtime.isoformat()}")
        resp = requests.get(url_json)
        if resp.status_code==404: #continue to next station if not found
            raise Exception(f'url 404: {resp.text}')    
        if resp.text == '[]':
            print("[NODATA]", end="")
            continue
        resp_json = resp.json()
        if 'error' in resp.json()[0].keys():
            raise Exception(resp.text)
        data_pd_one = pd.DataFrame.from_dict(resp_json)
        results_list.append(data_pd_one)
    print(" ", end="")
    
    if len(results_list)==0:
        # continue with next station if no data present in entire period
        return
    
    # convert to xarray
    data_pd_all = pd.concat(results_list)
    data_pd_all = data_pd_all.rename({"stime":"time"},axis=1)
    data_pd_all = pd.DataFrame({'slevel':data_pd_all['slevel'].values},
                            index=pd.to_datetime(data_pd_all['time']))
    data_pd_all = data_pd_all[~data_pd_all.index.duplicated(keep='last')]
    if subset_hourly:
        data_pd_all = data_pd_all.loc[data_pd_all.index.minute==0]
    ds = data_pd_all.to_xarray()
    ds = ds.rename_vars(slevel="waterlevel")
    ds["waterlevel"] = ds["waterlevel"].assign_attrs({"units":"m"})
    return ds


def psmsl_gnssir_ssh_retrieve_data(row, dir_output, time_min=None, time_max=None):
    
    # https://psmsl.org/data/gnssir/gnssir_daily_means.html
    # https://psmsl.org/data/gnssir/gnssir_example.html (also contains IOC retrieval example)
    
    station_id = row.name
    resp = urlopen(rf"https://psmsl.org/data/gnssir/data/main/{station_id}.zip")
    myzip = ZipFile(BytesIO(resp.read()))
    with myzip.open(f"{station_id}.csv") as f:
        data = pd.read_csv(f, comment="#", parse_dates=["time"])
    
    url = f"https://psmsl.org/data/gnssir/data/sites/{station_id}.json"
    station_meta = pd.read_json(url)["properties"]
    #TODO: the below equation is a guess to get the reference level about right
    data['slev'] = data['adjusted_height'] - station_meta['ellipsoidalHeight'] + station_meta['reflectorHeight']
    data = data.set_index("time")
    
    ds = data.to_xarray()
    ds['slev'] = ds['slev'].assign_attrs(units="m")
    ds = ds.rename_vars(slev="waterlevel")

    ds = ds.sel(time=slice(time_min, time_max))
    if len(ds.time) == 0:
        return
    return ds


def rwsddl_ssh_retrieve_data(row, dir_output, time_min, time_max, meta_dict=None):
    
    if time_min is None or time_max is None:
        raise ValueError("cannot supply None for 'time_min' or 'time_max' to 'rwsddl_ssh_retrieve_data()'")

    if meta_dict is None:
        meta_dict = rwsddl_ssh_meta_dict()
    
    # if we pass one row to the measurements function you can get all the measurements
    measurements = ddlpy.measurements(row, time_min, time_max)
    
    if measurements.empty:
        # no output so this station is skipped
        return
    
    # drop alfanumeriek if duplicate of numeriek
    if "Meetwaarde.Waarde_Alfanumeriek" in measurements.columns and 'Meetwaarde.Waarde_Numeriek' in measurements.columns:
        measurements = measurements.drop("Meetwaarde.Waarde_Alfanumeriek", axis=1)
    
    rename_dict = {'Meetwaarde.Waarde_Numeriek':'waterlevel',
                   'WaarnemingMetadata.KwaliteitswaardecodeLijst':'QC',
                   'WaarnemingMetadata.StatuswaardeLijst':'Status'}
    measurements = measurements.rename(columns=rename_dict)
    
    # simplify dataframe
    simplified = ddlpy.simplify_dataframe(measurements)
    
    # get dataframe attrs
    ds_attrs = simplified.attrs
    # drop irrelevant attrs
    ds_attrs = {k:v for k,v in ds_attrs.items() if not k.startswith("Bemonstering") and not k.startswith("BioTaxon")}
    
    # dropping timezone is required to get proper encoding in time variable (in netcdf file)
    simplified.index = simplified.index.tz_convert(None)
    ds = simplified.to_xarray()
    ds = ds.assign_attrs(ds_attrs)
    
    # convert meters to cm
    if ds_attrs['Eenheid.Code'] != 'cm':
        raise Exception("unexpected unit")
    ds['waterlevel'] = ds['waterlevel'].assign_attrs(units="m")
    ds['waterlevel'] /= 100 #convert from cm to m
    ds.attrs.pop('Eenheid.Code')
    ds.attrs.pop('Eenheid.Omschrijving')
    return ds


def ssh_catalog_subset(source=None,
                       lon_min=-180, lon_max=180, 
                       lat_min=-90, lat_max=90, 
                       time_min=None, time_max=None,
                       **kwargs):
    # TODO: check if min<max
    # TODO: accept None but replace with min/max value from dict. Time min/max are pd.min() and pd.max()
    # TODO: accept partial None (now one None is same as all None)
    
    ssh_sources = {"ssc": ssc_ssh_read_catalog,
                   "gesla3": gesla3_ssh_read_catalog,
                   "ioc": ioc_ssh_read_catalog,
                   "cmems": cmems_my_ssh_read_catalog,
                   "cmems-nrt": cmems_nrt_ssh_read_catalog,
                   "uhslc-fast": uhslc_fast_ssh_read_catalog,
                   "uhslc-rqds": uhslc_rqds_ssh_read_catalog,
                   "psmsl-gnssir": psmsl_gnssir_ssh_read_catalog,
                   "rwsddl": rwsddl_ssh_read_catalog,
                   }
    
    if source not in ssh_sources.keys():
        raise ValueError(f"source for ssh_catalog_subset should be one of {list(ssh_sources.keys())}, recieved '{source}'")
    
    catalog_read_func = ssh_sources[source]
    
    ssh_catalog_gpd = catalog_read_func(**kwargs)
    ssh_catalog_gpd["source"] = source

    # spatial subsetting and sort again: https://github.com/geopandas/geopandas/issues/2937
    ssh_catalog_gpd = ssh_catalog_gpd.clip((lon_min, lat_min, lon_max, lat_max))
    ssh_catalog_gpd = ssh_catalog_gpd.sort_index()
    
    if None not in [time_min, time_max]:
        if source=="psmsl-gnssir":
            ssh_catalog_gpd = psmsl_gnssir_ssh_read_catalog_gettimes(ssh_catalog_gpd)
        if "time_min" not in ssh_catalog_gpd.columns:
            raise KeyError(f"ssh_catalog_gpd for source='{source}' does not contain time_min and time_max, no time subsetting possible.")
        intime_bool = ((ssh_catalog_gpd['time_min']<time_max) &
                       (ssh_catalog_gpd['time_max']>time_min)
                       )
        ssh_catalog_gpd = ssh_catalog_gpd.loc[intime_bool].copy()
    return ssh_catalog_gpd


def ssh_catalog_toxynfile(ssc_catalog_gpd, file_xyn):
    lon = ssc_catalog_gpd.geometry.x
    lat = ssc_catalog_gpd.geometry.y
    name = ssc_catalog_gpd['station_name_unique']
    data = np.c_[lon, lat, name]
    np.savetxt(file_xyn, data, fmt='%11.6f %11.6f %-s')


def ssh_retrieve_data(ssh_catalog_gpd, dir_output, time_min=None, time_max=None,
                      **kwargs):
    
    if ssh_catalog_gpd.empty:
        raise ValueError("empty ssh_catalog_gpd provided to ssh_retrieve_data")
    
    source_list = ssh_catalog_gpd["source"].unique()
    if len(source_list) > 1:
        raise Exception("A ssh_catalog_gpd with multiple unique 'source' values was passed to ssh_retrieve_data(), this is not supported.")
    source = source_list[0]
    
    ssh_sources = {"gesla3": gesla3_ssh_retrieve_data,
                   "ioc": ioc_ssh_retrieve_data,
                   "cmems": cmems_ssh_retrieve_data,
                   "cmems-nrt": cmems_ssh_retrieve_data,
                   "uhslc-fast": uhslc_ssh_retrieve_data,
                   "uhslc-rqds": uhslc_ssh_retrieve_data,
                   "psmsl-gnssir": psmsl_gnssir_ssh_retrieve_data,
                   "rwsddl": rwsddl_ssh_retrieve_data,
                   }
    
    if source not in ssh_sources.keys():
        raise ValueError(f"source for ssh_retrieve_data should be one of {list(ssh_sources.keys())}, recieved '{source}'")
    
    retrieve_data_func = ssh_sources[source]
    
    # retrieve
    print(f"retrieving data for {len(ssh_catalog_gpd)} {source} stations:", end=" ")
    for idx_arbitrary, row in ssh_catalog_gpd.iterrows():
        irow = ssh_catalog_gpd.index.tolist().index(idx_arbitrary)
        print(irow+1, end=" ")
        ds = retrieve_data_func(row, dir_output, time_min=time_min, time_max=time_max, **kwargs)
        if ds is None:
            print("[NODATA] ",end="")
            continue
        
        # assign attrs from station catalog row
        ds = ds.assign_attrs(station_name=row["station_name"],
                             station_id=row["station_id"],
                             station_name_unique=row["station_name_unique"],
                             longitude=row.geometry.x,
                             latitude=row.geometry.x,
                             country=row["country"])
        
        ds["waterlevel"] = ds["waterlevel"].astype("float32")
        _make_hydrotools_consistent(ds)
        
        stat_name = ds.attrs["station_name_unique"]
        file_out = os.path.join(dir_output, f"{stat_name}.nc")
        ds.to_netcdf(file_out)
        del ds
    print()

