# -*- coding: utf-8 -*-
"""
Created on Sat Dec  9 06:52:54 2023

@author: veenstra
"""


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
from erddapy import ERDDAP #pip install erddapy
import requests
import matplotlib.pyplot as plt
plt.close('all')
from zipfile import ZipFile
from io import BytesIO


# example scripts in https://repos.deltares.nl/repos/global_tide_surge_model/trunk/scripts_gtsm5/observationdata/
__all__ = ["ssc_ssh_subset",
           "get_sscid_fromSSCcatalog",
           "cmems_ssh_subset",
           "cmems_ssh_retrieve_data",
           "uhslc_ssh_subset",
           "uhslc_ssh_retrieve_data",
           "ioc_ssh_subset",
           "ioc_ssh_retrieve_data",
           "psmsl_gnssir_ssh_subset",
           "psmsl_gnssir_retrieve_data",
           "gesla3_ssh_subset",
           "gesla3_ssh_retrieve_data",
       ]


def _remove_accents(input_str):
    nfkd_form = unicodedata.normalize('NFKD', input_str)
    only_ascii = nfkd_form.encode('ASCII', 'ignore').decode('ASCII')
    return only_ascii


def _check_groups_valid(groups):
    list_validgroups = ['psmsl','ioc','ptwc','gloss','uhslc']
    if not isinstance(groups, list):
        groups = [groups]
    for groupname in groups:
        if groupname not in list_validgroups:
            raise Exception(f"groupname should be one of {list_validgroups}, '{groupname}' is not valid")


def get_sscid_fromSSCcatalog(group_id, groupname):
    """
    sscid_fromcatalog = get_sscid_fromIOCcatalog(group_id=347, groupname='uhslc')
    fname_fromcatalog = IOC_catalog_pd.loc[sscid_fromcatalog,'station_name_fname']
    """
    group_id = str(group_id)
    SSC_catalog_pd = ssc_ssh_read_catalog()
    _check_groups_valid(groupname)
        
    bool_strinseries = SSC_catalog_pd[groupname].apply(lambda x: group_id in x)
    if bool_strinseries.sum() < 1:
        raise Exception('sscid not found for id %s in group %s'%(group_id, groupname))
    if bool_strinseries.sum() > 1:
        raise Exception('More than 1 sscid found for id %s in group %s:\n%s'%(group_id, groupname, SSC_catalog_pd.loc[bool_strinseries,['name','country', 'geo:lat', 'geo:lon',groupname]]))
    
    sscid = SSC_catalog_pd.loc[bool_strinseries].index[0]
    return sscid


def ssc_ssh_read_catalog():
    """
    The SSC catalog contains e.g. UHSLC and GLOSS ids, this can be used
    for unique station naming across several observation datasets

    info: https://www.ioc-sealevelmonitoring.org/ssc.php
    station list: https://www.ioc-sealevelmonitoring.org/ssc/
    """
    
    url_json = 'https://www.ioc-sealevelmonitoring.org/ssc/service.php?format=json'
    SSC_catalog_pd = pd.read_json(url_json)
    
    #convert all cells with ids to list of strings or NaN
    for colname in ['psmsl','ioc','ptwc','gloss','uhslc','sonel_gps','sonel_tg']:
        SSC_catalog_pd[colname] = SSC_catalog_pd[colname].apply(lambda x: x if isinstance(x,list) else [] if x is np.nan else [x])
    
    #generate station_name_fname (remove all non numeric/letter characters and replace by -, strip '-' from begin/end, set ssc_id as prefix)
    SSC_catalog_pd['station_name_fname'] = SSC_catalog_pd['name'].str.replace('ø','o') # necessary since otherwise these letters are dropped with remove_accents()
    SSC_catalog_pd['station_name_fname'] = SSC_catalog_pd['station_name_fname'].apply(lambda x: _remove_accents(x)) #remove accents from letters
    bool_somethingdropped = (SSC_catalog_pd['station_name_fname'].str.len() != SSC_catalog_pd['name'].str.len())
    if bool_somethingdropped.any():
        raise Exception('lengths mismatch, characters were dropped:\n%s'%(SSC_catalog_pd.loc[bool_somethingdropped,['name','station_name_fname']]))
    SSC_catalog_pd['station_name_fname'] = SSC_catalog_pd['station_name_fname'].apply(lambda x: re.sub("[^0-9a-zA-Z]+", "-", x)) #replace comma and brackets with dash
    SSC_catalog_pd['station_name_fname'] = SSC_catalog_pd['station_name_fname'].str.strip('-') #remove first/last dash from name if present
    col = SSC_catalog_pd.pop('station_name_fname')
    SSC_catalog_pd.insert(3, col.name, col)
    SSC_catalog_pd['station_name_unique'] = SSC_catalog_pd['ssc_id']+'_'+SSC_catalog_pd['station_name_fname']
    
    #set ssc_id as index
    SSC_catalog_pd = SSC_catalog_pd.set_index('ssc_id',drop=False)
    
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
    geom = [Point(x["geo:lon"], x["geo:lat"]) for irow, x in SSC_catalog_pd.iterrows()]
    SSC_catalog_gpd = gpd.GeoDataFrame(data=SSC_catalog_pd, geometry=geom)

    # compare coordinates of station metadata with coordinates of IOC/UHSLC linked stations
    if 0: 
        for station_ssc_id, row in SSC_catalog_gpd.iterrows():
            idx = SSC_catalog_gpd.index.tolist().index(station_ssc_id)
            print(f'station {idx+1} of {len(SSC_catalog_gpd)}: {station_ssc_id}')
            SSC_catalog_pd_stat_IOCUHSLC = SSC_catalog_gpd.loc[station_ssc_id,['ioc','uhslc']]
    
            url_station = f'https://www.ioc-sealevelmonitoring.org/ssc/stationdetails.php?id={station_ssc_id}'
            url_response = urlopen(url_station)
            url_response_read = url_response.read()
            
            station_meta_lon = SSC_catalog_gpd.loc[station_ssc_id].geometry.x
            station_meta_lat = SSC_catalog_gpd.loc[station_ssc_id].geometry.y
    
            if (SSC_catalog_pd_stat_IOCUHSLC.str.len()==0).all(): #skip station if no IOC/UHSLC id present (after retrieval of precise coordinate)
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
            for iR, row in tab_linked_tocheck.iterrows():
                station_check_lat = float(row['Latitude'])
                station_check_lon = float(row['Longitude'])
                station_check_dist = np.sqrt((station_meta_lat-station_check_lat)**2+(station_meta_lon-station_check_lon)**2)
                station_check_dist_all.append(station_check_dist)
                station_check_dict[row['Codes']] = [station_check_dist,station_check_lat,station_check_lon]
            SSC_catalog_gpd.loc[station_ssc_id,'dist_dict'] = [station_check_dict]
            SSC_catalog_gpd.loc[station_ssc_id,'dist_min'] = np.min(station_check_dist_all)
            SSC_catalog_gpd.loc[station_ssc_id,'dist_max'] = np.max(station_check_dist_all)

    return SSC_catalog_gpd


def ssc_ssh_subset(lon_min=-180, lon_max=180, 
                   lat_min=-90, lat_max=90, 
                   # time_min=None, time_max=None
                   groups=None
                   ):
    # TODO: check if min<max
    # TODO: accept None but replace with min/max value from dict. Time min/max are pd.min() and pd.max()
    # TODO: accept partial None (now one None is same as all None)
    
    SSC_catalog_gpd = ssc_ssh_read_catalog()
    
    # spatial subsetting
    SSC_catalog_gpd = SSC_catalog_gpd.clip((lon_min, lat_min, lon_max, lat_max))
    
    # if None not in [time_min, time_max]:
    #     intime_bool = ((SSC_catalog_gpd['time_coverage_start']<time_max) &
    #                    (SSC_catalog_gpd['time_coverage_end']>time_min)
    #                    )
    #     SSC_catalog_gpd = SSC_catalog_gpd.loc[intime_bool].copy()
    
    if groups is not None:
        if not isinstance(groups, list):
            groups = [groups]
        _check_groups_valid(groups)
        bool_ingroup = SSC_catalog_gpd[groups].apply(lambda x: x.str.len()).sum(axis=1)!=0
        SSC_catalog_gpd = SSC_catalog_gpd[bool_ingroup]

    return SSC_catalog_gpd


def cmems_ssh_read_catalog():
    # setup ftp connection
    host = 'my.cmems-du.eu'
    ftp = FTP(host=host)
    username, password = copernicusmarine_credentials()
    ftp.login(user=username, passwd=password)
    # TODO: maybe get path from catalogue?
    ftp.cwd('Core/INSITU_GLO_PHY_SSH_DISCRETE_MY_013_053/cmems_obs-ins_glo_phy-ssh_my_na_PT1H')
    
    # read index
    fname = 'index_history.txt'
    fname_out = os.path.join('.', fname) #TODO: write in cachedir
    with open(fname_out, 'wb') as fp:
        ftp.retrbinary(f'RETR {fname}', fp.write)
    
    with open(fname_out, 'r') as f:
        for line in f:
            if line.startswith('#'):
                header = line
            else:
                break #stop when there are no more #
    colnames = header.strip('#').strip().split(',')
    index_history_pd = pd.read_csv(fname_out,comment='#',names=colnames)
    
    # check equality before we ignore max
    assert (index_history_pd["geospatial_lon_min"] == index_history_pd["geospatial_lon_max"]).all()
    assert (index_history_pd["geospatial_lat_min"] == index_history_pd["geospatial_lat_max"]).all()
    
    # generate geom and geodataframe
    geom = [Point(x["geospatial_lon_min"], x["geospatial_lat_min"]) for irow, x in index_history_pd.iterrows()]
    index_history_gpd = gpd.GeoDataFrame(data=index_history_pd, geometry=geom)
    return index_history_gpd


def cmems_ssh_subset(lon_min=-180, lon_max=180, 
                     lat_min=-90, lat_max=90, 
                     time_min=None, time_max=None):
    # TODO: check if min<max
    # TODO: accept None but replace with min/max value from dict. Time min/max are pd.min() and pd.max()
    # TODO: accept partial None (now one None is same as all None)
    
    index_history_gpd = cmems_ssh_read_catalog()

    # spatial subsetting
    index_history_gpd = index_history_gpd.clip((lon_min, lat_min, lon_max, lat_max))
    
    if None not in [time_min, time_max]:
        intime_bool = ((index_history_gpd['time_coverage_start']<time_max) &
                       (index_history_gpd['time_coverage_end']>time_min)
                       )
        index_history_gpd = index_history_gpd.loc[intime_bool].copy()
    return index_history_gpd


def cmems_ssh_retrieve_data(index_history_gpd, dir_output):
    """
    Retrieve data from FTP
    Can only retrieve entire files, no temporal subsetting possible
    
    # # TODO: the below does not work anymore
    # import copernicus_marine_client as cmc
    # catalogue = cmc.fetch_catalogue()
    # dataset_id = 'cmems_mod_ibi_bgc_anfc_0.027deg-3D_P1D-m'
    # dataset_id = 'cmems_obs-ins_glo_phy-ssh_my_na_PT1H'
    # assert(dataset_id in cmc.get_all_dataset_ids())
    # dataset = catalogue.filter([dataset_id]).products[0].datasets[0]
    # # Object "dataset" can be used to display all the metadata necessary to build a SubsetRequest
    # variable_short = [variable for variable in dataset.variables]
    # variable = [variable for variable in dataset.variables if variable.short_name in ['zooc']][0]
    # coordinates = {coordinate.coordinates_id: (coordinate.minimum_value, coordinate.maximum_value) for coordinate in variable.coordinates}
    # #working example on https://pypi.org/project/copernicus-marine-client/
    
    """
    # setup ftp connection
    host = 'my.cmems-du.eu'
    ftp = FTP(host=host)
    username, password = copernicusmarine_credentials()
    ftp.login(user=username, passwd=password)
    # ftp.pwd()
    
    dir_data = os.path.dirname(index_history_gpd['file_name'].iloc[0]).split(host)[1]
    ftp.cwd(dir_data)
    # retrieve
    print(f"retrieving data for {len(index_history_gpd)} stations:", end=" ")
    ftp.cwd("")
    for idx_arbitrary, row in index_history_gpd.iterrows():
        irow = index_history_gpd.index.tolist().index(idx_arbitrary)
        print(irow+1, end=" ")
        fname = os.path.basename(row.loc["file_name"])
        fname_out = os.path.join(dir_output, fname) #TODO: write in cachedir
        with open(fname_out, 'wb') as fp:
            ftp.retrbinary(f'RETR {fname}', fp.write)
        # TODO: add fname to netcdf attrs
        # ds = xr.open_dataset(fname_out, mode="a")
    print()


def uhslc_ssh_read_catalog():
    uhslc_json = gpd.read_file("https://uhslc.soest.hawaii.edu/data/meta.geojson")
    # col_list = ['name', 'uhslc_id', 'ssc_id', 'gloss_id', 'country', 'country_code',
    #             'fd_span', 'rq_span', 'rq_basin', 'rq_versions', 'geometry']
    
    for drop_col in ["rq_basin", "rq_versions"]:
        if drop_col in uhslc_json.columns:
            uhslc_json = uhslc_json.drop(drop_col, axis=1)
    
    uhslc_json = uhslc_json.set_index('uhslc_id', drop=False)
    
    # shift from 0to360 to -180to180
    from shapely import Point
    geom_shift = [Point(((pnt.x + 180)%360 - 180), pnt.y) for pnt in uhslc_json.geometry]
    uhslc_json.geometry = geom_shift
    return uhslc_json


def uhslc_ssh_subset(lon_min=-180, lon_max=180, 
                     lat_min=-90, lat_max=90, 
                     time_min=None, time_max=None,
                     ):
    # TODO: check if min<max
    # TODO: accept None but replace with min/max value from dict. Time min/max are pd.min() and pd.max()
    # TODO: maybe use min of rqds and max of fast for time subsetting
    # TODO: accept partial None (now one None is same as all None)
    
    uhslc_json = uhslc_ssh_read_catalog()

    # spatial subsetting
    uhslc_json = uhslc_json.clip((lon_min, lat_min, lon_max, lat_max))
    
    uhslc_json['coords_x'] = uhslc_json.geometry.x
    uhslc_json['coords_y'] = uhslc_json.geometry.y
    tstart_rq = uhslc_json["rq_span"].apply(lambda x: x["oldest"])
    tstop_rq = uhslc_json["rq_span"].apply(lambda x: x["latest"])
    tstart_fd = uhslc_json["fd_span"].apply(lambda x: x["oldest"])
    tstop_fd = uhslc_json["fd_span"].apply(lambda x: x["latest"])
    
    # split catalog for rq and fd
    uhslc_json_rq = uhslc_json.loc[~tstart_rq.isnull()].copy()
    uhslc_json_rq["tstart"] = tstart_rq
    uhslc_json_rq["tstop"] = tstop_rq
    uhslc_json_rq["dataset_id"] = "global_hourly_rqds"
    uhslc_json_fd = uhslc_json.loc[~tstart_fd.isnull()].copy()
    uhslc_json_fd["tstart"] = tstart_fd
    uhslc_json_fd["tstop"] = tstop_fd
    uhslc_json_fd["dataset_id"] = "global_hourly_fast"
    
    if None not in [time_min, time_max]:
        intime_bool_rq = ((uhslc_json_rq["tstart"]<time_max) &
                          (uhslc_json_rq["tstop"]>time_min)
                          )
        uhslc_json_rq = uhslc_json_rq.loc[intime_bool_rq].copy()
        intime_bool_fd = ((uhslc_json_fd["tstart"]<time_max) &
                          (uhslc_json_fd["tstop"]>time_min)
                          )
        uhslc_json_fd = uhslc_json_fd.loc[intime_bool_fd].copy()

    return uhslc_json_rq, uhslc_json_fd


def uhslc_ssh_retrieve_data(uhslc_json, dir_output, time_min=None, time_max=None):
    
    # TODO: maybe merge rqds/fast datasets automatically
    # TODO: 5 meter offset?
    
    # docs from https://ioos.github.io/erddapy/ and https://ioos.github.io/erddapy/02-extras-output.html#
    
    server = "https://uhslc.soest.hawaii.edu/erddap"
    e = ERDDAP(server=server, protocol="tabledap", response="nc") #opendap is way slower than nc/csv/html
    
    # search_for = "global_hourly_"
    # url = e.get_search_url(search_for=search_for, response="csv")
    # df = pd.read_csv(url)
    # df["Dataset ID"]
    
    # uhslc_id, station_name, station_country, station_country_code, record_id, gloss_id, ssc_id
    # rowSize, version, decimation_method, reference_code, reference_offset, sea_level
    # e.variables = [
    #     "time", # we cannot drop time, erddappy will raise error (might be from server)
    #     "latitude",
    #     "longitude",
    #     "uhslc_id",
    # ]
    
    uhslc_id_list = uhslc_json["uhslc_id"].tolist()
    for uhslc_id in uhslc_id_list:
        dataset_id = uhslc_json.loc[uhslc_id, "dataset_id"]
        dataset_id_short = dataset_id.split("_")[-1]
        e.dataset_id = dataset_id
    
        e.constraints = {
            # "time>=": pd.Timestamp("2013-01-01"),
            # "time<=": pd.Timestamp("2018-01-01"),
            "uhslc_id=": uhslc_id, #use the geojson to derive the uhslc_id for a specific station
        }
        if time_min is not None:
            e.constraints["time>="] = pd.Timestamp(time_min)
        if time_max is not None:
            e.constraints["time<="] = pd.Timestamp(time_max)
        
        from httpx import HTTPError
        try:
            ds = e.to_xarray()
        except HTTPError:
            print(f"station {uhslc_id} not found in {dataset_id} with "
                  f"timerange tstart={time_min} to tstop={time_max}")
            continue
        
        # change longitude from 0to360 to -180to180
        ds["longitude"] = (ds.longitude + 180)%360 - 180
        ds = ds.sortby("longitude")
        
        # change units to meters
        with xr.set_options(keep_attrs=True):
            assert ds["sea_level"].attrs["units"] == "millimeters"
            ds["sea_level"] = ds["sea_level"]/1000
            ds["sea_level"] = ds["sea_level"].assign_attrs(units="meters")
        
        # set time index
        ds = ds.set_index(obs="time").rename(obs="time")
        ds['time'] = ds.time.dt.round('S') #round to seconds
        
        # write to netcdf file
        file_out = os.path.join(dir_output, f"UHSLC-{uhslc_id}_ssh_{dataset_id_short}.nc")
        ds.to_netcdf(file_out)
        # del ds


def ioc_ssh_read_catalog(drop_uhslc=True, drop_nonutc=True):
    url_json = 'https://www.ioc-sealevelmonitoring.org/service.php?query=stationlist&showall=a' #showall=a means active. showall=all is all known and gives more results than no argument, but also returns nonexistent stations
    resp = requests.get(url_json)
    if resp.status_code==404: #continue to next station if not found
        raise Exception(f'url 404: {resp.text}')    
    resp_json = resp.json()
    IOC_catalog_pd = pd.DataFrame.from_dict(resp_json)
    
    # def remove_accents(input_str):
    #     import unicodedata
    #     nfkd_form = unicodedata.normalize('NFKD', input_str)
    #     only_ascii = nfkd_form.encode('ASCII', 'ignore').decode('ASCII')
    #     return only_ascii
    #generate station_name_fname (remove all non numeric/letter characters and replace by -, strip '-' from begin/end, set ssc_id as prefix)
    # IOC_catalog_pd['station_name_fname'] = IOC_catalog_pd['Location'].str.replace('ø','o').str.replace('’',"'") # necessary since otherwise these letters are dropped with remove_accents()
    # IOC_catalog_pd['station_name_fname'] = IOC_catalog_pd['station_name_fname'].apply(lambda x: remove_accents(x)) #remove accents from letters
    # bool_somethingdropped = (IOC_catalog_pd['station_name_fname'].str.len() != IOC_catalog_pd['Location'].str.len())
    # if bool_somethingdropped.any():
    #     raise Exception('lengths mismatch, characters were dropped:\n%s'%(IOC_catalog_pd.loc[bool_somethingdropped,['Location','station_name_fname']]))
    # IOC_catalog_pd['station_name_fname'] = IOC_catalog_pd['station_name_fname'].apply(lambda x: re.sub("[^0-9a-zA-Z]+", "-", x)) #replace comma and brackets with dash
    # IOC_catalog_pd['station_name_fname'] = IOC_catalog_pd['station_name_fname'].str.strip('-') #remove first/last dash from name if present
    
    # IOC_catalog_pd['station_lat'] = IOC_catalog_pd['Lat']
    # IOC_catalog_pd['station_lon'] = IOC_catalog_pd['Lon']
    # IOC_catalog_pd['station_name_unique'] = 'IOC-'+IOC_catalog_pd['Code']+'_'+IOC_catalog_pd['station_name_fname']
    
    #set ssc_id as index
    IOC_catalog_pd = IOC_catalog_pd.set_index('Code',drop=False)
    
    geom = [Point(x["lon"], x["lat"]) for irow, x in IOC_catalog_pd.iterrows()]
    ioc_catalog_gpd = gpd.GeoDataFrame(data=IOC_catalog_pd, geometry=geom)
    
    if drop_uhslc:
        # filter non-UHSLC stations
        # TODO: do this without pkl file
        SSC_catalog_pd = pd.read_pickle(r'p:\1230882-emodnet_hrsm\global_tide_surge_model\trunk\observation_locations_org\SSC-catalog.pkl')
        ioc_catalog_gpd['UHSLC'] = False
        ioc_catalog_gpd['inSSClist'] = False
        for ioc_code in ioc_catalog_gpd['Code']:
            iocstat_uhslcid = SSC_catalog_pd.loc[SSC_catalog_pd['ioc'].apply(lambda x: ioc_code in x),'uhslc']
            if len(iocstat_uhslcid)==0: #ioc_code not in the SSC list
                continue
            else:
                ioc_catalog_gpd.loc[ioc_code,'inSSClist'] = True
                if len(iocstat_uhslcid.iloc[0])==0: #station contains no UHSLC code
                    continue
                else:
                    ioc_catalog_gpd.loc[ioc_code,'UHSLC'] = True
        ioc_catalog_gpd = ioc_catalog_gpd.loc[~ioc_catalog_gpd['UHSLC']]
    
    if drop_nonutc:
        # filter out all non-UTC stations
        # TODO: maybe implement UTCOffset fix, but not important enough now
        ioc_catalog_gpd = ioc_catalog_gpd.loc[ioc_catalog_gpd['UTCOffset']==0]
        
    return ioc_catalog_gpd


def ioc_ssh_subset(lon_min=-180, lon_max=180, 
                   lat_min=-90, lat_max=90, 
                   #time_min=None, time_max=None,
                   drop_uhslc=True, drop_nonutc=True,
                   ):
    # TODO: check if min<max
    # TODO: accept None but replace with min/max value from dict. Time min/max are pd.min() and pd.max()
    # TODO: accept partial None (now one None is same as all None)
    
    ioc_catalog_gpd = ioc_ssh_read_catalog(drop_uhslc=drop_uhslc, drop_nonutc=drop_nonutc)

    # spatial subsetting
    ioc_catalog_gpd = ioc_catalog_gpd.clip((lon_min, lat_min, lon_max, lat_max))
    
    # if None not in [time_min, time_max]:
    #     intime_bool = ((ioc_catalog_gpd['time_coverage_start']<time_max) &
    #                    (ioc_catalog_gpd['time_coverage_end']>time_min)
    #                    )
    #     ioc_catalog_gpd = ioc_catalog_gpd.loc[intime_bool].copy()
    return ioc_catalog_gpd

# bool_countries = ((ioc_catalog_gpd['country']=='RUS') | #IOC Russia: 22
#                   (ioc_catalog_gpd['country']=='IND') | #IOC India :  8
#                   # (ioc_catalog_gpd['countryname']=='Sri Lanka') | #IOC Sri Lanka: 3
#                   # (ioc_catalog_gpd['countryname']=='Myanmar') | #IOC Myanmar: 1
#                   # (ioc_catalog_gpd['countryname']=='Indonesia') | #IOC Indonesia: 1
#                   # (ioc_catalog_gpd['countryname']=='Malaysia') | #IOC Malaysia: 5
#                   # (ioc_catalog_gpd['countryname']=='Vanuatu') | #IOC Vanuatu: 3
#                   # (ioc_catalog_gpd['countryname']=='Solomon Islands') | #IOC Solomon Islands: 2
#                   # (ioc_catalog_gpd['countryname']=='China') | #IOC China: 4
#                   # (ioc_catalog_gpd['countryname']=='Hong Kong - China') | #IOC Hong Kong: 1
#                   # (ioc_catalog_gpd['countryname']=='Madagascar') | #IOC Madagascar: 1
#                   # (ioc_catalog_gpd['countryname']=='Comores') | #IOC Comores: 1
#                   # (ioc_catalog_gpd['countryname']=='Mauritius') | #IOC Mauritius: 1
#                   # (ioc_catalog_gpd['countryname']=='Mauritania') |
#                   # (ioc_catalog_gpd['countryname']=='Gambia') |
#                   # (ioc_catalog_gpd['countryname']=='Ghana') |
#                   # (ioc_catalog_gpd['countryname']=='Antigua, Barbuda & Redonda') | #IOC Antigua, Barbuda & Redonda: 3
#                   # (ioc_catalog_gpd['countryname']=='St. Kitts & Nevis') | #IOC St. Kitts & Nevis: 1
#                   # (ioc_catalog_gpd['countryname']=='Colombia') | #IOC Colombia: 8
#                   # (ioc_catalog_gpd['countryname']=='Aruba; Nederland') | #IOC Aruba: 1
#                   # (ioc_catalog_gpd['countryname']=='Anguilla') | #IOC Anguilla: 1
#                   # (ioc_catalog_gpd['countryname']=='Trinidad & Tobago') | #IOC Trinidad & Tobago: 5
#                   # (ioc_catalog_gpd['countryname']=='Dominica island') | #IOC Dominica island: 3
#                   # (ioc_catalog_gpd['countryname']=='Barbados') | #IOC Barbados: 1
#                   # (ioc_catalog_gpd['countryname']=='Grenada') | #IOC Grenada: 1
#                   # (ioc_catalog_gpd['countryname']=='Saint Lucia') | #IOC Saint Lucia: 4
#                   # (ioc_catalog_gpd['countryname']=='Saint Vincent & Grenadines') | #IOC Saint Vincent & Grenadines: 1
#                   (ioc_catalog_gpd['country']=='TUR') | #IOC Turkey: 22
#                   (ioc_catalog_gpd['country']=='ISR') | #IOC Israel:  7
#                   (ioc_catalog_gpd['country']=='CYP') | #IOC Cyprus:  6
#                   (ioc_catalog_gpd['countryname']=='Greece') | #IOC Greece:  23
#                   # (ioc_catalog_gpd['countryname']=='Italy') | #IOC Italy: 42
#                   (ioc_catalog_gpd['country']=='KRS') | #IOC Korea :  1
#                   (ioc_catalog_gpd['country']=='HRV')#| #IOC Croatia: 3
#                   )


def ioc_ssh_retrieve_data(ioc_catalog_gpd, dir_output, date_min, date_max, subset_hourly=False):
    # https://www.ioc-sealevelmonitoring.org/service.php?query=help
    print(f"retrieving data for {len(ioc_catalog_gpd)} stations.")
    for station_code, row in ioc_catalog_gpd.iterrows():
        results_list = []
        period_range = pd.period_range(date_min,date_max,freq="1M")
        print(f'retrieving data for {station_code}: ', end="")
        for date in period_range:
            date_str = str(date)
            print(f'{date_str} ', end="")
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
                print("[NODATA] ", end="")
                continue
            resp_json = resp.json()
            if 'error' in resp.json()[0].keys():
                raise Exception(resp.text)
            data_pd_one = pd.DataFrame.from_dict(resp_json)
            results_list.append(data_pd_one)
        print()
        
        if len(results_list)==0:
            # continue with next station if no data present in entire period
            continue
        
        # convert to xarray
        data_pd_all = pd.concat(results_list)
        data_pd_all = data_pd_all.rename({"stime":"time"},axis=1)
        data_pd_all = pd.DataFrame({'slev':data_pd_all['slevel'].values},
                                index=pd.to_datetime(data_pd_all['time']))
        data_pd_all = data_pd_all[~data_pd_all.index.duplicated(keep='last')]
        if subset_hourly:
            data_pd_all = data_pd_all.loc[data_pd_all.index.minute==0]
        ds = data_pd_all.to_xarray()
        ds = ds.assign_attrs(station_id=row["Location"],
                             station_code=row["Code"],
                             longitude=row["Lon"],
                             latitude=row["Lat"],
                             country_code=row["country"],
                             )
        
        # write to netcdf file
        file_out = os.path.join(dir_output, f"IOC-{station_code}_ssh.nc")
        ds.to_netcdf(file_out)


def psmsl_gnssir_ssh_read_catalog():
    # https://psmsl.org/data/gnssir/metadatainfo.php
    # https://psmsl.org/data/gnssir/useful_files.php
    # url = "https://psmsl.org/data/gnssir/data/maplayers/good_sites.json"
    #TODO: use only good_sites instead of all
    #TODO: filter on time
    url = "https://psmsl.org/data/gnssir/data/sites.json"
    station_list = pd.read_json(url).T
    
    # generate geom and geodataframe
    geom = [Point(x["Longitude"], x["Latitude"]) for irow, x in station_list.iterrows()]
    station_list_gpd = gpd.GeoDataFrame(data=station_list, geometry=geom)
    return station_list_gpd


def psmsl_gnssir_ssh_read_catalog_gettimes(station_list_gpd):
    # the catalog json does not contain time ranges so we derive it from daily csv files
    station_list_gpd["tstart"] = pd.NaT
    station_list_gpd["tstop"] = pd.NaT
    print(f"retrieving time extents for {len(station_list_gpd)} stations:", end=" ")
    for station_id in station_list_gpd.index:
        irow = station_list_gpd.index.tolist().index(station_id)
        print(irow+1, end=" ")
        url = f"https://psmsl.org/data/gnssir/data/daily/{station_id}_daily.csv"
        data_daily = pd.read_csv(url)
        station_list_gpd.loc[station_id, "tstart"] = pd.Timestamp(data_daily["time"].iloc[0])
        station_list_gpd.loc[station_id, "tstop"] = pd.Timestamp(data_daily["time"].iloc[-1])
    print()

    return station_list_gpd


def psmsl_gnssir_ssh_subset(lon_min=-180, lon_max=180, 
                            lat_min=-90, lat_max=90, 
                            time_min=None, time_max=None):
    # TODO: check if min<max
    # TODO: accept None but replace with min/max value from dict. Time min/max are pd.min() and pd.max()
    # TODO: accept partial None (now one None is same as all None)
    
    station_list_gpd = psmsl_gnssir_ssh_read_catalog()

    # spatial subsetting
    station_list_gpd = station_list_gpd.clip((lon_min, lat_min, lon_max, lat_max))
    
    if None not in [time_min, time_max]:
        station_list_gpd = psmsl_gnssir_ssh_read_catalog_gettimes(station_list_gpd)
        intime_bool = ((station_list_gpd['tstart']<time_max) &
                       (station_list_gpd['tstop']>time_min)
                       )
        station_list_gpd = station_list_gpd.loc[intime_bool].copy()
    return station_list_gpd


def psmsl_gnssir_retrieve_data(station_list_gpd, dir_output):
    # https://psmsl.org/data/gnssir/gnssir_daily_means.html
    # https://psmsl.org/data/gnssir/gnssir_example.html (also contains IOC retrieval example)
    
    print(f"retrieving data for {len(station_list_gpd)} stations:", end=" ")
    for station_id, row in station_list_gpd.iterrows():
        irow = station_list_gpd.index.tolist().index(station_id)
        print(irow+1, end=" ")
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
        ds = ds.assign_attrs(station_id=row["Name"],
                             station_code=row["Code"],
                             longitude=row["Longitude"],
                             latitude=row["Latitude"],
                             country_code=row["CountryCode"],
                             )
        file_out = os.path.join(dir_output, f"PSMSL-GNSSIR-{station_id}_ssh.nc")
        ds.to_netcdf(file_out)
        del ds
    print()


def gesla3_ssh_read_catalog(file_gesla3_meta=None):
    if file_gesla3_meta is None:
        # from https://www.icloud.com/iclouddrive/01a8u37HiumNKbg6CpQUEA7-A#GESLA3_ALL_2
        # linked on https://gesla787883612.wordpress.com/downloads/
        file_gesla3_meta = r"p:\1230882-emodnet_hrsm\data\GESLA3\GESLA3_ALL 2.csv"
    
    station_list_pd = pd.read_csv(file_gesla3_meta)
    station_list_pd.columns = [c.replace(" ", "_").replace("(", "").replace(")", "").replace("/", "_").lower() for c in station_list_pd.columns]
    station_list_pd = station_list_pd.set_index('file_name', drop=False)
    station_list_pd["start_date_time"] = pd.to_datetime(station_list_pd["start_date_time"])
    station_list_pd["end_date_time"] = pd.to_datetime(station_list_pd["end_date_time"])
    
    # generate geom and geodataframe
    geom = [Point(x["longitude"], x["latitude"]) for irow, x in station_list_pd.iterrows()]
    station_list_gpd = gpd.GeoDataFrame(data=station_list_pd, geometry=geom)
    return station_list_gpd


def gesla3_ssh_subset(lon_min=-180, lon_max=180, 
                      lat_min=-90, lat_max=90, 
                      time_min=None, time_max=None,
                      file_gesla3_meta=None):
    # TODO: check if min<max
    # TODO: accept None but replace with min/max value from dict. Time min/max are pd.min() and pd.max()
    # TODO: accept partial None (now one None is same as all None)
    
    station_list_gpd = gesla3_ssh_read_catalog(file_gesla3_meta=file_gesla3_meta)

    # spatial subsetting
    station_list_gpd = station_list_gpd.clip((lon_min, lat_min, lon_max, lat_max))
    
    if None not in [time_min, time_max]:
        intime_bool = ((station_list_gpd['start_date_time']<time_max) &
                       (station_list_gpd['end_date_time']>time_min)
                       )
        station_list_gpd = station_list_gpd.loc[intime_bool].copy()
    return station_list_gpd


def gesla3_ssh_retrieve_data(gesla3_stations_gpd, dir_output,
                             time_min=None, time_max=None,
                             dir_gesla3_data=None, file_gesla3_meta=None,
                             ):
    
    if dir_gesla3_data is None:
        # https://www.icloud.com/iclouddrive/0tHXOLCgBBjgmpHecFsfBXLag#GESLA3
        # linked on https://gesla787883612.wordpress.com/downloads/
        dir_gesla3_data = r"p:\1230882-emodnet_hrsm\data\GESLA3\GESLA3.0_ALL"
    
    gesla_stations_gpd = gesla3_ssh_subset(file_gesla3_meta=file_gesla3_meta)
    
    for file_gesla, row in gesla3_stations_gpd.iterrows():
        irow = gesla3_stations_gpd.index.tolist().index(file_gesla)
        print(f'processing {irow+1} of {len(gesla3_stations_gpd)}: {file_gesla}')
        
        file_gesla_org = os.path.join(dir_gesla3_data, file_gesla)
        file_gesla_nc = os.path.join(dir_output, f'{file_gesla}.nc')
        
        with open(file_gesla_org, "r") as f:
            data = pd.read_csv(
                f,
                comment='#',
                names=["date", "time", "sea_level", "qc_flag", "use_flag"],
                delim_whitespace = True,
                parse_dates=[[0, 1]],
                index_col=0)
        
        # clean up time duplicates
        data.index.name = 'time'
        bool_duplicate = data.index.duplicated()
        if bool_duplicate.sum() > 0:
            data = data.loc[~bool_duplicate]
            print(f"[{bool_duplicate.sum()} duplicate timestamps removed]")
        
        # convert to xarray and add metadata
        meta_sel = gesla_stations_gpd.loc[file_gesla].copy()
        meta_sel["start_date_time"] = str(meta_sel["start_date_time"])
        meta_sel["end_date_time"] = str(meta_sel["end_date_time"])
        meta_sel = meta_sel.drop('geometry')
        ds = data.to_xarray().assign_attrs(meta_sel.to_dict())
        ds['site_name'] = xr.DataArray([meta_sel.site_name], dims=('stations'))
        ds['latitude'] = xr.DataArray([meta_sel.latitude], dims=('stations'))
        ds['longitude'] = xr.DataArray([meta_sel.longitude], dims=('stations'))
        
        # subset time
        ds = ds.sel(time=slice(time_min, time_max))
        if len(ds.time) == 0:
            print("[NODATA]")
            continue
        
        # filter bad quality data
        ds = ds.where(ds.qc_flag==1)
        
        # write to file
        ds.to_netcdf(file_gesla_nc)
        del ds

