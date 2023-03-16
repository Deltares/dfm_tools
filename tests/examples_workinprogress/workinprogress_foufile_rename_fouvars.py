# -*- coding: utf-8 -*-
"""
Created on Thu Mar 16 14:02:47 2023

@author: veenstra
"""

import dfm_tools as dfmt


#file_nc_fou = r'p:\archivedprojects\11203379-005-mwra-updated-bem\03_model\02_final\A72_ntsu0_kzlb2\DFM_OUTPUT_MB_02_fou\MB_02_0000_fou.nc'
file_nc_fou = r'p:\1230882-emodnet_hrsm\GTSMv3.0EMODnet\EMOD_MichaelTUM_yearcomponents\GTSMv4.1_yeartide_2014_2.20.06\output\gtsm_model_0000_fou.nc'

uds = dfmt.open_partitioned_dataset(file_nc_fou)
vars_pd = dfmt.get_ncvarproperties(uds)

import xarray as xr
import xugrid as xu
import pandas as pd
def rename_fouvars(ds:(xr.Dataset,xu.UgridDataset)):
    """
    Rename all fourier variables in a dataset (like mesh2d_fourier033_amp) to a unique name containing quantity/analysistype/tstart/tstop
    
    Parameters
    ----------
    ds : (xr.Dataset,xu.UgridDataset)
        DESCRIPTION.

    Returns
    -------
    ds : TYPE
        DESCRIPTION.

    """
    
    file_freqs = 'https://raw.githubusercontent.com/Deltares/hatyan/main/hatyan/data/data_foreman_frequencies.txt'
    freqs_pd = pd.read_csv(file_freqs,names=['freq','dependents'],delim_whitespace=True,comment='#')
    freqs_pd['angfreq'] = freqs_pd['freq'] * 360 #deg/hr
    
    gridname = ds.grid.name
    list_fouvars = [i for i in ds.data_vars if '_fourier' in i] #water_quality_output and water_quality_stat
    
    rename_dict = {}
    for fouvar in list_fouvars:
        
        #get tstart/tstop
        refdate = pd.Timestamp(str(ds[fouvar].attrs['reference_date_in_yyyymmdd']))
        if hasattr(ds[fouvar],'starttime_fourier_analysis_in_minutes_since_reference_date'):
            tstart_min = ds[fouvar].attrs['starttime_fourier_analysis_in_minutes_since_reference_date']
            tstop_min = ds[fouvar].attrs['stoptime_fourier_analysis_in_minutes_since_reference_date']
        elif hasattr(ds[fouvar],'starttime_min_max_analysis_in_minutes_since_reference_date'):
            tstart_min = ds[fouvar].attrs['starttime_min_max_analysis_in_minutes_since_reference_date']
            tstop_min = ds[fouvar].attrs['stoptime_min_max_analysis_in_minutes_since_reference_date']
        else:
            raise Exception(f'starttime/stoptime attribute not found in fouvar:\n{ds[fouvar].attrs}')
        tstart_str = (refdate + pd.Timedelta(minutes=tstart_min)).strftime('%Y%m%d%H%M%S')
        tstop_str = (refdate + pd.Timedelta(minutes=tstop_min)).strftime('%Y%m%d%H%M%S')
        
        
        if hasattr(ds[fouvar],'frequency_degrees_per_hour'): #for tidal components
            quantity = 'wl' #TODO: where to put this
            quantity = fouvar.split('_')[-1] # amp/phs
            freq = ds[fouvar].attrs['frequency_degrees_per_hour']
            compidx_closestfreq = (freqs_pd['angfreq'] - freq).abs().argmin()
            analysistype = freqs_pd.index[compidx_closestfreq] #M2/NU2
        else: #for all other quantities
            long_name = ds[fouvar].attrs['long_name']
            long_name_noprefix = long_name.split(': ')[1]
            quantity_long = long_name_noprefix.split(',')[0]
            quantity_dict = {'water level':'wl'} #TODO: extend this dictionary
            if not quantity_long in quantity_dict.keys():
                raise Exception(f'quantity_dict does not yet contain quantity for: {quantity_long}')
            quantity = quantity_dict[quantity_long]
            
            analysistype = fouvar.split('_')[-1] #min/max/mean
            if analysistype == 'depth': #ends with min_depth or max_depth
                analysistype = ''.join(fouvar.split('_')[-2:]) #mindepth/maxdepth

        rename_dict[fouvar] = f'{gridname}_{quantity}_{analysistype}_{tstart_str}_{tstop_str}'
    
    ds = ds.rename(rename_dict)
    return ds

uds_renamed = rename_fouvars(uds)
print(uds.data_vars)
print(uds_renamed.data_vars)


