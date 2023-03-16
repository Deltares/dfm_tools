# -*- coding: utf-8 -*-
"""
Created on Thu Mar 16 14:02:47 2023

@author: veenstra
"""

import dfm_tools as dfmt


#file_nc_fou = r'p:\archivedprojects\11203379-005-mwra-updated-bem\03_model\02_final\A72_ntsu0_kzlb2\DFM_OUTPUT_MB_02_fou\MB_02_0000_fou.nc'
file_nc_fou = r'p:\1230882-emodnet_hrsm\GTSMv3.0EMODnet\EMOD_MichaelTUM_yearcomponents\GTSMv4.1_yeartide_2014_2.20.06\output\gtsm_model_0000_fou.nc'
#file_nc_fou = r'c:\DATA\dfm_tools_testdata\DFM_fou_RMM\RMM_dflowfm_0000_fou.nc'

uds = dfmt.open_partitioned_dataset(file_nc_fou)
vars_pd = dfmt.get_ncvarproperties(uds)

import xarray as xr
import xugrid as xu
import pandas as pd
import warnings
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
        fouvar_attrs_lower = {k.lower():v for k,v in ds[fouvar].attrs.items()}
        fouvar_lowerattrs = ds[fouvar].assign_attrs(fouvar_attrs_lower) #to avoid case issues
        
        #quantity 
        long_name = fouvar_lowerattrs.attrs['long_name']
        long_name_noprefix = long_name.split(': ')[1]
        quantity_long = long_name_noprefix.split(',')[0]
        quantity_dict = {'water level':'s1', #dict based on https://svn.oss.deltares.nl/repos/delft3d/trunk/src/engines_gpl/dflowfm/packages/dflowfm_kernel/src/dflowfm_kernel/prepost/fourier_analysis.f90
                         #'energy head':'s1', #TODO: duplicate namfun/dictvalue is not convenient
                         'wind speed':'ws',
                         'U-component of cell-centre velocity':'ux',
                         'V-component of cell-centre velocity':'uy',
                         'U-component velocity, column average':'uxa',
                         'V-component velocity, column average':'uya',
                         'velocity magnitude':'uc',
                         #'':'r1', #TODO: unclear which namfun/dictvalue corresponds (trim(namcon(gdfourier%fconno(ifou))))
                         'velocity':'u1',
                         'unit discharge':'qx',
                         'bed stress':'ta',
                         'freeboard':'fb',
                         'waterdepth_on_ground':'wdog',
                         'volume_on_ground':'vog',
                         'discharge through flow link':'q1',
                         'water level at flow link':'su1',
                         'temperature':'tem', #not clear from fourier_analysis.f90
                         'salt':'sal', #not clear from fourier_analysis.f90
                         }
        if not quantity_long in quantity_dict.keys():
            raise Exception(f'quantity_dict does not yet contain quantity for: {quantity_long}')
        quantity = quantity_dict[quantity_long]
        
        #analysistype
        istidal = False
        if hasattr(fouvar_lowerattrs,'frequency_degrees_per_hour'):
            if fouvar_lowerattrs.attrs['frequency_degrees_per_hour'] > 0: #wl mean with numcyc=0 has frequency attribute (wl min with numcyc=0 does not)
                istidal = True #for tidal components with frequency >0
        
        if istidal: #for tidal analysistype
            tidepart = fouvar.split('_')[-1] # amp/phs
            freq = fouvar_lowerattrs.attrs['frequency_degrees_per_hour']
            compidx_closestfreq = (freqs_pd['angfreq'] - freq).abs().argmin()
            compname = freqs_pd.index[compidx_closestfreq] #M2/NU2
            analysistype = tidepart+compname
            warnings.warn(UserWarning('tidal components found in foufile, matching frequency with online list to get component names, which might go wrong. Also, be aware that v0 and knfac columns from fourier inputfile are not available in fourier output. Recommended is to set them to 0 and compute them in postprocessing.'))
        else: #for all other quantities
            analysistype = fouvar.split('_')[-1] #min/max/mean
            if analysistype == 'depth': #ends with min_depth or max_depth
                analysistype = ''.join(fouvar.split('_')[-2:]) #mindepth/maxdepth
        
        #tstart/tstop
        refdate = pd.Timestamp(str(fouvar_lowerattrs.attrs['reference_date_in_yyyymmdd']))
        if hasattr(fouvar_lowerattrs,'starttime_fourier_analysis_in_minutes_since_reference_date'):
            tstart_min = fouvar_lowerattrs.attrs['starttime_fourier_analysis_in_minutes_since_reference_date']
            tstop_min = fouvar_lowerattrs.attrs['stoptime_fourier_analysis_in_minutes_since_reference_date']
        elif hasattr(fouvar_lowerattrs,'starttime_min_max_analysis_in_minutes_since_reference_date'):
            tstart_min = fouvar_lowerattrs.attrs['starttime_min_max_analysis_in_minutes_since_reference_date']
            tstop_min = fouvar_lowerattrs.attrs['stoptime_min_max_analysis_in_minutes_since_reference_date']
        else:
            raise Exception(f'starttime/stoptime attribute not found in fouvar:\n{fouvar_lowerattrs.attrs}')
        tstart_str = (refdate + pd.Timedelta(minutes=tstart_min)).strftime('%Y%m%d%H%M%S')
        tstop_str = (refdate + pd.Timedelta(minutes=tstop_min)).strftime('%Y%m%d%H%M%S')
        
        rename_dict[fouvar] = f'{gridname}_{quantity}_{analysistype}_{tstart_str}_{tstop_str}'
    
    ds = ds.rename(rename_dict)
    return ds

uds_renamed = rename_fouvars(uds)
print(uds.data_vars)
print(uds_renamed.data_vars)


