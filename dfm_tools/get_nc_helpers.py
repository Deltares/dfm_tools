# -*- coding: utf-8 -*-
"""
dfm_tools are post-processing tools for Delft3D FM
Copyright (C) 2020 Deltares. All rights reserved.

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  if not, see <http://www.gnu.org/licenses/>.

All names, logos, and references to "Deltares" are registered trademarks of
Stichting Deltares and remain full property of Stichting Deltares at all times.
All rights reserved.


INFORMATION
This script is part of dfm_tools: https://github.com/openearth/dfm_tools
Check the README.rst on github for other available functions
Check the tests folder on github for example scripts (this is the dfm_tools pytest testbank)
Check the pptx and example figures in (created by the testbank): N:/Deltabox/Bulletin/veenstra/info dfm_tools

Created on Fri Feb 14 12:43:19 2020

@author: veenstra

helper functions for functions in get_nc.py
"""

import xarray as xr
import xugrid as xu
import pandas as pd
import warnings


def get_ncvarproperties(data_xr):
    if not isinstance(data_xr,(xr.Dataset,xu.UgridDataset)):
        raise TypeError('data_xr should be of type xr.Dataset or xu.UgridDataset')
    
    nc_varkeys = data_xr.variables.mapping.keys()
    
    list_varattrs_pd = []
    for varkey in nc_varkeys:
        varattrs_pd = pd.DataFrame({varkey:data_xr.variables.mapping[varkey].attrs}).T
        varattrs_pd[['shape','dimensions','dtype']] = 3*[''] #set dtype as str (float will raise an error when putting tuple in there)
        varattrs_pd.at[varkey,'shape'] = data_xr[varkey].shape
        varattrs_pd.at[varkey,'dimensions'] = data_xr.variables[varkey].dims
        varattrs_pd.at[varkey,'dtype'] = data_xr.variables[varkey].dtype
        list_varattrs_pd.append(varattrs_pd)
    
    vars_pd = pd.concat(list_varattrs_pd,axis=0)
    vars_pd[vars_pd.isnull()] = '' #avoid nan values
    
    data_xr.close()

    return vars_pd


def rename_waqvars(ds:(xr.Dataset,xu.UgridDataset)):
    """
    Rename all water quality variables in a dataset (like mesh2d_water_quality_output_24) to their long_name attribute (like mesh2d_DOscore)
    
    Parameters
    ----------
    ds : (xr.Dataset,xu.UgridDataset)
        DESCRIPTION.

    Returns
    -------
    ds : TYPE
        DESCRIPTION.

    """
    #TODO: results also in variable "mesh2d_Water quality mass balance areas" (with spaces), report in FM issue (remove spaces from long_name attr)
    
    if hasattr(ds,'grid'): #append gridname (e.g. mesh2d) in case of mapfile
        varn_prepend = f'{ds.grid.name}_'
    else:
        varn_prepend = ''
    list_waqvars = [i for i in ds.data_vars if 'water_quality_' in i] #water_quality_output and water_quality_stat
    rename_dict = {waqvar:varn_prepend+ds[waqvar].attrs['long_name'] for waqvar in list_waqvars}
    
    if len(rename_dict) == 0: #early return to silence "FutureWarning: The default dtype for empty Series will be 'object' instead of 'float64' in a future version. Specify a dtype explicitly to silence this warning."
        return ds
    
    #prevent renaming duplicate long_names
    rename_pd = pd.Series(rename_dict)
    if rename_pd.duplicated().sum():
        duplicated_pd = rename_pd.loc[rename_pd.duplicated(keep=False)]
        print(UserWarning(f'duplicate long_name attributes found with dfmt.rename_waqvars(), renaming only first variable:\n{duplicated_pd}'))
        rename_dict = rename_pd.loc[~rename_pd.duplicated()].to_dict()
    
    ds = ds.rename(rename_dict)
    return ds


def rename_fouvars(ds:(xr.Dataset,xu.UgridDataset), drop_tidal_times:bool = True):
    """
    Rename all fourier variables in a dataset (like mesh2d_fourier033_amp) to a unique name containing gridname/quantity/analysistype/tstart/tstop
    
    Parameters
    ----------
    ds : (xr.Dataset,xu.UgridDataset)
        DESCRIPTION.

    Returns
    -------
    ds : TYPE
        DESCRIPTION.

    """
    
    file_freqs = 'https://raw.githubusercontent.com/Deltares/hatyan/main/hatyan/data/data_foreman_frequencies.txt' #TODO: fix hatyan dependency (MSQM and M1 were also added, but file is not used by hatyan, so might disappear one day)
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
                         'temperature':'tem', #not clear from fourier_analysis.f90, ct in user manual C.13
                         'salt':'sal', #not clear from fourier_analysis.f90, cs in user manual C.13
                         }
        if quantity_long not in quantity_dict.keys():
            raise KeyError(f'quantity_dict does not yet contain quantity for: {quantity_long}')
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
            warnings.warn(UserWarning('tidal components found in foufile, matching frequency with online list to get component names, which might go wrong. Also, be aware that v0 and knfac columns from fourier inputfile are not available in fourier output, so it is not clear whether to correct for these.'))
        else: #for all other quantities
            fouvar_splitted = fouvar.split('_')
            analysistype = ''.join(fouvar_splitted[2:]) #min/max/mean. min_depth/max_depth etc are converted to mindepth/maxdepth
            
        #tstart/tstop
        refdate = pd.Timestamp(str(fouvar_lowerattrs.attrs['reference_date_in_yyyymmdd']))
        if hasattr(fouvar_lowerattrs,'starttime_fourier_analysis_in_minutes_since_reference_date'):
            tstart_min = fouvar_lowerattrs.attrs['starttime_fourier_analysis_in_minutes_since_reference_date']
            tstop_min = fouvar_lowerattrs.attrs['stoptime_fourier_analysis_in_minutes_since_reference_date']
        elif hasattr(fouvar_lowerattrs,'starttime_min_max_analysis_in_minutes_since_reference_date'):
            tstart_min = fouvar_lowerattrs.attrs['starttime_min_max_analysis_in_minutes_since_reference_date']
            tstop_min = fouvar_lowerattrs.attrs['stoptime_min_max_analysis_in_minutes_since_reference_date']
        else:
            raise AttributeError(f'starttime/stoptime attribute not found in fouvar:\n{fouvar_lowerattrs.attrs}')
        tstart_str = (refdate + pd.Timedelta(minutes=tstart_min)).strftime('%Y%m%d%H%M%S')
        tstop_str = (refdate + pd.Timedelta(minutes=tstop_min)).strftime('%Y%m%d%H%M%S')
        
        rename_dict[fouvar] = f'{gridname}_{quantity}_{analysistype}_{tstart_str}_{tstop_str}'
        if istidal and drop_tidal_times:
            rename_dict[fouvar] = f'{gridname}_{quantity}_{analysistype}' #TODO: might cause conflicting variable names if one component is analysed for multiple periods or if component is not defined in frequency list. Add duplicate check like rename_waqvars() that provides some basic info for debugging.
    
    ds = ds.rename(rename_dict)
    return ds

