import os
import glob
import datetime as dt
import numpy as np
import pandas as pd
import xarray as xr
import xugrid as xu
from pathlib import Path
from scipy.spatial import KDTree
import logging
import hydrolib.core.dflowfm as hcdfm
import geopandas

from dfm_tools.hydrolib_helpers import (Dataset_to_TimeSeries, 
                                        Dataset_to_T3D, 
                                        Dataset_to_Astronomic, 
                                        PolyFile_to_geodataframe_points, 
                                        get_ncbnd_construct,
                                        da_from_gdf_points,
                                        maybe_convert_fews_to_dfmt,
                                        validate_polyline_names)
from dfm_tools.errors import OutOfRangeError

__all__ = ["get_conversion_dict",
           "interpolate_tide_to_bc",
           "interpolate_tide_to_plipoints",
           "interp_regularnc_to_plipointsDataset",
           "interp_uds_to_plipoints",
           "interp_hisnc_to_plipoints",
           "plipointsDataset_to_ForcingModel",
    ]

logger = logging.getLogger(__name__)

if os.name == "nt":
    # windows drive letter should include trailing slash
    # https://github.com/Deltares/dfm_tools/issues/1084
    PDRIVE = "p:/"
else:
    PDRIVE = "/p"


def get_conversion_dict(ncvarname_updates={}):
    
    """
    get the conversion_dict, optionally with updated ncvarnames
    conversion_dict.keys() are the dflowfm quantities and the ncvarname the corresponding netcdf variable name/key
    alternative ncvarnames can be supplied via ncvarname_updates, e.g. get_conversion_dict(ncvarname_updates={'temperaturebnd':'tos'})
    
    interpolate_nc_to_bc() renames netcdf variable like this:
    data_xr = data_xr.rename({ncvarname:quantity})
    
    for CMCC:
    conversion_dict = {
    # conversion is phyc in mol/m3 to newvar in g/m3
    'tracerbndOXY'        : {'ncvarname': 'o2',          'unit': 'g/m3', 'conversion': 32.0 },
    'tracerbndNO3'        : {'ncvarname': 'no3',         'unit': 'g/m3', 'conversion': 14.0 },
    'tracerbndPO4'        : {'ncvarname': 'po4',         'unit': 'g/m3', 'conversion': 30.97 },
    'tracerbndSi'         : {'ncvarname': 'si',          'unit': 'g/m3', 'conversion': 28.08},
    'tracerbndPON1'       : {'ncvarname': 'phyc',        'unit': 'g/m3', 'conversion': 14.0}, # Caution: this empirical relation might not be applicable to your use case
    'tracerbndPOP1'       : {'ncvarname': 'phyc',        'unit': 'g/m3', 'conversion': 30.97}, # Caution: this empirical relation might not be applicable to your use case
    'tracerbndPOC1'       : {'ncvarname': 'phyc',        'unit': 'g/m3', 'conversion': 14.0 * (106 / 16)}, # Caution: this empirical relation might not be applicable to your use case
    'salinitybnd'         : {'ncvarname': 'sos'},         #'1e-3'
    'temperaturebnd'      : {'ncvarname': 'tos'},         #'degC'
    'ux'                  : {'ncvarname': 'uo'},          #'m/s'
    'uy'                  : {'ncvarname': 'vo'},          #'m/s'
    'waterlevelbnd'       : {'ncvarname': 'zos'},         #'m' #steric
    'tide'                : {'ncvarname': ''},            #'m' #tide (dummy entry)
    }
    
    
    The below clarification about the CMEMS conversion factors in this function is by Jos van Gils.
    The use of boundary conditions from CMEMS requires some conversion (scaling). 
    In the first place, there is the unit conversion from mmol/m3 in CMEMS to g/m3 in DFM. 
    This involves a factor M/1000, where M is the molar weight of the constituent in question: 
    M = 12 for constituents expressed in gC/m3, 14 for gN/m3, 30.97 for gP/m3, 28.08 for gSi/m3 and 32 for dissolved oxygen. 
    
    Next there is the issue that CMEMS does not provide all state variables in the water 
    quality model. Dead organic matter, particulate and dissolved, are missing. 
    In practice, we use phytoplankton biomass from CMEMS (phyc) to estimate these missing variables. 
    For particulate organic matter (POM), often a carbon-based ratio of 2.0 to phytoplankton biomass is used.
    PON and POP are derived from POC by classical Redfield ratios (C:N:P = 106:16:1) or revised Redfield ratios (C:N:P = 117:14:1). 
    For Opal (biogenic silica from diatoms), the example here is based on a factor 0.5 expressing the share of 
    diatoms in phytoplankton and a 0.13 Si:C stoichiometric ratio by Brzezinski (1985).
    For dissolved organic matter (DOM) there is limited experience. This
    conversion_dict uses the approach used for the “Bays Eutrophication Model” (Massachusetts Bay, 
    https://www.mwra.com/harbor/enquad/pdf/2021-02.pdf). This relies on local field data.
    
    It is noted that such scale factors are anyhow inaccurate and typically water-system-specific. 
    Therefore, the values included in this conversion_dict should be considered indicative only and 
    should by no means be seen as a default approach. Ideally, the scale factors to link missing 
    states to available CMEMS variables are derived from local field data. 
    A validation is strongly recommended, for example by a comparison of simulated and 
    measured concentrations of total N and total P in the study area. 
    Similar considerations would hold for other data sources than CMEMS.
    """
    
    # conversion_dict for CMEMS
    conversion_dict = { # conversion is phyc in mmol/m3 to newvar in g/m3
                        'tracerbndOXY'        : {'ncvarname': 'o2',          'unit': 'g/m3', 'conversion': 32/1000},
                        'tracerbndNO3'        : {'ncvarname': 'no3',         'unit': 'g/m3', 'conversion': 14/1000},
                        'tracerbndPO4'        : {'ncvarname': 'po4',         'unit': 'g/m3', 'conversion': 30.97/1000},
                        'tracerbndSi'         : {'ncvarname': 'si',          'unit': 'g/m3', 'conversion': 28.08/1000},
                        'tracerbndGreen'      : {'ncvarname': 'phyc',        'unit': 'g/m3', 'conversion': 12/1000},
                        'tracerbndPON1'       : {'ncvarname': 'phyc',        'unit': 'g/m3', 'conversion': 2 * 16/106 * 14/1000}, # Caution: this empirical relation might not be applicable to your use case
                        'tracerbndPOP1'       : {'ncvarname': 'phyc',        'unit': 'g/m3', 'conversion': 2 * 30.97 / (106 * 1000)}, # Caution: this empirical relation might not be applicable to your use case
                        'tracerbndPOC1'       : {'ncvarname': 'phyc',        'unit': 'g/m3', 'conversion': 2 * 12/1000}, # Caution: this empirical relation might not be applicable to your use case
                        'tracerbndDON'        : {'ncvarname': 'phyc',        'unit': 'g/m3', 'conversion': 3.24 * 2 * 16/106 * 14/1000}, # Caution: this empirical relation might not be applicable to your use case
                        'tracerbndDOP'        : {'ncvarname': 'phyc',        'unit': 'g/m3', 'conversion': 2 * 30.97 / (106 * 1000)}, # Caution: this empirical relation might not be applicable to your use case
                        'tracerbndDOC'        : {'ncvarname': 'phyc',        'unit': 'g/m3', 'conversion': (199 / 20) * 3.24 * 2 * 16 * 12 / (106 * 1000)}, # Caution: this empirical relation might not be applicable to your use case
                        'tracerbndOpal'       : {'ncvarname': 'phyc',        'unit': 'g/m3', 'conversion': 0.5 * 0.13 * 28.08/1000}, # Caution: this empirical relation might not be applicable to your use case
                        'salinitybnd'         : {'ncvarname': 'so'},          #'1e-3'
                        'temperaturebnd'      : {'ncvarname': 'thetao'},      #'degC'
                        'ux'                  : {'ncvarname': 'uo'},          #'m/s'
                        'uy'                  : {'ncvarname': 'vo'},          #'m/s'
                        'waterlevelbnd'       : {'ncvarname': 'zos'},         #'m' #steric
                        'tide'                : {'ncvarname': ''},            #'m' #tide (dummy entry)
                        }
    #do updates
    for k,v in ncvarname_updates.items():
        conversion_dict[k]['ncvarname'] = v
    return conversion_dict


def ext_add_boundary_object_per_polyline(ext_new:hcdfm.ExtModel, boundary_object:hcdfm.Boundary):
    """
    If the polyfile contains multiple polylines, only the first 
    one will be used by DFLOW-FM for the boundary conditions.
    This function will save each polyline as a separate polyfile 
    and add duplicate boundary entries for each polyline to the extfile.
    """
    file_pli_main = boundary_object.locationfile.filepath
    polyfile_obj = hcdfm.PolyFile(file_pli_main)
    validate_polyline_names(polyfile_obj)
    dir_output = os.path.dirname(file_pli_main)
    for polyline_obj in polyfile_obj.objects:
        if len(polyfile_obj.objects) > 1:
            # copy to avoid location file to be the same for all duplicated boundary objects
            boundary_object = boundary_object.copy()
            # create a polyfile with a single plifile
            polyfile_oneline_obj = hcdfm.PolyFile(objects=[polyline_obj])
            file_pli_oneline = os.path.join(dir_output, f"{polyline_obj.metadata.name}.pli")
            if str(file_pli_oneline) == str(file_pli_main):
                raise ValueError("The names of one of the polylines in the polyfile is the "
                                 "same as the polyfilename, resulting in the same filename")
            polyfile_oneline_obj.save(file_pli_oneline)
            # update locationfile to the plifile with a single polyline
            boundary_object.locationfile = file_pli_oneline
        # add boundary to extmodel
        ext_new.boundary.append(boundary_object)


def interpolate_tide_to_bc(ext_new: hcdfm.ExtModel, tidemodel, file_pli, component_list=None, load=True):
    # read polyfile as geodataframe
    polyfile_obj = hcdfm.PolyFile(file_pli)
    gdf_points = PolyFile_to_geodataframe_points(polyfile_obj)
    
    # interpolate tidal components to plipoints
    data_interp = interpolate_tide_to_plipoints(tidemodel=tidemodel, gdf_points=gdf_points, 
                                                component_list=component_list, load=load)
    
    # convert interpolated xarray.Dataset to hydrolib ForcingModel and save as bc file
    ForcingModel_object = plipointsDataset_to_ForcingModel(plipointsDataset=data_interp)
    dir_output = os.path.dirname(file_pli)
    file_pli_name = polyfile_obj.filepath.stem
    file_bc_out = os.path.join(dir_output,f'tide_{tidemodel}_{file_pli_name}.bc')
    ForcingModel_object.save(filepath=file_bc_out)
    
    # generate hydrolib-core Boundary object to be appended to the ext file
    boundary_object = hcdfm.Boundary(quantity='waterlevelbnd', #the FM quantity for tide is also waterlevelbnd
                                      locationfile=file_pli,
                                      forcingfile=ForcingModel_object)
    
    # add the boundary object to the ext file for each polyline in the polyfile
    ext_add_boundary_object_per_polyline(ext_new=ext_new, boundary_object=boundary_object)


def tidemodel_componentlist(tidemodel:str, convention:bool):
    
    comp_list_dict = {'FES2014':['2n2','eps2','j1','k1','k2','l2','la2','m2','m3','m4','m6','m8','mf','mks2','mm','mn4','ms4','msf','msqm','mtm','mu2','n2','n4','nu2','o1','p1','q1','r2','s1','s2','s4','sa','ssa','t2'],
                      'FES2012':['2N2','E2',  'J1','K1','K2','L2','La2','M2','M3','M4','M6','M8','Mf','MKS2','Mm','MN4','MS4','MSf',       'Mtm','Mu2','N2','N4','Nu2','O1','P1','Q1','R2','S1','S2','S4',     'Ssa','T2','Z0'],
                      'EOT20':['2N2','J1','K1','K2','M2','M4','MF','MM','N2','O1','P1','Q1','S1','S2','SA','SSA','T2'],
                      'GTSMv4.1':['2N2','2Q1','CHI1','EPS2','ETA2','J1','K1','K2','L2','LABDA2','M1','M2','M3','M4','MF','MM','MSF','MSQM','MU2','N2','NU2','O1','OO1','P1','PI1','Q1','R2','S1','S2','SA','SSA','T2'],
                      'tpxo80_opendap':['M2','S2','N2','K2','K1','O1','P1','Q1','MF','MM','M4','MS4','MN4'],
                      }
    comp_list_dict['GTSMv4.1_opendap'] = comp_list_dict['GTSMv4.1']
    
    component_list = comp_list_dict[tidemodel]
    
    if convention:
        return components_translate_upper(component_list)
    else:
        return component_list


def components_translate_upper(component_list):
    # translate dict from .\hydro_tools\FES\PreProcessing_FES_TideModel_imaginary.m
    translate_dict = {'LA2':'LABDA2', #TODO: use value instead of key in bc file? Support using value instead of key in const_list also (like line above)
                      'EPS2':'EPSILON2', 
                      'Z0':'A0',
                      'MTM':'MFM', #Needs to be verified
                      #added later
                      'E2':'EPSILON2',
                      }
    
    component_list = pd.Series([x.upper() for x in component_list]).replace(translate_dict, regex=True).to_list()
    return component_list


def interpolate_tide_to_plipoints(tidemodel, gdf_points, component_list=None, load=True):
    """
    empty docstring
    """
    
    dir_pattern_dict = {'FES2014': Path(PDRIVE, 'metocean-data', 'licensed', 'FES2014', '*.nc'), #ocean_tide_extrapolated
                        'FES2012': Path(PDRIVE, 'metocean-data', 'open', 'FES2012', 'data','*_FES2012_SLEV.nc'), #is eigenlijk ook licensed
                        'EOT20': Path(PDRIVE, 'metocean-data', 'open', 'EOT20', 'ocean_tides','*_ocean_eot20.nc'),
                        'GTSMv4.1': Path(PDRIVE, '1230882-emodnet_hrsm', 'GTSMv3.0EMODnet', 'EMOD_MichaelTUM_yearcomponents', 'GTSMv4.1_yeartide_2014_2.20.06', 'compare_fouhis_fouxyz_v4', 'GTSMv4.1_tide_2014_*_rasterized.nc'),
                        'GTSMv4.1_opendap': 'https://opendap.deltares.nl/thredds/dodsC/opendap/deltares/GTSM/GTSMv4.1_tide/GTSMv4.1_tide_2014_*_rasterized.nc',
                        'tpxo80_opendap':'https://opendap.deltares.nl/thredds/dodsC/opendap/deltares/delftdashboard/tidemodels/tpxo80/tpxo80.nc',
                        }
    if tidemodel not in dir_pattern_dict.keys():
        raise KeyError(f"Invalid tidemodel '{tidemodel}', options are: {str(list(dir_pattern_dict.keys()))}")
        
    dir_pattern = dir_pattern_dict[tidemodel]

    component_list_tidemodel = tidemodel_componentlist(tidemodel, convention=False)
    component_list_tidemodel_convention = tidemodel_componentlist(tidemodel, convention=True)
    component_tidemodel_translate_dict = {k:v for k,v in zip(component_list_tidemodel_convention,component_list_tidemodel)}
    if component_list is None:
        component_list = component_list_tidemodel_convention
    else:
        component_list = components_translate_upper(component_list)
        
    def extract_component(ds):
        if tidemodel in ['FES2012']:
            ds = ds.sel(lon=ds.lon<360) #drop last instance, since 0 and 360 are both present
            ds = ds.rename({'Ha':'amplitude','Hg':'phase'})
            ds['phase'] = ds['phase'].assign_attrs({'units':'degrees'}) # degree in original file
        elif tidemodel in ['EOT20']:
            ds = ds.rename({'imag':'wl_imag','real':'wl_real'})
            ds = ds.sel(lon=ds.lon<360) #drop last instance, since 0 and 360 are both present
        elif tidemodel in ['GTSMv4.1','GTSMv4.1_opendap']:
            ds = ds.rename({'wl_amp':'amplitude','wl_phs':'phase'})
        
        convert_360to180 = (ds['lon'].to_numpy()>180).any()
        if convert_360to180: # results in large chunks if it is done after concatenation, so do for each file before concatenation
            ds.coords['lon'] = (ds.coords['lon'] + 180) % 360 - 180
            ds = ds.sortby('lon')
        return ds
    
    if tidemodel=='tpxo80_opendap':
        ds = xr.open_dataset(dir_pattern)
        bool_land = ds.depth==0
        ds['tidal_amplitude_h'] = ds['tidal_amplitude_h'].where(~bool_land).assign_attrs({'units':'m'})
        ds['tidal_phase_h'] = ds['tidal_phase_h'].where(~bool_land).assign_attrs({'units':'degrees'})
        convert_360to180 = (ds['lon'].to_numpy()>180).any()
        if convert_360to180:
            ds.coords['lon'] = (ds.coords['lon'] + 180) % 360 - 180
            ds = ds.sortby('lon')
        ds = ds.rename({'tidal_amplitude_h':'amplitude','tidal_phase_h':'phase','constituents':'compno'})
        ds = ds.drop_dims(['lon_u','lon_v','lat_u','lat_v'])
        components_infile = ds['tidal_constituents'].load().str.decode('utf-8',errors='ignore').str.strip().str.upper()
        ds['compnames'] = components_infile
        ds = ds.set_index({'compno':'compnames'})
        if component_list is not None:
            ds = ds.sel(compno=component_list)
        data_xrsel = ds
    else:
        #use open_mfdataset() with preprocess argument to open all requested tide files into one Dataset
        file_list_nc = [str(dir_pattern).replace('*',component_tidemodel_translate_dict[comp]) for comp in component_list]
        data_xrsel = xr.open_mfdataset(file_list_nc, combine='nested', concat_dim='compno', preprocess=extract_component)
    
    data_xrsel = data_xrsel.rename({'lon':'longitude','lat':'latitude'})

    #convert cm to m
    if data_xrsel['amplitude'].attrs['units'] == 'cm':
        data_xrsel['amplitude'] /= 100
        data_xrsel['amplitude'].attrs['units'] = 'm'
    
    # compute imag/real components to avoid zero-crossing interpolation issues with phase
    if 'wl_real' not in data_xrsel.data_vars:
        data_xrsel_phs_rad = np.deg2rad(data_xrsel['phase'])
        data_xrsel['wl_real'] = data_xrsel['amplitude'] * np.cos(data_xrsel_phs_rad)
        data_xrsel['wl_imag'] = data_xrsel['amplitude'] * np.sin(data_xrsel_phs_rad)
    
    data_xrsel['compnames'] = xr.DataArray(component_list,dims=('compno'))
    data_xrsel = data_xrsel.set_index({'compno':'compnames'})
    
    data_interp = interp_regularnc_to_plipointsDataset(data_xr_reg=data_xrsel, gdf_points=gdf_points, load=load)
    data_interp['phase_new'] = np.rad2deg(np.arctan2(data_interp['wl_imag'],data_interp['wl_real']))
    return data_interp


def check_time_extent(data_xr, tstart, tstop):
    #convert tstart/tstop from str/dt.datetime/pd.Timestamp to pd.Timestamp. WARNING: when supplying '05-04-2016', monthfirst is assumed, so 2016-05-04 will be the result (may instead of april). So always supply yyyy-mm-dd datestrings
    tstart = pd.Timestamp(tstart)
    tstop = pd.Timestamp(tstop)
    
    #get timevar and compare requested dates
    timevar = data_xr['time']
    xr_tstartstop = pd.to_datetime(timevar.isel(time=[0,-1]).to_series())
    nc_tstart = xr_tstartstop.index[0]
    nc_tstop = xr_tstartstop.index[-1]
    if tstart < nc_tstart:
        raise OutOfRangeError(f'requested tstart {tstart} outside of available range {nc_tstart} to {nc_tstop}.')
    if tstop > nc_tstop:
        raise OutOfRangeError(f'requested tstop {tstop} outside of available range {nc_tstart} to {nc_tstop}.')


def ds_apply_conventions(data_xr):
    ncbnd_construct = get_ncbnd_construct()
    dimn_depth = ncbnd_construct['dimn_depth']
    varn_depth = ncbnd_construct['varn_depth']
    
    #rename dims and vars time/depth/lat/lon/x/y #TODO: this has to be phased out some time, or made as an argument or merged with conversion_dict?
    rename_dims_dict = {'depth':dimn_depth, # depth for CMEMS and many others
                        'lev':dimn_depth, # lev for GFDL
                        'lon':'longitude','lat':'latitude'}
    for k,v in rename_dims_dict.items():
        if k in data_xr.dims and v not in data_xr.dims:
            data_xr = data_xr.rename_dims({k:v})
            print(f'dimension {k} renamed to {v}')
        if k in data_xr.variables and v not in data_xr.variables:
            data_xr = data_xr.rename_vars({k:v})
            print(f'varname {k} renamed to {v}')
    
    # make depth positive up if not yet the case
    if varn_depth in data_xr.coords:
        if 'positive' in data_xr[varn_depth].attrs.keys():
            if data_xr[varn_depth].attrs['positive'] == 'down': #attribute appears in CMEMS, GFDL and CMCC, save to assume presence?
                data_xr[varn_depth] = -data_xr[varn_depth]
                data_xr[varn_depth].attrs['positive'] = 'up'
    
    # get calendar and maybe convert_calendar, makes sure that nc_tstart/nc_tstop are of type pd._libs.tslibs.timestamps.Timestamp
    data_xr_calendar = data_xr['time'].dt.calendar
    if data_xr_calendar != 'proleptic_gregorian': #this is for instance the case in case of noleap (or 365_days) calendars from GFDL and CMCC
        units_copy = data_xr['time'].encoding['units']
        print(f'WARNING: calendar different than proleptic_gregorian found ({data_xr_calendar}), convert_calendar is called so check output carefully. It should be no issue for datasets with a monthly interval.')
        data_xr = data_xr.convert_calendar('standard') #TODO: does this not result in 29feb nan values in e.g. GFDL model? Check missing argument at https://docs.xarray.dev/en/stable/generated/xarray.Dataset.convert_calendar.html
        data_xr['time'].encoding['units'] = units_copy #put back dropped units
    
    # 360 to 180 conversion
    if 'longitude' in data_xr.variables:
        convert_360to180 = (data_xr['longitude'].to_numpy()>180).any() #TODO: replace to_numpy() with load()
        if convert_360to180: #TODO: make more flexible for models that eg pass -180/+180 crossing (add overlap at lon edges).
            data_xr.coords['longitude'] = (data_xr.coords['longitude'] + 180) % 360 - 180
            data_xr = data_xr.sortby(data_xr['longitude'])
    
    return data_xr


def ds_apply_conversion_dict(data_xr, conversion_dict, quantity):
    # rename variables from conversion_dict
    for k,v in conversion_dict.items():
        ncvarn = v['ncvarname']
        if ncvarn in data_xr.variables.mapping.keys() and k == quantity: #k in quantity_list so phyc is not always renamed to tracerbndPON1 (first in conversion_dict)
            data_xr = data_xr.rename({ncvarn:k})
            print(f'variable {ncvarn} renamed to {k}')

    # optional conversion of units. Multiplications or other simple operatiors do not affect performance (dask.array(getitem) becomes dask.array(mul). With more complex operation it is better do do it on the interpolated array.
    # TODO: maybe do unit conversion before interp or is that computationally heavy?
    if 'conversion' in conversion_dict[quantity].keys(): #if conversion is present, unit key must also be in conversion_dict
        print(f'> converting units from [{data_xr[quantity].attrs["units"]}] to [{conversion_dict[quantity]["unit"]}]')
        # conversion drops all attributes of which units (which are changed anyway)
        data_xr[quantity] = data_xr[quantity] * conversion_dict[quantity]['conversion']
        # add unit attribute with resulting unit
        quantity_attrs = {'units': conversion_dict[quantity]['unit']}
        data_xr[quantity] = data_xr[quantity].assign_attrs(quantity_attrs)
    return data_xr


def open_prepare_dataset(dir_pattern, quantity, tstart, tstop, conversion_dict=None, refdate_str=None, chunks=None):
    """
    empty docstring
    """
    
    if conversion_dict is None:
        conversion_dict = get_conversion_dict()
    
    file_list_nc = glob.glob(str(dir_pattern))
    print(f'loading mfdataset of {len(file_list_nc)} files with pattern(s) {dir_pattern}')
    
    data_xr = xr.open_mfdataset(file_list_nc, chunks=chunks, join="exact") #TODO: does chunks argument solve "PerformanceWarning: Slicing is producing a large chunk."? {'time':1} is not a convenient chunking to use for timeseries extraction
    
    data_xr = ds_apply_conventions(data_xr=data_xr)
    data_xr = ds_apply_conversion_dict(data_xr=data_xr, conversion_dict=conversion_dict, quantity=quantity)

    # retrieve var(s) (after potential longitude conversion)
    data_xr_vars = data_xr[[quantity]]
    
    # slice dataset to times outside of requested time range
    data_xr_vars = _ds_sel_time_outside(
        ds=data_xr_vars,
        tstart=tstart,
        tstop=tstop,
        )
        
    #optional refdate changing
    if refdate_str is not None:
        if 'long_name' in data_xr_vars.time.attrs: #for CMEMS it is 'hours since 1950-01-01', which would be wrong now #TODO: consider also removing attrs for depth/varname, since we would like to have salinitybnd/waterlevel instead of Salinity/sea_surface_height in xr plots?
            del data_xr_vars.time.attrs['long_name']
        data_xr_vars.time.encoding['units'] = refdate_str
    
    return data_xr_vars


def _ds_sel_time_outside(ds: xr.Dataset, tstart, tstop) -> xr.Dataset:
    """
    Subset the dataset on time, making sure the requested times are always
    included. If there is no exact match for start/end times, the previous/next
    timestamp is taken as an extreme.
    Inspired by copernicusmarine.download_functions.subset_xarray.py
    """    
    check_time_extent(ds, tstart, tstop)
    external_minimum = ds.sel(time=tstart, method="pad")
    external_maximum = ds.sel(time=tstop, method="backfill")
    time_slice = slice(external_minimum.time.values, external_maximum.time.values)
    ds_sel = ds.sel(time=time_slice) 
    return ds_sel


def interp_regularnc_to_plipointsDataset(data_xr_reg, gdf_points, load=True):
    
    ncbnd_construct = get_ncbnd_construct()
    varn_pointx = ncbnd_construct['varn_pointx']
    varn_pointy = ncbnd_construct['varn_pointy']
    
    da_plipoints = da_from_gdf_points(gdf_points)
    
    #interpolation to lat/lon combinations
    print('> interp mfdataset to all PolyFile points (lat/lon coordinates)')
    
    # linear, nearest, combine_first
    data_interp_lin = data_xr_reg.interp(longitude=da_plipoints[varn_pointx], latitude=da_plipoints[varn_pointy],
                                         method='linear', 
                                         kwargs={'bounds_error':True}) #error is only raised upon load(), so when the actual value retrieval happens
    data_interp_near = data_xr_reg.interp(longitude=da_plipoints[varn_pointx], latitude=da_plipoints[varn_pointy],
                                          method='nearest', 
                                          kwargs={'bounds_error':True}) #error is only raised upon load(), so when the actual value retrieval happens
    data_interp = data_interp_lin.combine_first(data_interp_near)
    
    # drop original latitude/longitude vars (lat/lon in da_plipoints)
    data_interp = data_interp.drop_vars(['latitude','longitude'])
    
    if not load:
        return data_interp
    
    print(f'> actual extraction of data from netcdf with .load() (for {len(gdf_points)} plipoints at once, this might take a while)')
    dtstart = dt.datetime.now()
    try:
        data_interp_loaded = data_interp.load() #loading data for all points at once is more efficient compared to loading data per point in loop 
    except ValueError: #generate a proper error with outofbounds requested coordinates, default is "ValueError: One of the requested xi is out of bounds in dimension 0" #TODO: improve error in xarray
        lonvar_vals = data_xr_reg['longitude'].to_numpy()
        latvar_vals = data_xr_reg['latitude'].to_numpy()
        data_pol_pd = data_interp[[varn_pointx,varn_pointy]].to_dataframe()
        bool_reqlon_outbounds = (data_pol_pd[varn_pointx] <= lonvar_vals.min()) | (data_pol_pd[varn_pointx] >= lonvar_vals.max())
        bool_reqlat_outbounds = (data_pol_pd[varn_pointy] <= latvar_vals.min()) | (data_pol_pd[varn_pointy] >= latvar_vals.max())
        reqlatlon_pd = pd.DataFrame({'longitude':data_pol_pd[varn_pointx],'latitude':data_pol_pd[varn_pointy],'lon outbounds':bool_reqlon_outbounds,'lat outbounds':bool_reqlat_outbounds})
        reqlatlon_pd_outbounds = reqlatlon_pd.loc[bool_reqlon_outbounds | bool_reqlat_outbounds]
        raise ValueError(f'{len(reqlatlon_pd_outbounds)} of requested pli points are out of bounds (valid longitude range {lonvar_vals.min()} to {lonvar_vals.max()}, valid latitude range {latvar_vals.min()} to {latvar_vals.max()}):\n{reqlatlon_pd_outbounds}')
    time_passed = (dt.datetime.now()-dtstart).total_seconds()
    print(f'>>time passed: {time_passed:.2f} sec')

    return data_interp_loaded


def interp_uds_to_plipoints(uds:xu.UgridDataset, gdf:geopandas.GeoDataFrame) -> xr.Dataset:
    """
    To interpolate an unstructured dataset (like a _map.nc file) read with xugrid to plipoint locations
    
    Parameters
    ----------
    uds : xu.UgridDataset
        dfm model output read using dfm_tools. Dims: mesh2d_nLayers, mesh2d_nInterfaces, time, mesh2d_nNodes, mesh2d_nFaces, mesh2d_nMax_face_nodes, mesh2d_nEdges.
    gdf : geopandas.GeoDataFrame
        gdf with location, geometry (Point) and crs corresponding to model crs.

    Raises
    ------
    Exception
        DESCRIPTION.

    Returns
    -------
    ds : TYPE
        xr.Dataset with dims: node, time, z.

    """
    # TODO: this function requires gdf with points, but the interp_regularnc_to_plipoints requires paths to plifiles (and others also)
    
    facedim = uds.grid.face_dimension
    edgedim = uds.grid.edge_dimension
    nodedim = uds.grid.node_dimension
    ncbnd_construct = get_ncbnd_construct()
    dimn_point = ncbnd_construct['dimn_point']
    varn_pointname = ncbnd_construct['varn_pointname']
    
    # drop node/edge variables since they are not interpolated but 
    # get an additional face dimension if some points are out of bounds
    # TODO: revert after fixing https://github.com/Deltares/xugrid/issues/274
    vars_without_facedim = []
    for varn in uds.variables:
        if facedim not in uds.variables[varn].dims:
            vars_without_facedim.append(varn)
    uds_face = uds.drop_vars(vars_without_facedim)
    
    # interpolate to provided points
    ds = uds_face.ugrid.sel_points(x=gdf.geometry.x, y=gdf.geometry.y)
    
    # re-add removed variables again, sometimes important for e.g. depth
    # TODO: remove after fixing https://github.com/Deltares/xugrid/issues/274
    for varn in vars_without_facedim:
        vardims = uds.variables[varn].dims
        if edgedim not in vardims and nodedim not in vardims:
            ds[varn] = uds.variables[varn]
    
    # rename station dimname and varname (is index, are both mesh2d_nFaces to start with)
    ds = ds.rename({facedim:dimn_point}) # rename mesh2d_nFaces to plipoints
    ds = ds.rename_vars({dimn_point:varn_pointname})
    # change name of plipoint (node to gdf name)
    ds[varn_pointname] = xr.DataArray(gdf[varn_pointname].tolist(), dims=dimn_point)
    return ds


def interp_hisnc_to_plipoints(data_xr_his, file_pli, kdtree_k=3, load=True):
    """
    interpolate stations in hisfile to points of polyfile
    """
    #KDTree, finds minimal eucledian distance between points (haversine would be better). Alternatively sklearn.neighbors.BallTree: tree = BallTree(data_lonlat_pd)
    
    ncbnd_construct = get_ncbnd_construct()
    dimn_point = ncbnd_construct['dimn_point']
    
    datavars_list = list(data_xr_his.data_vars)
    if len(datavars_list)>5:
        print('WARNING: more than 5 data_vars, you might want to subset your data_xr_his')
    #read hisfile and make KDTree of xy of hisstations
    hisstations_pd = data_xr_his.stations.to_dataframe() #TODO: add check if stations are strings? (whether preprocess_hisnc was used)
    tree_nest2 = KDTree(hisstations_pd[['station_x_coordinate','station_y_coordinate']])
    
    #read polyfile
    polyfile_object = hcdfm.PolyFile(file_pli)
    gdf_points = PolyFile_to_geodataframe_points(polyfile_object)
    gdf_coords = gdf_points.get_coordinates()
    
    # query k nearest hisstations (names)
    plicoords_distance2, plicoords_nestpointidx = tree_nest2.query(gdf_coords, k=kdtree_k)
    da_plicoords_nestpointidx = xr.DataArray(plicoords_nestpointidx, dims=(dimn_point,'nearestkpoints'))
    da_plicoords_nestpointnames = data_xr_his.stations.isel(stations=da_plicoords_nestpointidx)
    
    # convert gdf to xarray dataset
    data_interp = da_from_gdf_points(gdf_points)
    
    #interpolate hisfile variables to plipoints
    for varone in datavars_list:
        da_dist = xr.DataArray(plicoords_distance2, dims=(dimn_point,'nearestkpoints'))
        da_invdistweight = (1/da_dist)/(1/da_dist).sum(dim='nearestkpoints') #TODO: set weigths for invalid points to 0 (and increase others weights so sum remains 1)
        data_xr_his_pli_knearest = data_xr_his[varone].sel(stations=da_plicoords_nestpointnames)
        data_interp[varone] = (data_xr_his_pli_knearest * da_invdistweight).sum(dim='nearestkpoints')
        data_interp[varone].attrs = data_xr_his[varone].attrs #copy units and other attributes
        
    if not load:
        return data_interp
    
    print('loading dataset (might take a while, but increases performance of point loop significantly)')
    data_interp = data_interp.load()
    print('done')

    return data_interp
    

def plipointsDataset_to_ForcingModel(plipointsDataset):
    """
    empty docstring
    """
    
    ncbnd_construct = get_ncbnd_construct()
    dimn_point = ncbnd_construct['dimn_point']
    dimn_depth = ncbnd_construct['dimn_depth']
    varn_pointname = ncbnd_construct['varn_pointname']
    
    plipointsDataset = maybe_convert_fews_to_dfmt(plipointsDataset)
    
    quantity_list = list(plipointsDataset.data_vars)
    npoints = len(plipointsDataset[dimn_point])
    
    #start conversion to Forcingmodel object
    print(f'Converting {npoints} plipoints to hcdfm.ForcingModel():',end='')
    dtstart = dt.datetime.now()
    ForcingModel_object = hcdfm.ForcingModel()
    for iP in range(npoints):
        print(f' {iP+1}',end='')
        
        #select data for this point, ffill nans, concatenating time column, constructing T3D/TimeSeries and append to hcdfm.ForcingModel()
        datablock_xr_onepoint = plipointsDataset.isel({dimn_point:iP})
        plipoint_name = str(datablock_xr_onepoint[varn_pointname].to_numpy())
        plipoint_onlynan = False
        for quan in quantity_list:
            datablock_xr_onepoint[quan].attrs['locationname'] = plipoint_name #TODO: is there a nicer way of passing this data?
            if datablock_xr_onepoint[quan].isnull().all(): # check if all values of plipoint are nan (on land)
                plipoint_onlynan = True
                logger.warning(f'Plipoint "{plipoint_name}" might be on land since it only contain nan values. '
                               'This point is skipped to avoid bc-writing errors. Consider altering your PolyFile or extrapolate the data.')
        
        # skip this point if quantity has only-nan values
        if plipoint_onlynan:
            continue
        
        if dimn_depth in plipointsDataset.dims:
            ts_one = Dataset_to_T3D(datablock_xr_onepoint)
        elif 'amplitude' in quantity_list:
            ts_one = Dataset_to_Astronomic(datablock_xr_onepoint)
        else:
            ts_one = Dataset_to_TimeSeries(datablock_xr_onepoint)
        ForcingModel_object.forcing.append(ts_one)
    print(f'. >> done in {(dt.datetime.now()-dtstart).total_seconds():.2f} sec')
    
    return ForcingModel_object

