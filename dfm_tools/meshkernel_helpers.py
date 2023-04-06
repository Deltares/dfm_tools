# -*- coding: utf-8 -*-
"""
Created on Thu Apr  6 18:46:35 2023

@author: veenstra
"""

import xugrid as xu
import meshkernel
import xarray as xr
import datetime as dt
import hydrolib.core.dflowfm as hcdfm
import pandas as pd
from dfm_tools.hydrolib_helpers import pointlike_to_DataFrame
import getpass
from importlib.metadata import version
dfmt_version = version('dfm_tools') #TODO: cleaner way to do this?


def meshkernel_delete_withpol(mk, file_ldb, minpoints=None):
    
    print('>> reading+converting ldb: ',end='')
    dtstart = dt.datetime.now()
    pol_ldb = hcdfm.PolyFile(file_ldb)
    pol_ldb_list = [pointlike_to_DataFrame(x) for x in pol_ldb.objects] #TODO, this is quite slow, speed up possible?
    if minpoints is not None:
        pol_ldb_list = [x for x in pol_ldb_list if len(x)>minpoints] #filter only large polygons for performance
    for iP, pol_ldb in enumerate(pol_ldb_list):
        if not (pol_ldb.iloc[0] == pol_ldb.iloc[-1]).all(): #close the polygon if it is not yet closed
            pol_ldb_list[iP] = pd.concat([pol_ldb,pol_ldb.iloc[[0]]],axis=0)
    print(f'{(dt.datetime.now()-dtstart).total_seconds():.2f} sec')
    
    for iP, pol_del in enumerate(pol_ldb_list): #TODO: also possible without loop? >> geometry_separator=-999.9 so that value can be used to concat polygons. >> use hydrolib poly as input? https://github.com/Deltares/MeshKernelPy/issues/35
        delete_pol_geom = meshkernel.GeometryList(x_coordinates=pol_del['x'].to_numpy(), y_coordinates=pol_del['y'].to_numpy()) #TODO: .copy()/to_numpy() makes the array contiguous in memory, which is necessary for meshkernel.mesh2d_delete()
        mk.mesh2d_delete(geometry_list=delete_pol_geom, 
                         delete_option=meshkernel.DeleteMeshOption(2), #ALL_COMPLETE_FACES/2: Delete all faces of which the complete face is inside the polygon
                         invert_deletion=False) #TODO: cuts away link that is neccesary, so results in non-orthogonal grid
    return mk


def meshkernel_to_UgridDataset(mk:meshkernel.meshkernel.MeshKernel) -> xr.Dataset:
    mesh2d_grid3 = mk.mesh2d_get()

    xu_grid = xu.Ugrid2d.from_meshkernel(mesh2d_grid3)
    #xu_grid_ds = xu_grid.to_dataset()
    xu_grid_ds = xu_grid.assign_face_coords(xu_grid.to_dataset()) #TODO: face_coords are necessary for hydrolib-core mdu.geometry.network validation, report this
    
    #convert 0-based to 1-based grid #TODO: FM kernel needs 1-based grid, but it should read the attributes instead. Report this (#ug_get_meshgeom, #12, ierr=0. ** WARNING: Could not read mesh face x-coordinates)
    varn_fnc = xu_grid_ds.mesh2d.attrs['face_node_connectivity'] #TODO: maybe also do for edge_node_connectivity
    if xu_grid_ds[varn_fnc].attrs['start_index']==0:
        xu_grid_ds[varn_fnc] += 1
        xu_grid_ds[varn_fnc].attrs["_FillValue"] += 1
        xu_grid_ds[varn_fnc].attrs["start_index"] += 1
        
    xu_grid_ds = xu_grid_ds.assign_attrs({#'Conventions': 'CF-1.8 UGRID-1.0 Deltares-0.10', #add Deltares convention (was CF-1.8 UGRID-1.0)
                                          'institution': 'Deltares',
                                          'references': 'https://www.deltares.nl',
                                          'source': f'Created with meshkernel {meshkernel.__version__}, xugrid {xu.__version__} and dfm_tools {dfmt_version}',
                                          'history': 'Created on %s, %s'%(dt.datetime.now().strftime('%Y-%m-%dT%H:%M:%S%z'),getpass.getuser()), #TODO: add timezone
                                          })
    
    # add attrs (to projected_coordinate_system/wgs84 empty int variable): #TODO: should depend on is_geographic flag in make_basegrid()
    # attribute_dict = {
    #     'name': 'WGS84',
    #     'epsg': np.array([4326], dtype=int),
    #     'grid_mapping_name': 'Unknown projected',
    #     'longitude_of_prime_meridian': np.array([0.0], dtype=float),
    #     'semi_major_axis': np.array([6378137.0], dtype=float),
    #     'semi_minor_axis': np.array([6356752.314245], dtype=float),
    #     'inverse_flattening': np.array([6356752.314245], dtype=float),
    #     'EPSG_code': 'EPSG:4326',
    #     'value': 'value is equal to EPSG code'}
    # xu_grid_ds['wgs84'] = xr.DataArray(np.array(0,dtype=int),dims=(),attrs=attribute_dict)
    
    xu_grid_uds = xu.UgridDataset(xu_grid_ds)
    return xu_grid_uds