# -*- coding: utf-8 -*-
"""
Created on Tue Jun 20 15:41:08 2023

@author: veenstra
"""

from pathlib import Path
import dfm_tools as dfmt
import matplotlib.pyplot as plt
plt.close('all')

import geopandas
import xarray as xr
import xugrid as xu
import pandas as pd
from shapely.geometry import Point


def dutch_crs_plipoint_list_to_decimaldeg_gdf(list_plifiles):
    '''
    input:
    - list of all the plifiles as WindowsPath. Format = from offshore boundary pli files. Dutch crs!
    
    output:
    - gdf of all the pli files with columns: location, geometry (Points). Crs = 4326 (decimal degrees)
    
    '''
    all_boundaries = []
    for i,file in list(enumerate(list_plifiles)):  # start with zuid (as did for testing)
        print(str(file).split('\\')[-1].split('.')[0])

        # Read points within this offshore boundary:
        pli_file = pd.read_csv(list_plifiles[i], skiprows=1, sep=' ')   # read pli file
        pli_file = pli_file.dropna(axis=1)                              # drop empty columns
        # Convert dutch crs to decimal degrees (for DCSM)
        points = [] 
        for i in range(0,len(pli_file)):
            points.append(Point(pli_file[pli_file.columns[0]].to_list()[i], pli_file[pli_file.columns[1]].to_list()[i]))
        d = {'location': pli_file[pli_file.columns[2]].to_list(), 'geometry': points}
        gdf = geopandas.GeoDataFrame(d, crs='epsg:28992')   # set Dutch crs
        gdf = gdf.to_crs(4326)                              # convert to decimal degrees
        all_boundaries.append(gdf)
    all_boundaries = pd.concat(all_boundaries, ignore_index=True)
    return all_boundaries



# Read the pli points:
list_plifiles = [Path(r'p:\11208479-sequestration-seaweed\Oosterschelde_DFM_hydro_waq\dflowfm3d-oosterschelde-wq\boundary_conditions\zee-noord.pli'),
                Path(r'p:\11208479-sequestration-seaweed\Oosterschelde_DFM_hydro_waq\dflowfm3d-oosterschelde-wq\boundary_conditions\zee-west.pli'),
                Path(r'p:\11208479-sequestration-seaweed\Oosterschelde_DFM_hydro_waq\dflowfm3d-oosterschelde-wq\boundary_conditions\zee-zuid.pli')]
# Use function
gdf = dutch_crs_plipoint_list_to_decimaldeg_gdf(list_plifiles)
gdf


def interpolate_uds_to_plipoints(gdf, uds, nPoints=None):
    '''
    input:
    - gdf with location, geometry (Point) and crs=4326
    - uds: dfm model output read using dfm_tools. 
           Dims: mesh2d_nLayers, mesh2d_nInterfaces, time, mesh2d_nNodes, mesh2d_nFaces, mesh2d_nMax_face_nodes, mesh2d_nEdges
    - nPoints: amount of points (None gives all)
    
    output:
    - xr.Dataset with dims: plipoints, time, depth
    
    '''
    #TODO: align plipoints/gdf with other functions
    
    gdf_sel = gdf.iloc[:nPoints]
    ds = uds.ugrid.sel_points(x=gdf_sel.geometry.x, y=gdf_sel.geometry.y)
    
    if len(gdf_sel)!=uds.dims['mesh2d_nFaces']: #TODO: check this until https://github.com/Deltares/xugrid/issues/100 is solved
        raise Exception(f'requested {len(gdf_sel)} points but resulted in ds with {ds.dims["mesh2d_nFaces"]} points, some points are outside of the model domain')
    
    ds = ds.rename({'depth_from_z0':'depth',
                    'mesh2d_nFaces':'plipoints'}) # rename mesh2d_nFaces to plipoints
    
    ds['plipoint_name'] = xr.DataArray(gdf_sel.location.tolist(), dims='plipoints') # change name of plipoint (node to gdf name)
    ds = ds.set_index({'plipoints':'plipoint_name'})

    return ds


# read model (multi-years)
chunks = {'time':1}
files_p1 = [r'p:\\archivedprojects\\11206044-002ospareutrophication\\final_results_Snellius_latestversion_20220210\\A16b_Numtopsiguniform1_2015_pH\\DFM_OUTPUT_DCSM-FM_0_5nm_waq\\DCSM-FM_0_5nm_waq_0146_map.nc',r'p:\\archivedprojects\\11206044-002ospareutrophication\\final_results_Snellius_latestversion_20220210\\A16b_Numtopsiguniform1_2016_pH\\DFM_OUTPUT_DCSM-FM_0_5nm_waq\\DCSM-FM_0_5nm_waq_0146_map.nc']
files_p2 = [r'p:\\archivedprojects\\11206044-002ospareutrophication\\final_results_Snellius_latestversion_20220210\\A16b_Numtopsiguniform1_2015_pH\\DFM_OUTPUT_DCSM-FM_0_5nm_waq\\DCSM-FM_0_5nm_waq_0147_map.nc',r'p:\\archivedprojects\\11206044-002ospareutrophication\\final_results_Snellius_latestversion_20220210\\A16b_Numtopsiguniform1_2016_pH\\DFM_OUTPUT_DCSM-FM_0_5nm_waq\\DCSM-FM_0_5nm_waq_0147_map.nc']
files_p3 = [r'p:\\archivedprojects\\11206044-002ospareutrophication\\final_results_Snellius_latestversion_20220210\\A16b_Numtopsiguniform1_2015_pH\\DFM_OUTPUT_DCSM-FM_0_5nm_waq\\DCSM-FM_0_5nm_waq_0148_map.nc',r'p:\\archivedprojects\\11206044-002ospareutrophication\\final_results_Snellius_latestversion_20220210\\A16b_Numtopsiguniform1_2016_pH\\DFM_OUTPUT_DCSM-FM_0_5nm_waq\\DCSM-FM_0_5nm_waq_0148_map.nc']
all_p1 = xu.open_mfdataset(files_p1,chunks=chunks)
all_p2 = xu.open_mfdataset(files_p2,chunks=chunks)
all_p3 = xu.open_mfdataset(files_p3,chunks=chunks)
uds = xu.merge_partitions([all_p1,all_p2, all_p3])

#conversion and depth extraction
uds = dfmt.rename_waqvars(uds)

depths_slice = uds.mesh2d_layer_z.values
# convert depths (sigma to meter):
ds_atdepths = dfmt.get_Dataset_atdepths(data_xr=uds, depths=depths_slice)

# use definition:
nPoints = None

uds_onplipoints = interpolate_uds_to_plipoints(gdf=gdf, uds=ds_atdepths, nPoints=nPoints)

fig,ax = plt.subplots()
#uds_onplipoints.mesh2d_sa1.isel(time=-1).plot()
uds_onplipoints.mesh2d_sa1.isel(plipoints=1).T.plot(ax=ax)
ax.set_ylim(-50,1)

