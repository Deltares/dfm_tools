# -*- coding: utf-8 -*-
"""
Created on Tue Jun 20 15:41:08 2023

@author: veenstra
"""

import dfm_tools as dfmt
import matplotlib.pyplot as plt
plt.close('all')
import xugrid as xu
import pandas as pd
import hydrolib.core.dflowfm as hcdfm

# Read the pli points:
list_plifiles = [r'p:\11208479-sequestration-seaweed\Oosterschelde_DFM_hydro_waq\dflowfm3d-oosterschelde-wq\boundary_conditions\zee-noord.pli',
                 r'p:\11208479-sequestration-seaweed\Oosterschelde_DFM_hydro_waq\dflowfm3d-oosterschelde-wq\boundary_conditions\zee-west.pli',
                 r'p:\11208479-sequestration-seaweed\Oosterschelde_DFM_hydro_waq\dflowfm3d-oosterschelde-wq\boundary_conditions\zee-zuid.pli']

# read plifiles into single gdf and convert coordinates
gdf_list = [dfmt.PolyFile_to_geodataframe_points(hcdfm.PolyFile(x), crs='epsg:28992') for x in list_plifiles]
gdf = pd.concat(gdf_list, ignore_index=True)
gdf = gdf.to_crs('epsg:4326')

# read model (multi-years) #TODO: simplify reading of model files after issue is solved: https://github.com/Deltares/xugrid/issues/80
chunks = {'time':1}
files_p1 = [r'p:\\archivedprojects\\11206044-002ospareutrophication\\final_results_Snellius_latestversion_20220210\\A16b_Numtopsiguniform1_2015_pH\\DFM_OUTPUT_DCSM-FM_0_5nm_waq\\DCSM-FM_0_5nm_waq_0146_map.nc',r'p:\\archivedprojects\\11206044-002ospareutrophication\\final_results_Snellius_latestversion_20220210\\A16b_Numtopsiguniform1_2016_pH\\DFM_OUTPUT_DCSM-FM_0_5nm_waq\\DCSM-FM_0_5nm_waq_0146_map.nc']
files_p2 = [r'p:\\archivedprojects\\11206044-002ospareutrophication\\final_results_Snellius_latestversion_20220210\\A16b_Numtopsiguniform1_2015_pH\\DFM_OUTPUT_DCSM-FM_0_5nm_waq\\DCSM-FM_0_5nm_waq_0147_map.nc',r'p:\\archivedprojects\\11206044-002ospareutrophication\\final_results_Snellius_latestversion_20220210\\A16b_Numtopsiguniform1_2016_pH\\DFM_OUTPUT_DCSM-FM_0_5nm_waq\\DCSM-FM_0_5nm_waq_0147_map.nc']
files_p3 = [r'p:\\archivedprojects\\11206044-002ospareutrophication\\final_results_Snellius_latestversion_20220210\\A16b_Numtopsiguniform1_2015_pH\\DFM_OUTPUT_DCSM-FM_0_5nm_waq\\DCSM-FM_0_5nm_waq_0148_map.nc',r'p:\\archivedprojects\\11206044-002ospareutrophication\\final_results_Snellius_latestversion_20220210\\A16b_Numtopsiguniform1_2016_pH\\DFM_OUTPUT_DCSM-FM_0_5nm_waq\\DCSM-FM_0_5nm_waq_0148_map.nc']
all_p1 = xu.open_mfdataset(files_p1,chunks=chunks)
all_p2 = xu.open_mfdataset(files_p2,chunks=chunks)
all_p3 = xu.open_mfdataset(files_p3,chunks=chunks)
uds = xu.merge_partitions([all_p1,all_p2, all_p3])
uds = uds.ugrid.sel(y=slice(None,54)) #to cut off some rogue cells

#conversion and depth extraction
uds = dfmt.rename_waqvars(uds)

# interpolate depths (zsigma to z):
depths_slice = range(-25,0,1)# uds.mesh2d_layer_z.values
ds_atdepths = dfmt.get_Dataset_atdepths(data_xr=uds, depths=depths_slice)
ds_atdepths = ds_atdepths.rename({'depth_from_z0':'depth'})

#interpolate to plipoints
ds_plipoints = dfmt.interp_uds_to_plipoints(uds=ds_atdepths, gdf=gdf)
fig,ax = plt.subplots()
uds.mesh2d_flowelem_bl.ugrid.plot(ax=ax,center=False)
gdf.plot(ax=ax)

fig,ax = plt.subplots()
ds_plipoints.mesh2d_sa1.isel(time=-1).T.plot()

for i in range(1,24,8):
    fig,ax = plt.subplots()
    ds_plipoints.mesh2d_sa1.isel(node=i).T.plot(ax=ax)
    ax.set_ylim(-30,1)

