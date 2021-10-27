# -*- coding: utf-8 -*-
"""
Created on Wed Oct 27 14:40:34 2021

@author: veenstra

https://contextily.readthedocs.io/en/latest/reference.html
https://contextily.readthedocs.io/en/latest/intro_guide.html
ctx.add_basemap() defaults:
    source: None defaults to ctx.providers.Stamen.Terrain
    crs: coordinate reference system (CRS). If None (default), no warping is performed and the original Spherical Mercator (EPSG:3857) is used.
"""

import os
import matplotlib.pyplot as plt
plt.close('all')

from dfm_tools.testutils import try_importmodule
try_importmodule(modulename='contextily') #check if contextily was installed since it is an optional module, also happens in plot_cartopybasemap()
import contextily as ctx

from dfm_tools.get_nc import get_netdata, get_ncmodeldata, plot_netmapdata
dir_testinput = r'c:\DATA\dfm_tools_testdata'
dir_output = '.'

file_nc_map = os.path.join(dir_testinput,'DFM_3D_z_Grevelingen\\computations\\run01\\DFM_OUTPUT_Grevelingen-FM\\Grevelingen-FM_0000_map.nc')
ugrid = get_netdata(file_nc=file_nc_map)
data_frommap_bl = get_ncmodeldata(file_nc=file_nc_map, varname='mesh2d_flowelem_bl')

source_list = [ctx.providers.Stamen.Terrain, #default source
               ctx.providers.Esri.WorldImagery,
               ctx.providers.CartoDB.Voyager,
               #ctx.providers.NASAGIBS.ViirsEarthAtNight2012,
               ctx.providers.Stamen.Watercolor]

for source_ctx in source_list:
    source_name = source_ctx['name'].replace('.','_')
    fig, ax = plt.subplots(1,1,figsize=(10,6))
    pc = plot_netmapdata(ugrid.verts, values=data_frommap_bl, ax=ax, linewidth=0.5, cmap='jet')
    fig.colorbar(pc, ax=ax)
    fig.tight_layout()
    ctx.add_basemap(ax, source=source_ctx, crs="EPSG:28992", attribution_size=5)
    fig.savefig(os.path.join(dir_output,'contextily_grevelingen_RD_%s'%(source_name)))
