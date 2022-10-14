# -*- coding: utf-8 -*-
"""
Created on Wed Oct 27 14:21:28 2021

@author: veenstra
"""

import os
import matplotlib.pyplot as plt
plt.close('all')
import numpy as np
import contextily as ctx
import xarray as xr

from dfm_tools.get_nc import get_netdata, get_ncmodeldata, get_xzcoords_onintersection, plot_netmapdata, plot_ztdata
from dfm_tools.get_nc_helpers import get_ncvarproperties, get_stationid_fromstationlist
from dfm_tools.linebuilder import LineBuilder

dir_testinput = r'c:\DATA\dfm_tools_testdata'

file_nc = os.path.join(dir_testinput,'DFM_3D_z_Grevelingen','computations','run01','DFM_OUTPUT_Grevelingen-FM','Grevelingen-FM_0000_map.nc')

source = ctx.providers.Esri.WorldImagery # ctx.providers.Stamen.Terrain (default), ctx.providers.CartoDB.Voyager, ctx.providers.NASAGIBS.ViirsEarthAtNight2012, ctx.providers.Stamen.Watercolor
line_array = np.array([[ 53181.96942503, 424270.83361629],
                       [ 55160.15232593, 416913.77136685]])
#line_array = np.array([[ 49324.74390138, 423723.47046544],
#       [ 50751.23158142, 417109.75485795],
#       [ 62163.13302179, 423139.9073236 ]])
val_ylim = [-25,5]

#get variables and ugrid
vars_pd = get_ncvarproperties(file_nc=file_nc)
vars_pd_sel = vars_pd[['standard_name','long_name','shape','dimensions','dtype']]
ugrid = get_netdata(file_nc=file_nc)


print('get grid and satellite background')
fig, ax = plt.subplots(figsize=(8,5.5))
pc = plot_netmapdata(ugrid.verts, values=None, ax=None, linewidth=0.5, color="crimson", facecolor="None")
ctx.add_basemap(ax, source=source, crs="EPSG:28992", attribution=False)
#ax.set_aspect('equal')
fig.tight_layout()


print('get bedlevel and create plot with ugrid and cross section line')
data_frommap_bl = get_ncmodeldata(file_nc=file_nc, varname='mesh2d_flowelem_bl')
fig, ax_input = plt.subplots(figsize=(9,5.5))
pc = plot_netmapdata(ugrid.verts, values=data_frommap_bl, ax=ax_input, linewidth=0.5, edgecolors='face', cmap='jet')#, color='crimson', facecolor="None")
fig.colorbar(pc, ax=ax_input)
if 0: #click interactive polygon #TODO: this is useful but should work also without killing the code
    line, = ax_input.plot([], [],'o-')  # empty line
    linebuilder = LineBuilder(line) #after this click your line and then run the line below
    #breakit
    line_array = linebuilder.line_array
ax_input.plot(line_array[0,0],line_array[0,1],'bx',linewidth=3,markersize=10)
ax_input.plot(line_array[:,0],line_array[:,1],'b',linewidth=3)
ctx.add_basemap(ax=ax_input, source=source, crs="EPSG:28992", attribution=False)
#ax_input.set_aspect('equal')
fig.tight_layout()


print('plot cross section')
intersect_pd = ugrid.polygon_intersect(line_array, optimize_dist=False, calcdist_fromlatlon=None)
crs_verts, crs_plotdata = get_xzcoords_onintersection(file_nc=file_nc, varname='mesh2d_sa1', intersect_pd=intersect_pd, timestep=3)
fig, ax = plt.subplots(figsize=(8,5.5))
pc = plot_netmapdata(crs_verts, values=crs_plotdata, ax=ax, cmap='jet')
fig.colorbar(pc, ax=ax)
ax.set_ylim(val_ylim)
fig.tight_layout()


print('plot grid and values from mapdata (salinity on layer, 3dim, on cell centers)')
for iL in [28,33,34]:
    data_frommap = get_ncmodeldata(file_nc=file_nc, varname='mesh2d_sa1', timestep=3, layer=iL)
    fig, ax = plt.subplots(figsize=(9,5.5))
    pc = plot_netmapdata(ugrid.verts, values=data_frommap, ax=None, linewidth=0.5, cmap="jet", edgecolor='face')
    cbar = fig.colorbar(pc, ax=ax)
    cbar.set_label('%s [%s]'%(data_frommap.var_ncattrs['long_name'], data_frommap.var_ncattrs['units']))
    ctx.add_basemap(ax, source=source, crs="EPSG:28992", attribution=False)
    #ax.set_aspect('equal')
    fig.tight_layout()








file_nc_his = os.path.join(dir_testinput,'DFM_3D_z_Grevelingen\\computations\\run01\\DFM_OUTPUT_Grevelingen-FM\\Grevelingen-FM_0000_his.nc')

station = ['GTSO-01','GTSO-02','GTSO-03','GTSO-04','GTSO-05','GTSO-06','GTSO-07',
           #'GTSO-08','GTSO-09','GTSO-10','GTSO-11','GTSO-12','GTSO-13','GTSO-14',
           #'GTSO-15','GTSO-16','GTSO-17','GTSO-18','GTSO-19','GTSO-20',
           'Bommenede','Grevelingen hevel West','Brouwerssluis binnen','Brouwerssluis binnen-hand']
station_zt = ['GTSO-02']
    
data_xr = xr.open_dataset(file_nc_his)
data_xr['station_name_str'] = data_xr['station_name'].str.decode('utf-8',errors='ignore').str.strip()
data_xr = data_xr.set_coords('station_name_str')

idx_stations = get_stationid_fromstationlist(data_xr, stationlist=station)
idx_stations_zt = get_stationid_fromstationlist(data_xr, stationlist=station_zt)[0] #if provide single station (string, no list), the shape of the resulting xarray is correct

print('plot waterlevel from his')
data_fromhis_xr = data_xr.waterlevel.isel(stations=idx_stations)#.sel(time=slice('2007-11-01','2007-11-04'))
fig, ax = plt.subplots()
data_fromhis_xr.plot.line('-',ax=ax,x='time')
ax.legend(data_fromhis_xr['station_name_str'].to_numpy())

print('plot salinity from his')
data_fromhis_xr = data_xr.salinity.isel(stations=idx_stations,laydim=20)
fig, ax = plt.subplots()
data_fromhis_xr.plot.line('-',ax=ax,x='time')
ax.legend(data_fromhis_xr['station_name_str'].to_numpy())

print('zt temperature plot and wl')
data_xr_selzt = data_xr.isel(stations=idx_stations_zt,time=slice(40,100))
data_fromhis_wl_xr = data_xr_selzt['waterlevel']
fig, (axwl,ax1) = plt.subplots(2,1,figsize=(12,7),gridspec_kw={'height_ratios':[1,2]},sharex=True,sharey=True)
axwl.plot(data_xr_selzt.time[[0,-1]],[0,0],'k-',linewidth=0.5)
ax1.plot(data_xr_selzt.time[[0,-1]],[0,0],'k-',linewidth=0.5)
axwl.plot(data_xr_selzt.time,data_fromhis_wl_xr,'-',label=f'wl {station_zt}')
c = plot_ztdata(data_xr_sel=data_xr_selzt, varname='temperature', ax=ax1, cmap='jet')
fig.colorbar(c,ax=axwl)
fig.colorbar(c,ax=ax1)
#contour
CS = plot_ztdata(data_xr_sel=data_xr_selzt, varname='temperature', ax=ax1, only_contour=True, levels=6, colors='k', linewidths=0.8, linestyles='solid')
ax1.clabel(CS, fontsize=10)
fig.tight_layout()


