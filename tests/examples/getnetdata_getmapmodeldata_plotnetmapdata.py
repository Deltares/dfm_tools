# -*- coding: utf-8 -*-
"""
Created on Wed Oct 27 14:42:52 2021

@author: veenstra
"""

import os
import matplotlib.pyplot as plt
plt.close('all')
import numpy as np

from dfm_tools.get_nc import get_netdata, get_ncmodeldata, plot_netmapdata
from dfm_tools.regulargrid import scatter_to_regulargrid
from dfm_tools.get_nc_helpers import get_ncvardimlist
from dfm_tools.testutils import try_importmodule
try_importmodule(modulename='contextily') #check if contextily was installed since it is an optional module, also happens in plot_cartopybasemap()
import contextily as ctx

dir_testinput = r'c:\DATA\dfm_tools_testdata'
dir_output = '.'

file_nc_list = [os.path.join(dir_testinput,r'DFM_sigma_curved_bend\DFM_OUTPUT_cb_3d\cb_3d_map.nc'),
                os.path.join(dir_testinput,r'DFM_3D_z_Grevelingen\computations\run01\DFM_OUTPUT_Grevelingen-FM\Grevelingen-FM_0000_map.nc'),
                r'p:\11206813-006-kpp2021_rmm-2d\C_Work\31_RMM_FMmodel\computations\model_setup\run_207\results\RMM_dflowfm_0000_map.nc',
                ]

for file_nc in file_nc_list:
    print('processing %s'%(os.path.basename(file_nc)))
    
    if 'cb_3d_map' in file_nc:
        timestep = 3
        layer = 5
        clim_bl = None
        clim_wl = [-0.5,1]
        clim_sal = None
        clim_tem = None
        crs = None
        file_nc_fou = None
        file_nc_rst = None
    elif 'Grevelingen-FM_0000_map' in file_nc:
        timestep = 3
        layer = 33
        clim_bl = None
        clim_wl = [-0.5,1]
        clim_sal = [28,30.2]
        clim_tem = [4,10]
        crs = "EPSG:28992"
        file_nc_fou = None
        file_nc_rst = None
    elif 'RMM_dflowfm_0000_map' in file_nc:
        timestep = 50
        layer = None
        clim_bl = [-10,10]
        clim_wl = [-2,2]
        clim_sal = None
        clim_tem = None
        crs = "EPSG:28992"
        file_nc_fou = os.path.join(dir_testinput,r'DFM_fou_RMM\RMM_dflowfm_0000_fou.nc')
        file_nc_rst = None # os.path.join(dir_testinput,r'DFM_fou_RMM\RMM_dflowfm_0006_20131127_000000_rst.nc') #not corresponding to current modelrun/grid anymore
    else:
        raise Exception('ERROR: no settings provided for this mapfile')
    

    #get ugrid data, vars informatin and grid units (latter from bedlevel coordinates)
    vars_pd, dims_pd = get_ncvardimlist(file_nc=file_nc)
    ugrid_all = get_netdata(file_nc=file_nc)#,multipart=False)
    data_frommap_dummy = get_ncmodeldata(file_nc=file_nc, varname='mesh2d_flowelem_bl', multipart=False)
    coordnames_xy = data_frommap_dummy.var_ncattrs['coordinates'].split()
    varcoords_x = get_ncmodeldata(file_nc=file_nc, varname=coordnames_xy[0],multipart=False)
    gridunits = varcoords_x.var_ncattrs['units']

    
    print('plot grid from mapdata')
    fig, ax = plt.subplots()
    pc = plot_netmapdata(ugrid_all.verts, values=None, ax=None, linewidth=0.5, color="crimson", facecolor="None")
    ax.set_xlabel('x [%s]'%(gridunits))
    ax.set_ylabel('y [%s]'%(gridunits))
    ax.set_aspect('equal')
    fig.tight_layout()
    fig.savefig(os.path.join(dir_output,'%s_grid'%(os.path.basename(file_nc).replace('.',''))))


    print('plot grid and bedlevel (constantvalue, 1 dim)')
    #iT = 3 #for iT in range(10):
    data_frommap = get_ncmodeldata(file_nc=file_nc, varname='mesh2d_flowelem_bl')#, multipart=False)
    fig, ax = plt.subplots()
    pc = plot_netmapdata(ugrid_all.verts, values=data_frommap, ax=None, linewidth=0.5, cmap="jet", edgecolor='face')
    pc.set_clim(clim_bl)
    cbar = fig.colorbar(pc, ax=ax)
    cbar.set_label('%s [%s]'%(data_frommap.var_ncattrs['long_name'], data_frommap.var_ncattrs['units']))
    ax.set_xlabel('x [%s]'%(gridunits))
    ax.set_ylabel('y [%s]'%(gridunits))
    ax.set_aspect('equal')
    fig.tight_layout()
    fig.savefig(os.path.join(dir_output,'%s_mesh2d_flowelem_bl'%(os.path.basename(file_nc).replace('.',''))))
    if crs is not None:
        """
        https://contextily.readthedocs.io/en/latest/reference.html
        https://contextily.readthedocs.io/en/latest/intro_guide.html
        ctx.add_basemap() defaults:
        source: None defaults to ctx.providers.Stamen.Terrain
        crs: coordinate reference system (CRS). If None (default), no warping is performed and the original Spherical Mercator (EPSG:3857) is used.
        More complex basemaps/coastlines are available in dfm_tools.net_nc.plot_background()
        """
        source = ctx.providers.Esri.WorldImagery # ctx.providers.Stamen.Terrain (default), ctx.providers.CartoDB.Voyager, ctx.providers.NASAGIBS.ViirsEarthAtNight2012, ctx.providers.Stamen.Watercolor
        ctx.add_basemap(ax, source=source, crs=crs, attribution=False)
        fig.savefig(os.path.join(dir_output,'%s_mesh2d_flowelem_bl_withbasemap'%(os.path.basename(file_nc).replace('.',''))))
        
    
    print('plot grid and values from mapdata (salinity on layer, 3dim, on cell centers)')
    data_frommap = get_ncmodeldata(file_nc=file_nc, varname='mesh2d_sa1', timestep=timestep, layer=layer)#, multipart=False)
    fig, ax = plt.subplots()
    pc = plot_netmapdata(ugrid_all.verts, values=data_frommap, ax=None, linewidth=0.5, cmap="jet", edgecolor='face')
    pc.set_clim(clim_sal)
    cbar = fig.colorbar(pc, ax=ax)
    cbar.set_label('%s [%s]'%(data_frommap.var_ncattrs['long_name'], data_frommap.var_ncattrs['units']))
    ax.set_xlabel('x [%s]'%(gridunits))
    ax.set_ylabel('y [%s]'%(gridunits))
    ax.set_aspect('equal')
    fig.tight_layout()
    fig.savefig(os.path.join(dir_output,'%s_mesh2d_sa1'%(os.path.basename(file_nc).replace('.',''))))


    print('plot grid and values from mapdata on net links (water/wind velocity on cell edges)')
    if 'mesh2d_u1' in vars_pd['nc_varkeys'].tolist():
        varname_edge = 'mesh2d_u1'
    else: #RMM does not contain mesh2d_u1 variable, so alternative is used
        varname_edge = 'mesh2d_windxu'
    data_frommap = get_ncmodeldata(file_nc=file_nc, varname=varname_edge, timestep=timestep, layer=layer)#, multipart=False)
    fig, ax = plt.subplots()
    pc = plot_netmapdata(ugrid_all.edge_verts, values=data_frommap, ax=None, linewidth=0.5, cmap="jet", edgecolor='face')
    #pc.set_clim(0,5)
    cbar = fig.colorbar(pc, ax=ax)
    cbar.set_label('%s [%s]'%(data_frommap.var_ncattrs['long_name'], data_frommap.var_ncattrs['units']))
    ax.set_xlabel('x [%s]'%(gridunits))
    ax.set_ylabel('y [%s]'%(gridunits))
    ax.set_aspect('equal')
    fig.tight_layout()
    fig.savefig(os.path.join(dir_output,'%s_%s_edges'%(os.path.basename(file_nc).replace('.',''),varname_edge)))
    
        
    if file_nc_fou is not None:
        #RMM foufile met quivers
        vars_pd, dims_pd = get_ncvardimlist(file_nc=file_nc_fou)
        #stations_pd = get_hisstationlist(file_nc,varname='waterlevel')
        
        ugrid_all_fou = get_netdata(file_nc=file_nc_fou)
        ux_mean = get_ncmodeldata(file_nc=file_nc_fou, varname='mesh2d_fourier001_mean')
        uy_mean = get_ncmodeldata(file_nc=file_nc_fou, varname='mesh2d_fourier002_mean')
        magn_mean = np.sqrt(ux_mean**2+uy_mean**2)
        #uc_mean = get_ncmodeldata(file_nc=file_nc, varname='mesh2d_fourier003_mean')
        #uc_max = get_ncmodeldata(file_nc=file_nc, varname='mesh2d_fourier004_max')
        facex = get_ncmodeldata(file_nc=file_nc_fou, varname='mesh2d_face_x')
        facey = get_ncmodeldata(file_nc=file_nc_fou, varname='mesh2d_face_y')
        
        X,Y,U = scatter_to_regulargrid(xcoords=facex, ycoords=facey, ncellx=60, ncelly=35, values=ux_mean)
        X,Y,V = scatter_to_regulargrid(xcoords=facex, ycoords=facey, ncellx=60, ncelly=35, values=uy_mean)
        
        #thinning = 3
        fig1,ax1 = plt.subplots(figsize=(9,5))
        pc1 = plot_netmapdata(ugrid_all_fou.verts, magn_mean, edgecolor='face')
        #ax1.quiver(facex[::thinning], facey[::thinning], ux_mean[::thinning], uy_mean[::thinning], color='w',scale=20)#,width=0.005)#, edgecolor='face', cmap='jet')
        ax1.quiver(X,Y,U,V, color='w',scale=5)#,width=0.005)#, edgecolor='face', cmap='jet')
        pc1.set_clim([0,0.10])
        #ax1.set_title('sqrt(x^2+y^2)\nx=%s\ny=%s'%(ux_mean.var_ncvarobject.long_name,uy_mean.var_ncvarobject.long_name))
        ax1.set_aspect('equal')
        ax1.set_xlabel('x [%s]'%(gridunits))
        ax1.set_ylabel('y [%s]'%(gridunits))
        cbar = fig1.colorbar(pc1)
        cbar.set_label('residuele stroming [m/s]')
        fig1.tight_layout()
        fig1.savefig(os.path.join(dir_output,os.path.basename(file_nc).replace('.','')))
        
        
    if file_nc_rst is not None:
        #RMM rst file
        vars_pd, dims_pd = get_ncvardimlist(file_nc=file_nc_rst)
        #ugrid = get_netdata(file_nc=file_nc_rst)
        ugrid_0006 = get_netdata(file_nc=file_nc.replace('0000','0006'),multipart=False) #cannot be retrieved from rstfile ("mesh2d_node_x" and "mesh2d_node_y" missing), so get from partitioned mapfile or partitioned network file
        data_s1 = get_ncmodeldata(file_nc=file_nc_rst, varname='s1',timestep=0, multipart=False)
        
        fig, ax = plt.subplots(figsize=(10,4))
        pc = plot_netmapdata(ugrid_0006.verts, values=data_s1, ax=None, linewidth=0.5, cmap="jet",edgecolor='face')
        pc.set_clim([0,2])
        fig.colorbar(pc)
        ax.set_xlabel('x [%s]'%(gridunits))
        ax.set_ylabel('y [%s]'%(gridunits))
        ax.set_aspect('equal')
        fig.tight_layout()
        plt.savefig(os.path.join(dir_output,os.path.basename(file_nc_rst).replace('.','')))
        
        #keep this test to check flow link handling by get_nc.py (does not have domains)
        #ugrid_FlowLink_xu = get_ncmodeldata(file_nc=file_nc_rst, varname='FlowLink_xu', multipart=False)
        #ugrid_FlowLink_yu = get_ncmodeldata(file_nc=file_nc_rst, varname='FlowLink_yu', multipart=False)
        data_q1 = get_ncmodeldata(file_nc=file_nc_rst, varname='q1',timestep=0, multipart=False)
    