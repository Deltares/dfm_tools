# -*- coding: utf-8 -*-
"""
Created on Wed Oct 27 14:42:52 2021

@author: veenstra
"""

import os
import matplotlib.pyplot as plt
plt.close('all')
import numpy as np
import contextily as ctx
import datetime as dt
import dfm_tools as dfmt


dir_output = '.'

file_nc_list = [dfmt.data.fm_curvedbend_map(return_filepath=True), # sigmalayer
                dfmt.data.fm_grevelingen_map(return_filepath=True), # zlayer
                r'p:\1204257-dcsmzuno\2006-2012\3D-DCSM-FM\A18b_ntsu1\DFM_OUTPUT_DCSM-FM_0_5nm\DCSM-FM_0_5nm_0*_map.nc', # szigma fullgrid
                dfmt.data.fm_westernscheldt_map(return_filepath=True), # zsigma model without fullgrid output but with new ocean_sigma_z_coordinate variable
                # r'p:\archivedprojects\11206813-006-kpp2021_rmm-2d\C_Work\31_RMM_FMmodel\computations\model_setup\run_207\results\RMM_dflowfm_0*_map.nc', # 2D model
                # r'p:\archivedprojects\11203379-005-mwra-updated-bem\03_model\02_final\A72_ntsu0_kzlb2\DFM_OUTPUT_MB_02\MB_02_0*_map.nc', # zlayer
                ]

for file_nc in file_nc_list:
    plt.close('all')

    # defaults (can be overwitten by file specific settings)
    remove_edges = False
    rename_dict = None
    
    if 'cb_3d_map' in file_nc:
        timestep = 72
        layno = 5
        sel_slice_x, sel_slice_y = slice(1500,3500), slice(1000,3500)
        section_y = 2000
        line_array = np.array([[ 185.08667065, 2461.11775254],
                               [2934.63837418, 1134.16019127]])
        line_array = np.array([[ 104.15421399, 2042.7077107 ],
                               [2913.47878063, 2102.48057382]])
        #line_array = np.array([[2084.67741935, 3353.02419355], #with linebend in cell en with line crossing same cell twice
        #   [2255.79637097, 3307.15725806],
        #   [2222.27822581, 3206.60282258],
        #   [2128.78024194, 3266.58266129]])
        val_ylim = None
        clim_bl = None
        clim_sal = None
        crs = None
        raster_res = 200
        umag_clim = None
    elif 'Grevelingen' in file_nc:
        timestep = 3
        layno = 33 #35 is top
        sel_slice_x, sel_slice_y = slice(50000,55000), slice(None,424000)
        section_y = 418000
        line_array = np.array([[ 56267.59146475, 415644.67447155],
                               [ 64053.73427496, 419407.58239502]])
        line_array = np.array([[ 53181.96942503, 424270.83361629],
                               [ 55160.15232593, 416913.77136685]])
        #line_array = np.array([[ 53181.96942503, 424270.83361629],
        #                       [ 55160.15232593, 416913.77136685],
        #                       [ 65288.15232593, 419360.77136685]])
        val_ylim = [-25,5]
        clim_bl = [-40,10]
        clim_sal = [28,30.2]
        crs = "EPSG:28992"
        raster_res = 1000
        umag_clim = (None,0.1)
    elif 'DCSM-FM_0_5nm' in file_nc:
        remove_edges = True
        
        timestep = 365
        layno = 45
        sel_slice_x, sel_slice_y = slice(0,5), slice(50,55)
        section_y = 52.5
        #provide xy order, so lonlat
        line_array = np.array([[ 0.97452229, 51.13407643],
                               [ 1.89808917, 50.75191083]])
        line_array = np.array([[10.17702481, 57.03663877], #dummy for partition 0000
                               [12.38583134, 57.61284917]])
        line_array = np.array([[ 8.92659074, 56.91538014],
                               [ 8.58447136, 58.66874192]])
        val_ylim = [-600,1]
        clim_bl = [-500,0]
        clim_sal = [25,36]
        crs = "EPSG:4326"
        raster_res = 0.3
        umag_clim = (None,1)
    elif 'westerscheldt01_0subst_map' in file_nc:
        remove_edges = True
        rename_dict = {'mesh2d_ucmag':'mesh2d_sa1'}
        
        timestep = 1
        layno = -2
        sel_slice_x, sel_slice_y = slice(None,None), slice(None,None)
        section_y = 380000
        line_array = np.array([[19108.74,386404.0],
                               [40255.92,377797.6],
                               [48739.38,376322.2],
                               [52796.69,377428.7],
                               [55255.67,382838.5]])
        clim_bl = None
        clim_sal = None
        crs = "EPSG:28992"
        raster_res = 2500
        umag_clim = None
    elif 'RMM_dflowfm' in file_nc:
        timestep = 365 #50
        layno = None
        sel_slice_x, sel_slice_y = slice(None,None), slice(None,None)
        section_y = 440000
        #provide xy order, so lonlat
        line_array = np.array([[ 65655.72699961, 444092.54776465],
                               [ 78880.42720631, 435019.78832052]])
        line_array = np.array([[ 52444.56849912, 434039.27970214], #HVSL
                               [ 61304.25484967, 430703.86837017],
                               [ 62164.16558369, 428619.23628769]])
        line_array = np.array([[ 61013.8966525 , 446291.69129373], #NWW
                               [ 67151.68543524, 444096.96681991],
                               [ 69011.62143001, 442981.00522304],
                               [ 72210.71134101, 440302.69739058],
                               [ 74405.43581484, 438889.14603455],
                               [ 75632.99357138, 437401.19723874],
                               [ 79018.07708186, 435169.27404501],
                               [ 81324.39771538, 434536.89580679],
                               [ 82923.94267088, 434611.29324658],
                               [ 84449.09018659, 435132.07532512],
                               [ 86606.61594052, 434685.69068637],
                               [ 88689.74425466, 435355.26764449],
                               [ 90772.8725688 , 434983.28044554],
                               [ 91926.03288556, 435132.07532512]])
        val_ylim = None
        clim_bl = [-10,10]
        clim_sal = None
        crs = "EPSG:28992"
        raster_res = 2500
        umag_clim = (None,0.5)
    elif 'MB_02_' in file_nc:
        timestep = 10
        layno = 45
        sel_slice_x, sel_slice_y = slice(None,None), slice(None,None)
        section_y = 41
        #provide xy order, so lonlat
        line_array = np.array([[-71.81578813,  42.68460697],
                               [-65.2535983 ,  41.8699903 ]])
        val_ylim = [-600,1]
        clim_bl = [-500,0]
        clim_sal = [25,36]
        crs = "EPSG:4326"
        raster_res = 0.3
        umag_clim = (None,0.8)
    else:
        raise KeyError('ERROR: no settings provided for this mapfile')
    
    
    print('processing %s'%(os.path.basename(file_nc)))
    basename = os.path.basename(file_nc).replace('.','').replace('_0*_','_0000_')
    
    uds = dfmt.open_partitioned_dataset(file_nc, remove_edges=remove_edges)
    
    # optionally rename variables to allow for hardcoded plotting
    if rename_dict is not None:
        uds = uds.rename(rename_dict)
    
    # get dataframe of variables
    vars_pd = dfmt.get_ncvarproperties(uds)
    
    
    print('plot grid from mapdata')
    fig, ax = plt.subplots()
    pc = uds.grid.plot(edgecolor='crimson', linewidth=0.5)
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_aspect('equal')
    fig.tight_layout()
    fig.savefig(os.path.join(dir_output,f'{basename}_grid'))
    
    
    print('plot bedlevel')
    #get bedlevel and create plot with ugrid and cross section line
    fig, ax_input = plt.subplots()
    pc = uds['mesh2d_flowelem_bl'].ugrid.plot(cmap='jet') #TODO: default is edgecolor='face', should work even better with edgecolor='none', but that results in seethrough edges anyway, report to matplotlib?
    pc.set_clim(clim_bl)
    ax_input.set_aspect('equal')
    # line_array is defined above, alternatively click a cross-section line_array in the figure interactively with dfmt.LineBuilder
    # line_array = dfmt.LineBuilder(ax=ax_input).line_array
    ax_input.plot(line_array[0,0],line_array[0,1],'bx',linewidth=3,markersize=10)
    ax_input.plot(line_array[:,0],line_array[:,1],'b',linewidth=3)
    fig.tight_layout()
    fig.savefig(os.path.join(dir_output,f'{basename}_mesh2d_flowelem_bl'))
    if crs is not None:
        """
        https://contextily.readthedocs.io/en/latest/reference.html
        https://contextily.readthedocs.io/en/latest/intro_guide.html
        ctx.add_basemap() defaults:
        source: None defaults to ctx.providers.Stamen.Terrain, suggestions: ctx.providers.Stamen.Terrain (default), ctx.providers.Esri.WorldStreetMap, ctx.providers.Esri.WorldImagery, ctx.providers.CartoDB.Voyager, ctx.providers.NASAGIBS.ViirsEarthAtNight2012, ctx.providers.Stamen.Watercolor
        crs: coordinate reference system (CRS). If None (default), no warping is performed and the original Spherical Mercator (EPSG:3857) is used.
        """
        source = ctx.providers.Esri.WorldImagery
        ctx.add_basemap(ax=ax_input, source=source, crs=crs, attribution=False)
        fig.savefig(os.path.join(dir_output,f'{basename}_mesh2d_flowelem_bl_withbasemap'))
    
    
    print('plot bedlevel in different coordinate systems')
    if crs == 'EPSG:28992':
        to_crs = 'EPSG:4326'
        uds.ugrid.set_crs(crs)
        uds_wgs84 = uds.ugrid.to_crs(to_crs)
        fig, (ax1,ax2) = plt.subplots(2,1,figsize=(7,8))
        uds["mesh2d_waterdepth"].isel(time=0).ugrid.plot(ax=ax1)
        ctx.add_basemap(ax=ax1, source=None, crs=crs, attribution=False)
        uds_wgs84["mesh2d_waterdepth"].isel(time=0).ugrid.plot(ax=ax2)
        ctx.add_basemap(ax=ax2, source=None, crs=to_crs, attribution=False)
        fig.tight_layout()
        fig.savefig(os.path.join(dir_output,f'{basename}_convertedcoords'))
    
    
    #ugrid sel via x/y
    uds_sel = uds.ugrid.sel(x=sel_slice_x,y=sel_slice_y)
    fig, ax = plt.subplots()
    pc = uds_sel['mesh2d_flowelem_bl'].ugrid.plot(ax=ax, linewidth=0.5, cmap='jet')
    pc.set_clim(clim_bl)
    fig.tight_layout()
    fig.savefig(os.path.join(dir_output,f'{basename}_selxyslice'))
    
    
    print('plot bedlevel as polycollection, contourf, contour, rasterized')
    #create fancy plots, more options at https://deltares.github.io/xugrid/examples/plotting.html
    if clim_bl is None:
        vmin = vmax = None
    else:
        vmin, vmax = clim_bl # vmin/vmax are necessary upon plot initialization (instead of pc.set_clim(clim_bl)) for proper colorbar, this is also matplotlib behaviour
    fig, ((ax1,ax2),(ax3,ax4)) = plt.subplots(2,2,figsize=(12,7),sharex=True,sharey=True)
    pc = uds['mesh2d_flowelem_bl'].ugrid.plot(ax=ax1, linewidth=0.5, cmap='jet', vmin=vmin, vmax=vmax)
    pc = uds['mesh2d_flowelem_bl'].ugrid.plot.contourf(ax=ax2, levels=11, cmap='jet', vmin=vmin, vmax=vmax)
    pc = uds['mesh2d_flowelem_bl'].ugrid.plot.contour(ax=ax3, levels=11, cmap='jet', vmin=vmin, vmax=vmax, add_colorbar=True)
    bl_raster = dfmt.rasterize_ugrid(uds['mesh2d_flowelem_bl'],resolution=raster_res) #rasterize ugrid uds/uda
    pc = bl_raster.plot(ax=ax4, cmap='jet', vmin=vmin, vmax=vmax) #plot with non-ugrid method
    fig.tight_layout()
    fig.savefig(os.path.join(dir_output,f'{basename}_gridbedcontour'))
    
    
    #filter for dry cells
    bool_drycells = uds['mesh2d_s1']==uds['mesh2d_flowelem_bl']
    uds['mesh2d_s1_filt'] = uds['mesh2d_s1'].where(~bool_drycells)
    print('plot grid and values from mapdata (waterlevel on layer, 2dim, on cell centers)')
    fig, ax = plt.subplots()
    pc = uds['mesh2d_s1_filt'].isel(time=timestep).ugrid.plot(cmap='jet')
    ax.set_aspect('equal')
    fig.tight_layout()
    fig.savefig(os.path.join(dir_output,f'{basename}_mesh2d_s1_filt'))
    
    
    print('calculating and plotting cross section')
    crs_tstart = dt.datetime.now() #start timer
    xr_crs_ugrid = dfmt.polyline_mapslice(uds.isel(time=timestep), line_array)
    fig, ax = plt.subplots()
    xr_crs_ugrid['mesh2d_sa1'].ugrid.plot(cmap='jet')
    fig.tight_layout()
    plt.savefig(os.path.join(dir_output,f'{basename}_crossect'))
    print(f'calculating and plotting cross section finished in {dt.datetime.now()-crs_tstart}')
    
    
    print('plot grid and values from mapdata (salinity on layer, 3dim, on cell centers), on layer')
    fig, ax = plt.subplots()
    pc = uds['mesh2d_sa1'].isel(time=timestep, mesh2d_nLayers=layno, nmesh2d_layer=layno, missing_dims='ignore').ugrid.plot(cmap='jet') #missing_dims='ignore' ignores .isel() on mesh2d_nLayers/nmesh2d_layer if that dimension is not present
    pc.set_clim(clim_sal)
    ax.set_aspect('equal')
    fig.tight_layout()
    fig.savefig(os.path.join(dir_output,f'{basename}_mesh2d_sa1'))
    
    
    print('plot grid and values from mapdata (salinity on layer, 3dim, on cell centers), on fixed depth(s)')
    data_frommap_timesel = uds.isel(time=timestep)
    data_frommap_timesel_atdepths = dfmt.get_Dataset_atdepths(data_xr=data_frommap_timesel, depths=-4, reference='z0') #depth w.r.t. z0/waterlevel/bedlevel (also possible to provide list of floats)
    fig, ax = plt.subplots()
    pc = data_frommap_timesel_atdepths['mesh2d_sa1'].ugrid.plot(cmap='jet') #TODO: dask\array\reductions.py:640: RuntimeWarning: All-NaN slice encountered
    pc.set_clim(clim_sal)
    ax.set_aspect('equal')
    fig.tight_layout()
    fig.savefig(os.path.join(dir_output,f'{basename}_mesh2d_sa1_onfixeddepth'))
    
    
    print('plot grid and values from mapdata on net links (water/wind velocity on cell edges)')
    if 'mesh2d_u1' in uds.data_vars: #for cb_3d_map and Grevelingen
        fig, ax = plt.subplots()
        pc = uds['mesh2d_u1'].isel(time=timestep, mesh2d_nLayers=layno, nmesh2d_layer=layno, missing_dims='ignore').ugrid.plot(cmap='jet') #missing_dims='ignore' ignores .isel() on mesh2d_nLayers/nmesh2d_layer if that dimension is not present
        ax.set_aspect('equal')
        fig.tight_layout()
        fig.savefig(os.path.join(dir_output,f'{basename}_edges'))
    
    
    print('plot velocity magnitude and quiver')
    uds_quiv = uds.isel(time=-1, mesh2d_nLayers=-2, nmesh2d_layer=-2, missing_dims='ignore')
    varn_ucx, varn_ucy = 'mesh2d_ucx', 'mesh2d_ucy'
    magn_attrs = {'long_name':'velocity magnitude', 'units':'m/s'}
    uds_quiv['magn'] = np.sqrt(uds_quiv[varn_ucx]**2+uds_quiv[varn_ucy]**2).assign_attrs(magn_attrs)
    raster_quiv = dfmt.rasterize_ugrid(uds_quiv[[varn_ucx,varn_ucy]], resolution=raster_res)
    fig,ax = plt.subplots(figsize=(9,5))
    pc = uds_quiv['magn'].ugrid.plot()
    raster_quiv.plot.quiver(x='mesh2d_face_x',y='mesh2d_face_y',u=varn_ucx,v=varn_ucy,color='w',scale=25,add_guide=False)
    ax.set_aspect('equal')
    pc.set_clim(umag_clim)
    fig.tight_layout()
    fig.savefig(os.path.join(dir_output,f'{basename}_quiver'))
    
    
    #TODO: add hovmoller to notebook. x='x' does not work for spherical models, since it is sorted by 's'
    print('hovmoller plot: mean salinity over depth along section_y over time')
    #ax.axhline(y=section_y, color="red")
    uds_sel = uds.isel(time=slice(-30,None)).ugrid.sel(y=section_y)
    fig, ax = plt.subplots(figsize=(10,5.5))
    if 'nmesh2d_layer' in uds_sel.dims:
        slice_sa1 = uds_sel.mesh2d_sa1.mean(dim='nmesh2d_layer')
    elif 'mesh2d_nLayers' in uds_sel.dims:
        slice_sa1 = uds_sel.mesh2d_sa1.mean(dim='mesh2d_nLayers')
    else:
        slice_sa1 = uds_sel.mesh2d_sa1
    slice_sa1.plot(x='mesh2d_s',y='time')
    fig.tight_layout()
    fig.savefig(os.path.join(dir_output,f'{basename}_hovmoller'))
        
    
