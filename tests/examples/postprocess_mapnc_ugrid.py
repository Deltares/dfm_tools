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


dir_testinput = r'c:\DATA\dfm_tools_testdata'
dir_output = '.'

file_nc_list = [os.path.join(dir_testinput,'DFM_sigma_curved_bend\\DFM_OUTPUT_cb_3d\\cb_3d_map.nc'), #sigmalayer
                os.path.join(dir_testinput,'DFM_3D_z_Grevelingen','computations','run01','DFM_OUTPUT_Grevelingen-FM','Grevelingen-FM_0*_map.nc'), #zlayer
                r'p:\1204257-dcsmzuno\2006-2012\3D-DCSM-FM\A18b_ntsu1\DFM_OUTPUT_DCSM-FM_0_5nm\DCSM-FM_0_5nm_0*_map.nc', #fullgrid
                r'p:\11206813-006-kpp2021_rmm-2d\C_Work\31_RMM_FMmodel\computations\model_setup\run_207\results\RMM_dflowfm_0*_map.nc', #2D model
                r'p:\archivedprojects\11203379-005-mwra-updated-bem\03_model\02_final\A72_ntsu0_kzlb2\DFM_OUTPUT_MB_02\MB_02_0*_map.nc',
                ]


for file_nc in file_nc_list:
    print('processing %s'%(os.path.basename(file_nc)))
    basename = os.path.basename(file_nc).replace('.','').replace('_0*_','_0000_')
    
    if 'cb_3d_map' in file_nc:
        timestep = 72
        layno = 5
        sel_slice_x, sel_slice_y = slice(1500,3500), slice(1000,3500)
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
        file_nc_fou = None
    elif 'Grevelingen' in file_nc:
        timestep = 3
        layno = 33 #35 is top
        sel_slice_x, sel_slice_y = slice(50000,55000), slice(None,424000)
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
        file_nc_fou = None
    elif 'DCSM-FM_0_5nm' in file_nc:
        timestep = 365
        layno = 45
        sel_slice_x, sel_slice_y = slice(0,5), slice(50,55)
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
        file_nc_fou = None
    elif 'RMM_dflowfm' in file_nc:
        timestep = 365 #50
        layno = None
        sel_slice_x, sel_slice_y = slice(None,None), slice(None,None)
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
        file_nc_fou = os.path.join(dir_testinput,r'DFM_fou_RMM\RMM_dflowfm_0*_fou.nc')
        fou_varname_u, fou_varname_v = 'mesh2d_fourier001_mean', 'mesh2d_fourier002_mean'
    elif 'MB_02_' in file_nc:
        timestep = 10
        layno = 45
        sel_slice_x, sel_slice_y = slice(None,None), slice(None,None)
        #provide xy order, so lonlat
        line_array = np.array([[-71.81578813,  42.68460697],
                               [-65.2535983 ,  41.8699903 ]])
        val_ylim = [-600,1]
        clim_bl = [-500,0]
        clim_sal = [25,36]
        crs = "EPSG:4326"
        file_nc_fou = r'p:\archivedprojects\11203379-005-mwra-updated-bem\03_model\02_final\A72_ntsu0_kzlb2\DFM_OUTPUT_MB_02_fou\MB_02_0*_fou.nc'
        fou_varname_u, fou_varname_v = 'mesh2d_fourier027_mean', 'mesh2d_fourier040_mean'
    else:
        raise Exception('ERROR: no settings provided for this mapfile')
    
    
    data_frommap_merged = dfmt.open_partitioned_dataset(file_nc)
    vars_pd = dfmt.get_ncvarproperties(data_frommap_merged)
    
    
    print('plot grid from mapdata')
    fig, ax = plt.subplots()
    pc = data_frommap_merged.ugrid.grid.plot(edgecolor='crimson', linewidth=0.5)
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_aspect('equal')
    fig.tight_layout()
    fig.savefig(os.path.join(dir_output,f'{basename}_grid'))
    
    
    print('plot grid and bedlevel (constantvalue, 1 dim)')
    #get bedlevel and create plot with ugrid and cross section line
    fig, ax_input = plt.subplots()
    pc = data_frommap_merged['mesh2d_flowelem_bl'].ugrid.plot(edgecolor='face',cmap='jet') #TODO: should work even better with edgecolor='none', but that results in seethrough edges anyway, report to matplotlib?
    pc.set_clim(clim_bl)
    ax_input.set_aspect('equal')
    line, = ax_input.plot([], [],'o-') # empty line
    linebuilder = dfmt.LineBuilder(line) #this makes it possible to interactively click a line in the bedlevel figure. Use linebuilder.line_array as alternative line_array
    ax_input.plot(line_array[0,0],line_array[0,1],'bx',linewidth=3,markersize=10)
    ax_input.plot(line_array[:,0],line_array[:,1],'b',linewidth=3)
    fig.tight_layout()
    fig.savefig(os.path.join(dir_output,f'{basename}_mesh2d_flowelem_bl'))
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
        ctx.add_basemap(ax=ax_input, source=source, crs=crs, attribution=False)
        fig.savefig(os.path.join(dir_output,f'{basename}_mesh2d_flowelem_bl_withbasemap'))
    
    
    #ugrid sel via x/y
    data_frommap_merged_sel = data_frommap_merged.ugrid.sel(x=sel_slice_x,y=sel_slice_y)
    fig, ax = plt.subplots()
    pc = data_frommap_merged_sel['mesh2d_flowelem_bl'].ugrid.plot(ax=ax, linewidth=0.5, edgecolors='face', cmap='jet')
    pc.set_clim(clim_bl)
    
    
    print('plot bedlevel as polycollection, contourf, contour')
    #create fancy plots, more options at https://deltares.github.io/xugrid/examples/plotting.html
    if clim_bl is None:
        vmin = vmax = None
    else:
        vmin, vmax = clim_bl #TODO: vmin/vmax are necessary upon plot initialization (instead of pc.set_clim(clim_bl)) for proper colorbar, is this an xugrid or matplotlib issue?
    fig, (ax1,ax2,ax3) = plt.subplots(3,1,figsize=(6,9))
    pc = data_frommap_merged['mesh2d_flowelem_bl'].ugrid.plot(ax=ax1, linewidth=0.5, edgecolors='face', cmap='jet', vmin=vmin, vmax=vmax)
    pc = data_frommap_merged['mesh2d_flowelem_bl'].ugrid.plot.contourf(ax=ax2, levels=11, cmap='jet', vmin=vmin, vmax=vmax)
    if 'cb_3d_map' not in file_nc: #TODO: cb_3d_map fails on contour with "UserWarning: No contour levels were found within the data range." (because all bedlevels are -5m) >> colorbar gives error, is this an xugrid or matplotlib issue?
        pc = data_frommap_merged['mesh2d_flowelem_bl'].ugrid.plot.contour(ax=ax3, levels=11, cmap='jet', vmin=vmin, vmax=vmax, add_colorbar=True)
    fig.tight_layout()
    fig.savefig(os.path.join(dir_output,f'{basename}_gridbedcontour'))
    
    
    #filter for dry cells
    bool_drycells = data_frommap_merged['mesh2d_s1']==data_frommap_merged['mesh2d_flowelem_bl']
    data_frommap_merged['mesh2d_s1_filt'] = data_frommap_merged['mesh2d_s1'].where(~bool_drycells)
    print('plot grid and values from mapdata (waterlevel on layer, 2dim, on cell centers)')
    fig, ax = plt.subplots()
    pc = data_frommap_merged['mesh2d_s1_filt'].isel(time=timestep).ugrid.plot(edgecolor='face',cmap='jet')
    ax.set_aspect('equal')
    fig.tight_layout()
    fig.savefig(os.path.join(dir_output,f'{basename}_mesh2d_s1_filt'))
    
    
    print('calculating and plotting cross section')
    crs_tstart = dt.datetime.now() #start timer
    xr_crs_ugrid = dfmt.polyline_mapslice(data_frommap_merged, line_array, timestep=timestep)
    fig, ax = plt.subplots()
    xr_crs_ugrid['mesh2d_sa1'].ugrid.plot(cmap='jet')
    fig.tight_layout()
    plt.savefig(os.path.join(dir_output,f'{basename}_crossect'))
    print(f'calculating and plotting cross section finished in {dt.datetime.now()-crs_tstart}')
    
    
    print('plot grid and values from mapdata (salinity on layer, 3dim, on cell centers), on layer')
    fig, ax = plt.subplots()
    pc = data_frommap_merged['mesh2d_sa1'].isel(time=timestep, mesh2d_nLayers=layno, nmesh2d_layer=layno, missing_dims='ignore').ugrid.plot(edgecolor='face',cmap='jet') #missing_dims='ignore' ignores .isel() on mesh2d_nLayers/nmesh2d_layer if that dimension is not present
    pc.set_clim(clim_sal)
    ax.set_aspect('equal')
    fig.tight_layout()
    fig.savefig(os.path.join(dir_output,f'{basename}_mesh2d_sa1'))
    
    
    print('plot grid and values from mapdata (salinity on layer, 3dim, on cell centers), on fixed depth(s)')
    data_frommap_timesel = data_frommap_merged.isel(time=timestep)
    data_frommap_timesel_atdepths = dfmt.get_Dataset_atdepths(data_xr=data_frommap_timesel, depths=-4, reference='z0') #depth w.r.t. z0/waterlevel/bedlevel (also possible to provide list of floats)
    fig, ax = plt.subplots()
    pc = data_frommap_timesel_atdepths['mesh2d_sa1'].ugrid.plot(edgecolor='face',cmap='jet') #TODO: dask\array\reductions.py:640: RuntimeWarning: All-NaN slice encountered
    pc.set_clim(clim_sal)
    ax.set_aspect('equal')
    fig.tight_layout()
    fig.savefig(os.path.join(dir_output,f'{basename}_mesh2d_sa1_onfixeddepth'))
    
    
    print('plot grid and values from mapdata on net links (water/wind velocity on cell edges)')
    if 'mesh2d_u1' in data_frommap_merged.data_vars: #for cb_3d_map and Grevelingen
        fig, ax = plt.subplots()
        pc = data_frommap_merged['mesh2d_u1'].isel(time=timestep, mesh2d_nLayers=layno, nmesh2d_layer=layno, missing_dims='ignore').ugrid.plot(cmap='jet') #missing_dims='ignore' ignores .isel() on mesh2d_nLayers/nmesh2d_layer if that dimension is not present
        ax.set_aspect('equal')
        fig.tight_layout()
        fig.savefig(os.path.join(dir_output,f'{basename}_edges'))
    
    
    if file_nc_fou is not None:
        #RMM/MBAY foufile met quivers #TODO: maybe fancy xugridplotting can help out here? (imshow regrids to 500x500 dataset also)
        #pc = data_frommap_fou_atdepth[['mesh2d_ucx','mesh2d_ucy']].ugrid.plot.quiver(ax=ax,x='mesh2d_face_x',y='mesh2d_face_y',u='mesh2d_ucx',v='mesh2d_ucy') #TODO: quiver is now not possible: "AttributeError: 'UgridDatasetAccessor' object has no attribute 'plot'"
        #xugrid issue: https://github.com/Deltares/xugrid/issues/31 (plotting quiver on regridded dataset). If it works, also add to notebook (for mapfile)
        
        data_frommap_fou = dfmt.open_partitioned_dataset(file_nc_fou)
        vars_pd_fou = dfmt.get_ncvarproperties(data_frommap_fou)
        if 'mesh2d_nLayers' in data_frommap_fou.dims: #reduce layer dimension via isel/sel/interp. TODO: slicing over depth is not possible with dfmt.get_Dataset_atdepths(), since waterlevel is missing from file. (does it work for rstfiles?)
            data_frommap_fou = data_frommap_fou.set_index(mesh2d_nLayers='mesh2d_layer_z') #TODO: not supported for sigmalayers, zlayers is for some reason in foufile of this zsigma model (or not the case with a rerun?) TODO: should these not be coordinate variables to begin with? (zw/zcc are also coordinates)
            if 1:
                data_frommap_fou_atdepth = data_frommap_fou.isel(mesh2d_nLayers=-2) #second to last layer
            elif 0: #nearest
                data_frommap_fou_atdepth = data_frommap_fou.sel(mesh2d_nLayers=-4, method='nearest') #layer closest to z==-4m
            else: #interp
                data_frommap_fou_atdepth = data_frommap_fou.interp(mesh2d_nLayers=-4) #interp to -4m depth
        else:
            data_frommap_fou_atdepth = data_frommap_fou
        
        facex = data_frommap_fou_atdepth['mesh2d_face_x'].to_numpy()
        facey = data_frommap_fou_atdepth['mesh2d_face_y'].to_numpy()
        ux_mean = data_frommap_fou_atdepth[fou_varname_u]
        uy_mean = data_frommap_fou_atdepth[fou_varname_v]
        magn_mean_attrs = {'long_name':'residuele stroming', 'units':'m/s'}
        data_frommap_fou_atdepth['magn_mean'] = np.sqrt(ux_mean**2+uy_mean**2).assign_attrs(magn_mean_attrs)
        X,Y,U = dfmt.scatter_to_regulargrid(xcoords=facex, ycoords=facey, ncellx=60, ncelly=35, values=ux_mean.to_numpy())
        X,Y,V = dfmt.scatter_to_regulargrid(xcoords=facex, ycoords=facey, ncellx=60, ncelly=35, values=uy_mean.to_numpy())
        
        fig,ax = plt.subplots(figsize=(9,5))
        pc = data_frommap_fou_atdepth['magn_mean'].ugrid.plot(edgecolor='face')
        ax.quiver(X,Y,U,V, color='w',scale=5)#,width=0.005)#, edgecolor='face', cmap='jet')
        pc.set_clim(0,0.10)
        ax.set_aspect('equal')
        fig.tight_layout()
        fig.savefig(os.path.join(dir_output,f'{basename}_fou'))
    
        