# -*- coding: utf-8 -*-
"""
Created on Fri Nov 18 15:43:34 2022

@author: veenstra
"""


import os
import matplotlib.pyplot as plt
plt.close('all')
import numpy as np
import contextily as ctx
import datetime as dt
import xarray as xr
import dfm_tools as dfmt
import warnings

dir_testinput = r'c:\DATA\dfm_tools_testdata'
dir_output = '.'

file_nc_list = [#os.path.join(dir_testinput,'DFM_sigma_curved_bend\\DFM_OUTPUT_cb_3d\\cb_3d_map.nc'), #sigmalayer
                os.path.join(dir_testinput,'DFM_3D_z_Grevelingen','computations','run01','DFM_OUTPUT_Grevelingen-FM','Grevelingen-FM_0000_map.nc'), #zlayer
                #r'p:\1204257-dcsmzuno\2006-2012\3D-DCSM-FM\A18b_ntsu1\DFM_OUTPUT_DCSM-FM_0_5nm\DCSM-FM_0_5nm_0000_map.nc', #fullgrid
                #r'p:\11206813-006-kpp2021_rmm-2d\C_Work\31_RMM_FMmodel\computations\model_setup\run_207\results\RMM_dflowfm_0000_map.nc', #2D model
                ]

for file_nc in file_nc_list:
    print('processing %s'%(os.path.basename(file_nc)))
    basename = os.path.basename(file_nc).replace('.','')
    
    if 'cb_3d_map' in file_nc:
        timestep = 72
        layno = 5
        calcdist_fromlatlon = None
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
        calcdist_fromlatlon = None
        line_array = np.array([[ 56267.59146475, 415644.67447155],
                               [ 64053.73427496, 419407.58239502]])
        line_array = np.array([[ 53181.96942503, 424270.83361629],
                               [ 55160.15232593, 416913.77136685]])
        #line_array = np.array([[ 53181.96942503, 424270.83361629],
        #                       [ 55160.15232593, 416913.77136685],
        #                       [ 65288.15232593, 419360.77136685]])
        val_ylim = [-25,5]
        clim_bl = None
        clim_sal = [28,30.2]
        crs = "EPSG:28992"
        file_nc_fou = None
    elif 'DCSM-FM_0_5nm' in file_nc:
        timestep = 365
        layno = 45
        calcdist_fromlatlon = True
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
        calcdist_fromlatlon = None
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
        file_nc_fou = os.path.join(dir_testinput,r'DFM_fou_RMM\RMM_dflowfm_0000_fou.nc')
    else:
        raise Exception('ERROR: no settings provided for this mapfile')
    
    file_nc_star = file_nc.replace('_0000_','_0*_') #TODO: make this the default
    data_xr = xr.open_dataset(file_nc)
    data_frommap_merged = dfmt.open_partitioned_dataset(file_nc_star)
    
    #get ugrid data, vars informatin and grid units (latter from bedlevel coordinates)
    vars_pd = dfmt.get_ncvarproperties(file_nc=file_nc)
    ugrid_all = dfmt.get_netdata(file_nc=file_nc) 
    
    reproduce = False
        
    def polygon_intersect(ugrid_all, data_frommap_merged, line_array, optimize_dist=False, calcdist_fromlatlon=False): #TODO: copy of ugrid function, remove ugrid when xugrid works for all stuf
        #data_frommap_merged: xugrid dataset (contains ds and grid)
        
        import numpy as np
        from matplotlib.path import Path
        import shapely #separate import, since sometimes this works, while import shapely.geometry fails
        from shapely.geometry import LineString, Polygon, MultiLineString, Point
        from dfm_tools.get_nc import calc_dist_pythagoras, calc_dist_haversine
    
        print('defining celinlinebox')
        
        line_section = LineString(line_array)
        
        face_nos = data_frommap_merged.ugrid.grid.to_dataset().mesh2d_face_nodes
        bool_nonemptyfacenode = face_nos!=-1
        #face_nos = face_nos.where(face_nos!=-1) #replace nans
        face_nnodecoords_x = data_frommap_merged.ugrid.grid.to_dataset().mesh2d_node_x.isel(mesh2d_nNodes=face_nos).where(bool_nonemptyfacenode)
        face_nnodecoords_y = data_frommap_merged.ugrid.grid.to_dataset().mesh2d_node_y.isel(mesh2d_nNodes=face_nos).where(bool_nonemptyfacenode)
        if reproduce:
            ugrid_all_verts = ugrid_all.verts
        else:
            ugrid_all_verts = np.c_[face_nnodecoords_x.to_numpy()[...,np.newaxis],face_nnodecoords_y.to_numpy()[...,np.newaxis]]
        
        # verts_xmax = np.nanmax(ugrid_all.verts[:,:,0].data,axis=1)
        # verts_xmin = np.nanmin(ugrid_all.verts[:,:,0].data,axis=1)
        # verts_ymax = np.nanmax(ugrid_all.verts[:,:,1].data,axis=1)
        # verts_ymin = np.nanmin(ugrid_all.verts[:,:,1].data,axis=1)
        verts_xmax = np.nanmax(face_nnodecoords_x.data,axis=1)
        verts_xmin = np.nanmin(face_nnodecoords_x.data,axis=1)
        verts_ymax = np.nanmax(face_nnodecoords_y.data,axis=1)
        verts_ymin = np.nanmin(face_nnodecoords_y.data,axis=1)
        
        if not optimize_dist: #TODO: replace this with sel once it works
            cellinlinebox_all_bool = (((np.min(line_array[:,0]) <= verts_xmax) &
                                       (np.max(line_array[:,0]) >= verts_xmin)) &
                                      ((np.min(line_array[:,1]) <= verts_ymax) & 
                                       (np.max(line_array[:,1]) >= verts_ymin))
                                      )
        elif type(optimize_dist) in [int,float]: #not properly tested and documented
            #calculate angles wrt x axis
            angles_wrtx = []
            nlinecoords = line_array.shape[0]
            for iL in range(nlinecoords-1):
                dx = line_array[iL+1,0] - line_array[iL,0]
                dy = line_array[iL+1,1] - line_array[iL,1]
                angles_wrtx.append(np.rad2deg(np.arctan2(dy,dx)))
            angles_toprev = np.concatenate([[90],np.diff(angles_wrtx),[90]])
            angles_wrtx_ext = np.concatenate([[angles_wrtx[0]-90],np.array(angles_wrtx),[angles_wrtx[-1]+90]])
            angtot_wrtx = angles_wrtx_ext[:-1] + 0.5*(180+angles_toprev)
            #distance over xy-axis from original points
            dxynewpoints = optimize_dist * np.array([np.cos(np.deg2rad(angtot_wrtx)),np.sin(np.deg2rad(angtot_wrtx))]).T
            newpoints1 = line_array+dxynewpoints
            newpoints2 = line_array-dxynewpoints
            pol_inpol = np.concatenate([newpoints1, np.flip(newpoints2,axis=0)])
            pol_inpol_path = Path(pol_inpol)
            bool_all = []
            #for iC in range(ugrid_all.verts.shape[1]):
            for iC in range(face_nnodecoords_x.shape[1]):
                #data_arr_ic = ugrid_all.verts[:,iC,:]
                data_arr_ic = ugrid_all_verts[:,iC,:]
                #data_arr_ic = np.c_[face_nnodecoords_x.isel(mesh2d_nMax_face_nodes=iC).to_numpy(),face_nnodecoords_y.isel(mesh2d_nMax_face_nodes=iC).to_numpy()]
                test = pol_inpol_path.contains_points(data_arr_ic)
                bool_all.append(test)
            test_all = np.array(bool_all)
            cellinlinebox_all_bool = (test_all==True).any(axis=0)
        else:
            raise Exception('ERROR: invalid type for optimize_dist argument')
        
        #intersect_coords = np.empty((0,2,2))
        intersect_coords = np.empty((0,4))
        intersect_gridnos = np.empty((0),dtype=int) #has to be numbers, since a boolean is differently ordered
        #verts_inlinebox = ugrid_all.verts[cellinlinebox_all_bool,:,:]
        verts_inlinebox = ugrid_all_verts[cellinlinebox_all_bool,:,:]
        verts_inlinebox_nos = np.where(cellinlinebox_all_bool)[0]
        print('finding crossing flow links (can take a while if linebox over xy covers a lot of cells, %i of %i cells are being processed)'%(cellinlinebox_all_bool.sum(),len(cellinlinebox_all_bool)))
        
        for iP, pol_data in enumerate(verts_inlinebox):
            pol_shp = Polygon(pol_data[~np.isnan(pol_data).all(axis=1)])
            intersect_result = pol_shp.intersection(line_section)
            if isinstance(intersect_result,shapely.geometry.multilinestring.MultiLineString): #in the rare case that a cell (pol_shp) is crossed by multiple parts of the line
                intersect_result_multi = intersect_result
            elif isinstance(intersect_result,shapely.geometry.linestring.LineString): #if one linepart trough cell (ex/including node), make multilinestring anyway
                if intersect_result.coords == []: #when the line does not cross this cell, intersect_results.coords is an empty linestring and this cell can be skipped (continue makes forloop continue with next in line without finishing the rest of the steps for this instance)
                    continue
                elif len(intersect_result.coords.xy[0]) == 0: #for newer cartopy versions, when line does not cross this cell, intersect_result.coords.xy is (array('d'), array('d')), and both arrays in tuple have len 0.
                    continue
                intersect_result_multi = MultiLineString([intersect_result])
                
            for iLL, intesect_result_one in enumerate(intersect_result_multi.geoms): #loop over multilinestrings, will mostly only contain one linestring. Will be two if the line crosses a cell more than once.
                intersection_line = intesect_result_one.coords
                intline_xyshape = np.array(intersection_line.xy).shape
                #print('len(intersection_line.xy): %s'%([intline_xyshape]))
                for numlinepart_incell in range(1,intline_xyshape[1]): #is mostly 1, but more if there is a linebreakpoint in this cell (then there are two or more lineparts)
                    intersect_gridnos = np.append(intersect_gridnos,verts_inlinebox_nos[iP])
                    #intersect_coords = np.concatenate([intersect_coords,np.array(intersection_line.xy)[np.newaxis,:,numlinepart_incell-1:numlinepart_incell+1]],axis=0)
                    intersect_coords = np.concatenate([intersect_coords,np.array(intersection_line.xy).T[numlinepart_incell-1:numlinepart_incell+1].flatten()[np.newaxis]])
        
        if intersect_coords.shape[0] != len(intersect_gridnos):
            raise Exception('something went wrong, intersect_coords.shape[0] and len(intersect_gridnos) are not equal')
        
        import pandas as pd
        intersect_pd = pd.DataFrame(intersect_coords,index=intersect_gridnos,columns=['x1','y1','x2','y2'])
        intersect_pd.index.name = 'gridnumber'
        
        #TODO up to here could come from meshkernelpy
        
        print('calculating distance for all crossed cells, from first point of line (should not take long, but if it does, optimisation is needed)')
        nlinecoords = line_array.shape[0]
        nlinedims = len(line_array.shape)
        ncrosscellparts = len(intersect_pd)
        if nlinecoords<2 or nlinedims != 2:
            raise Exception('ERROR: line_array should at least contain two xy points [[x,y],[x,y]]')
        
        #calculate distance between celledge-linepart crossing (is zero when line iL crosses cell)
        distperline_tostart = np.zeros((ncrosscellparts,nlinecoords-1))
        distperline_tostop = np.zeros((ncrosscellparts,nlinecoords-1))
        linepart_length = np.zeros((nlinecoords))
        for iL in range(nlinecoords-1):
            #calculate length of lineparts
            line_section_part = LineString(line_array[iL:iL+2,:])
            if calcdist_fromlatlon:
                linepart_length[iL+1] = calc_dist_haversine(line_array[iL,0],line_array[iL+1,0],line_array[iL,1],line_array[iL+1,1])
            else:
                linepart_length[iL+1] = line_section_part.length
        
            #get distance between all lineparts and point (later used to calculate distance from beginpoint of closest linepart)
            for iP in range(ncrosscellparts):
                distperline_tostart[iP,iL] = line_section_part.distance(Point(intersect_coords[:,0][iP],intersect_coords[:,1][iP]))
                distperline_tostop[iP,iL] = line_section_part.distance(Point(intersect_coords[:,2][iP],intersect_coords[:,3][iP]))
        linepart_lengthcum = np.cumsum(linepart_length)
        cross_points_closestlineid = np.argmin(np.maximum(distperline_tostart,distperline_tostop),axis=1)
        intersect_pd['closestlineid'] = cross_points_closestlineid
        print('finished calculating distance for all crossed cells, from first point of line')
        
        if not calcdist_fromlatlon:
            crs_dist_starts = calc_dist_pythagoras(line_array[cross_points_closestlineid,0], intersect_coords[:,0], line_array[cross_points_closestlineid,1], intersect_coords[:,1]) + linepart_lengthcum[cross_points_closestlineid]
            crs_dist_stops = calc_dist_pythagoras(line_array[cross_points_closestlineid,0], intersect_coords[:,2], line_array[cross_points_closestlineid,1], intersect_coords[:,3]) + linepart_lengthcum[cross_points_closestlineid]
        else:
            crs_dist_starts = calc_dist_haversine(line_array[cross_points_closestlineid,0], intersect_coords[:,0], line_array[cross_points_closestlineid,1], intersect_coords[:,1]) + linepart_lengthcum[cross_points_closestlineid]
            crs_dist_stops = calc_dist_haversine(line_array[cross_points_closestlineid,0], intersect_coords[:,2], line_array[cross_points_closestlineid,1], intersect_coords[:,3]) + linepart_lengthcum[cross_points_closestlineid]
        intersect_pd['crs_dist_starts'] = crs_dist_starts
        intersect_pd['crs_dist_stops'] = crs_dist_stops
        intersect_pd['linepartlen'] = crs_dist_stops-crs_dist_starts
        intersect_pd = intersect_pd.sort_values('crs_dist_starts')
        
        #dimensions (gridnos, xy, firstsecond)
        print('done finding crossing flow links: %i of %i'%(len(intersect_gridnos),len(cellinlinebox_all_bool)))
        return intersect_pd
    
    
    
    def get_xzcoords_onintersection(data_frommap_merged, intersect_pd, timestep=None, multipart=None, varname=None):
        #check if all necessary arguments are provided
        if timestep is None:
            raise Exception('ERROR: argument timestep not provided, this is necessary to retrieve correct waterlevel or fullgrid output')
        
        varkeys_list = list(data_frommap_merged.variables.keys())
        dimn_layer = 'nmesh2d_layer'#dfmt.get_varname_fromnc(data_frommap_merged,'nmesh2d_layer',vardim='dim')
        if dimn_layer is None: #no layers, 2D model
            nlay = 1
        else:
            nlay = data_frommap_merged.dims[dimn_layer]
            
        intersect_gridnos = intersect_pd.index
        if 'mesh2d_flowelem_zw' in varkeys_list:
            print('layertype: fullgrid output') #TODO: still to be fixed
            zvals_interface_allfaces = get_ncmodeldata(file_nc, varname='mesh2d_flowelem_zw', timestep=timestep, multipart=multipart)
            zvals_interface = zvals_interface_allfaces[0,intersect_gridnos,:].T #timestep=0 since correct timestep was already retrieved with get_ncmodeldata. Transpose to make in line with 2D sigma dataset
        else: #no full grid output, so reconstruct
            data_frommap_merged_sel = data_frommap_merged.isel(time=timestep,mesh2d_nFaces=intersect_gridnos)
            #varn_mesh2d_s1 = get_varname_fromnc(data_nc,'mesh2d_s1',vardim='var')
            #data_frommap_wl3 = get_ncmodeldata(file_nc, varname=varn_mesh2d_s1, timestep=timestep, multipart=multipart)
            #data_frommap_wl3_sel = data_frommap_wl3[0,intersect_gridnos] #timestep=0 since correct timestep was already retrieved with get_ncmodeldata
            #varn_mesh2d_flowelem_bl = get_varname_fromnc(data_nc,'mesh2d_flowelem_bl',vardim='var')
            #data_frommap_bl = get_ncmodeldata(file_nc, varname=varn_mesh2d_flowelem_bl, multipart=multipart)
            #data_frommap_bl_sel = data_frommap_bl[intersect_gridnos]
            data_frommap_wl3_sel = data_frommap_merged_sel['mesh2d_s1'].to_numpy() #TODO: no hardcoding?
            data_frommap_bl_sel = data_frommap_merged_sel['mesh2d_flowelem_bl'].to_numpy()
            if 'mesh2d_layer_z' in varkeys_list or 'LayCoord_cc' in varkeys_list:
                print('layertype: zlayer')
                warnings.warn('WARNING: your model seems to contain only z-layers. if the modeloutput is generated with an older version of dflowfm, the coordinates can be incorrect. if your model contains z-sigma-layers, use the fulloutput option in the mdu and rerun (happens automatically in newer dflowfm versions).')
                zvals_interface_vec = data_frommap_merged_sel['mesh2d_interface_z'].to_numpy()[:,np.newaxis]
                zvals_interface = np.repeat(zvals_interface_vec,len(data_frommap_wl3_sel),axis=1)
                # zvalues lower than bedlevel should be overwritten with bedlevel
                for iL in range(nlay):
                    zvalbot_belowbl_bool = zvals_interface[iL,:]<data_frommap_bl_sel
                    zvals_interface[iL,zvalbot_belowbl_bool] = data_frommap_bl_sel[zvalbot_belowbl_bool]
                #top z-layer is extended to water level, if wl is higher than zval_lay_top
                zvals_interface[-1,:] = np.maximum(zvals_interface[-1,:],data_frommap_wl3_sel)
            elif 'mesh2d_layer_sigma' in varkeys_list:
                print('layertype: sigmalayer')
                zvals_interface_percentage = data_frommap_merged_sel['mesh2d_interface_sigma'].to_numpy()[:,np.newaxis]
                zvals_interface = data_frommap_wl3_sel+(data_frommap_wl3_sel-data_frommap_bl_sel)[np.newaxis]*zvals_interface_percentage
            else: # 2D model
                print('layertype: 2D model')
                if nlay!=1:
                    raise Exception('recheck this')
                #zvals_cen = np.linspace(data_frommap_bl_sel,data_frommap_wl3_sel,nlay)
                zvals_interface = np.linspace(data_frommap_bl_sel,data_frommap_wl3_sel,nlay+1)
        
        #convert to output for plot_netmapdata
        crs_dist_starts_matrix = np.repeat(intersect_pd['crs_dist_starts'].values[np.newaxis],nlay,axis=0)
        crs_dist_stops_matrix = np.repeat(intersect_pd['crs_dist_stops'].values[np.newaxis],nlay,axis=0)
        crs_verts_x_all = np.array([[crs_dist_starts_matrix.ravel(),crs_dist_stops_matrix.ravel(),crs_dist_stops_matrix.ravel(),crs_dist_starts_matrix.ravel()]]).T
        crs_verts_z_all = np.ma.array([zvals_interface[1:,:].ravel(),zvals_interface[1:,:].ravel(),zvals_interface[:-1,:].ravel(),zvals_interface[:-1,:].ravel()]).T[:,:,np.newaxis]
        crs_verts = np.ma.concatenate([crs_verts_x_all, crs_verts_z_all], axis=2)
        
        if varname is not None: #retrieve data for varname and return
            data_frommap_selvar = data_frommap_merged_sel[varname]
            if 'nmesh2d_layer' in data_frommap_selvar.dims:
                crs_plotdata = data_frommap_selvar.T.to_numpy().flatten() #TODO: flatten necessary?
            else: #for 2D models, no layers
                crs_plotdata = data_frommap_selvar
            return crs_verts, crs_plotdata
        else:
            return crs_verts
    
    
    
    
    print('calculating and plotting cross section')
    runtime_tstart = dt.datetime.now() #start timer
    #intersect function, find crossed cell numbers (gridnos) and coordinates of intersection (2 per crossed cell)
    intersect_pd = polygon_intersect(ugrid_all, data_frommap_merged, line_array, optimize_dist=False, calcdist_fromlatlon=calcdist_fromlatlon) #TODO: move to xarray
    #derive vertices from cross section (distance from first point)
    if reproduce:
        crs_verts, crs_plotdata = dfmt.get_xzcoords_onintersection(file_nc=file_nc, varname='mesh2d_sa1', intersect_pd=intersect_pd, timestep=timestep) #TODO: move to xarray
    else:
        crs_verts, crs_plotdata = get_xzcoords_onintersection(data_frommap_merged, varname='mesh2d_sa1', intersect_pd=intersect_pd, timestep=timestep) #TODO: move to xarray
    fig, ax = plt.subplots()
    pc = dfmt.plot_netmapdata(crs_verts, values=crs_plotdata, ax=ax, cmap='jet')#, linewidth=0.5, edgecolor='k')
    fig.colorbar(pc, ax=ax)
    ax.set_ylim(val_ylim)
    plt.savefig(os.path.join(dir_output,f'{basename}_crossect'))
    runtime_timedelta = (dt.datetime.now()-runtime_tstart)
    print(f'calculating and plotting cross section finished in {runtime_timedelta}')