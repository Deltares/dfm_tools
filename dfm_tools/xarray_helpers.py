# -*- coding: utf-8 -*-
"""
Created on Fri Oct 14 19:58:36 2022

@author: veenstra
"""

from netCDF4 import Dataset
import xarray as xr
import numpy as np
import xugrid as xu
import matplotlib.pyplot as plt
plt.close('all')
import datetime as dt
import glob
import pandas as pd


def preprocess_hisnc(ds):
    """
    Look for dim/coord combination and use this for Dataset.set_index(), to enable station/gs/crs/laterals label based indexing. If duplicate labels are found (like duplicate stations), these are dropped to avoid indexing issues.
    
    Parameters
    ----------
    ds : xarray.Dataset
        DESCRIPTION.
    drop_duplicate_stations : TYPE, optional
        DESCRIPTION. The default is True.

    Returns
    -------
    ds : TYPE
        DESCRIPTION.

    """
    
    #generate dim_coord_dict to set indexes, this will be something like {'stations':'station_name','cross_section':'cross_section_name'} after loop
    dim_coord_dict = {}
    for ds_coord in ds.coords.keys():
        ds_coord_dtype = ds[ds_coord].dtype
        ds_coord_dim = ds[ds_coord].dims[0] #these vars always have only one dim
        if ds_coord_dtype.str.startswith('|S'): #these are station/crs/laterals/gs names/ids
            dim_coord_dict[ds_coord_dim] = ds_coord
    
    #loop over dimensions and set corresponding coordinates/variables from dim_coord_dict as their index
    for dim in dim_coord_dict.keys():
        coord = dim_coord_dict[dim]
        coord_str = f'{coord}'#_str' #avoid losing the original variable by creating a new name
        ds[coord_str] = ds[coord].load().str.decode('utf-8',errors='ignore').str.strip() #.load() is essential to convert not only first letter of string.
        ds = ds.set_index({dim:coord_str})
        
        #drop duplicate indices (stations/crs/gs), this avoids "InvalidIndexError: Reindexing only valid with uniquely valued Index objects"
        duplicated_keepfirst = ds[dim].to_series().duplicated(keep='first')
        if duplicated_keepfirst.sum()>0:
            print(f'dropping {duplicated_keepfirst.sum()} duplicate "{coord}" labels to avoid InvalidIndexError')
            ds = ds[{dim:~duplicated_keepfirst}]
    return ds


def preprocess_hirlam(ds):
    """
    add xy variables as longitude/latitude to avoid duplicate var/dim names
    add xy as variables again with help of NetCDF4 
    #TODO: this part is hopefully temporary, necessary since variables cannot have the same name as dimensions in xarray
    # background and future solution: https://github.com/pydata/xarray/issues/6293
    """
    
    print('adding x/y variables again as lon/lat')
    file_nc_one = ds.encoding['source']
    with Dataset(file_nc_one) as data_nc:
        data_nc_x = data_nc['x']
        data_nc_y = data_nc['y']
        ds['longitude'] = xr.DataArray(data_nc_x,dims=data_nc_x.dimensions,attrs=data_nc_x.__dict__)
        ds['latitude'] = xr.DataArray(data_nc_y,dims=data_nc_y.dimensions,attrs=data_nc_y.__dict__)
    ds = ds.set_coords(['latitude','longitude'])
    for varkey in ds.data_vars:
        del ds[varkey].encoding['coordinates'] #remove {'coordinates':'y x'} from encoding (otherwise set twice)
    return ds


def Dataset_varswithdim(ds,dimname):
    if dimname not in ds.dims:
        raise Exception(f'dimension {dimname} not in dataset, available are: {list(ds.dims)}')
    
    varlist_keep = []
    for varname in ds.variables.keys():
        if dimname in ds[varname].dims:
            varlist_keep.append(varname)
    ds = ds[varlist_keep]
    
    return ds


def open_partitioned_dataset(file_nc, chunks={'time':1}): #chunks={'time':1} increases performance significantly
    """
    
    Opmerkingen HB:
        - Dit werkt nu alleen voor data op de faces (die je nu expliciet aangeeft in bovenstaande functie). Voor edge data zal het ook wel kunnen werken, is uiteraard meer werk, moet even nadenken hoe dat qua nummering samenhangt.
        - Een nogwat suffe limitatie in xugrid: je kunt nog niet allerhande namen opgeven aan het grid. Dus ik genereer nu een nieuw grid (die gaat nu automatisch uit van een dimensie naam van "{naam_mesh}_nFaces". Beter zou zijn om alle nemen bij de initialisatie van het grid op te geven, dan kun je alle kanten uit. Ga ik even issue van maken.
        - Voor data op de edges zou het ook werken, maar dan is nog een andere isel noodzakelijk, specifiek voor de edge data.
        - Dit werkt nu ook alleen als je enkel grid in je dataset hebt. Bij meerdere grids zouden we een keyword moeten toevoegen dat je aangeeft welke je gemerged wilt zien.
    #TODO: maybe optimize by parallellization?
    """
    
    dtstart_all = dt.datetime.now()
    if isinstance(file_nc,list):
        file_nc_list = file_nc
    else:
        file_nc_list = glob.glob(file_nc)
    if len(file_nc_list)==0:
        raise Exception('file(s) not found, empty file_nc_list')
    
    dtstart = dt.datetime.now()
    print(f'>> xu.open_dataset() with {len(file_nc_list)} partition(s): ',end='')
    partitions = [xu.open_dataset(file_nc_one,chunks=chunks) for file_nc_one in file_nc_list]
    print(f'{(dt.datetime.now()-dtstart).total_seconds():.2f} sec')
    
    #rename old dimension and some variable names
    gridname = 'mesh2d' #partitions[0].ugrid.grid.name #'mesh2d' #TODO: works if xugrid accepts arbitrary grid names
    rename_dict = {}
    varn_maxfnodes = f'max_n{gridname}_face_nodes' #TODO: replace mesh2d with grid.name
    maxfnodes_opts = [f'{gridname}_nMax_face_nodes','nNetElemMaxNode'] #options for old domain variable name
    for opt in maxfnodes_opts:
        if opt in partitions[0].dims:
            rename_dict.update({opt:varn_maxfnodes})
    layer_opts = ['mesh2d_nLayers','laydim'] # options for old layer dimension name #TODO: others from get_varname_fromnc: ['nmesh2d_layer_dlwq']
    for opt in layer_opts:
        if opt in partitions[0].dims:
            #print(f'hardcoded replacing {opt} with nmesh2d_layer. Auto would replace "{partitions[0].ugrid.grid.to_dataset().mesh2d.vertical_dimensions}"')
            rename_dict.update({opt:'nmesh2d_layer'})
    #TODO: below works if xugrid handles arbitrary grid names
    # gridspecs = partitions[0].ugrid.grid.to_dataset()[gridname]
    # if hasattr(gridspecs,'vertical_dimensions'):
    #     layer_dimn = gridspecs.vertical_dimensions.split(':')[0]
    #     rename_dict.update({layer_dimn:f'n{gridname}_layer'}) #old varnames for mesh2d_layer: 'mesh2d_nLayers','laydim','nmesh2d_layer_dlwq'
    # else:
    #     print(f'no layer dimension found in gridspecs ds.{gridname}')
    varn_domain = f'{gridname}_flowelem_domain' #TODO: replace mesh2d with grid.name
    domain_opts = ['idomain','FlowElemDomain'] #options for old domain variable name
    for opt in domain_opts:
        if opt in partitions[0].data_vars:
            rename_dict.update({opt:varn_domain})
    varn_globalnr = f'{gridname}_flowelem_globalnr'
    globalnr_opts = ['iglobal_s'] #options for old globalnr variable name
    for opt in globalnr_opts:
        if opt in partitions[0].data_vars:
            rename_dict.update({opt:varn_globalnr})
    partitions = [part.rename(rename_dict) for part in partitions]
    
    varlist_onepart = list(partitions[0].variables.keys())
    
    all_indices = []
    all_faces = []
    all_nodes_x = []
    all_nodes_y = []
    accumulator = 0
    domainno_all = []
    #dtstart = dt.datetime.now()
    #print('>> process partitions facenumbers/ghostcells: ',end='')
    for i, part in enumerate(partitions):
        # For ghost nodes, keep the values of the domain number that occurs most.
        grid = part.ugrid.grid
        if varn_domain not in varlist_onepart:
            if len(partitions)==1:
                return partitions[0] #escape for non-partitioned files (domainno not found and one file provided). skipp rest of function
            else:
                raise Exception('no domain variable found, while there are multiple partition files supplied, this is not expected')
        da_domainno = part[varn_domain]
        try: #derive domainno from filename #TODO: this fails for restarts since it is _0000_20200101_120000_rst.nc (still the case?)
            part_domainno = int(part.encoding['source'][-11:-7])
        except: #derive domainno via domainno variable
            print('getting domainno from filename failed, now trying with bincount (might be costly)')
            part_domainno = np.bincount(da_domainno).argmax() 
        finally:
            if part_domainno in domainno_all:
                raise Exception(f'something went wrong, domainno {part_domainno} already occured: {domainno_all}')
            domainno_all.append(part_domainno)
        idx = np.flatnonzero(da_domainno == part_domainno) #something like >=i is applicable to edges/nodes
        faces = grid.face_node_connectivity[idx]
        faces[faces != grid.fill_value] += accumulator
        accumulator += grid.n_node
        all_indices.append(idx)
        all_faces.append(faces)
        all_nodes_x.append(grid.node_x)
        all_nodes_y.append(grid.node_y)
    #print(f'{(dt.datetime.now()-dtstart).total_seconds():.2f} sec')
    node_x = np.concatenate(all_nodes_x)
    node_y = np.concatenate(all_nodes_y)
    node_xy = np.column_stack([node_x, node_y])
    merged_nodes, inverse = np.unique(node_xy, return_inverse=True, axis=0)
    n_face_total = sum(len(faces) for faces in all_faces)
    n_max_node = max(faces.shape[1] for faces in all_faces)
    merged_faces = np.full((n_face_total, n_max_node), -1, dtype=np.intp)
    start = 0
    for faces in all_faces:
        n_face, n_max_node = faces.shape
        end = start + n_face
        merged_faces[start:end, :n_max_node] = faces
        start = end
    isnode = merged_faces != -1
    faces_flat = merged_faces[isnode]
    renumbered = inverse[faces_flat]
    merged_faces[isnode] = renumbered
    merged_grid = xu.Ugrid2d(
        node_x=merged_nodes[:, 0],
        node_y=merged_nodes[:, 1],
        fill_value=-1,
        face_node_connectivity=merged_faces,
    )
    facedim = partitions[0].ugrid.grid.face_dimension
    nodedim = partitions[0].ugrid.grid.node_dimension
    edgedim = partitions[0].ugrid.grid.edge_dimension
    #print(facedim,nodedim,edgedim)
    
    #define list of variables per dimension

    ds_face_list = []
    ds_node_list = []
    ds_edge_list = []
    #ds_rest_list = []
    dtstart = dt.datetime.now()
    print('>> ds.isel()/xr.append(): ',end='')
    for idx, uds in zip(all_indices, partitions):
        face_variables = []
        node_variables = []
        edge_variables = []
        for varname in uds.variables.keys():
            if varn_maxfnodes in uds[varname].dims: # not possible to concatenate this dim (size varies per partition) #therefore, vars mesh2d_face_x_bnd and mesh2d_face_y_bnd cannot be included currently. Maybe drop topology_dimension?: partitions[0].ugrid.grid.to_dataset().mesh2d.topology_dimension
                continue
            if facedim in uds[varname].dims:
                face_variables.append(varname)
            if nodedim in uds[varname].dims:
                node_variables.append(varname)
            if edgedim in uds[varname].dims:
                edge_variables.append(varname)
        ds_face = uds.ugrid.obj[face_variables]
        ds_node = uds.ugrid.obj[node_variables]
        ds_edge = uds.ugrid.obj[edge_variables]
        ds_rest = uds.ugrid.obj.drop_dims([facedim,nodedim,edgedim])
        ds_face_list.append(ds_face.isel({facedim: idx}))
        ds_node_list.append(ds_node)#.isel({nodedim: idx})) #TODO: add ghostcell removal for nodes and edges?
        ds_edge_list.append(ds_edge)#.isel({edgedim: idx}))
        #ds_rest_list.append(ds_rest)
    print(f'{(dt.datetime.now()-dtstart).total_seconds():.2f} sec')
    dtstart = dt.datetime.now()
    
    """
    # get a subset of the data
    da.isel(dim=[0, 1, 2])  # returns dimensions, etc. #isel is factor 10 langzamer dan rechtstreeks indexeren op dask/numpy array: https://github.com/pydata/xarray/issues/2227
    # Do the work yourself
    da.data[[0, 1, 2], ...]  # returns dask
    # Concatenate the data
    xr.concat([da1, da2, da3]) #is ook trager dan dask.array.stack()
    # Do it yourself
    dask.array.stack([da1.data, da2.data, da3.data])
    """
    
    print('>> xr.concat(): ',end='')
    ds_face_concat = xr.concat(ds_face_list, dim=facedim) #TODO: replace this with dask.stack() but that requires more book keeping?
    ds_node_concat = xr.concat(ds_node_list, dim=nodedim) #TODO: evt compat="override" proberen
    ds_edge_concat = xr.concat(ds_edge_list, dim=edgedim)
    #ds_rest_concat = xr.concat(ds_rest_list, dim=None)
    print(f'{(dt.datetime.now()-dtstart).total_seconds():.2f} sec')
    
    ds_merged = xr.merge([ds_face_concat,ds_node_concat,ds_edge_concat,ds_rest])
        
    varlist_merged = list(ds_merged.variables.keys())
    varlist_dropped_bool = ~pd.Series(varlist_onepart).isin(varlist_merged)
    varlist_dropped = pd.Series(varlist_onepart).loc[varlist_dropped_bool]
    if varlist_dropped_bool.any():
        print(f'WARNING: some variables dropped with merging of partitions:\n{varlist_dropped}')
    
    ds_merged = ds_merged.rename({facedim: merged_grid.face_dimension,
                                  nodedim: merged_grid.node_dimension,
                                  edgedim: merged_grid.edge_dimension}) #TODO: xugrid does not support other dimnames, xugrid issue is created: https://github.com/Deltares/xugrid/issues/25
    ds_merged_xu = xu.UgridDataset(ds_merged, grids=[merged_grid])
    print(f'>> open_partitioned_dataset total: {(dt.datetime.now()-dtstart_all).total_seconds():.2f} sec')
    return ds_merged_xu
