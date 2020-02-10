# -*- coding: utf-8 -*-
"""
Created on Thu Feb  6 16:06:54 2020

@author: veenstra
"""

def get_varname_mapnc(data_nc,varname_requested):
    import pandas as pd
    
    #VARIABLE names used within different versions of Delft3D-Flexible Mesh
    varnames_list = pd.DataFrame()
    varnames_list['time'] = ['time','nmesh2d_dlwq_time',''] # time
    
    varnames_list['mesh2d_node_x'] = ['mesh2d_node_x','NetNode_x','mesh2d_agg_node_x'] # x-coordinate of nodes
    varnames_list['mesh2d_node_y'] = ['mesh2d_node_y','NetNode_y','mesh2d_agg_node_y'] # y-coordinate of nodes
    varnames_list['mesh2d_node_z'] = ['mesh2d_node_z','NetNode_z',''] # z-coordinate of nodes
    
    varnames_list['mesh2d_face_x'] = ['mesh2d_face_x','FlowElem_xzw','mesh2d_agg_face_x'] # x-coordinate of faces
    varnames_list['mesh2d_face_y'] = ['mesh2d_face_y','FlowElem_yzw','mesh2d_agg_face_y'] # y-coordinate of faces
    
    varnames_list['mesh2d_edge_x'] = ['mesh2d_edge_x','',''] # x-coordinate of velocity-points
    varnames_list['mesh2d_edge_y'] = ['mesh2d_edge_y','',''] # y-coordinate of velocity-points
    
    varnames_list['mesh2d_edge_nodes'] = ['mesh2d_edge_nodes','NetLink',''] # 'link between two netnodes' / 'Mapping from every edge to the two nodes that it connects'
    varnames_list['mesh2d_face_nodes'] = ['mesh2d_face_nodes','NetElemNode','mesh2d_agg_face_nodes'] # 
    
    varnames_list['mesh2d_face_x_bnd'] = ['mesh2d_face_x_bnd','FlowElemContour_x','mesh2d_agg_face_x_bnd'] # x-coordinates of flow element contours
    varnames_list['mesh2d_face_y_bnd'] = ['mesh2d_face_y_bnd','FlowElemContour_y','mesh2d_agg_face_y_bnd'] # y-coordinates of flow element contours
    
    varnames_list['mesh2d_flowelem_domain'] = ['mesh2d_flowelem_domain','FlowElemDomain',''] # flow element domain
    varnames_list['mesh2d_flowelem_bl'] = ['mesh2d_flowelem_bl','FlowElem_bl',''] # bed level
    varnames_list['mesh2d_flowelem_ba'] = ['mesh2d_flowelem_ba','FlowElem_bac',''] # area (m2) of cell faces
    
    varnames_list['mesh2d_ucx'] = ['mesh2d_ucx','ucx',''] # 
    varnames_list['mesh2d_ucy'] = ['mesh2d_ucy','ucy',''] # 
    varnames_list['mesh2d_layer_z'] = ['mesh2d_layer_z','LayCoord_cc',''] # 
    varnames_list['mesh2d_sa1'] = ['mesh2d_sa1','sa1',''] # 
    varnames_list['mesh2d_tem1'] = ['mesh2d_tem1','tem1',''] # 
    
    
    ### DIMENSION names used within different versions of Delft3D-Flexible Mesh
    #dimnames_list = pd.DataFrame()
    varnames_list['nmesh2d_node'] = ['nmesh2d_node','mesh2d_nNodes','']#,'nNetNode','NetElemNode'] # number of nodes
    varnames_list['nmesh2d_face'] = ['nmesh2d_face','mesh2d_nFaces','']#,'nNetElem','nFlowElem'] # number of faces
    varnames_list['nmesh2d_edge'] = ['nmesh2d_edge','','']#,'nNetLink'] # number of velocity-points
    
    varnames_list['nmesh2d_layer'] = ['nmesh2d_layer','mesh2d_nLayers','']#,'laydim'] # layer
    
    #look for correct pd column
    pdcol_bool = varnames_list.eq(varname_requested).any()
    varname_pdcol = pdcol_bool.index[pdcol_bool].tolist()
    if len(varname_pdcol) == 0:
        raise Exception('varname %s not found in internal database'%(varname_requested))
    elif len(varname_pdcol)>1:
        raise Exception('varname %s not found but multiple equivalents found in internal database: %s'%(varname_requested,varname_pdcol))
    else:
        varname_pdcol = varname_pdcol[0]
    
    data_nc_varnames_list = list(data_nc.variables.keys())
    data_nc_dimnames_list = list(data_nc.dimensions.keys())
    
    def get_vardimname(data_nc_names_list):
        #check what is in netcdf file
        if varname_requested in data_nc_names_list:
            varname = varname_requested
        elif varname_pdcol in data_nc_names_list:
            varname = varname_pdcol
        else:
            var_options = list(varnames_list[varname_pdcol])
            varname = [var for var in var_options if var in data_nc_names_list]
            if varname == []:
                varname = None
            else:
                varname = varname[0]
        return varname
    
    varname = get_vardimname(data_nc_varnames_list)
    if varname is None:
        varname = get_vardimname(data_nc_dimnames_list)
    if varname is None:
        print('WARNING: var/dim name %s or equivalent not found in netCDF file with variables:\n%s \nand dimensions:\n%s'%(varname_requested, data_nc_varnames_list, data_nc_dimnames_list))
    
    return varname











