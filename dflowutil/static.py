import netCDF4

def nc_format(grd):
    """
    returns grid variables depending on net type

    Arguments:
        grd {str} -- path to a *_net.nc file
    """

    ds = netCDF4.Dataset(grd)
    map1 = {'xnode': 'NetNode_x',
            'ynode': 'NetNode_y',
            'xvelocity': 'ucx',
            'yvelocity': 'ucy',
            'layerz': 'LayCoord_cc',
            'cellnodes': 'NetElemNode',
            'face_x' : 'FlowElem_xzw',
            'face_y' : 'FlowElem_yzw',
            'domain_number': 'FlowElemDomain',
            'salinity': 'sa1',
            'temperature': 'tem1'}

    map4 = {'xnode': 'mesh2d_node_x',
            'ynode': 'mesh2d_node_y',
            'xvelocity': 'mesh2d_ucx',
            'yvelocity': 'mesh2d_ucy',
            'layerz': 'mesh2d_layer_z',
            'cellnodes': 'mesh2d_face_nodes',
            'face_x': 'mesh2d_face_x',
            'face_y': 'mesh2d_face_y',
            'domain_number': 'mesh2d_flowelem_domain',
            'salinity': 'mesh2d_sa1',
            'temperature': 'mesh2d_tem1'}

    agg = {'xnode': 'mesh2d_agg_node_x',
           'ynode': 'mesh2d_agg_node_y',
           'cellnodes': 'mesh2d_agg_face_nodes',
           'face_x': 'mesh2d_agg_face_x',
           'face_y': 'mesh2d_agg_face_y'}

    try:
        ds.variables['NetNode_x']
        varnames = map1
        print('Map type = 1')
    except:
        try:
            ds.variables['mesh2d_node_x']
            varnames = map4
            print('Map type = 4')
        except:
            try:
                ds.variables['mesh2d_agg_node_x']
                varnames = agg
            except:
                print('ERROR: file is broken!')
                raise
            
    return varnames
