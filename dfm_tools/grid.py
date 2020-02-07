import numpy as np
class UGrid:
    """Unstructured grid"""
    def __init__(self, mesh2d_node_x, mesh2d_node_y, mesh2d_face_nodes, verts, mesh2d_node_z=None, *args, **kwargs):
        self.mesh2d_node_x = mesh2d_node_x
        self.mesh2d_node_y = mesh2d_node_y
        self.mesh2d_face_nodes = mesh2d_face_nodes
        self.verts = verts
        if mesh2d_node_z is not None:
            self.mesh2d_node_z = mesh2d_node_z
        else:
            self.mesh2d_node_z = np.zeros(self.mesh2d_node_x.shape)
    @staticmethod
    def fromfile(filename_nc):
        import netCDF4
        from dfm_tools.get_varname_mapnc import get_varname_mapnc
        
        def nodexyfaces2verts(node_x,node_y, faces):
            quatrangles = faces-1 #convert 1-based indexing to 0-based indexing
            #https://stackoverflow.com/questions/49640311/matplotlib-unstructered-quadrilaterals-instead-of-triangles
            #https://stackoverflow.com/questions/52202014/how-can-i-plot-2d-fem-results-using-matplotlib
            yz = np.c_[node_x,node_y]
            verts= yz[quatrangles]
            verts[quatrangles.mask==True,:] = np.nan #remove all masked values by making them nan
            return verts
        
        data_nc = netCDF4.Dataset(filename_nc)

        mesh2d_node_x = data_nc.variables[get_varname_mapnc(data_nc,'mesh2d_node_x')][:]
        mesh2d_node_y = data_nc.variables[get_varname_mapnc(data_nc,'mesh2d_node_y')][:]
        #mesh2d_node_z = data_nc.variables[get_varname_mapnc(data_nc,'mesh2d_node_z')][:]
        #TODO: add node_z again, necessary for slicing?
        mesh2d_face_nodes = data_nc.variables[get_varname_mapnc(data_nc,'mesh2d_face_nodes')][:, :]
        verts = nodexyfaces2verts(mesh2d_node_x, mesh2d_node_y, mesh2d_face_nodes) #xy coordinates of face nodes

        #remove ghost cells
        varn = get_varname_mapnc(data_nc,'mesh2d_flowelem_domain')
        if varn != []: # domain variable is present, so there are multiple domains
            domain = data_nc.variables[varn][:]
            domain_no = np.bincount(domain).argmax() #meest voorkomende domeinnummer
            verts = verts[domain==domain_no]
            mesh2d_face_nodes = mesh2d_face_nodes[domain==domain_no]
        
        
        data_nc.close()
        ugrid = UGrid(mesh2d_node_x, mesh2d_node_y, mesh2d_face_nodes, verts, mesh2d_node_z=None)
        return ugrid


def get_mapfilelist(file_nc, multipart=None):
    #get list of mapfiles
    import re
    import glob
    lastpart = file_nc.split('_')[-2]
    if multipart != False and len(lastpart) == 4 and lastpart.isdigit(): #if part before '_map.nc' is eg '0000'
        filename_start = re.compile('(.*)_([0-9]+)_map.nc').search(file_nc).group(1)
        #filename_number = re.compile('(.*)_([0-9]+)_map.nc').search(file_nc).group(2)
        #file_ncs = [file_nc.replace('_%s_map.nc','_%04d_map.nc'%(filename_number, domain_id)) for domain_id in range(ndomains)]
        file_ncs = glob.glob('%s*_map.nc'%(filename_start))
    else:
        file_ncs = [file_nc]
    return file_ncs



def get_mapmodeldata(file_nc, var_values=None, multipart=None, timestep=None, lay=None):
    from netCDF4 import Dataset
    import numpy as np
    
    from dfm_tools.get_varname_mapnc import get_varname_mapnc
    from dfm_tools.grid import get_mapfilelist
    
    if timestep == None:
        raise Exception('ERROR: paramter timestep not privided')
    if lay == None:
        raise Exception('ERROR: paramter lay not privided')
    
    #convert timestep to list if it is not already
    if type(timestep)==int:
        list_timestep = [timestep]
    elif type(timestep)==list:
        list_timestep = timestep
    elif type(timestep)==range:
        list_timestep = list(timestep)
    #TODO: add check for list of ints
    #TODO: add option for all timesteps
    #TODO: add option for timesteps by timestamp instead of integer
    #TODO: add check for minmax timesteps (not lower than 0 and not higher than ntimesteps)
    #TODO: 0/1-based indexing for timesteps?
    
    #convert layer to list if it is not already
    if type(lay)==int:
        list_lay = [lay]
    elif type(lay)==list:
        list_lay = lay
    elif type(lay)==range:
        list_lay = list(lay)
    #TODO: add check for list of ints
    #TODO: add option for all layers
    #TODO: add check for minmax layers (not lower than 0 and not higher than nlayers)
    #TODO: 0/1-based indexing for layers?
    
    #TODO: add retrieval via depth (then dflowutil.mesh can be removed?)
    
    
    file_ncs = get_mapfilelist(file_nc, multipart)
    for iF, file_nc_sel in enumerate(file_ncs):
        print('processing mapdata from domain %04d of %04d'%(iF, len(file_ncs)-1))
        data_nc = Dataset(file_nc_sel)
        #list(data_nc.variables.keys())
        
        #values = get_mapmodeldata(data_nc, var_values=var_values, multipart=multipart, timestep=timestep, lay=lay)
        nc_values = data_nc.variables[var_values]
        #nc_values_shape = nc_values.shape
        nc_values_dims = nc_values.dimensions
        nc_values_ndims = len(nc_values_dims)
    

        #select values
        if nc_values_ndims == 2:
            values = nc_values[list_timestep,:]
        elif nc_values_ndims == 3:
            values = nc_values[list_timestep,:,list_lay]
        else:
            raise Exception('incorrect number of dimensions: %s'%(nc_values_ndims))
        
        if iF == 0:
            #setup initial array
            if nc_values_ndims == 2:
                values_all = np.ma.empty((len(list_timestep),0))
                concat_axis = 1
            elif nc_values_ndims == 3:
                values_all = np.ma.empty((len(list_timestep),0,len(list_lay)))
                concat_axis = 1
            else:
                raise Exception('incorrect number of dimensions: %s'%(nc_values_ndims))

        #filter ghost cells
        varn = get_varname_mapnc(data_nc,'mesh2d_flowelem_domain')
        if varn != []: # domain variable is present, so there are multiple domains
            domain = data_nc.variables['mesh2d_flowelem_domain'][:]
            domain_no = np.bincount(domain).argmax() #meest voorkomende domeinnummer
            nonghost_ids = domain==domain_no
        
            if nc_values_ndims == 2:
                values_all = np.ma.concatenate([values_all,values[:,nonghost_ids]],axis=concat_axis)
            if nc_values_ndims == 3:
                values_all = np.ma.concatenate([values_all,values[:,nonghost_ids,:]],axis=concat_axis)
        else: #1 domain
            #TODO: 1 domain with multiple layers would now go wrong, fix it
            values_all = np.ma.concatenate([values_all,values],axis=concat_axis)
    
    #TODO: add requested times and layers to outputdata        
    return values_all






def get_mapnetdata(file_nc, multipart=None):
    import numpy as np
    
    from dfm_tools.grid import get_mapfilelist, UGrid

    file_ncs = get_mapfilelist(file_nc, multipart)
    #get all data
    num_nodes = [0]
    for iF, file_nc_sel in enumerate(file_ncs):
        print('processing netdata from domain %04d of %04d'%(iF, len(file_ncs)-1))
        #data_nc = Dataset(file_nc_sel)
        #list(data_nc.variables.keys())
        
        ugrid = UGrid.fromfile(file_nc_sel)
        node_x = ugrid.mesh2d_node_x
        node_y = ugrid.mesh2d_node_y
        faces = ugrid.mesh2d_face_nodes
        verts = ugrid.verts

        #setup initial array
        if iF == 0:
            node_x_all = np.ma.empty((0,))
            node_y_all = np.ma.empty((0,))
            verts_all = np.ma.empty((0,verts.shape[1],verts.shape[2]))
            faces_all = np.ma.empty((0,faces.shape[1]),dtype='int32')
        
        #TODO: simplify the two blocks below
        #increase size of verts_all if too small for verts
        if verts_all.shape[1] < verts.shape[1]:
            tofew_cols = verts_all.shape[1] - verts.shape[1]
            verts_all_cordimsize = np.ma.zeros((verts_all.shape[0],verts.shape[1],verts_all.shape[2]))
            verts_all_cordimsize[:,:tofew_cols,:] = verts_all
            verts_all_cordimsize.mask = True
            verts_all_cordimsize.mask[:,:tofew_cols,:] = verts_all.mask
            faces_all_cordimsize = np.ma.zeros((faces_all.shape[0],faces.shape[1]),dtype='int32')
            faces_all_cordimsize[:,:tofew_cols] = faces_all
            faces_all_cordimsize.mask = True
            faces_all_cordimsize.mask[:,:tofew_cols] = faces_all.mask
        else:
            verts_all_cordimsize = verts_all
            faces_all_cordimsize = faces_all
            
        #increase size of verts if too small for verts_all
        if verts.shape[1] < verts_all.shape[1]:
            tofew_cols = verts.shape[1] - verts_all.shape[1]
            verts_cordimsize = np.ma.zeros((verts.shape[0],verts_all.shape[1],verts.shape[2]))
            verts_cordimsize[:,:tofew_cols,:] = verts
            verts_cordimsize.mask = True
            verts_cordimsize.mask[:,:tofew_cols,:] = verts.mask
            faces_cordimsize = np.ma.zeros((faces.shape[0],faces_all.shape[1]),dtype='int32')
            faces_cordimsize[:,:tofew_cols] = faces
            faces_cordimsize.mask = True
            faces_cordimsize.mask[:,:tofew_cols] = faces.mask
        else:
            verts_cordimsize = verts
            faces_cordimsize = faces

        #merge all
        node_x_all = np.ma.concatenate([node_x_all,node_x])
        node_y_all = np.ma.concatenate([node_y_all,node_y])
        verts_all = np.ma.concatenate([verts_all_cordimsize,verts_cordimsize])
        faces_all = np.ma.concatenate([faces_all_cordimsize,faces_cordimsize+np.sum(num_nodes)])
        num_nodes.append(node_x.shape[0])

        
    #set all invalid values to the same value (tends to differ between partitions)
    faces_all.data[faces_all.mask] = -999
    faces_all.fill_value = -999
    
    ugrid_all = UGrid(node_x_all, node_y_all, faces_all, verts_all)
    
    return ugrid_all


def plot_netmapdata(verts, values=None, ax=None, **kwargs):
    #https://stackoverflow.com/questions/52202014/how-can-i-plot-2d-fem-results-using-matplotlib
    #https://stackoverflow.com/questions/49640311/matplotlib-unstructered-quadrilaterals-instead-of-triangles
    import matplotlib.pyplot as plt
    import matplotlib.collections

    if not ax: ax=plt.gca()
    pc = matplotlib.collections.PolyCollection(verts, **kwargs)
    pc.set_array(values)
    ax.add_collection(pc)
    ax.autoscale()
    return pc





