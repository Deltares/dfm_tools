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
        mesh2d_node_z = data_nc.variables[get_varname_mapnc(data_nc,'mesh2d_node_z')][:]
        mesh2d_face_nodes = data_nc.variables[get_varname_mapnc(data_nc,'mesh2d_face_nodes')][:, :]
        verts = nodexyfaces2verts(mesh2d_node_x, mesh2d_node_y, mesh2d_face_nodes) #xy coordinates of face nodes
        #TODO: couple with get_varname_mapnc function to make more generic for older grids/mapfiles
        
        #remove ghost cells
        varn = get_varname_mapnc(data_nc,'mesh2d_flowelem_domain')
        if varn != []: # domain variable is present, so there are multiple domains
            domain = data_nc.variables[varn][:]
            domain_no = np.bincount(domain).argmax() #meest voorkomende domeinnummer
            verts = verts[domain==domain_no]
            mesh2d_face_nodes = mesh2d_face_nodes[domain==domain_no]
        
        
        data_nc.close()
        ugrid = UGrid(mesh2d_node_x, mesh2d_node_y, mesh2d_face_nodes, verts, mesh2d_node_z=mesh2d_node_z)
        return ugrid


def get_mapfilelist(file_nc, multipart=None):
    #get list of mapfiles
    import re
    import glob
    lastpart = file_nc.split('_')[-2]
    if file_nc.endswith('_map.nc') and multipart != False and len(lastpart) == 4 and lastpart.isdigit(): #if part before '_map.nc' is eg '0000'
        filename_start = re.compile('(.*)_([0-9]+)_map.nc').search(file_nc).group(1)
        #filename_number = re.compile('(.*)_([0-9]+)_map.nc').search(file_nc).group(2)
        #file_ncs = [file_nc.replace('_%s_map.nc','_%04d_map.nc'%(filename_number, domain_id)) for domain_id in range(ndomains)]
        file_ncs = glob.glob('%s*_map.nc'%(filename_start))
    else:
        file_ncs = [file_nc]
    return file_ncs

def get_timesfromnc(file_nc):
    from netCDF4 import Dataset,num2date#,date2num
    import numpy as np
    
    data_nc = Dataset(file_nc)
    data_nc_timevar = data_nc.variables['time']    
    
    time0 = data_nc_timevar[0]
    time1 = data_nc_timevar[1]
    timeend = data_nc_timevar[-1]
    timeinc = time1-time0
    
    data_nc_times = np.arange(time0,timeend+timeinc,timeinc)
    data_nc_datetimes = num2date(data_nc_times, units = data_nc_timevar.units)
    
    
    return data_nc_times, data_nc_datetimes


def get_hismapmodeldata(file_nc, var_values=None, multipart=None, timestep=None, lay=None):
    from netCDF4 import Dataset
    import numpy as np
    
    from dfm_tools.get_varname_mapnc import get_varname_mapnc
    from dfm_tools.grid import get_mapfilelist, get_timesfromnc
    
    if timestep is None:
        raise Exception('ERROR: paramter timestep not privided')
    if lay == None:
        raise Exception('ERROR: paramter lay not privided')
    
    #get times
    data_nc_times, data_nc_datetimes = get_timesfromnc(file_nc)
    
    #convert timestep to list if it is not already
    if type(timestep)==int:
        list_timestep = [timestep]
    elif type(timestep)==list:
        list_timestep = timestep
    elif type(timestep)==range:
        list_timestep = timestep
    elif timestep == 'all':
        list_timestep = range(len(data_nc_times))
    else:
        print('WARNING: timestep variable type not predifined (%s), take as it comes'%(type(timestep)))
        list_timestep = timestep
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
        list_lay = lay
    else:
        print('WARNING: lay variable type not predifined (%s), take as it comes'%(type(lay)))
        list_lay = lay
    #TODO: add check for list of ints
    #TODO: add option for all layers
    #TODO: add check for minmax layers (not lower than 0 and not higher than nlayers)
    #TODO: 0/1-based indexing for layers?
    
    #TODO: add retrieval via depth (then dflowutil.mesh can be removed?)
    
    
    file_ncs = get_mapfilelist(file_nc, multipart)
    
    for iF, file_nc_sel in enumerate(file_ncs):
        print('processing mapdata from domain %04d of %04d'%(iF, len(file_ncs)-1))
        data_nc = Dataset(file_nc_sel)
        
        # check if requested variable is in netcdf
        nc_varkeys = list(data_nc.variables.keys())
        if var_values not in nc_varkeys:
            raise Exception('ERROR: requested variable %s not in netcdf, available are: %s'%(var_values, nc_varkeys))
        #TODO: check ndims of requested variable, eg no timestep/lay arguments should be provided for bedlevel)
        
        nc_values = data_nc.variables[var_values]
        #nc_values_shape = nc_values.shape
        nc_values_dims = nc_values.dimensions
        nc_values_ndims = len(nc_values_dims)
    
        #filter ghost cells
        varn = get_varname_mapnc(data_nc,'mesh2d_flowelem_domain')
        if varn != []: # domain variable is present, so there are multiple domains
            domain = data_nc.variables['mesh2d_flowelem_domain'][:]
            domain_no = np.bincount(domain).argmax() #meest voorkomende domeinnummer
            nonghost_ids = domain==domain_no

            
        # 1 dimension (faces)
        if nc_values_ndims == 1:
            #print('WARNING: untested number of dimensions: %s'%(nc_values_ndims))
            if iF == 0: #setup initial array
                values_all = np.ma.empty((0))    
            values = nc_values[:]
            concat_axis = 0
            if varn != []: # domain variable is present, so there are multiple domains
                values_all = np.ma.concatenate([values_all,values[nonghost_ids]],axis=concat_axis)
            else:
                values_all = np.ma.concatenate([values_all,values],axis=concat_axis)

        # 2 dimensions (time, faces)
        elif nc_values_ndims == 2:
            if iF == 0: #setup initial array
                values_all = np.ma.empty((len(list_timestep),0))    
            #select values
            values = nc_values[list_timestep,:]
            concat_axis = 1
            if varn != []: # domain variable is present, so there are multiple domains
                values_all = np.ma.concatenate([values_all,values[:,nonghost_ids]],axis=concat_axis)
            else:
                values_all = np.ma.concatenate([values_all,values],axis=concat_axis)
        # 3 dimensions (time, faces, layers)?
        elif nc_values_ndims == 3:
            if iF == 0: #setup initial array
                values_all = np.ma.empty((len(list_timestep),0,len(list_lay)))
            #select values
            values = nc_values[list_timestep,:,list_lay]
            concat_axis = 1
            if varn != []: # domain variable is present, so there are multiple domains
                values_all = np.ma.concatenate([values_all,values[:,nonghost_ids,:]],axis=concat_axis)
            else:
                values_all = np.ma.concatenate([values_all,values],axis=concat_axis)

        else:
            raise Exception('unanticipated number of dimensions: %s'%(nc_values_ndims))
    
    #TODO: add requested times and layers to outputdata, also multiple variables with different dimensions? (class?)     
    return values_all






def get_netdata(file_nc, multipart=None):
    import numpy as np
    
    from dfm_tools.grid import get_mapfilelist, UGrid

    file_ncs = get_mapfilelist(file_nc, multipart)
    #get all data
    num_nodes = [0]
    verts_shape2_all = []
    for iF, file_nc_sel in enumerate(file_ncs):
        print('analyzing netdata from domain %04d of %04d'%(iF, len(file_ncs)-1))
        ugrid = UGrid.fromfile(file_nc_sel)
        verts_shape2_all.append(ugrid.verts.shape[1])
    verts_shape2_max = np.max(verts_shape2_all)
        
    for iF, file_nc_sel in enumerate(file_ncs):
        print('processing netdata from domain %04d of %04d'%(iF, len(file_ncs)-1))
        #data_nc = Dataset(file_nc_sel)
        #list(data_nc.variables.keys())
        
        ugrid = UGrid.fromfile(file_nc_sel)
        node_x = ugrid.mesh2d_node_x
        node_y = ugrid.mesh2d_node_y
        node_z = ugrid.mesh2d_node_z
        faces = ugrid.mesh2d_face_nodes
        verts = ugrid.verts

        #setup initial array
        if iF == 0:
            node_x_all = np.ma.empty((0,))
            node_y_all = np.ma.empty((0,))
            node_z_all = np.ma.empty((0,))
            verts_all = np.ma.empty((0,verts_shape2_max,verts.shape[2]))
            faces_all = np.ma.empty((0,verts_shape2_max),dtype='int32')
        
        #if necessary, add masked column(s) to increase size to max in domains
        if verts.shape[1] < verts_shape2_max:
            tofew_cols = -(verts.shape[1] - verts_shape2_max)
            vcol_extra = verts[:,[0],:]
            vcol_extra.mask = True
            fcol_extra = faces[:,[0]]
            fcol_extra.mask = True
            for iC in range(tofew_cols):
                verts = np.hstack([verts,vcol_extra])
                faces = np.hstack([faces,fcol_extra])
        """
        #increase size of verts if too small for verts_all
        if verts.shape[1] < verts_shape2_max:
            tofew_cols = verts.shape[1] - verts_shape2_max
            verts_cordimsize = np.ma.zeros((verts.shape[0],verts_shape2_max,verts.shape[2]))
            verts_cordimsize[:,:tofew_cols,:] = verts
            verts_cordimsize.mask = True
            verts_cordimsize.mask[:,:tofew_cols,:] = verts.mask
            faces_cordimsize = np.ma.zeros((faces.shape[0],verts_shape2_max),dtype='int32')
            faces_cordimsize[:,:tofew_cols] = faces
            faces_cordimsize.mask = True
            faces_cordimsize.mask[:,:tofew_cols] = faces.mask
        else:
            verts_cordimsize = verts
            faces_cordimsize = faces
        """
        #merge all
        node_x_all = np.ma.concatenate([node_x_all,node_x])
        node_y_all = np.ma.concatenate([node_y_all,node_y])
        node_z_all = np.ma.concatenate([node_z_all,node_z])
        verts_all = np.ma.concatenate([verts_all,verts])
        faces_all = np.ma.concatenate([faces_all,faces+np.sum(num_nodes)])
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
    
    #check if data size is equal
    if not values is None:
        if verts.shape[0] != values.shape[0]:
            raise Exception('ERROR: size of grid and values is not equal, cannot plot')
    
    if not ax: ax=plt.gca()
    pc = matplotlib.collections.PolyCollection(verts, **kwargs)
    pc.set_array(values)
    ax.add_collection(pc)
    ax.autoscale()
    return pc





