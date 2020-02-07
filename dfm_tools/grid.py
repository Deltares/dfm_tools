import numpy as np
class UGrid:
    """Unstructured grid"""
    def __init__(self, netnodex, netnodey, netnodez=None, netelemnode=None, *args, **kwargs):
        self.netnodex = netnodex
        self.netnodey = netnodey
        if netnodez is not None:
            self.netnodez = netnodez
        else:
            self.netnodez = np.zeros(self.netnodex.shape)
        self.netelemnode = netelemnode
    @staticmethod
    def fromfile(filename):
        import netCDF4
        from dfm_tools.get_varname_mapnc import get_varname_mapnc
        ds = netCDF4.Dataset(filename)
        varname_netnodex = get_varname_mapnc(ds,'NetNode_x')
        varname_netnodey = get_varname_mapnc(ds,'NetNode_y')
        varname_netnodez = get_varname_mapnc(ds,'NetNode_z')
        varname_netelemnode = get_varname_mapnc(ds,'NetElemNode')
        #varname_netlink = get_varname_mapnc(ds,'NetLink')
        netnodex = ds.variables[varname_netnodex][:]
        netnodey = ds.variables[varname_netnodey][:]
        netnodez = ds.variables[varname_netnodez][:]
        netelemnode = ds.variables[varname_netelemnode][:]
        #netlink = ds.variables[varname_netlink][:]
        ds.close()
        ugrid = UGrid(netnodex, netnodey, netnodez, netelemnode)
        return ugrid
    @staticmethod
    def frombmi(bmi):
        raise NotImplemented('todo')
    def celltypes(self):
        return (~self.netelemnode.mask).sum(1)
    def cellcoords(self):
        """Create a list of coordinates (xy) per cell"""
        # Create split locations
        splitidx = np.cumsum(np.r_[self.celltypes()][:-1])
        # Convert to 1d filled idx
        idx = self.netelemnode[(~self.netelemnode.mask)]-1
        xpoly = np.array(np.split(self.netnodex[idx],splitidx)) # x vector per poly
        ypoly = np.array(np.split(self.netnodey[idx],splitidx))
        zpoly = np.array(np.split(self.netnodez[idx],splitidx))
        cellxycoords = np.concatenate([xpoly[:,np.newaxis],ypoly[:,np.newaxis],zpoly[:,np.newaxis]],1)
        polycoords = [np.c_[xyz[0],xyz[1], xyz[2]] for xyz in cellxycoords]
        return polycoords
    def export(self, drivername='Memory', filename='', epsg=None):
        """create a file using ogr"""
        import ogr
        import osr
        driver = ogr.GetDriverByName(drivername)
        datasource = driver.CreateDataSource(filename)
        wgs = osr.SpatialReference()
        if epsg is None:
            def project(*args):
                return args
        else:
            wgs.ImportFromEPSG(4326)
            rd = osr.SpatialReference()
            rd.ImportFromEPSG(epsg)
            transform = osr.CoordinateTransformation(rd, wgs)
            project = transform.TransformPoint

        # Add fields here, if required...
        fielddefn = ogr.FieldDefn()
        fielddefn.SetName('cellid')
        fielddefn.SetType(ogr.OFTInteger)
        layer = datasource.CreateLayer('cells', wgs, ogr.wkbPolygon25D)
        layer.CreateField(fielddefn)

        featuredefn = ogr.FeatureDefn()
        featuredefn.AddFieldDefn(fielddefn)
        for i, cell in enumerate(self.cellcoords()):
            feature = ogr.Feature(featuredefn)
            linearring = ogr.Geometry(ogr.wkbLinearRing)
            for point in cell:
                linearring.AddPoint(*project(point[0], point[1], point[2]))
            # close ring
            linearring.AddPoint(*project(point[0], point[1], point[2]))
            geometry = ogr.Geometry(ogr.wkbPolygon25D)
            geometry.AddGeometry(linearring)
            feature.SetGeometry(geometry)
            feature.SetField(0, i)
            layer.CreateFeature(feature)
            feature.Destroy()
        datasource.SyncToDisk()
        if drivername == 'Memory':
            featurestr = ",\n".join([layer.GetFeature(i).ExportToJson() for i in range(layer.GetFeatureCount())])
            return '{ "type": "FeatureCollection",\n"features": [\n%s\n] }' % (featurestr, )
        datasource= None



def plot_net(file_nc, var_values=None):
    import netCDF4
    import matplotlib.pyplot as plt
    #plt.close('all')
    import matplotlib.collections
    #import matplotlib.cm as cm
    import numpy as np

    data_nc = netCDF4.Dataset(file_nc)
    node_x = data_nc.variables['mesh2d_node_x'][:]
    node_y = data_nc.variables['mesh2d_node_y'][:]
    faces = data_nc.variables['mesh2d_face_nodes'][:, :]
    if var_values == None:
        values = None
    else:
        values = data_nc.variables[var_values][0, :]
    #list(data_nc.variables.keys())

    quatrangles = faces-1 #convert 1-based indexing to 0-based indexing
    def quatplot(node_x,node_y, quatrangles, values, ax=None, **kwargs):
        #https://stackoverflow.com/questions/49640311/matplotlib-unstructered-quadrilaterals-instead-of-triangles
        #https://stackoverflow.com/questions/52202014/how-can-i-plot-2d-fem-results-using-matplotlib
        if not ax: ax=plt.gca()
        yz = np.c_[node_x,node_y]
        verts= yz[quatrangles]
        verts[quatrangles.mask==True,:] = np.nan #remove all masked values by making them nan
        pc = matplotlib.collections.PolyCollection(verts, **kwargs)
        pc.set_array(values)
        ax.add_collection(pc)
        ax.autoscale()
        return pc
    
    
    fig, ax = plt.subplots()
    ax.set_aspect('equal')

    if var_values == None:
        pc = quatplot(node_x, node_y, quatrangles, values, ax=None, linewidth=0.5, color="crimson", facecolor="None")
    else:
        pc = quatplot(node_x, node_y, quatrangles, values, ax=None, linewidth=0.5)#, cmap="jet")
        fig.colorbar(pc, ax=ax)
    #plt.plot(node_x,node_y, marker="o", ls="", color="crimson")
    plt.title('This is the plot for: quad')
    plt.xlabel('Y Axis')
    plt.ylabel('Z Axis')

