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
        ds = netCDF4.Dataset(filename)
        netnodex = ds.variables['NetNode_x'][:]
        netnodey = ds.variables['NetNode_y'][:]
        netnodez = ds.variables['NetNode_z'][:]
        netelemnode = ds.variables['NetElemNode'][:]
        netlink = ds.variables['NetLink'][:]
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
