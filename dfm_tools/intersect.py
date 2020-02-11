# -*- coding: utf-8 -*-
"""
Created on Tue Feb 11 15:01:06 2020

@author: veenstra
"""
from netCDF4 import Dataset
from shapely.geometry import MultiPolygon, Polygon, LineString, MultiLineString, multipoint, shape
import matplotlib.pyplot as plt
plt.close('all')
import numpy as np

class LineBuilder:
    import numpy as np
    def __init__(self, line):
        self.line = line
        self.xs = list(line.get_xdata())
        self.ys = list(line.get_ydata())
        self.cid = line.figure.canvas.mpl_connect('button_press_event', self)

    def __call__(self, event):
        print('click', event)
        if event.inaxes!=self.line.axes: return
        self.xs.append(event.xdata)
        self.ys.append(event.ydata)
        self.line_array = np.c_[self.xs, self.ys]
        self.line.set_data(self.xs, self.ys)
        self.line.figure.canvas.draw()


def nodexyfaces2verts(node_x,node_y, faces):
    quatrangles = faces-1 #convert 1-based indexing of cell numbering in ugrid to 0-based indexing
    #https://stackoverflow.com/questions/49640311/matplotlib-unstructered-quadrilaterals-instead-of-triangles
    #https://stackoverflow.com/questions/52202014/how-can-i-plot-2d-fem-results-using-matplotlib
    yz = np.c_[node_x,node_y]
    verts= yz[quatrangles]
    verts[quatrangles.mask==True,:] = np.nan #remove all masked values by making them nan
    return verts




from dfm_tools.grid import get_netdata, get_hismapmodeldata, plot_netmapdata

file_map = r'c:\DATA\werkmap\vanJulien_shortmodelfiles\DFM_sigma_curved_bend\DFM_OUTPUT_cb_3d\cb_3d_map.nc'
file_map = r'c:\DATA\werkmap\vanJulien_shortmodelfiles\DFM_3D_z_Grevelingen\computations\run01\DFM_OUTPUT_Grevelingen-FM\Grevelingen-FM_0000_map.nc'

#CURVIBEND (datetime)
print('plot grid and values from mapdata (constantvalue, 1 dim)')
ugrid = get_netdata(file_nc=file_map)#,multipart=False)
fig, ax = plt.subplots()
pc = plot_netmapdata(ugrid.verts, values=None, ax=ax, linewidth=0.5, color='crimson', facecolor="None")
ax.set_aspect('equal')

if 0:
    line, = ax.plot([], [],'o-')  # empty line
    linebuilder = LineBuilder(line)
    line_array = linebuilder.line_array
else:
    if 'cb_3d_map' in file_map:
        line_array = np.array([[ 185.08667065, 2461.11775254],
                               [2934.63837418, 1134.16019127]])
        line_array = np.array([[ 104.15421399, 2042.7077107 ],
                               [2913.47878063, 2102.48057382]])
    elif 'Grevelingen' in file_map:
        line_array = np.array([[ 56267.59146475, 415644.67447155],
                               [ 64053.73427496, 419407.58239502]])
        line_array = np.array([[ 53181.96942503, 424270.83361629],
                               [ 55160.15232593, 416913.77136685]])
    
#allpol = []
all_gridnos = []
all_intersect = []
line_section = LineString(line_array)
for iP, pol_data in enumerate(ugrid.verts):
    pol_data_nonan = pol_data[~np.isnan(pol_data).all(axis=1)]
    pol_shp = Polygon(pol_data_nonan)
    #allpol.append(pol_shp)
    intersection_line = pol_shp.intersection(line_section).coords
    if intersection_line != []:
        all_gridnos.append(iP)
        all_intersect.append(list(intersection_line))
    #print(iP)
#allpol_multi = MultiPolygon(allpol)
#intersection_line = allpol_multi.intersection(line_section).coords #does not work


if 'cb_3d_map' in file_map:
    timestep = 72
    layno = 5
elif 'Grevelingen' in file_map:
    timestep = 3
    layno = 35
data_frommap = get_hismapmodeldata(file_nc=file_map, varname='mesh2d_sa1', timestep=timestep, lay='all')#, multipart=False)
data_frommap_flat = data_frommap[0,all_gridnos,layno]
#fig, ax = plt.subplots()
ax.plot(line_array[:,0],line_array[:,1])
pc = plot_netmapdata(ugrid.verts[all_gridnos,:,:], values=data_frommap_flat, ax=None, linewidth=0.5, cmap="jet")
#pc.set_clim([28,30.2])
fig.colorbar(pc, ax=ax)
ax.set_aspect('equal')



data_frommap_wl3 = get_hismapmodeldata(file_nc=file_map, varname='mesh2d_s1', timestep=timestep)#, multipart=False)
data_frommap_wl3_sel = data_frommap_wl3[0,all_gridnos]
data_frommap_bl = get_hismapmodeldata(file_nc=file_map, varname='mesh2d_flowelem_bl')#, multipart=False)
data_frommap_bl_sel = data_frommap_bl[all_gridnos]
data_frommap_sel = data_frommap[0,all_gridnos,:]

crs_xstart = [x[0][0] for x in all_intersect]
crs_xstop = [x[1][0] for x in all_intersect]
crs_ystart = [x[0][1] for x in all_intersect]
crs_ystop = [x[1][1] for x in all_intersect]

nlay = data_frommap.shape[2]

crs_dist_starts = np.sqrt((crs_xstart - line_array[0,0])**2 + (crs_ystart - line_array[0,1])**2)
crs_dist_stops = np.sqrt((crs_xstop - line_array[0,0])**2 + (crs_ystop - line_array[0,1])**2)

data_nc = Dataset(file_map)
if 'cb_3d_map' in file_map:
    zvals_cen = np.linspace(data_frommap_bl_sel,data_frommap_wl3_sel,nlay)
    zvals_interface = np.linspace(data_frommap_bl_sel,data_frommap_wl3_sel,nlay+1)
elif 'Grevelingen' in file_map:
    #zvals_cen = get_hismapmodeldata(file_nc=file_map, varname='mesh2d_layer_z', lay='all')#, multipart=False)
    #zvals_interface = get_hismapmodeldata(file_nc=file_map, varname='mesh2d_interface_z')#, multipart=False)
    zvals_interface = data_nc.variables['mesh2d_interface_z'][:]


fig, ax_crs = plt.subplots()
crs_verts_x_all = np.empty((0,4,1))
crs_verts_z_all = np.empty((0,4,1))
#data_frommap_sel_flat = np.empty((0))
for iL in range(nlay):
    zval_lay_bot = zvals_interface[iL]
    zval_lay_top = zvals_interface[iL+1]
    crs_verts_x = np.array([[crs_dist_starts,crs_dist_stops,crs_dist_stops,crs_dist_starts]]).T
    if 'cb_3d_map' in file_map:
        crs_verts_z = np.array([[zval_lay_bot,zval_lay_bot,zval_lay_top,zval_lay_top]]).T
    elif 'Grevelingen' in file_map:
        crs_verts_z = np.repeat(np.array([[zval_lay_bot,zval_lay_bot,zval_lay_top,zval_lay_top]]).T[np.newaxis],repeats=crs_verts_x.shape[0],axis=0)
        #top z-layer is extended to water level, if wl is higher than zval_lay_top
        if iL == nlay-1:
            crs_verts_z[:,2,0] = np.maximum(zval_lay_top,data_frommap_wl3_sel.data)
            crs_verts_z[:,3,0] = np.maximum(zval_lay_top,data_frommap_wl3_sel.data)
        # zval_lay_bot lower than bedlevel should be overwritten with bedlevel
        zvalbot_belowbl_bool = crs_verts_z[:,0,0]<data_frommap_bl_sel
        crs_verts_z[zvalbot_belowbl_bool,0,0] = data_frommap_bl_sel[zvalbot_belowbl_bool]
        crs_verts_z[zvalbot_belowbl_bool,1,0] = data_frommap_bl_sel[zvalbot_belowbl_bool]
    crs_verts_x_all = np.concatenate([crs_verts_x_all, crs_verts_x])
    crs_verts_z_all = np.concatenate([crs_verts_z_all, crs_verts_z])
    #data_frommap_sel_flat = np.concatenate([data_frommap_sel_flat,data_frommap_sel[:,iL]])

data_frommap_sel_flat = data_frommap_sel.T.flatten()

crs_verts = np.concatenate([crs_verts_x_all, crs_verts_z_all], axis=2)

pc = plot_netmapdata(crs_verts, values=data_frommap_sel_flat, ax=ax_crs, linewidth=0.5, cmap="jet")
fig.colorbar(pc, ax=ax_crs)
ax_crs.set_ylim([-25,5])

#data_frommap = get_hismapmodeldata(file_nc=file_map, varname='mesh2d_sa12', timestep=timestep, lay='all')#, multipart=False)








