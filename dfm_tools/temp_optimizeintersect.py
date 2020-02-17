# -*- coding: utf-8 -*-
"""
Created on Sun Feb 16 14:37:31 2020

@author: veenstra
"""

import os
import matplotlib.pyplot as plt
plt.close('all')
import numpy as np
import datetime as dt

from dfm_tools.get_nc import get_netdata, get_ncmodeldata, get_xzcoords_onintersection, plot_netmapdata
from dfm_tools.polygon import LineBuilder


#grevelingen variables
dir_testinput = os.path.join(r'c:/DATA/werkmap','dfm_tools_testdata')
file_nc = os.path.join(dir_testinput, r'DFM_3D_z_Grevelingen\computations\run01\DFM_OUTPUT_Grevelingen-FM\Grevelingen-FM_0000_map.nc')
timestep = 3
layno = 35
calcdist_fromlatlon = None
multipart = None
line_array = np.array([[ 56267.59146475, 415644.67447155],
                       [ 64053.73427496, 419407.58239502]])
line_array = np.array([[ 53181.96942503, 424270.83361629],
                       [ 55160.15232593, 416913.77136685]])
line_array = np.array([[ 52787.21854294, 424392.10414528],
                       [ 55017.72655174, 416403.77313703],
                       [ 65288.43784807, 419360.49305567]])
line_array = np.array([[ 49843.46669412, 423787.66241104],
                       [ 48287.29831589, 415877.13982169],
                       [ 54317.45078154, 425019.62904381],
                       [ 64173.18384367, 420545.64495639]])
val_ylim = [-25,5]
clim_bl = None


ugrid = get_netdata(file_nc=file_nc, multipart=multipart)
#get bed layer
data_frommap_bl = get_ncmodeldata(file_nc=file_nc, varname='mesh2d_flowelem_bl', multipart=multipart)

#create plot with ugrid and cross section line
fig, ax_input = plt.subplots()
pc = plot_netmapdata(ugrid.verts, values=data_frommap_bl, ax=ax_input, linewidth=0.5, edgecolors='face', cmap='jet')#, color='crimson', facecolor="None")
pc.set_clim(clim_bl)
fig.colorbar(pc, ax=ax_input)
ax_input.set_aspect('equal')
if 0: #click interactive polygon
    #pol_frominput = Polygon.frominteractive(ax)
    line, = ax_input.plot([], [],'o-')  # empty line
    linebuilder = LineBuilder(line)
    line_array = linebuilder.line_array
ax_input.plot(line_array[:,0],line_array[:,1],'b',linewidth=3)


runtime_tstart = dt.datetime.now() #start timer
#intersect function, find crossed cell numbers (gridnos) and coordinates of intersection (2 per crossed cell)
intersect_gridnos, intersect_coords = ugrid.polygon_intersect(line_array)
#derive vertices from cross section (distance from first point)
#crs_verts = get_xzcoords_onintersection(file_nc=file_nc, line_array=line_array, intersect_gridnos=intersect_gridnos, intersect_coords=intersect_coords, timestep=timestep, calcdist_fromlatlon=calcdist_fromlatlon, multipart=multipart)
"""
#get data to plot
data_frommap = get_ncmodeldata(file_nc=file_nc, varname='mesh2d_sa1', timestep=timestep, layer='all', multipart=multipart)

#plot crossed cells (gridnos) in first plot
print(layno)#data_frommap_flat = data_frommap[0,intersect_gridnos,layno]
#pc = plot_netmapdata(ugrid.verts[intersect_gridnos,:,:], values=data_frommap_flat, ax=ax_input, linewidth=0.5, cmap="jet")
#plt.savefig(os.path.join(dir_output,'%s_gridbed'%(os.path.basename(file_nc).replace('.nc',''))))

#plot cross section
if len(data_frommap.shape) == 3:
    data_frommap_sel = data_frommap[0,intersect_gridnos,:]
    data_frommap_sel_flat = data_frommap_sel.T.flatten()
elif len(data_frommap.shape) == 2: #for 2D models, no layers 
    data_frommap_sel = data_frommap[0,intersect_gridnos]
    data_frommap_sel_flat = data_frommap_sel
fig, ax = plt.subplots()
pc = plot_netmapdata(crs_verts, values=data_frommap_sel_flat, ax=ax, linewidth=0.5, cmap='jet')
fig.colorbar(pc, ax=ax)
ax.set_ylim(val_ylim)
#plt.savefig(os.path.join(dir_output,'%s_crossect'%(os.path.basename(file_nc).replace('.nc',''))))

runtime_tstop = dt.datetime.now()
runtime_timedelta = (runtime_tstop-runtime_tstart).total_seconds()
print('caculating and plotting cross section finished in %.1f seconds'%(runtime_timedelta))

"""


def pairs(lst):
    """Iterate over a list in overlapping pairs.

    Args:
        lst: an iterable/list

    Returns:
        Yields a pair of consecutive elements (lst[k], lst[k+1]) of lst. Last
        call yields (lst[-2], lst[-1]).

    Example:
        lst = [4, 7, 11, 2]
        pairs(lst) yields (4, 7), (7, 11), (11, 2)

    Source:
        http://stackoverflow.com/questions/1257413/1257446#1257446
    """
    i = iter(lst)
    prev = next(i)
    for item in i:
        yield prev, item
        prev = item

def angle(v1, v2):
    """return angle between vector v1 and vector v2"""
    return np.arctan2(np.abs(np.cross(v1, v2)), np.dot(v1, v2))/np.pi*180


def segment_angles(line_array_ls):
    """calculate angles between all segments of a line

    :param line:    a shapely LineStrings
    :return:        a list angles and points
    #https://svn.oss.deltares.nl/repos/openearthtools/trunk/python/applications/hydrotools/hydrotools/gis/shapely_tools.py
    """
    # make pairs of segment coordinates
    segmentpairs = list(pairs([s for s in pairs(line_array_ls.coords)]))
    # find coordinates of vertices between two segments
    points = [np.array(cs[0][1]) for cs in segmentpairs]
    # calculate vectors by transformation of middle point to (0,0) coordinates
    vectors = [(np.array(cs[0][0])-p0, np.array(cs[1][1])-p0)
               for cs, p0 in zip(segmentpairs, points)]
    # calculate angle between vectors
    angles = [angle(v[0], v[1]) for v in vectors]
    angles = np.array(angles)

    
    return angles, points


from shapely.geometry import LineString, Point

#crs_xstart = intersect_coords[:,0,0]
#crs_xstop = intersect_coords[:,0,1]
#crs_ystart = intersect_coords[:,1,0]
#crs_ystop = intersect_coords[:,1,1]

#nlinedims = len(line_array.shape)
#ncrosscells = intersect_coords.shape[0]


#line_array_ls = LineString(line_array)
#angles, points = segment_angles(line_array_ls)

#calculate angles wrt x axis
angles_wrtx = []
#line_part_list = []
nlinecoords = line_array.shape[0]
for iL in range(nlinecoords-1):
    #line_part = line_array[iL:iL+2,:]
    #line_part_list.append(line_part)
    dx = line_array[iL+1,0] - line_array[iL,0]
    dy = line_array[iL+1,1] - line_array[iL,1]
    angles_wrtx.append(np.rad2deg(np.arctan2(dy,dx)))
angles_wrtx = np.array(angles_wrtx)
angles = 180-np.sign(np.diff(angles_wrtx))*np.diff(angles_wrtx)
angles_toprev = np.diff(angles_wrtx)


dist = 2000
angtot = angles_wrtx[:-1] + np.sign(angles_wrtx[:-1])*0.5*angles
angtot_wrtx = angles_wrtx[:-1] + 0.5*(180+angles_toprev)
#angtot_wrtx_otherside = 180-(angles_wrtx[:-1] + 0.5*(180+angles_toprev))
#print(angtot)

if 0:
    dxynewpoints = dist * np.array([np.cos(np.deg2rad(angtot_wrtx)),np.sin(np.deg2rad(angtot_wrtx))]).T
    newpoints1 = line_array[1:-1]+dxynewpoints
    newpoints2 = line_array[1:-1]-dxynewpoints
else:
    angtot_wrtx_ext = np.insert(angtot_wrtx,[0,-1],[angtot_wrtx[0],angtot_wrtx[-1]])
    #angtot_wrtx_ext = np.array([0]+list(angtot_wrtx)+[0])
    dxynewpoints = dist * np.array([np.cos(np.deg2rad(angtot_wrtx_ext)),np.sin(np.deg2rad(angtot_wrtx_ext))]).T
    newpoints1 = line_array+dxynewpoints
    newpoints2 = line_array-dxynewpoints


plt.plot(newpoints1[:,0],newpoints1[:,1], 'o-')
plt.plot(newpoints2[:,0],newpoints2[:,1], 'o-')
ax_input.grid()




"""
#perpendicular line
def perpendicular_line(l1, length):
    #Create a new Line perpendicular to this linear entity which passes
    through the point `p`.


    
    dx = l1.coords[1][0] - l1.coords[0][0]
    dy = l1.coords[1][1] - l1.coords[0][1]

    p = Point(l1.coords[0][0] + 0.5*dx, l1.coords[0][1] + 0.5*dy)
    x, y = p.coords[0][0],  p.coords[0][1]

    if (dy == 0) or (dx == 0):
        a = length / l1.length
        l2 = LineString([(x - 0.5*a*dy, y - 0.5*a*dx),
                         (x + 0.5*a*dy, y + 0.5*a*dx)])

    else:
        s = -dx/dy
        a = ((length * 0.5)**2 / (1 + s**2))**0.5
        l2 = LineString([(x + a, y + s*a),
                         (x - a, y - s*a)])

    return l2


polyline_half1 = []
polyline_half2 = []


xcoords = [x[0] for x in line2_list[0].coords]
ycoords = [x[1] for x in line2_list[0].coords]

plt.plot(xcoords, ycoords, 'o-r')

line_section_list = []
for iL in range(nlinecoords-1):
    #calculate length of lineparts
    line_section_list.append(LineString(line_array[iL:iL+2,:]))

    
iP = 0
testpoint = Point(crs_xstart[iP],crs_ystart[iP])
#linepart_lengthcum = np.cumsum(linepart_length)
"""





