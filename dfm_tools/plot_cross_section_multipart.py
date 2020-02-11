import os
import sys
import math
import glob
import numpy as np
import datetime as dt
from netCDF4 import Dataset, num2date
from geopy import distance
from scipy.interpolate import interp2d, griddata
from shapely.geometry import LineString, Point
from shapely.ops import nearest_points
import matplotlib.pyplot as plt
import pickle as pkl
import pandas as pd
from pyproj import Proj, transform

# Convert lat/lon corrdinates to x/y mercator projection
def merc(lat, lon):
    r_major = 6378137.000
    x = r_major * math.radians(lon)
    scale = x/lon
    y = 180.0/math.pi * math.log(math.tan(math.pi/4.0 + lat * (math.pi/180.0)/2.0)) * scale
    #x = lon
    #y = lat
    return (x, y)
merc = np.vectorize(merc)

def intersect(seg1_x, seg1_y, seg2_x, seg2_y):

    A1 = (seg1_y[0]-seg1_y[1]) / (seg1_x[0]-seg1_x[1])
    b1 = seg1_y[0]-A1*seg1_x[0]
    A2 = (seg2_y[0]-seg2_y[1]) / (seg2_x[0]-seg2_x[1])
    b2 = seg2_y[0]-A2*seg2_x[0]
    
    if (A1 == A2):
        return False # parallel segments
    else:
        Xint = (b2-b1) / (A1-A2)
        Yint = A1*Xint + b1
                
        if (Xint<=max(seg1_x)) and (Xint<=max(seg2_x)) and\
           (Xint>=min(seg1_x)) and (Xint>=min(seg2_x)) and\
           (Yint<=max(seg1_y)) and (Yint<=max(seg2_y)) and\
           (Yint>=min(seg1_y)) and (Yint>=min(seg2_y)):
            print(seg1_x, seg1_y, seg2_x, seg2_y, Xint, Yint)
            return True
        else:
            return False
###
class Pli:
    def __init__(self, lats, lons, x=None, y=None):
        self.lats = lats
        self.lons = lons
        self.x, self.y = merc(self.lats, self.lons)

class GridData:
    def __init___(self, mapfile, pli):
        self.mapfile = mapfile
        self.pli = pli
        self.xcor = []
        self.zcor = []
        self.scor = []
        self.xcen = []
        self.ycen = []
        self.scen = []
        self.zint = []
        self.zcen = []
        self.wl = []
        self.bed = []
        self.val = []
        self.times = []
        
###
#mapfile = r'P:\11203379-mwra-new-bem-model\waq_model\simulations\A31_1year_20191219\DFM_OUTPUT_MB_02_waq\MB_02_waq_0000_map.nc'
mapfiles = glob.glob(r'p:\11203379-mwra-new-bem-model\waq_model\simulations\A31_1year_20191219\DFM_OUTPUT_MB_02_waq\MB_02_waq_*_map.nc')
varname = 'mesh2d_OXY'
t0 = dt.datetime(2012,9,20)
mergepartitions = 1
#plifile = r'C:\Users\vilmin\OneDrive - Stichting Deltares\Documents\MassBay\BEM_plots\test.pli'
pliline = Pli([42.38, 42.38], [-70.98, -70.25])
print(pliline.x, pliline.y)

ifaces = []
xfaces = []
yfaces = []
sfaces = []
bfaces = []
zfaces = []
vfaces = []

ntot = 0

# Convert pli line to x, y
for imap, mapfile in enumerate(mapfiles):
#for imap, mapfile in enumerate([mapfiles[0]]):

    istart = len(sfaces)

    map_part = Dataset(mapfile)
    basename = os.path.basename(mapfile)
    npart = int(os.path.splitext(basename)[0][-8:-4])
    print('Reading mapfile for partition '+str(npart)+'...')
    maptimes = map_part['time'][:]
    mapdates = num2date(maptimes, map_part['time'].units)
    domain = map_part.variables['mesh2d_flowelem_domain'][:]
    lons_nodes = map_part.variables['mesh2d_node_x'][:]
    lats_nodes = map_part.variables['mesh2d_node_y'][:]
    lons_faces = map_part.variables['mesh2d_face_x'][:]
    lats_faces = map_part.variables['mesh2d_face_y'][:]
    global_ifaces = map_part.variables['mesh2d_flowelem_globalnr'][:]
    print(max(global_ifaces))
    nodes_edges = map_part.variables['mesh2d_edge_nodes'][:]
    face_edges = map_part.variables['mesh2d_edge_faces'][:]
    print(np.amax(face_edges))
    depth = map_part.variables['mesh2d_waterdepth'][:,:]
    sigma = map_part.variables['mesh2d_layer_sigma'][:]
    val = map_part.variables[varname][:,:,:]
    nedges = len(nodes_edges[:,0])
    nnodes = len(lons_nodes)
    nfaces = len(lons_faces)
    print(nfaces)
    lats_edges = []
    lons_edges = []
    for ie in range(nedges):
        #print(nodes_edges[ie,:])
        if (max(nodes_edges[ie,:]) < nnodes):
            lats_edges.append([lats_nodes[nodes_edges[ie,0]], lats_nodes[nodes_edges[ie,1]]])
            lons_edges.append([lons_nodes[nodes_edges[ie,0]], lons_nodes[nodes_edges[ie,1]]])
    xedges, yedges = merc(lats_edges, lons_edges)
    #print(xedges, yedges)
    nedges = len(xedges)

    # Check if part of the pliline is contained in this partition
    if (np.amax(xedges) < min(pliline.x)) or (np.amin(xedges) > max(pliline.x)) or\
       (np.amax(yedges) < min(pliline.y)) or (np.amin(yedges) > max(pliline.y)):
       print('The map '+mapfile+' file does not contain any point of the transect.')
       continue

    # Check which edges cross the pliline
    for ie in range(nedges):
        if (intersect(xedges[ie], yedges[ie], pliline.x, pliline.y)):
            #print(xedges[ie], yedges[ie], pliline.x, pliline.y)
            ntot += 1
            iface1 = face_edges[ie,0]
            iface2 = face_edges[ie,1]
            #print(face_edges[ie,:])
            #print(lons_faces[iface1], lons_faces[iface2])
            #print(lats_faces[iface1], lats_faces[iface2])
            #print(iface1, nfaces, domain[iface1], npart, ifaces)
            #print(iface2, nfaces, domain[iface2], npart, ifaces)
            
            # Get rid of ghost cells+cells that where already adjacent to another edge
            if (iface1>0) and (iface1<nfaces):
                if (domain[iface1] == npart) and not(iface1 in ifaces[istart:]):
                    print('Yes!')
                    ifaces.append(iface1)
                    xfaces.append(lons_faces[iface1])
                    yfaces.append(lats_faces[iface1])
                    sfaces.append(distance.distance((yfaces[-1],xfaces[-1]), (pliline.lats[0],pliline.lons[0])).km)
                
            if (iface2>0) and (iface2<nfaces):
                if (domain[iface2] == npart) and not(iface2 in ifaces[istart:]):
                    print('Yes!')
                    ifaces.append(iface2)
                    xfaces.append(lons_faces[iface2])
                    yfaces.append(lats_faces[iface2])
                    sfaces.append(distance.distance((yfaces[-1],xfaces[-1]), (pliline.lats[0],pliline.lons[0])).km)
    print('Number of faces: '+str(len(sfaces)))
    
    # Get depths and WQ variable values
    #for idate in range(len(mapdates)):
    print(istart)
    for idate in range(10,12): #just trying out on a couple of time steps
        if (imap == 0):
            bfaces.append([])
            zfaces.append([])
            vfaces.append([])
        for item in range(istart, len(sfaces)):
            iface = ifaces[item]
            bfaces[idate-10].append(-1.*depth[idate,iface])
            zfaces[idate-10].append([])
            vfaces[idate-10].append([])
            for ilay in range(len(sigma)):
                zfaces[idate-10][-1].append(depth[idate,iface]*sigma[ilay])
                vfaces[idate-10][-1].append(val[idate,iface,ilay])

plt_x = []
plt_y = []
plt_z = []
for item in range(len(sfaces)):
    for ilay in range(len(sigma)):
        plt_x.append(sfaces[item])
        plt_y.append(zfaces[1][item][ilay]) #2nd time step only for now
        plt_z.append(vfaces[1][item][ilay])
        
sorted_sfaces = sorted(sfaces)
sorted_bfaces = [x for _,x in sorted(zip(sfaces, bfaces[1]))]

print(plt_z)
fig = plt.figure(figsize=(12,6))
ax = fig.add_subplot(1,1,1)       
ax.tricontour(plt_x, plt_y, plt_z, levels=np.arange(6,12,0.5), linewidths=0, colors='k')
cntr2 = ax.tricontourf(plt_x, plt_y, plt_z, levels=np.arange(6,12,0.5), cmap="RdBu_r")
print(sorted_bfaces)
ax.fill_between(sorted_sfaces, sorted_bfaces, -100., facecolor='#e6ccb3', linewidth=0.)
# set the limits of the plot to the limits of the data
ax.axis([0., 60., -100., 0.])
fig.colorbar(cntr2, ax=ax)
plt.savefig('test_cross_section_allparts.png')
plt.close(fig)

print(ntot)

fp = open('truc2.csv', 'w')
for item in range(len(ifaces)):
    fp.write(str(lons_faces[item])+';'+str(lats_faces[item])+';'+str(sfaces[item])+';'+str(bfaces[1][item])+'\n')
fp.close()