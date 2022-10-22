# -*- coding: utf-8 -*-
"""
Created on Wed Oct 27 22:05:55 2021

@author: veenstra
"""

import os
import numpy as np
from netCDF4 import Dataset
import matplotlib.pyplot as plt
plt.close('all')

from dfm_tools.get_nc import get_ncmodeldata#, get_netdata
from dfm_tools.get_nc_helpers import get_varname_fromnc

dir_testinput = r'c:\DATA\dfm_tools_testdata'
dir_output = '.'

#code from test_get_nc test d
file_nc = os.path.join(dir_testinput,r'DFM_sigma_curved_bend\DFM_OUTPUT_cb_3d\cb_3d_map.nc')

timestep = 72
#layno = 5
#calcdist_fromlatlon = None
multipart = None
#line_array = np.array([[ 104.15421399, 2042.7077107 ],
#                       [2913.47878063, 2102.48057382]])
#val_ylim = None
#clim_bl = None
#optimize_dist = None
#ugrid = get_netdata(file_nc=file_nc, multipart=multipart)
#intersect_gridnos, intersect_coords = ugrid.polygon_intersect(line_array, optimize_dist=None)

#code from get_xzcoords_onintersection
data_nc = Dataset(file_nc)

varn_mesh2d_s1 = get_varname_fromnc(data_nc,'mesh2d_s1', vardim='var')
data_frommap_wl3 = get_ncmodeldata(file_nc, varname=varn_mesh2d_s1, timestep=timestep, multipart=multipart)
data_frommap_wl3 = data_frommap_wl3[0,:]
#data_frommap_wl3_sel = data_frommap_wl3[0,intersect_gridnos]
varn_mesh2d_flowelem_bl = get_varname_fromnc(data_nc,'mesh2d_flowelem_bl', vardim='var')
data_frommap_bl = get_ncmodeldata(file_nc, varname=varn_mesh2d_flowelem_bl, multipart=multipart)
#data_frommap_bl_sel = data_frommap_bl[intersect_gridnos]

dimn_layer = get_varname_fromnc(data_nc,'nmesh2d_layer', vardim='dim')
if dimn_layer is None: #no layers, 2D model
    nlay = 1
else:
    nlay = data_nc.dimensions[dimn_layer].size

if 'mesh2d_layer_z' in data_nc.variables.keys(): #TODO: check get_nc.get_xzcoords_onintersection() for extra checks for layertypes (incl fullgrid output)
    laytyp = 'zlayer'
    #zvals_cen = get_ncmodeldata(file_nc=file_map, varname='mesh2d_layer_z', lay='all')#, multipart=False)
    #zvals_interface = get_ncmodeldata(file_nc=file_map, varname='mesh2d_interface_z')#, multipart=False)
    zvals_interface = data_nc.variables['mesh2d_interface_z'][:]
else:
    laytyp = 'sigmalayer'
    #zvals_cen = np.linspace(data_frommap_bl_sel,data_frommap_wl3_sel,nlay)
    #zvals_interface = np.linspace(data_frommap_bl_sel,data_frommap_wl3_sel,nlay+1)
    zvals_interface = np.linspace(data_frommap_bl,data_frommap_wl3,nlay+1)


print(laytyp)
depth = -1
z_test_higher = np.argmax((zvals_interface > depth),axis=0)
z_test_lower = np.argmin((zvals_interface < depth),axis=0)
z_test_all = z_test_higher==z_test_lower

