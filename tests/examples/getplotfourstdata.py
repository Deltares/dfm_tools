# -*- coding: utf-8 -*-
"""
Created on Wed Oct 27 14:49:14 2021

@author: veenstra
"""

import os
import numpy as np
import matplotlib.pyplot as plt
plt.close('all')

from dfm_tools.get_nc import get_netdata, get_ncmodeldata, plot_netmapdata
from dfm_tools.get_nc_helpers import get_ncvardimlist
from dfm_tools.regulargrid import scatter_to_regulargrid

dir_testinput = r'c:\DATA\dfm_tools_testdata'
dir_output = '.'

#RMM foufile met quivers
file_nc = os.path.join(dir_testinput,r'DFM_fou_RMM\RMM_dflowfm_0000_fou.nc')
vars_pd, dims_pd = get_ncvardimlist(file_nc=file_nc)
#stations_pd = get_hisstationlist(file_nc,varname='waterlevel')

ugrid_all = get_netdata(file_nc=file_nc)
ux_mean = get_ncmodeldata(file_nc=file_nc, varname='mesh2d_fourier001_mean')
uy_mean = get_ncmodeldata(file_nc=file_nc, varname='mesh2d_fourier002_mean')
magn_mean = np.sqrt(ux_mean**2+uy_mean**2)
#uc_mean = get_ncmodeldata(file_nc=file_nc, varname='mesh2d_fourier003_mean')
#uc_max = get_ncmodeldata(file_nc=file_nc, varname='mesh2d_fourier004_max')
facex = get_ncmodeldata(file_nc=file_nc, varname='mesh2d_face_x')
facey = get_ncmodeldata(file_nc=file_nc, varname='mesh2d_face_y')

X,Y,U = scatter_to_regulargrid(xcoords=facex, ycoords=facey, ncellx=60, ncelly=35, values=ux_mean)
X,Y,V = scatter_to_regulargrid(xcoords=facex, ycoords=facey, ncellx=60, ncelly=35, values=uy_mean)

#thinning = 3
fig1,ax1 = plt.subplots(figsize=(9,5))
pc1 = plot_netmapdata(ugrid_all.verts, magn_mean, edgecolor='face')
#ax1.quiver(facex[::thinning], facey[::thinning], ux_mean[::thinning], uy_mean[::thinning], color='w',scale=20)#,width=0.005)#, edgecolor='face', cmap='jet')
ax1.quiver(X,Y,U,V, color='w',scale=5)#,width=0.005)#, edgecolor='face', cmap='jet')
pc1.set_clim([0,0.10])
#ax1.set_title('sqrt(x^2+y^2)\nx=%s\ny=%s'%(ux_mean.var_ncvarobject.long_name,uy_mean.var_ncvarobject.long_name))
ax1.set_aspect('equal')
ax1.set_xlabel('RD x [m]')
ax1.set_ylabel('RD y [m]')
cbar = fig1.colorbar(pc1)
cbar.set_label('residuele stroming [m/s]')
fig1.tight_layout()
fig1.savefig(os.path.join(dir_output,os.path.basename(file_nc).replace('.','')))


#RMM rst file
file_nc = os.path.join(dir_testinput,r'DFM_fou_RMM\RMM_dflowfm_0006_20131127_000000_rst.nc')
vars_pd, dims_pd = get_ncvardimlist(file_nc=file_nc)
#ugrid = get_netdata(file_nc=file_nc) #does not work, so scatter has to be used
ugrid_FlowElem_xzw = get_ncmodeldata(file_nc=file_nc, varname='FlowElem_xzw', multipart=False)
ugrid_FlowElem_yzw = get_ncmodeldata(file_nc=file_nc, varname='FlowElem_yzw', multipart=False)
data_s1 = get_ncmodeldata(file_nc=file_nc, varname='s1',timestep=0, multipart=False)

fig, ax = plt.subplots(figsize=(10,4))
pc = plt.scatter(ugrid_FlowElem_xzw,ugrid_FlowElem_yzw,1, data_s1[0,:],cmap='jet')
pc.set_clim([0,2])
fig.colorbar(pc)
ax.set_aspect('equal')
fig.tight_layout()
plt.savefig(os.path.join(dir_output,os.path.basename(file_nc).replace('.','')))

#keep this test to check flow link handling by get_nc.py (does not have domains)
#ugrid_FlowLink_xu = get_ncmodeldata(file_nc=file_nc, varname='FlowLink_xu', multipart=False)
#ugrid_FlowLink_yu = get_ncmodeldata(file_nc=file_nc, varname='FlowLink_yu', multipart=False)
data_q1 = get_ncmodeldata(file_nc=file_nc, varname='q1',timestep=0, multipart=False)

"""
fig, ax = plt.subplots(figsize=(10,4))
pc = plt.scatter(ugrid_FlowLink_xu,ugrid_FlowLink_yu,1, data_q1[0,:],cmap='jet')
pc.set_clim(None)
fig.colorbar(pc)
ax.set_aspect('equal')
fig.tight_layout()
plt.savefig(os.path.join(dir_output,os.path.basename(file_nc).replace('.','')))
"""

