# -*- coding: utf-8 -*-
"""
Created on Wed Oct 27 22:04:02 2021

@author: veenstra

WARNING: THIS TEST IS NOT YET FINISHED, WILL BE IMPROVED AND LINKED TO INTERNAL FUNCTIONS ASAP
"""

import os
import numpy as np
import matplotlib.pyplot as plt
plt.close('all')
import geopandas as gpd #conda install --channel conda-forge geopandas (breaks dfm_tools environment because of Qt issue)
from shapely.geometry import Polygon
import dfm_tools as dfmt

dir_testinput = r'c:\DATA\dfm_tools_testdata'
dir_output = '.'


varlist = ['Chlfa']#,'mesh2d_s1']
dir_shp = dir_output
if not os.path.exists(dir_shp):
    os.makedirs(dir_shp)
file_nc = os.path.join(r'p:\archivedprojects\11203850-coastserv\06-Model\waq_model\simulations\run0_20200319\DFM_OUTPUT_kzn_waq', 'kzn_waq_0000_map.nc')

vars_pd = dfmt.get_ncvarproperties(file_nc=file_nc)
vars_pd_matching = vars_pd[vars_pd.loc[:,'long_name'].str.match('.*Chl.*')]
#vars_pd_matching = vars_pd[vars_pd.loc[:,'long_name'].str.startswith('') & vars_pd.loc[:,'long_name'].str.endswith('Chlo')]
varns_Chl = vars_pd_matching.index.tolist()
varns_Chl_long = vars_pd_matching['long_name'].tolist()

ugrid = dfmt.get_netdata(file_nc=file_nc)#, multipart=False)

pol_shp_list = []
#partly from dfm_tools.polygon_intersect()
for iP, pol_data in enumerate(ugrid.verts): #[range(5000),:,:]
    pol_data_nonan = pol_data[~np.isnan(pol_data).all(axis=1)]
    pol_shp = Polygon(pol_data_nonan)
    pol_shp_list.append(pol_shp)

print('creating geodataframe with cells')
newdata = gpd.GeoDataFrame({'geometry': pol_shp_list},crs="EPSG:4326") #way more time efficient than doing it the loop

for iV, varname in enumerate(varlist):
    newdata[varname] = None

for timestep in [6]:#[0,10,20,30]:
    for iV, varname in enumerate(varlist):
        try:
            data_fromnc_all = dfmt.get_ncmodeldata(file_nc=file_nc, varname=varname, timestep=timestep, layer='all')
            data_fromnc_bot = dfmt.get_ncmodeldata(file_nc=file_nc, varname=varname, timestep=timestep, layer='bottom')
            data_fromnc_top = dfmt.get_ncmodeldata(file_nc=file_nc, varname=varname, timestep=timestep, layer='top')
        except:
            data_fromnc_top = dfmt.get_ncmodeldata(file_nc=file_nc, varname=varname, timestep=timestep)

        data_fromnc_nonan = data_fromnc_top[:]
        data_fromnc_nonan[data_fromnc_nonan.mask] = np.nan
        newdata[varname] = data_fromnc_nonan.data.flatten()
    file_shp = os.path.join(dir_shp,'shp_%s_%s'%(varname,data_fromnc_top.var_times.iloc[0].strftime('%Y%m%d')))
    newdata.to_file(file_shp)
    """
    fig, ax = plt.subplots(figsize=(6,7))
    pc = plot_netmapdata(ugrid.verts, values=data_fromnc_top.data.flatten(), ax=None, linewidth=0.5, cmap='jet')
    #pc.set_clim([-1,0])
    fig.colorbar(pc)
    ax.set_aspect(1./np.cos(np.mean(ax.get_ylim())/180*np.pi),adjustable='box')
    fig.tight_layout()
    """
