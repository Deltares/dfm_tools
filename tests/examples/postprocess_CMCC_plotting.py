# -*- coding: utf-8 -*-
"""
Created on Mon Oct 24 22:20:28 2022

@author: veenstra
"""

import matplotlib.pyplot as plt
plt.close('all')
import dfm_tools as dfmt


for varn in ['uo', 'vo', 'thetao']:
    file_nc = f'p:\\archivedprojects\\11206304-futuremares-rawdata-preps\\data\\CMIP6_BC\\CMCC-ESM2\\{varn}_Omon_CMCC-ESM2_ssp126_r1i1p1f1_gn_201501-203412.nc'
    
    uds = dfmt.open_dataset_curvilinear(
        file_nc,
        x_dim='i',
        y_dim='j',
        x_bounds='vertices_longitude',
        y_bounds='vertices_latitude',
        convert_360to180=True,
        )
    uds = uds.ugrid.to_nonperiodic(xmax=180)
    # slice to europe (arbitrary, but to visualy compare grids)
    uds = uds.ugrid.sel(x=slice(-15,30), y=slice(30,70))
    
    fig, ax = plt.subplots()
    uds[varn].isel(time=0,lev=0).ugrid.plot(center=False)
    dfmt.plot_coastlines(ax=ax,res='l')
    ax.set_xlim(-10,30)
    ax.set_ylim(30,70)
    
    #this can be used to interpolate points of a polyline to a bcfile
    ds_onpoints = uds[varn].ugrid.sel_points(x=[-12,-12,-12,-8], y=[40,45,50,50])
