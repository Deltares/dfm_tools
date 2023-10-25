# -*- coding: utf-8 -*-
"""
Created on Tue Mar 14 17:38:58 2023

@author: veenstra
"""

#TODO: A SIMPLE WORKFLOW WITH MATPLOTLIB SEEMS TO WORK MORE EFFICIENT THAN SHAPELY?

import matplotlib as mpl
import numpy as np

delete_pol = np.array([[ 1.91741935, 49.76580645],
                        [ 0.20387097, 49.9       ],
                        [-0.25032258, 48.71290323],
                        [ 1.92774194, 48.59935484],
                        ])

# xPO = delete_pol[:,0]
# yPO = delete_pol[:,1]
# xyPO = [(x,y) for x,y in zip(xPO,yPO)]
# poly = mpl.path.Path(xyPO)
poly = mpl.path.Path(delete_pol)

lonvals,latvals = [1,2],[49,49]
gridPoints = np.vstack((lonvals,latvals)).T
gridMask   = poly.contains_points(gridPoints)
print(gridMask)



#BELOW IS A WORKFLOW WITH SHAPELY

import os
import numpy as np
import xarray as xr
import pandas as pd
import geopandas as gpd
import matplotlib.pyplot as plt
plt.close('all')
import dfm_tools as dfmt
from shapely.geometry import Point

# modelrun settings
crs = 4326
plot_shp = False

## All the shapefiles in the Netherlands
shpdir = os.path.join('p:\\archivedprojects\\11208800-kdvoedsel-noordzee\\3.data_analysis\\inputData','windFarmAreas')
shpfile = ['Operational and tendered wind farms.shp', 'Search areas offshore wind.shp']  # 2x folders
idcol =  'wndfrm_'

area_shp = pd.DataFrame()
areaids = []
for i, shp_id in enumerate(shpfile):
    print(i, shp_id)
    # import area shapefile
    area_shp_single = gpd.read_file(os.path.join(shpdir, shp_id))
    area_shp_single = area_shp_single[area_shp_single.area>5e8] #TODO: temporarily using only the large shapefiles
    
    if shp_id == 'Operational and tendered wind farms.shp':
        area_shp_single = area_shp_single[area_shp_single.country == 'Netherlands']
    
    if shp_id == 'Search areas offshore wind.shp':
        # rename column to fit other shape file
        area_shp_single = area_shp_single.rename(columns={'OBJECTID':'wndfrm_', 'Name':'windfrm'})
        # add extra label because there were already numbers in the other shapefile ids
        area_shp_single ['wndfrm_'] = 'Id' + area_shp_single ['wndfrm_'].astype(str)
        # remove .0 because they were somehow added in the previous step 
        area_shp_single ['wndfrm_'] = area_shp_single ['wndfrm_'].str[:-2]
    
    # change crs to match mod_gdf
    area_shp_single  = area_shp_single.to_crs(crs)
    # repair polygons by adding a zero buffer
    area_shp_single.geometry = area_shp_single.geometry.buffer(0)
    # create list with area ID
    areaids_single = list(area_shp_single[idcol].unique())
    # add area shapes to total data frame
    area_shp = pd.concat([area_shp, area_shp_single])
    # add areaids to total list
    areaids = areaids + areaids_single

## Plotting the shapefiles
fig, ax = plt.subplots()
ax.set_xlim(0.0, 10)
ax.set_ylim(50, 56)
dfmt.plot_coastlines(ax=ax,res='i')
ax.grid()
area_shp.plot(ax=ax, linewidth = 1, facecolor = (1, 1, 1, 0), edgecolor = (0.5, 0.5, 0.5, 1))


## Loop over years
years = ['2017'] #['2009', '2010', '2011', '2012', '2013', '2014', '2015', '2016', '2017']
basedir = r"p:\archivedprojects\11208800-kdvoedsel-noordzee\3.data_analysis\inputData\extractedModelData" 
for year in years: # last few years - from 2015
    print(year)
    
    file_nc = os.path.join(basedir, year, f'DCSM-FM_0_5nm_waq_statistics_{year}_01*_map.nc') # read filepath #TODO: revert back from 01* to 0* again
    part_nos = np.array([100, 101, 102, 103, 121, 128, 129, 130, 131, 132, 133, 134, 135,
                         136, 137, 138, 139, 140, 141, 142, 143, 144, 145, 146, 147, 148,
                         149, 150, 151, 152, 153, 154, 155, 156, 157, 158, 159, 160, 161,
                         162, 163, 164, 165, 166, 167, 168, 169, 170, 171, 172, 173, 174, 175])
    file_nc = [os.path.join(basedir, year, f'DCSM-FM_0_5nm_waq_statistics_{year}_{partno:04d}_map.nc') for partno in part_nos]
    waq_xr = dfmt.open_partitioned_dataset(file_nc) # read data and combine parts
    
    # 1. Grid conversion
    faces = [Point(face_xy) for face_xy in waq_xr.grid.face_coordinates] # create points to compare to polygon
    
    # 2. Crop and average the output per polygon
    variables = ['SPM_gro']#['DIN_win', 'DIP_win', 'ucmag_ann', 'Chfla_gro', 'SPM_gro']
    for var in variables:
        print(var)
        average = []
        for i, row in area_shp.iterrows():  # iterate over all polygons
            poly = row.geometry
            filter_boolean = np.array([poly.contains(face) for face in faces])
            filter_boolean_xr = xr.DataArray(filter_boolean,dims=('mesh2d_nFaces'))
            # Select model points within the polygon:
            crop = waq_xr[var].isel(mesh2d_nLayers=-1,missing_dims='ignore').where(filter_boolean_xr, drop=True) # surface data (e.g. DIN_win) or Depth data (e.g. ucmag_ann)
            #plot cropped data and polygon
            if plot_shp:
                fig, ax = plt.subplots()
                crop.ugrid.plot()
                area_shp.loc[[i]].plot(ax=ax,facecolor='none')
            # Average over polygon
            wl_average = crop.mean(dim='mesh2d_nFaces') 
            print('i = ', wl_average.values) # To see progress 
            average.append(np.around(wl_average.values, 4))  # Rounding to 4 decimals

        ## Add column of average for that year and varibale for each polygon to the df:
        area_shp[f'{year}_{var}'] = average
        

### Calc average over years, per polygon and var and add as columns to the df:
area_shp['average_SPM_gro'] = area_shp.filter(regex='SPM_gro$',axis=1).mean(numeric_only=True, axis=1) # SPM_gro

## Plot the values onto each polygon, per variable:
fig, ax = plt.subplots()
ax.set_xlim(0.0, 10)
ax.set_ylim(50, 56)
dfmt.plot_coastlines(ax=ax,res='i')
ax.grid()
## SPM:
area_shp.plot(column='average_SPM_gro', ax=ax, legend=True, vmin=0, vmax=8)

