# -*- coding: utf-8 -*-
"""
dfm_tools are post-processing tools for Delft3D FM
Copyright (C) 2020 Deltares. All rights reserved.

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  if not, see <http://www.gnu.org/licenses/>.

All names, logos, and references to "Deltares" are registered trademarks of
Stichting Deltares and remain full property of Stichting Deltares at all times.
All rights reserved.


INFORMATION
This script is part of dfm_tools: https://github.com/openearth/dfm_tools
Check the README.rst on github for other available functions
Check the tests folder on github for example scripts (this is the dfm_tools pytest testbank)
Check the pptx and example figures in (created by the testbank): N:/Deltabox/Bulletin/veenstra/info dfm_tools

Created on Sun Mar 22 08:41:00 2020

@author: veenstra
"""




def meshgridxy2verts(x_coords, y_coords):
    import numpy as np
    
    if len(x_coords.shape) != len(y_coords.shape):
        raise Exception('number of dimensions of input arrays should be equal')
    
    #from here equal amount of dimensions, so only one array is tested for n dimensions
    if len(x_coords.shape) == 1:
        x_mesh, y_mesh = np.meshgrid(x_coords, y_coords)
    elif len(x_coords.shape) == 2:
        x_mesh = x_coords
        y_mesh = y_coords
    else:
        raise Exception('input arrays can only be 1D or 2D')
    
    regulargrid_verts = np.empty(shape=((x_coords.shape[0]-1),(y_coords.shape[1]-1),4,2))
    regulargrid_verts[:] = np.nan
    regulargrid_verts[:,:,0,0] = x_mesh[ :-1,1:  ]
    regulargrid_verts[:,:,1,0] = x_mesh[ :-1, :-1]
    regulargrid_verts[:,:,2,0] = x_mesh[1:  , :-1]
    regulargrid_verts[:,:,3,0] = x_mesh[1:  ,1:  ]
    regulargrid_verts[:,:,0,1] = y_mesh[ :-1,1:  ]
    regulargrid_verts[:,:,1,1] = y_mesh[ :-1, :-1]
    regulargrid_verts[:,:,2,1] = y_mesh[1:  , :-1]
    regulargrid_verts[:,:,3,1] = y_mesh[1:  ,1:  ]
    
    return regulargrid_verts



def scatter_to_regulargrid(xcoords=None, ycoords=None, ncellx=None, ncelly=None, values=None, method='nearest', maskland_dist=None):
    """
    interpolates scatter values (x,y,z) or meshgrids to regular grid

    Parameters
    ----------
    xcoords : TYPE, optional
        DESCRIPTION. The default is None.
    ycoords : TYPE, optional
        DESCRIPTION. The default is None.
    ncellx : TYPE, optional
        DESCRIPTION. The default is None.
    ncelly : TYPE, optional
        DESCRIPTION. The default is None.
    values : TYPE, optional
        DESCRIPTION. The default is None.
    method : TYPE, optional
        DESCRIPTION. The default is 'nearest'.
    maskland_dist : TYPE, optional
        Distance used to mask out land cells. Value differs per model resolution. For sperical use maskland_dist of eg 0.01, for cartesian use maskland_dist of eg 100. The default is None.

    Returns
    -------
    x_grid : TYPE
        DESCRIPTION.
    y_grid : TYPE
        DESCRIPTION.
    value_grid : TYPE
        DESCRIPTION.

    """
    import numpy as np
    from scipy.interpolate import griddata
    from scipy.spatial import KDTree
    
    reg_x_vec = np.linspace(np.min(xcoords),np.max(xcoords),ncellx)
    reg_y_vec = np.linspace(np.min(ycoords),np.max(ycoords),ncelly)
    x_grid,y_grid = np.meshgrid(reg_x_vec,reg_y_vec)
    
    #first replace masked value with nan (mask is not used in griddata)
    try:
        values[values.mask] = np.nan
    except:
        pass
    value_grid = griddata((xcoords,ycoords),values,(x_grid,y_grid),method=method)
    
    if maskland_dist is not None:
        #Mask out land values with maskland_dist
        #https://stackoverflow.com/questions/10456143/get-only-valid-points-in-2d-interpolation-of-cloud-point-using-scipy-numpy
        tree = KDTree(np.c_[xcoords, ycoords])
        dist, _ = tree.query(np.c_[x_grid.ravel(), y_grid.ravel()], k=1)
        dist = dist.reshape(x_grid.shape)
        value_grid[dist > maskland_dist] = np.nan

    return x_grid, y_grid, value_grid



def center2corner(cen):
    import numpy as np
    
    if len(cen.shape) != 2:
        raise Exception('input array should have 2 dimensions')
    
    cen_nobnd = corner2center(cen)
    cen_nobnd_diff_ax0 = np.diff(cen_nobnd, axis=0)
    add_top = cen_nobnd[0,:] - cen_nobnd_diff_ax0[0,:]
    add_bot = cen_nobnd[-1,:] + cen_nobnd_diff_ax0[-1,:]
    cen_vertbnd = np.vstack([add_top, cen_nobnd, add_bot])
    cen_nobnd_diff_ax1 = np.diff(cen_vertbnd, axis=1)
    add_left = cen_vertbnd[:,0] - cen_nobnd_diff_ax1[:,0]
    add_right = cen_vertbnd[:,-1] + cen_nobnd_diff_ax1[:,-1]
    cen_withbnd = np.vstack([add_left.T, cen_vertbnd.T, add_right.T]).T
    cor = cen_withbnd
        
    return cor





def corner2center(cor):
    """
    from OET but edited, original author is Gerben de Boer

    cen = corner2center(cor) calculates the value of the center
    of the pixels by avareging the surrounding 4 corner values
    for arrays, or the surrounding 2 corner for vectors or [1,n]
    or [n,1 2D arrays. corner2center works both for coordinate meshes
    as well as data values defined on those meshes.
    
    >>> corner2center([1,3,5])
    array([ 2.,  4.])
    
    >>> corner2center([[1,3,5]])
    array([[ 2.,  4.]])
    
    >>> corner2center([[1],[3],[5]])
    array([[ 2.],
           [ 4.]])
    
    >>> corner2center([[1,3,5],[2,6,10]])
    array([[ 3.,  6.]])
    """
    import numpy as np
    cor = np.asarray(cor)
    shp = cor.shape
    
    if len(shp)==1 and len(cor)<2:
        raise ValueError('at least 2 elements required')
    elif len(shp)==1:
        cen = np.zeros(tuple([shp[0]-1,1]))
        cen = (cor[ :-1] + cor[1:  ])/2
    elif len(shp)==2 and shp[0]==1:
        cen = np.zeros(tuple([1,shp[1]-1]))
        cen[:,:] = (cor[:, :-1] + cor[:,1:  ])/2
    elif len(shp)==2 and shp[1]==1:
        cen = np.zeros(tuple([shp[0]-1,1]))
        cen[:,:] = (cor[ :-1,:] + cor[1:  ,:])/2
    elif len(shp)==2: 
        cen = np.zeros(tuple(np.array(shp)-1))
        cen[:,:] = (cor[ :-1, :-1] + 
                    cor[ :-1,1:  ] + 
                    cor[1:  ,1:  ] + 
                    cor[1:  , :-1])/4
    elif len(shp)>3:
        raise NotImplementedError('only 1D and 2D arrays implemented, only intervals and pixels, no voxels')
    return cen
      





def uva2xymagdeg(U1, V1, ALFAS, KCU=None, KCV=None, inactivewhen4x0=True):
    """
    this function converts velocities in m,n-direction (defined mathematically, so 0 on x-axis and increasing counter-clockwise)
    alpha is a matrix with orientations of cells, with respect to the north (varname='ALFAS') in D3D output
    output:
        vec_x - velocity in x-direction (east)
        vec_y - velocity in y-direction (east)
        vec_x - velocity in x-direction (east)
        vec_x - velocity in x-direction (east)

    De snelheden U1 zijn gedefinieerd op de U-punten. Die moet je eerst middelen naar de celcentra. Laat Uc de celcentre snelheid zijn:
    Uc(m,n) = (U1(m,n) + U1(m-1,n))/2
    Vc(m,n) = (V1(m,n) + V1(m,n-1))/2
    Vervolgens moeten die componenten gedraaid worden daarvoor is de hoek ALFAS nodig.
    Ux (m,n) = Uc(m,n)*cos(ALFAS) - Vc(m,n)*sin(ALFAS)
    Uy (m,n) = Uc(m,n)*sin(ALFAS) + Vc(m,n)*cos(ALFAS)
    Alles zou goed moeten gaan als je bij de middeling ervoor zorgt dat Uc(m,n) weer op de positie (m,n) in de array komt.
    Uc = (U1(1:end-1,:) + U1(2:end-1,:))/2 gaat dus net mis.
    Uc(2:end-1,: ) = (U1(1:end-1,: )) + U1(2:end,: ))/2 klopt wel, maar werkt alleen als je eerst Uc op de juiste afmeting hebt geïnitialiseerd.
    """
    import numpy as np
    """
    vel_magn = np.sqrt(U1**2 + V1**2)
    direction_math_deg = np.rad2deg(np.arctan2(V1, U1))+ALFAS
    direction_naut_deg = (90-direction_math_deg)%360
    vel_x = vel_magn*np.cos(np.deg2rad(direction_math_deg))
    vel_y = vel_magn*np.sin(np.deg2rad(direction_math_deg))
    """
    #replace -999. by 0
    U1[(U1==-999.)] = np.nan
    V1[(V1==-999.)] = np.nan

    #calculate UV velocities on centers
    Uc = np.empty(shape=U1.shape)
    Uc[:] = np.nan
    Vc = np.empty(shape=U1.shape)
    Vc[:] = np.nan
    Uc[1:,:] = np.nansum([U1[:-1,:], U1[1:,:]], axis=0)/2
    Vc[:,1:] = np.nansum([V1[:,:-1], V1[:,1:]], axis=0)/2
    
    #remove cells when all four velocity points (U1 and V1) were -999.
    mask_u999 = np.zeros(shape=U1.shape)
    mask_v999 = np.zeros(shape=U1.shape)
    mask_u999[1:,:] = ~(np.isnan(U1[:-1,:]) & np.isnan(U1[1:,:]))
    mask_v999[:,1:] = ~(np.isnan(V1[:,:-1]) & np.isnan(V1[:,1:]))
    mask_UV999 = mask_u999 + mask_v999
    Uc[mask_UV999==0] = np.nan
    Vc[mask_UV999==0] = np.nan
    
    if inactivewhen4x0:
        #mask cells when all four velocity points (U1 and V1) are 0.
        mask_u999 = np.zeros(shape=U1.shape)
        mask_v999 = np.zeros(shape=U1.shape)
        mask_u999[1:,:] = ~((U1[:-1,:]==0) & (U1[1:,:]==0))
        mask_v999[:,1:] = ~((V1[:,:-1]==0) & (V1[:,1:]==0))
        mask_UV999 = mask_u999 + mask_v999
        Uc[mask_UV999==0] = np.nan
        Vc[mask_UV999==0] = np.nan
    
    #remove cell with uv mask (inactive when all four velocity points are masked, so sum is zero)
    if KCU is not None and KCV is not None:
        mask_u = np.zeros(shape=U1.shape)
        mask_v = np.zeros(shape=U1.shape)
        mask_u[1:,:] = KCU[:-1,:]+KCU[1:,:]
        mask_v[:,1:] = KCV[:,:-1]+KCV[:,1:]
        mask_UV = mask_u + mask_v
        Uc[mask_UV==0] = np.nan
        Vc[mask_UV==0] = np.nan
    
    vel_x = Uc*np.cos(np.deg2rad(ALFAS)) - Vc*np.sin(np.deg2rad(ALFAS))
    vel_y = Uc*np.sin(np.deg2rad(ALFAS)) + Vc*np.cos(np.deg2rad(ALFAS))
    #vel_x = vel_x[1:,1:]
    #vel_y = vel_y[1:,1:]
    vel_magn = np.sqrt(vel_x**2 + vel_y**2)
    direction_naut_deg = np.rad2deg(np.arctan2(vel_y, vel_x))%360
   
    return vel_x, vel_y, vel_magn, direction_naut_deg







