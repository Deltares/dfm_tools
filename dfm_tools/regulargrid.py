# -*- coding: utf-8 -*-
"""
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
      
