# -*- coding: utf-8 -*-
"""
Created on Sat Mar 21 22:07:17 2020

@author: veenstra
"""

def center2corner(array_x, array_y):
    import numpy as np
    
    from dfm_tools.corner2center import corner2center
    
    x_cen_nobnd = corner2center(array_x)
    x_cen_nobnd_diff_ax0 = np.diff(x_cen_nobnd, axis=0)
    add_top = x_cen_nobnd[0,:] - x_cen_nobnd_diff_ax0[0,:]
    add_bot = x_cen_nobnd[-1,:] + x_cen_nobnd_diff_ax0[-1,:]
    x_cen_vertbnd = np.vstack([add_top, x_cen_nobnd, add_bot])
    x_cen_nobnd_diff_ax1 = np.diff(x_cen_vertbnd, axis=1)
    add_left = x_cen_vertbnd[:,0] - x_cen_nobnd_diff_ax1[:,0]
    add_right = x_cen_vertbnd[:,-1] + x_cen_nobnd_diff_ax1[:,-1]
    x_cen_withbnd = np.vstack([add_left.T, x_cen_vertbnd.T, add_right.T]).T
    
    y_cen_nobnd = corner2center(array_y)
    y_cen_nobnd_diff_ax0 = np.diff(y_cen_nobnd, axis=0)
    add_top = y_cen_nobnd[0,:] - y_cen_nobnd_diff_ax0[0,:]
    add_bot = y_cen_nobnd[-1,:] + y_cen_nobnd_diff_ax0[-1,:]
    y_cen_vertbnd = np.vstack([add_top, y_cen_nobnd, add_bot])
    y_cen_nobnd_diff_ax1 = np.diff(y_cen_vertbnd, axis=1)
    add_left = y_cen_vertbnd[:,0] - y_cen_nobnd_diff_ax1[:,0]
    add_right = y_cen_vertbnd[:,-1] + y_cen_nobnd_diff_ax1[:,-1]
    y_cen_withbnd = np.vstack([add_left.T, y_cen_vertbnd.T, add_right.T]).T
    
    return x_cen_withbnd, y_cen_withbnd
