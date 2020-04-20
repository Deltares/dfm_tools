# -*- coding: utf-8 -*-
"""
Created on Fri Apr 17 13:02:48 2020

@author: veenstra
"""
#WARNING, THIS IS WORK IN PROGRESS, NOT TESTED YET

def convert_xycoords(loc_in_xy,epsg_in,epsg_out): #default epsg_out value is RD
    from pyproj import Proj, transform
    import pandas as pd
    import numpy as np
    
    if epsg_in == epsg_out:
        loc_out_xy = loc_in_xy
    else:
        inProj = Proj(init='epsg:%i'%(epsg_in))
        outProj = Proj(init='epsg:%i'%(epsg_out))
        
        if isinstance(loc_in_xy, pd.DataFrame):
            x_in = loc_in_xy['RDx'].values
            y_in = loc_in_xy['RDy'].values
            x_out,y_out = transform(inProj,outProj,x_in,y_in)
            data_forpd = np.empty([len(x_out),2])
            data_forpd[:,0] = x_out
            data_forpd[:,1] = y_out
            loc_out_xy = pd.DataFrame(data=data_forpd, columns=['RDx','RDy'])
        else:
            x_in = loc_in_xy[0]
            y_in = loc_in_xy[1]
            x,y = transform(inProj,outProj,x_in,y_in)
            loc_out_xy = [x,y]
        
    return loc_out_xy







