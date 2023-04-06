# -*- coding: utf-8 -*-
"""
Created on Thu Apr  6 21:25:48 2023

@author: veenstra
"""

#TODO: this can probably be used to select all faces within a polygon

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

