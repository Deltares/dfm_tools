# -*- coding: utf-8 -*-
"""
Created on Wed Feb 12 01:55:57 2020

@author: veenstra
"""

class LineBuilder:
    def __init__(self, line):
        self.line = line
        self.xs = list(line.get_xdata())
        self.ys = list(line.get_ydata())
        self.cid = line.figure.canvas.mpl_connect('button_press_event', self)

    def __call__(self, event):
        import numpy as np
        print('click', event)
        if event.inaxes!=self.line.axes: return
        self.xs.append(event.xdata)
        self.ys.append(event.ydata)
        self.line_array = np.c_[self.xs, self.ys]
        self.line.set_data(self.xs, self.ys)
        self.line.figure.canvas.draw()


class Polygon:
    def __init__(self, x, y):
        import numpy as np
        self.x = x
        self.y = y
        self.line_array = np.c_[self.xs, self.ys]
        
    def fromfile(file_pol):
        raise Exception('ERROR: reading polygon from file is not yet possible')

    def frominteractive(fig, ax):
        raise Exception('ERROR: interactive polygon from file is not yet possible')

        
        #line_array = linebuilder.line_array
        pol_frominput = Polygon(linebuilder.xs, linebuilder.ys)
        
        return pol_frominput