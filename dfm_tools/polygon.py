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
    def __init__(self, data, name):
        #import numpy as np
        self.data = data
        self.name = name
        #self.line_array = np.c_[self.x, self.y]
        
    def fromfile(file_pol):
        import numpy as np
        print('WARNING: the function dfm_tools.polygon.Polygon.fromfile() will be improved, outputformat will change')
        '''
        reads a pli boundary into an array
        '''

        with open(file_pol) as plifile:
            lines = plifile.readlines()
            pol_data_list = []
            pol_name_list = []
            #pol_object_list = []
            iLine=0
            while iLine<=len(lines)-1:
                line_str = lines[iLine].split()
                iLine=iLine+1
                try:
                    float(line_str[0]) #this line should start with a string, so this should fail
                    raise Exception('ERROR: line_str[0] is not a string: (%s, somewhere around line %d), polygon file format might be wrong'%(line_str, iLine))
                except:
                    pass
                if not line_str[0].startswith('*'):
                    #the start of a new polygon
                    pol_name_list.append(line_str[0]) #get name from name line_str
                    
                    #convert shape numbers to int
                    line_shp = [int(x) for x in lines[iLine].split()]
                    iLine=iLine+1
                    len_pol=int(line_shp[0])
                    
                    #retrieve data numbers as floats
                    data_pol_str = lines[iLine:iLine+len_pol]
                    iLine = iLine+len_pol
                    data_pol = np.array([x.split() for x in data_pol_str],dtype=float)
                    data_pol_nanbool = data_pol==999.999
                    data_pol[data_pol_nanbool] = np.nan
                    pol_data_list.append(data_pol)
                    
                    #pol_object = Polygon(data=data_pol, name=line_str[0])
                    #pol_object_list.append(pol_object)
        return pol_data_list, pol_name_list#, pol_object_list
        

    def frominteractive(fig, ax):
        raise Exception('ERROR: interactive polygon from file is not yet possible')

        
        #line_array = linebuilder.line_array
        pol_frominput = Polygon(linebuilder.xs, linebuilder.ys)
        
        return pol_frominput