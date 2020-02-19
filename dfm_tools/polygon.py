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
    def __init__(self, x, y, name):
        import numpy as np
        self.x = x
        self.y = y
        self.name = name
        self.line_array = np.c_[self.x, self.y]
        
    def fromfile(file_pol):
        import numpy as np
        print('WARNING: this function is to be improved, outputformat will change')
        '''
        reads a pli boundary into an array
        '''
        def row2array(text):
            '''
            takes a string of space seperated floats and returns an array
            '''
            text = text.split()
            arr = []
            for ch in text:
                try:
                    val = float(ch)
                    arr.append(val)
                except:
                    val = str(ch)
                    arr.append(val)
            return np.array(arr)

        with open(file_pol) as plifile:
            lines = plifile.readlines()
            X = list([]) #list of pol x
            Y = list([])#list of pol y
            NAME=list([])
            countpol=0       
            count=0
            
            while count<=len(lines)-1:
                row=lines[count]
                count=count+1
                line = row2array(row)
                if type(line[0]) is np.str_ or type(line[0]) is np.string_ or type(line[0]) is str:
                    if not line[0].startswith('*'):
                        countpol=countpol+1
                        #add to new polygon
                        namei=str(line[0])
                        line=row2array(lines[count])
                        count=count+1
                        len_pol=int(line[0])
                        xpol=list([])
                        ypol=list([])
                        for ipoint in np.arange(len_pol):
                            xpol.append(row2array(lines[count])[0])
                            ypol.append(row2array(lines[count])[1])             
                            count=count+1                
                        #pol finished, add to total
                        X.append(xpol)
                        Y.append(ypol)
                        NAME.append(namei)
        data_pol = Polygon(x=X, y=Y, name=NAME)
        return data_pol
        

    def frominteractive(fig, ax):
        raise Exception('ERROR: interactive polygon from file is not yet possible')

        
        #line_array = linebuilder.line_array
        pol_frominput = Polygon(linebuilder.xs, linebuilder.ys)
        
        return pol_frominput