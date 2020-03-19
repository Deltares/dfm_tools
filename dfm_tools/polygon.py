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
        
    def fromfile(file_pol, pd_output=False):
        import numpy as np
        import pandas as pd
        print('WARNING: the function dfm_tools.polygon.Polygon.fromfile() will be improved, outputformat will change')
        '''
        reads a pli boundary into an array
        tekal_output==True: create pandas array of tekal file, including comments
        tekal_output==False (default): create list of numpy arrays and list of polygon names
        '''

        with open(file_pol) as plifile:
            lines = plifile.readlines()
            pol_data_list = []
            pol_name_list = []
            pol_comment_list = []
            pol_comment_temp = []
            pol_data_pd_list = []
            #pol_object_list = []
            iLine=0
            while iLine<=len(lines)-1:
                line_str = lines[iLine].split()
                try:
                    float(line_str[0]) #this line should start with a string, so this should fail
                    raise Exception('ERROR: line_str[0] is not a string: (%s, on line %d), polygon file format might be wrong'%(line_str, iLine))
                except:
                    pass
                iLine=iLine+1
                if line_str[0].startswith('*'): #comments
                    line_str_joined = ' '.join(line_str)
                    pol_comment_temp.append(line_str_joined)
                else:
                    #the start of a new polygon
                    pol_name_list.append(line_str[0]) #get name from name line_str
                    pol_comment_list.append(pol_comment_temp)
                    
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
                    
                    if len(pol_comment_temp)==data_pol.shape[1]: #expected format
                        pol_data_pd = pd.DataFrame(data_pol,columns=pol_comment_temp)
                        if 'date' in pol_comment_temp[0].lower():
                            try:
                                pol_data_datetime = pd.to_datetime(pol_data_pd.iloc[:,0]*1000000+pol_data_pd.iloc[:,1],format='%Y%m%d%H%M%S')
                                pol_data_pd.insert(0, 'datetime', pol_data_datetime)
                            except:
                                print('WARNING: conversion from date/time column to datetime failed, incorrect format of first two columns?')
                    else:
                        pol_data_pd = pd.DataFrame(data_pol)
                    pol_data_pd_list.append(pol_data_pd)
                    
                    # reset comment list for next block
                    pol_comment_temp = []
                    #pol_object = Polygon(data=data_pol, name=line_str[0])
                    #pol_object_list.append(pol_object)
        if pd_output==True:
            return pol_data_pd_list
        else:
            return pol_data_list, pol_name_list, pol_comment_list
            

    def frominteractive(fig, ax):
        raise Exception('ERROR: interactive polygon from file is not yet possible')

        
        #line_array = linebuilder.line_array
        pol_frominput = Polygon(linebuilder.xs, linebuilder.ys)
        
        return pol_frominput