# -*- coding: utf-8 -*-
"""
Created on Wed Feb 19 12:18:27 2020

@author: irazoqui
"""
import numpy as np


def read_pol_multiple(file_pol):
    '''
    reads a pli boundary into an array
    '''
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
                
    return X,Y,NAME

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


file_pol = r'c:\DATA\werkmap\dfm_tools_testdata\DFM_3D_z_Grevelingen\geometry\structures\Grevelingen-FM_BL_fxw.pliz'
file_pol = r'c:\DATA\werkmap\dfm_tools_testdata\world.ldb'
pol_grevX,pol_grevY,pol_grevNAME = read_pol_multiple(file_pol)
