# -*- coding: utf-8 -*-
"""
Created on Wed Dec 19
@author: schueder

"""

import numpy as np
import matplotlib.pyplot as plt
import glob
import netCDF4
import os
import pandas as pd
import shutil as sh
import datetime
from .static import nc_format
from .mesh import plot_net
from .SubFile import SubFile

####################################
#           UTILS
####################################


def find_last(var, ss):
    """    
    returns index of last instance of char 'ss' in string 'var'

    Arguments:
        var {str} -- a string
        ss {str} -- a character in var you want to find the last instance of
    """
    ind = 0
    lstInd = ind
    it = 0
    while ind >= 0:
        ind = var.find(ss,ind + it,len(var))
        it = 1
        if ind < 0:
            return lstInd + 1
        lstInd = ind


def change_os(var, enforce=None):
    """ 
    returns a linux path if fed windows and vice versa
    
    Arguments:
        var {str} -- a path
        enfore{str} -- makes sure returned path is equivalent to passes string
    """

    osys = []
    for ch in var:
        if ':' in ch:
            osys = 'windows'
        if ch == '\\':
            osys = 'windows'

    if len(osys) == 0:
        osys = 'linux'

    if enforce is not None:
        # regardless of path os, ensure it is returned as the enforced os
        if enforce == 'linux':
            osys = 'windows'
        elif enforce == 'windows':
            osys = 'linux'

    if '/p/' in var and osys == 'linux':
        return var.replace('/p/','p:\\').replace('/','\\')
    elif osys == 'linux':
        return var.replace('/','\\')
    elif ':\\' in var and osys == 'windows':
        return '/' + var.replace(':\\','/').replace('\\','/')

def nc_station(stations):
    """
    takes a netcdf array and returns a list

    Arguments:
        stations {[type]} -- netCDF4.Dataset.variables['stations][:,:]
    """

    nstations = []
    for line in stations:
        char = ''
        st = line.data
        st = str(list(st))
        st = st.replace('[','')
        st = st.replace(']','')
        for pos,letter in enumerate(st):
            if letter == 'b':
                if pos != 0:
                    if st[pos-1] == "'" and st[pos+1] == "'":
                        char = char + letter
                    else:
                        pass
            else:
                if letter != "'" and letter != ',' and letter != ' ' and letter != '\n':
                    char = char + letter
        nstations.append(char)
    return nstations   


def rst_to_xyz(mapdir, subfile, tind, out, rst = False):
    """
    makes a series of xyz files to be used as an initial condition in an ext file
    uses rst files found in the directory by default, but can also use map files
    based on flag rst = False
    
    Arguments:
        mapdir {str} -- path to map files from which you wish to form initial conditions
        sublist {list} -- list of substance names
        tind {int} -- time index, = -1 takes the last available time
        out {str} -- location of out files
    
    Keyword Arguments:
        rst {bool} -- use rst files and not map files (default: {False})


    example:
    mapdir = 'p:\\11200975-hongkongwaq\\WAQ\\03_baseCase\\A05\\DFM_OUTPUT_HK-FMWAQ\\'
    subfile = r'p:\11200975-hongkongwaq\WAQ\03_baseCase\01_substances\HATS_PCA_v3ep.sub'
    out = 'p:\\11200975-hongkongwaq\\WAQ\\03_baseCase\\A06\\'

    dflowutil.rst_to_xyz(mapdir, subfile, -1, out)
    """
    if not os.path.exists(out):
        os.makedirs(out)

    subs = SubFile(subfile)
    sublist = subs.substances


    files = list(glob.glob(mapdir + '*_map.nc'))
    if len(files) == 0:
        raise FileNotFoundError('No map files available')

    ds = netCDF4.Dataset(files[0])
    params = list(ds.variables.keys())
    physchem = []
    for par in params:
        if 'sa1' in par:
            physchem.append(par.replace('mesh2d_',''))
        elif 'tem' in par:
            physchem.append(par.replace('mesh2d_',''))
        elif 's1' in par:
            physchem.append(par.replace('mesh2d_',''))

    if len(physchem) > 0:
        sublist = sublist + physchem

    varnames = nc_format(files[0])

    print(sublist)

    with open(out + 'ini.ext', 'w') as ext:
        for sub in sublist:
            if sub not in subs.transportable.keys() or subs.transportable[sub] == 'active':
                ext.write('QUANTITY=initialtracer%s\n' % sub)
                ext.write('FILENAME=%s.xyz\n' % sub)
                ext.write('FILETYPE=7\n')
                ext.write('METHOD=5\n')

            else:
                ext.write('QUANTITY=initialwaqbot%s\n' % sub)
                ext.write('FILENAME=%s.xyz\n' % sub)
                ext.write('FILETYPE=7\n')
                ext.write('METHOD=5\n')

            ext.write('OPERAND=O\n')
            ext.write('AVERAGINGTYPE=2\n')
            ext.write('RELATIVESEARCHCELLSIZE=1\n')
            ext.write('\n')


            file_names = ['.xyz']
                
            for file_name in file_names:
                with open(out + '%s' % (sub) + file_name, 'w') as ini:
                    for imap, filei in enumerate(files):
                        mapid = filei[:-7]
                        ds = netCDF4.Dataset(filei)
                        x = ds.variables[varnames['face_x']][:]
                        y = ds.variables[varnames['face_y']][:]
                        # time, space, depth

                        if tind == -1:
                            tmp_times = ds.variables['time'][:]
                            tind = len(tmp_times) - 1
                        if sub not in subs.transportable.keys() or subs.transportable[sub] == 'active':
                            try:
                                s1 = ds.variables['mesh2d_' + sub][tind, :, :]
                                mn_s1 = np.mean(s1, axis = 1)
                            except ValueError:                            
                                mn_s1 = ds.variables['mesh2d_' + sub][tind, :]
                        else:
                            # 2d variable
                            mn_s1 = ds.variables['mesh2d_' + sub][tind, :]
                        
                        for pos, _ in enumerate(x):
                            ini.write('%.6f  %.6f  %.4e\n' % (x[pos], y[pos], mn_s1[pos]))

                    print('finished substance ' + sub)


def find_limit_cell(mapdir):
    """
    plots scatter of limiting cells

    Arguments:
        mapdir {str} --  location of the mapfiles, a directory (str)
    """

    for imap,filei in enumerate(glob.glob(mapdir + '*_map.nc')):
        mapid=filei[filei.index('_map.nc')-4:filei.index('_map.nc')]
        ds = netCDF4.Dataset(filei)
        if imap == 0:
           varnames = nc_format(filei)
                    
        print('processing domain ' + str(mapid))                        
        mesh2d_face_nodes=ds.variables[varnames['cellnodes']][:]
        try:
            domainno=ds.variables[varnames['domain_number']][:]
            ghost = True
        except:
            ghost = False
            print('missing extra information, ghost cells not cleaned')

        tridata= dflow_grid_2_tri(mesh2d_face_nodes)
        index = tridata['index']
        tri=tridata['triangles']
        xnode=ds.variables[varnames['xnode']][:]
        ynode=ds.variables[varnames['ynode']][:]
        name = 'numlimdt'
        var=ds.variables[name][len(ds.variables['time'])-1,:] 
            
        newvar=var[index.astype(np.int64)]
        tri2=tri-1

        if ghost:
            selectcells=(domainno==np.int(mapid))
            selectcellstri=selectcells[index.astype(np.int64)]
            totalselected=np.array([selectcellstri]).squeeze()
            tri2=tri2[totalselected,:]
            newvar=newvar[totalselected]

        plt.tripcolor(xnode,ynode,tri2,facecolors=np.nan*newvar,edgecolors='k',cmap='jet')
        ind = np.argmax(newvar)
        # return tri2,xnode,ynode
        tind = np.array([int(ii) for ii in tri2[ind]])
        plt.scatter(xnode[tind],ynode[tind],20,'r')
        plt.text(xnode[tind[0]],ynode[tind[0]],('%.2e' % newvar[ind]))
        plt.gca().set_aspect('equal', adjustable='box')


def read_polygon(file):
    '''
    imports a polygon as nx2 array
    '''
    with open(file) as ldbfile:
        lines = ldbfile.readlines()
        X = []
        Y = []
        for ind, row in enumerate(lines):
            if '*' not in row:
                if '.' in row and '999.' not in row and '-999' not in row:
                    line = row2array(row)
                    X.append(float(line[0]))
                    Y.append(float(line[1]))            
                elif len(row.replace('\n','')) > 0:
                    X.append(np.nan)
                    Y.append(np.nan)
                else:
                    pass
    return np.array([X, Y]).T


def read_pli(var):
    '''
    reads a pli boundary into an array
    '''
    with open(var) as plifile:
        lines = plifile.readlines()
        X = []
        Y = []
        for ind, row in enumerate(lines):
            if '.' in row:
                line = row2array(row)
                X.append(float(line[0]))
                Y.append(float(line[1]))
    return np.array([X, Y]).T


def read_pliz(var):
    '''
    reads a pli boundary into an array
    '''
    with open(var) as plifile:
        lines = plifile.readlines()
        X = []
        Y = []
        Z1 = []
        Z2 = []
        Z3 = []
        for ind, row in enumerate(lines):
            if '.' in row:
                line = row2array(row)
                X.append(float(line[0]))
                Y.append(float(line[1]))
                Z1.append(float(line[2]))
                Z2.append(float(line[3]))
                Z3.append(float(line[4]))
    return np.array([X, Y, Z1, Z2, Z3]).T


def check_data_path(boundaries, root, name, bnd_type):
    '''
    makes all paths os accessible via literal absolute string path
    acts on strings that have already been converted to windows
    does not touch if absolute path is found
    '''
    os.chdir(root)
    orig = os.getcwd()

    if '..' in boundaries[name][bnd_type]['pli_loc']:
        # is relative, need to navigate
        up = boundaries[name][bnd_type]['pli_loc'].count('..')
        for jump in range(0, up):
            try:
                os.chdir('..\\')
            except:
                os.chdir('../')
        new_root = os.getcwd()
        boundaries[name][bnd_type]['pli_loc']  = os.path.join(new_root, boundaries[name][bnd_type]['pli_loc'].replace('..\\','').replace('../',''))
        os.chdir(orig)

    elif '\\' not in boundaries[name][bnd_type]['pli_loc'] and '/' not in boundaries[name][bnd_type]['pli_loc']:
        # is local
        boundaries[name][bnd_type]['pli_loc']  = root + boundaries[name][bnd_type]['pli_loc']  
    
    if '..' in boundaries[name][bnd_type]['data_loc']:
        # is relative, need to navigate
        up = boundaries[name][bnd_type]['data_loc'].count('..')
        for jump in range(0, up):
            os.chdir('..\\')
        new_root = os.getcwd()
        boundaries[name][bnd_type]['data_loc']  = new_root + '\\' + boundaries[name][bnd_type]['data_loc'].replace('..\\','')
        os.chdir(orig)

    elif '\\' not in boundaries[name][bnd_type]['data_loc'] and '/' not in boundaries[name][bnd_type]['data_loc']:
        # is local
        boundaries[name][bnd_type]['data_loc']  = root + boundaries[name][bnd_type]['data_loc']  

    return boundaries


def show_waq_segment(grd,nolay,segments):
    """
    visualize the location of a delwaq segment in x,y given the segment number

    
    Arguments:
        grd {str} -- a path to a waq geom
        nolay {int} -- the number of layers (not known to the WAQ geom)
        segments {dict} -- a dictionary with name (key) segment (int) 
    """

    varnames = nc_format(grd)
    ds    = netCDF4.Dataset(grd)
    x     = ds.variables[varnames['face_x']][:]
    y     = ds.variables[varnames['face_y']][:]
    elem  = ds.variables[varnames['cellnodes']][:,:]

    segspl = np.size(elem,0)
    # make an array of all of the numbers of the first segment in each layer,
    # including the nth + 1 layer
    blay   = np.arange(1,segspl*nolay+2,segspl)

    plot_net(grd)

    plt.gca().set_aspect('equal', adjustable='box')

    for sind, ss in enumerate(segments.keys()):
        # find the first segment in this layer     
        first_in_next_layer = np.min([ii for ii,jj in enumerate(blay) if segments[ss] < jj])                           
        topseg = segments[ss] - blay[first_in_next_layer - 1]
        ind = topseg-1
        xi = x[ind]
        yi = y[ind]
        plt.scatter(xi, yi, 30, 'r')
        plt.text(xi, yi, ss)


def pdistf(X1, Y1, X2, Y2):
    '''
    returns array of euclidean distances between a point and an array
    '''
    return np.sqrt((X2 - X1)**2 + (Y2 - Y1)**2)


def row2array(line_orig):
    '''
    takes a string of space seperated floats and returns an array
    '''
    line = line_orig.split(' ')
    arr = []
    for ch in line:
        try:
            val = float(ch)
            arr.append(val)
        except:
            pass
    if len(arr) > 1:
        return np.array(arr)
    else:
        line = line_orig.split('\t')
        arr = []
        for ch in line:
            try:
                val = float(ch)
                arr.append(val)
            except:
                pass
        if len(arr) > 1:
            return np.array(arr)
