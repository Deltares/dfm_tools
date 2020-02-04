import numpy as np
from .utils import nc_format
import matplotlib.pyplot as plt
import netCDF4
import glob


def dflow_grid_2_tri(mesh2d_face_nodes):
    """
    returns the indicies of nodes in a mesh grid that constitute each facenode

    Arguments:
        mesh2d_face_nodes {np.array} -- mxn array where m is the number of face
        nodes and n is the maximum number of nodes constituting a face node
    """

    n = mesh2d_face_nodes.shape[0]
    count = np.sum(~np.isnan(mesh2d_face_nodes), axis=1)

    m4 = np.sum(count == 4)
    m5 = np.sum(count == 5)
    m6 = np.sum(count == 6)

    index3 = np.arange(0, n)
    index4 = np.where(count == 4)
    index5 = np.where(count == 5)
    index6 = np.where(count == 6)

    # initialize
    tri = np.zeros([n + m4 + 2 * m5 + 3 * m6, 3])
    tri[0:n, :] = mesh2d_face_nodes[0:n, 0:3]

    index = np.zeros(tri.shape[0], dtype=np.int64)
    index[0:n] = index3
    # 4 nodes
    offset = n
    tri[offset + np.arange(0, m4), 0:3] = mesh2d_face_nodes[np.ix_(index4[0], np.asarray([0, 2, 3]))]
    index[offset + np.arange(0, m4)] = index4[0]
    offset = offset + m4

    # 5 nodes
    if m5 > 0:
        tri[offset + np.arange(0, m5), 0:3] = mesh2d_face_nodes[np.ix_(index5[0], np.asarray([0, 2, 3]))]
        index[offset + np.arange(0, m5)] = index5[0]
        offset = offset + m5

        tri[offset + np.arange(0, m5), 0:3] = mesh2d_face_nodes[np.ix_(index5[0], np.asarray([0, 3, 4]))]
        index[offset + np.arange(0, m5)] = index5[0]
        offset = offset + m5
    # 6 nodes
    if m6 > 0:
        tri[offset + np.arange(0, m6), 0:3] = mesh2d_face_nodes[np.ix_(index6[0], np.asarray([0, 2, 3]))]
        index[offset + np.arange(0, m6)] = index6[0]
        offset = offset + m6

        tri[offset + np.arange(0, m6), 0:3] = mesh2d_face_nodes[np.ix_(index6[0], np.asarray([0, 3, 4]))]
        index[offset + np.arange(0, m6)] = index6[0]
        offset = offset + m6

        tri[offset + np.arange(0, m6), 0:3] = mesh2d_face_nodes[np.ix_(index6[0], np.asarray([0, 4, 5]))]
        index[offset + np.arange(0, m6)] = index6[0]

    return {'triangles': tri, 'index': index}


def plot_nc_map(mapdir, elem, time, depth=None, layer=None, lim=None, c_map='jet'):
    """
    plots a 2D patch plot of a variable in a certain layer or depth

    Arguments:
        mapdir {str} -- location of the mapfiles, a directory
        elem {str} -- name of constituent, such as 'salinity', must be in map1 and map4 dictionary
        time {int} -- time index

    Keyword Arguments:
        depth {float} -- depth, negative down (default: {None})
        layer {int} -- layer index (default: {None})
        lim {tuple} -- color limit to use (default: {None})
        c_map {str} -- colormap to use (default: {'jet'})
    """
    '''

    '''
    if depth is None and layer is None:
        print('Error: depth or layer must be specified')
    else:
        for imap, filei in enumerate(glob.glob(mapdir + '*_map.nc')):
            mapid = filei[filei.index('_map.nc') - 4:filei.index('_map.nc')]
            ds = netCDF4.Dataset(filei)
            if imap == 0:
                varnames = nc_format(filei)
            print('processing domain ' + str(mapid))
            mesh2d_face_nodes = ds.variables[varnames['cellnodes']][:]
            try:
                domainno = ds.variables[varnames['domain_number']][:]
                ghost = True
            except:
                ghost = False
                print('missing extra information, ghost cells not cleaned')

            tridata = dflow_grid_2_tri(mesh2d_face_nodes)
            index = tridata['index']
            tri = tridata['triangles']
            xnode = ds.variables[varnames['xnode']][:]
            ynode = ds.variables[varnames['ynode']][:]
            try:
                name = varnames[elem]
            except(KeyError):
                # not a standard constituent, try WAQ convention
                if ('mesh2d_' + elem) in ds.variables.keys():
                    name = 'mesh2d_' + elem
                else:
                    if (elem) in ds.variables.keys():
                        name = elem
                    else:
                        print('ERROR: cannot guess element name')
                        raise

            if layer is not None:
                if layer != 0:
                    var = ds.variables[name][time, :, layer - 1]
                else:
                    var = ds.variables[name][time, :]

            elif depth is not None:
                wd = ds.variables['mesh2d_waterdepth'][time, :]
                wd = wd.reshape((-1, 1))
                frac = ds.variables['mesh2d_interface_sigma'][:]
                frac = frac.reshape((-1, 1))
                depths = np.dot(frac, wd.T)
                # find the number of cell interfaces in each column that the
                # query is less than
                # the sum is always >= 1 since the min is ~ 0
                # is the sum == 1, we want to query the first index (0)
                # so we minus one from the sum to get the index

                ind = np.sum(depths[:, :] > depth, axis=0) - 1

                # for cells shallower than the query depth, clip to bottom cell

                too_deep = ind == len(frac) - 1
                too_deep = [ii for ii, jj in enumerate(too_deep) if jj]
                ind[too_deep] = len(frac) - 2

                # need to flip array because top cells are last aray elements
                # ex. 0->19, 1->18, 2->17 if no. layers = 20
                # -2 because frac is no. layers + 1
                ind = len(frac) - 2 - ind
                # ind contains the indicies for each xy position that
                # will give the desired depth

                # get indicies for each position's depth

                var = ds.variables[name][time, :, :]
                # this does not reduce the dimensions for some reason...
                # var = var[ind]
                # var[too_deep] = np.nan

                arr = np.zeros((var.shape[0], 1))
                for pos in range(0, var.shape[0]):
                    arr[pos] = var[pos, ind[pos]]
                var = arr
                var[too_deep] = np.nan

            newvar = var[index.astype(np.int64)]
            tri2 = tri - 1

            if ghost:
                selectcells = (domainno == np.int(mapid))
                selectcellstri = selectcells[index.astype(np.int64)]
                totalselected = np.array([selectcellstri]).squeeze()
                tri2 = tri2[totalselected, :]
                newvar = newvar[totalselected]

            # append partition to the image
            time_unit = ds.variables['time'].units
            mod_ref = pd.Timestamp(time_unit.replace('seconds since ', ''))
            times = ds.variables['time'][:]
            times = np.array([mod_ref + pd.Timedelta(seconds=tt) for tt in times])
            print('plotting time = %s ' % times[time])
            nmask = np.isnan(newvar.ravel())
            plt.tripcolor(xnode, ynode, tri2, facecolors=newvar.ravel(), edgecolors='none', cmap=c_map,
                          mask=nmask)  # , shading = 'gouraud')
            if lim is not None:
                plt.clim(lim)
            plt.gca().set_aspect('equal', adjustable='box')

    return plt.gcf()


def plot_net(net_file):
    '''
    plots a netCDF mesh
    is slower than matlab counterpart
    '''
    dat = netCDF4.Dataset(net_file)
    varnames = nc_format(net_file)

    node_x = dat.variables[varnames['xnode']][:]
    node_y = dat.variables[varnames['ynode']][:]
    face = dat.variables[varnames['cellnodes']][:, :]

    # written as transpose to eliminate transposes in loop
    # sometimes they are nan, and sometimes they are masked?
    xy = np.nan * np.zeros((2, (np.sum(np.sum(~np.isnan(face))) + len(face[:, 0]) - 1 + len(face[:, 0]))))
    inter = 0

    for ff in range(0, len(face)):
        face_nodes = face[ff, ~np.isnan(face[ff, :])]
        # from segment index (min == 1) to position index (min == 0)
        face_nodes = face_nodes[face_nodes.mask == False] - 1
        xy[:, inter:(inter) + len(face_nodes)] = np.array([node_x[face_nodes], node_y[face_nodes]])

        # leave a nan in between for plotting
        xy[:, inter + len(face_nodes)] = np.array([node_x[face_nodes[0]], node_y[face_nodes[0]]])
        inter = inter + len(face_nodes) + 2

    plt.plot(xy[0, :], xy[1, :], '-')
    plt.gca().set_aspect('equal', adjustable='box')