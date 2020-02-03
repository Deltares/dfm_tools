import numpy as np
import pandas as pd
from .utils import find_last, check_data_path, change_os, pdistf, row2array
from .SubFile import SubFile
import datetime

def boundary_from_ext(var):
    """
    returns a dictionary containing the boundary names, types, location files and data files
    from a boundary definition

    Arguments:
        var {str} -- path to ext file
    """

    root = var[:find_last(var, '\\')]
    boundaries = {}
    with open(var, 'r') as nmf:
        page = nmf.readlines()
        ext_type = 'old'
        for line, text in enumerate(page):
            if '[boundary]' in text or '.bc' in text:
                ext_type = 'new'
        if ext_type == 'new':
            for line, text in enumerate(page):
                if '*' not in text:
                    if '[boundary]' in text:
                        name = page[line + 2].replace('locationfile', '').replace('=', '').replace('.pliz', '').replace(
                            '.pli', '').replace('\n', '').strip()
                        if '/' in name:
                            # if it is a path
                            name = name[find_last(name, '/'):]
                        if name not in boundaries.keys():
                            # append if new boundary
                            boundaries[name] = {}

                        bnd_type = page[line + 1].replace('quantity', '').replace('=', '').replace('\n', '').strip()
                        boundaries[name][bnd_type] = {}

                        boundaries[name][bnd_type]['type'] = bnd_type
                        # change the os to make it windows readable for later processing
                        boundaries[name][bnd_type]['pli_loc_orig'] = page[line + 2].replace('locationfile', '').replace(
                            '=', '').replace('\n', '').strip()
                        boundaries[name][bnd_type]['data_loc_orig'] = page[line + 3].replace('forcingfile', '').replace(
                            '=', '').replace('\n', '').strip()

                        boundaries[name][bnd_type]['pli_loc'] = change_os(
                            page[line + 2].replace('locationfile', '').replace('=', '').replace('\n', '')).strip()
                        boundaries[name][bnd_type]['data_loc'] = change_os(
                            page[line + 3].replace('forcingfile', '').replace('=', '').replace('\n', '')).strip()

                        boundaries = check_data_path(boundaries, root, name, bnd_type)

        else:
            # old style
            for line, text in enumerate(page):
                if '*' not in text:
                    if 'QUANTITY' in text and '=' in text:
                        name = page[line + 1].replace('FILENAME', '').replace('=', '').replace('.pliz', '').replace(
                            '.pli', '').replace('\n', '').strip()
                        if '/' in name or '\\' in name:
                            # is a path and we need to extract the name
                            name = name[find_last(name, '/'):]
                        if name not in boundaries.keys():
                            boundaries[name] = {}

                        # NOTE: THIS AHRD CODED ORDERING IS INVALIDATED BY THE USE OF VAR NAME
                        bnd_type = text.replace('QUANTITY', '').replace('=', '').replace('\n', '').strip()
                        boundaries[name][bnd_type] = {}
                        boundaries[name][bnd_type]['type'] = bnd_type
                        boundaries[name][bnd_type]['pli_loc_orig'] = page[line + 1].replace('FILENAME', '').replace('=',
                                                                                                                    '').replace(
                            '\n', '').strip()
                        boundaries[name][bnd_type]['data_loc_orig'] = page[line + 1].replace('FILENAME', '').replace(
                            '=', '').replace('\n', '').replace('.pliz', '.pli').replace('.pli', '.tim').strip()

                        boundaries[name][bnd_type]['pli_loc'] = change_os(
                            page[line + 1].replace('FILENAME', '').replace('=', '').replace('\n', '')).strip()
                        boundaries[name][bnd_type]['data_loc'] = change_os(
                            page[line + 1].replace('FILENAME', '').replace('=', '').replace('\n', '').replace('.pliz',
                                                                                                              '.pli').replace(
                                '.pli', '.tim')).strip()

                        boundaries[name][bnd_type]['FILETYPE'] = page[line + 2].replace('FILETYPE', '').replace('=',
                                                                                                                '').replace(
                            '\n', '').strip()
                        boundaries[name][bnd_type]['METHOD'] = page[line + 3].replace('METHOD', '').replace('=',
                                                                                                            '').replace(
                            '\n', '').strip()
                        boundaries[name][bnd_type]['OPERAND'] = page[line + 4].replace('OPERAND', '').replace('=',
                                                                                                              '').replace(
                            '\n', '').strip()

                        boundaries = check_data_path(boundaries, root, name, bnd_type)

    return boundaries


def read_bc(pli_file, bc_file):
    '''
    reads a bc file into a format useful for plotting a cross sections
    usage restrictions:
    * does not work for uxuyadvectionboundaries
    * file must contain only one variable
    * distance does not check for projection, so will be wrong if spherical

    plotting is expected to look as follows:

    data = read_bc(pli_file, bc_file)

    meshX, meshY = np.meshgrid(data['distance'], data['zprofile'])
    C = np.squeeze(data['salinitybnd'][:,:,time])
    plt.pcolormesh(meshX, meshY, C)

    '''
    non_data = ['ame', 'orcing', 'unction', 'ertical', 'ime', 'uantity', 'nit', 'ince']
    pli = read_pli(pli_file)
    data = {}
    # create an array of distances
    dist = np.zeros((len(pli)))
    for position in range(1, len(pli)):
        dist[position] = dist[position - 1] + np.abs(
            pdistf(pli[position, 0], pli[position, 1], pli[position - 1, 0], pli[position - 1, 1]))

    with open(bc_file, 'r') as bc:
        ind = []
        page = bc.readlines()
        for row, line in enumerate(page):
            # first pass, obtain metadata
            if '[forcing]' in line:
                ind.append(row)
            if 'Vertical position type          = zdatum' in line:
                data['vertical position type'] = 'zdatum'
            if 'Vertical position specification' in line:
                line = line.replace('Vertical position specification =', '')
                arr = row2array(line)
                data['zprofile'] = arr
            if 'uantity' in line and 'time' not in line:
                data['quantity'] = line.replace('Quantity', '').replace('quantity', '').replace('=', '').strip()
            if 'nit' in line:
                data['unit'] = line.replace('Unit', '').replace('unit', '').replace('=', '').strip()
            if 'ince' in line or 'INCE' in line:
                if 'inutes' in line or 'INUTES' in line:
                    data['timeunit'] = 'minutes'
                elif 'econds' in line or 'ECONDS' in line:
                    data['timeunit'] = 'seconds'

                line = line.split(' ')
                time = line[-2] + ' ' + line[-1].replace('\n', '')
                data['reftime'] = pd.Timestamp(time)

    assert (len(ind) == len(pli))
    if 'vertical position type' not in data.keys():
        print(
            'ERROR: vertical specification is not zdatum, bc file zprofile is not self describing and not implemented')
        return None
    else:
        with open(bc_file, 'r') as bc:
            # second pass, load data into memory
            page = bc.readlines()
            # estimate number of times based on first position
            # purpose is for allocation of array
            times = []
            for row in np.arange(ind[0], ind[1]):
                line = page[row]
                chk = sum([word in line for word in non_data])
                if chk == 0 and '.' in line:
                    arr = row2array(line)
                    time_val = arr[0]
                    if data['timeunit'] == 'minutes':
                        times.append(data['reftime'] + pd.Timedelta(days=time_val / 1440.0))
                    elif data['timeunit'] == 'seconds':
                        times.append(data['reftime'] + pd.Timedelta(days=time_val / 86400.0))

            data['distance'] = dist
            data['coordinates'] = pli
            # data['times'] = np.array(times)
            data['times'] = times

            data[data['quantity']] = np.zeros((len(data['zprofile']), len(ind), len(times)))

            for position in range(0, len(ind)):
                tt = 0
                if position + 1 == len(ind):
                    curr_data = page[ind[position]:]
                else:
                    curr_data = page[ind[position]: ind[position + 1]]
                for line in curr_data:
                    chk = sum([word in line for word in non_data])
                    if chk == 0 and '.' in line:
                        arr = row2array(line)
                        data[data['quantity']][:, position, tt] = arr[1:]
                        tt += 1

        return data


def make_tim_ts(tref, var):
    return tref + datetime.timedelta(minutes=var)


def read_sours(file, tref, sal=True, temp=True, subfile=None):
    """
    reads a sours file to a pandas dataframe

    Arguments:
        file {str} -- path to .tim

    Keyword Arguments:
        sal {bool} -- salinity is modelled (default: {True})
        temp {bool} -- temperature is modelled (default: {True})
        subfile {str/list} -- path to sub file, or list (default: {None})
    """
    cols = list()
    cols.append('time')
    cols.append('flow')

    if sal:
        cols.append('salinity')
    if temp:
        cols.append('temperature')

    if subfile is not None:

        if isinstance(subfile, str):
            all_subs = SubFile(subfile).substances
        elif isinstance(subfile, list):
            all_subs = subfile
        else:
            print('ERROR: incorrect sub file')
            raise

        for sub in all_subs:
            if 'S1' in sub or 'S2' in sub or 'Det' in sub or sub == 'SOD':
                pass
            else:
                cols.append(sub)
    else:
        print('Error: substances not specified')
        raise FileNotFoundError

    df = pd.DataFrame(columns=cols)
    row = 0
    with open(file, 'r') as tim:
        lines = tim.readlines()
        for ind, line in enumerate(lines):
            if '*' in line:
                pass
            if len(line) == 0:
                pass
            else:
                tmp_df = pd.DataFrame(columns=cols, index=[row])
                vec = row2array(line)
                for col, val in enumerate(vec):
                    if cols[col] == 'time':
                        val = make_tim_ts(tref, val)
                    tmp_df[cols[col]] = val
                row += 1
                df = pd.concat([df, tmp_df], axis=0, ignore_index=False)

    return df