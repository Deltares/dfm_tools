import shutil as sh
from .utils import find_last, change_os, row2array
from .static import nc_format
from .boundary import boundary_from_ext, read_sours
from .SubFile import SubFile
import datetime
import numpy as np
import netCDF4
import pandas as pd
import os

class DFMWAQModel():

    def __init__(self, mdu, ext, subfile, new_dir, tref, ini=None, run_sys='linux', v=None, cores=None,
                 loads_data=None, bounds_data=None,
                 process_path=None):
        """
        DFMWAQ model initialized from DFM model inputs

        call DFMWAQModel.build() to build an initialized model
        Arguments:
            mdu {str} -- path to mdu
            ext {list} -- [ext1, ext2], if only one, still must be a list
            subfile {dflowutil.SubFile} -- a dflowutil.SubFile object
            new_dir {path} -- path where model will be built
            tref {datetime.datetime}
            ini {dict} -- sub name value pair for initial conditions. May be empty.
            v {str} -- version - i.e. 1.2.56.xxx.
            cores {list} -- [nodes, threads]
        """

        if bounds_data is None:
            bounds_data = dict()
            self.bounds_data = bounds_data
        if ini is None:
            ini = dict()
            self.ini = ini
        if loads_data is None:
            loads_data = dict()
            self.loads_data = loads_data

        if v is not None:
            self.version = v
        else:
            self.version = None
        if cores is not None:
            self.nodes = cores[0]
            self.threads = cores[1]
        else:
            self.nodes = None
            self.threads = None

        self.source_mdu = mdu
        self.mdu = os.path.join(new_dir, os.path.split(self.source_mdu)[1])
        self.ext = ext

        self.subfile = subfile
        self.substances = subfile.substances
        self.out = new_dir
        self.run_sys = run_sys

        if isinstance(loads_data, list):
            print('processing loads...')
            self.loads_data = self.process_loads_data(loads_data)
            print('loads processed')
        elif not isinstance(loads_data, dict):
            print('ERROR: loads data must be passed as list of strings')
            raise

        if isinstance(bounds_data, list):
            print('processing boundaries')
            self.bounds_data = self.process_bounds_data(bounds_data)
            print('boundaries processed')
        elif not isinstance(bounds_data, dict):
            print('ERROR: bounds data must be passed as list of strings')
            raise

        if not os.path.exists(self.out):
            os.makedirs(self.out)

        self.ini = ini
        self.tref = tref
        print('parsing mdu...')
        with open(self.mdu, 'w') as mdu_file:
            with open(self.source_mdu, 'r') as template:
                lines = template.readlines()
                for line in lines:
                    if 'ExtForceFileNew' in line:
                        mdu_file.write('ExtForceFileNew = FlowFM_DFMWAQ_new.ext \n')
                    elif 'ExtForceFile' in line:
                        mdu_file.write('ExtForceFile = FlowFM_DFMWAQ.ext \n')
                    elif 'NetFile' in line:
                        self.migrate_model_asset(line, 'NetFile', mdu_file)
                    elif 'DryPointsFile' in line:
                        self.migrate_model_asset(line, 'DryPointsFile', mdu_file)
                    elif 'LandBoundaryFile' in line:
                        self.migrate_model_asset(line, 'LandBoundaryFile', mdu_file)
                    elif 'ThinDamFile' in line:
                        self.migrate_model_asset(line, 'ThinDamFile', mdu_file)
                    elif 'ObsFile' in line:
                        self.migrate_model_asset(line, 'ObsFile', mdu_file)
                    else:
                        mdu_file.write(line)

            # DFMWAQ lines
            mdu_file.write('\n')
            mdu_file.write('[processes]\n')
            mdu_file.write('SubstanceFile = %s \n' % change_os(subfile.path, self.run_sys))  #
            mdu_file.write('AdditionalHistoryOutputFile       =                                         \n')  #
            mdu_file.write('ProcesDataBaseFile                = proc_def.dat\n')
            mdu_file.write('DtProcesses                       = 600.                  # waq processes time step\n')
            mdu_file.write('DtMassBalance                     = 86400.  # waq mass balance output time step\n')
            mdu_file.write('ProcessFluxIntegration            = 1       # integration option (1: WAQ, 2: D-Flow FM)\n')
            mdu_file.write('\n')

            # process library
            try:
                sh.copyfile(r'p:\h6\opt\delft3d\delwaq\5.08.00.64083\x64\dwaq\default\proc_def.dat',
                            os.path.join(new_dir, 'proc_def.dat'))
                sh.copyfile(r'p:\h6\opt\delft3d\delwaq\5.08.00.64083\x64\dwaq\default\proc_def.def',
                            os.path.join(new_dir, 'proc_def.def'))
                sh.copyfile(r'p:\h6\opt\delft3d\delwaq\5.08.00.64083\x64\dwaq\default\bloom.spe', os.path.join(new_dir, 'bloom.spe'))
                sh.copyfile(r'p:\h6\opt\delft3d\delwaq\5.08.00.64083\x64\dwaq\default\bloominp.d09',
                            os.path.join(new_dir, 'bloominp.d09'))

            except FileNotFoundError:
                if process_path is not None:
                    sh.copyfile(os.path.join(process_path, 'proc_def.dat'),
                                os.path.join(new_dir, 'proc_def.dat'))
                    sh.copyfile(os.path.join(process_path, 'proc_def.def'),
                                os.path.join(new_dir, 'proc_def.def'))
                    try:
                        sh.copyfile(os.path.join(process_path, 'bloom.spe'),
                                    os.path.join(new_dir, 'bloom.spe'))
                        sh.copyfile(os.path.join(process_path, 'bloominp.d09'),
                                    os.path.join(new_dir, 'bloominp.d09'))
                    except FileNotFoundError:
                        print('WARNING: No BLOOM files found, only DYNAMO can be used')
                else:
                    print('using tool outside Deltares, please specify proc_def path')
                    raise
        print('mdu created')
        print('registering boundaries...')
        # register boundaries
        if isinstance(ext, list):
            self.boundaries = []
            for file in self.ext:
                self.boundaries.append(boundary_from_ext(file))
        else:
            print('ERROR: *.ext must be passed as a list of files')
            raise
        print('boundaries registered')
        print('initialization complete')

    def build(self):
        '''
        builds model based solely on attributes from initialization
        '''
        print('building model')
        self.soursin()
        self.open_bnd()
        self.administrate_ext_files()
        self.merge_ext_files()
        self.run_set_file()
        print('model built')

    def soursin(self):
        '''
        concatenates additional substances to a new local copy of a soursin file
        can be done from data tables, whereby the data will be interpolated and the new
        concentration calculated based on the flow in the tim file
        '''
        for boundaries in self.boundaries:  # each file
            for ind, name in enumerate(boundaries.keys()):
                for bndtype in boundaries[name].keys():
                    if 'discharge_salinity_temperature_sorsin' in bndtype:  # source sink tim files
                        if name in self.loads_data.keys():
                            fill = True
                        else:
                            fill = False
                        data = boundaries[name][bndtype]['data_loc']
                        with open(os.path.join(self.out, '%s' % os.path.split(data)[1]), 'w') as tim:
                            # should be one matching tim for each sours pli
                            # parse the original data file for flows, substances will be concatenated
                            with open(data, 'r') as flow:
                                arr = []
                                head = []
                                flows = flow.readlines()
                                for line in flows:
                                    if '*' in line:
                                        head.append(line)
                                    elif len(line) == 0:
                                        pass
                                    else:
                                        line = line.strip()
                                        ar = row2array(line)
                                        arr.append(ar)

                            # concatenation
                            arr = np.array(arr)
                            times = np.array([self.tref + datetime.timedelta(minutes=int(tt)) for tt in arr[:, 0]])
                            times = pd.to_datetime(times)
                            flow_vals = arr[:, 1]
                            # does not matter if sal and temp or not, will be concatenated
                            for sub in self.substances:
                                if self.subfile.transportable[sub] == 'active':
                                    # not an inactive substance
                                    if fill:
                                        # find interpolated value
                                        if sub in self.loads_data[name].keys():
                                            loads_times = pd.to_datetime(np.array(self.loads_data[name]['time']))
                                            loads_val = np.array(self.loads_data[name][sub])

                                            # g/s
                                            val = np.interp(times, loads_times, loads_val)
                                            # g/m3
                                            val = val / flow_vals
                                        else:
                                            val = np.zeros(len(times))
                                    else:
                                        val = np.zeros(len(times))

                                    val = np.array([val, ]).T
                                    arr = np.concatenate((arr, val), axis=1)

                            # header
                            for hh in head:
                                tim.write(hh)
                            for sind, sub in enumerate(self.substances):
                                if self.subfile.transportable[sub] == 'active':
                                    tim.write('* COLUMN%i=%s\n' % (sind + len(head), sub))
                            # data
                            for tind, time in enumerate(times):
                                minute_time = (time - self.tref).seconds / 60.0 + (time - self.tref).days * 1440.0
                                tim.write('%.4e    ' % minute_time)
                                for val in arr[tind, 1:]:
                                    tim.write('%.4e    ' % val)
                                tim.write('\n')

                            # copy the pli, including internal support point references
                            pli = boundaries[name][bndtype]['pli_loc']
                            # could be pliz or pli, so use orig name
                            with open(os.path.join(self.out, '%s' % os.path.split(pli)[1]), 'w') as pli_file:
                                with open(pli, 'r') as bndFile:
                                    lines = bndFile.readlines()
                                    for line in lines:
                                        pli_file.write(line)

    def open_bnd(self):
        '''
        creates a unique pli and tim file for each substance at each unique open boundary
        '''
        for boundaries in self.boundaries:
            for ind, name in enumerate(boundaries.keys()):
                for bndtype in boundaries[name].keys():
                    if 'waterlevelbnd' in bndtype or 'dischargebnd' in bndtype:
                        dummy_times = np.array([self.tref, self.tref + datetime.timedelta(days=20*365)])
                        if name in self.bounds_data.keys():
                            fill = True
                        else:
                            fill = False
                        # create a pli and tim for every substance
                        for sub in self.substances:
                            if self.subfile.transportable[sub] == 'active':
                                if fill:
                                    if sub in self.bounds_data[name].keys():
                                        times = pd.to_datetime(np.array(self.bounds_data[name]['time']))
                                        vals = np.array(self.bounds_data[name][sub])
                                    else:
                                        times = dummy_times
                                        vals = np.zeros(len(dummy_times))
                                else:
                                    times = dummy_times
                                    vals = np.zeros(len(dummy_times))

                                # pli
                                pli = boundaries[name][bndtype]['pli_loc']
                                file_name = os.path.split(pli)[1]
                                point = file_name.find('.')
                                with open(os.path.join(self.out, '%s%s%s' % (file_name[:point], sub, file_name[point:])),
                                          'w') as pli_file:
                                    # copy the existing pli
                                    with open(pli, 'r') as bndFile:
                                        lines = bndFile.readlines()
                                        for line in lines:
                                            pli_file.write(line.replace(name, name + sub))
                                # tim
                                with open(os.path.join(self.out, '%s%s_0001.tim' % (name, sub)), 'w') as tim:
                                    tim.write('* column 1 = time in minutes since \n')
                                    tim.write('* column 2 = concentration of  ' + sub + '\n')
                                    for tind, time in enumerate(times):
                                        min_time = (time - self.tref).seconds / 60.0 + (time - self.tref).days * 1440.0
                                        tim.write('%.4e    %.4e\n' % (min_time, vals[tind]))

                        # copy the original boundary to make folder self-contained model
                        # sh.copyfile(pli, self.out + '%s.pli' % (name))
                        # data = boundaries[name][bndtype]['data_loc']
                        # sh.copyfile(data, self.out + data[find_last(data, '\\'):])


    def migrate_model_asset(self, line, keyword, mdu_file):
        '''
        for certain lines in the mdu, parse the file location and ensure the new mdu can find it
        '''
        orig = os.getcwd()
        line = line.replace(keyword, '').replace('=', '').strip()
        if '#' in line:
            asset = line[:line.find('#')].strip()
        else:
            asset = line.strip()

        if '..' in asset:
            # is relative to old mdu path, need to navigate to where it is pointing
            os.chdir(os.path.split(self.source_mdu)[0])
            up = asset.count('..')
            for jump in range(0, up):
                os.chdir('..\\')
            new_root = os.getcwd()
            # replace all 'up dir' commands
            asset = asset.replace('..\\', '').replace('../', '')
            asset = change_os(asset)
            # copy to new directory
            sh.copyfile(os.path.join(new_root, asset), os.path.join(self.out, os.path.split(asset)[1]))
            # refer to new directory
            if 'net.nc' in asset:
                self.grid = os.path.join(self.out, os.path.split(asset)[1])
            os.chdir(orig)

        elif '\\' not in asset and len(asset) > 0:
            # is local
            # must obtain from original mdu directory, copy, and reference new local
            sh.copyfile(os.path.join(os.path.split(self.source_mdu)[0], asset), os.path.join(self.out, asset))
            if 'net.nc' in asset:
                self.grid = os.path.join(self.out, asset)

        elif len(asset) > 0:
            # is absolute, linux
            # copy from -> to
            sh.copyfile(asset, os.path.join(self.out, os.path.split(asset)[1]))
            if 'net.nc' in asset:
                self.grid = os.path.join(self.out, os.path.split(asset)[1])
        else:
            pass # no file
        # must write local reference
        mdu_file.write('%s                           = %s \n' % (keyword, asset))

    def administrate_ext_files(self):
        '''
        writes the ext files that relate to the newly written data, which are local tim files
        this is challenging because there could be two ext files
        uses self.run_sys to decide path types
        '''
        with open(os.path.join(self.out, 'FlowFM_DFMWAQ_source.ext'), 'w') as nmfs:  # new model file sours ins
            with open(os.path.join(self.out, 'FlowFM_DFMWAQ_open.ext'), 'w') as nmfo:  # new model file open bnds
                with open(os.path.join(self.out, 'FlowFM_DFMWAQ_new.ext'), 'w') as nmfn:  # new model file new
                    with open(os.path.join(self.out, 'FlowFM_DFMWAQ_extra.ext'), 'w') as nmfe:  # new model file extra
                        # write pre-existing, this will account for new concatenated tim files, as the name is unchanged
                        for boundaries in self.boundaries:
                            # for all of the boundary files
                            for ind, name in enumerate(boundaries.keys()):
                                for bndtype in boundaries[name].keys():
                                    # re-write the original boundary data
                                    # water level and discharge data will be written to the old ext, even though
                                    # they came from the new ext

                                    if 'waterlevelbnd' in bndtype or 'dischargebnd' in bndtype:
                                        # inflow from open boundaries
                                        # not copied to return path to linux
                                        nmfn.write('[boundary]\n')
                                        nmfn.write('quantity=%s\n' % boundaries[name][bndtype]['type'])
                                        nmfn.write('locationfile=%s\n' %
                                                   change_os(boundaries[name][bndtype]['pli_loc'], self.run_sys))
                                        nmfn.write('forcingfile=%s\n' %
                                                   change_os(boundaries[name][bndtype]['data_loc'], self.run_sys))
                                        nmfn.write('\n')

                                        for sub in self.substances:
                                            if self.subfile.transportable[sub] == 'active':
                                                nmfo.write('QUANTITY=tracerbnd%s\n' % sub)
                                                nmfo.write('FILENAME=%s%s.pli\n' % (name, sub))
                                                nmfo.write('FILETYPE=9\n')
                                                nmfo.write('METHOD=3\n')
                                                nmfo.write('OPERAND=O\n')
                                                nmfo.write('\n')

                                    elif 'discharge_salinity_temperature_sorsin' in bndtype:
                                        # simply rewrite the boundaries, as they are either
                                        # not relevant or have been silently edited/copied
                                        # reference the local copy
                                        nmfs.write('QUANTITY=%s\n' % bndtype)
                                        # refer to new local copy
                                        nmfs.write('FILENAME=%s\n' % boundaries[name][bndtype]['pli_loc'][
                                                                     find_last(boundaries[name][bndtype]['pli_loc'],
                                                                               '\\'):])
                                        nmfs.write('FILETYPE=%s\n' % boundaries[name][bndtype]['FILETYPE'])
                                        nmfs.write('METHOD=%s\n' % boundaries[name][bndtype]['METHOD'])
                                        nmfs.write('OPERAND=%s\n' % boundaries[name][bndtype]['OPERAND'])
                                        nmfs.write('\n')
                                    else:
                                        if 'FILETYPE' in boundaries[name][bndtype].keys():
                                            nmfe.write('QUANTITY=%s\n' % bndtype)
                                            # refer to original, not copied
                                            nmfe.write('FILENAME=%s\n' %
                                                       change_os(boundaries[name][bndtype]['pli_loc'], self.run_sys))
                                            nmfe.write('FILETYPE=%s\n' % boundaries[name][bndtype]['FILETYPE'])
                                            nmfe.write('METHOD=%s\n' % boundaries[name][bndtype]['METHOD'])
                                            nmfe.write('OPERAND=%s\n' % boundaries[name][bndtype]['OPERAND'])
                                            nmfe.write('\n')
                                        else:
                                            nmfn.write('[boundary]\n')
                                            nmfn.write('quantity=%s\n' % boundaries[name][bndtype]['type'])
                                            # refer to original, not copied
                                            nmfn.write('locationfile=%s\n' %
                                                       change_os(boundaries[name][bndtype]['pli_loc'], self.run_sys))
                                            nmfn.write('forcingfile=%s\n'
                                                       % change_os(boundaries[name][bndtype]['data_loc'], self.run_sys))
                                            nmfn.write('\n')

                    # initials polygon
                    varnames = nc_format(self.grid)
                    grd = netCDF4.Dataset(self.grid)

                    x_min = np.min(grd.variables[varnames['xnode']][:])
                    x_max = np.max(grd.variables[varnames['xnode']][:])
                    y_min = np.min(grd.variables[varnames['ynode']][:])
                    y_max = np.max(grd.variables[varnames['ynode']][:])

                    with open(os.path.join(self.out, 'domain.pol'), 'w') as pol:
                        pol.write('domain\n')
                        pol.write('4   2\n')
                        pol.write('%.4e    %.4e\n' % (x_min, y_min))
                        pol.write('%.4e    %.4e\n' % (x_min, y_max))
                        pol.write('%.4e    %.4e\n' % (x_max, y_max))
                        pol.write('%.4e    %.4e\n' % (x_max, y_min))

                    for sub in self.substances:
                        if self.subfile.transportable[sub] == 'active':
                            nmfo.write('QUANTITY=initialtracer%s\n' % sub)
                        else:
                            nmfo.write('QUANTITY=initialwaqbot%s\n' % sub)
                        nmfo.write('FILENAME=domain.pol\n')
                        nmfo.write('FILETYPE=10\n')
                        nmfo.write('METHOD=4\n')
                        nmfo.write('OPERAND=O\n')
                        if sub in self.ini.keys():
                            nmfo.write('VALUE=%.4e\n' % self.ini[sub])
                        else:
                            nmfo.write('VALUE=0.0\n')
                        nmfo.write('\n')


    def merge_ext_files(self):
        '''
        merge temporary partial ext files into master file
        '''
        with open(os.path.join(self.out, 'FlowFM_DFMWAQ.ext'), 'w') as nmff:  # new model file final
            with open(os.path.join(self.out, 'FlowFM_DFMWAQ_extra.ext'), 'r') as cpfile:
                lines = cpfile.readlines()
                for line in lines:
                    nmff.write(line)
            with open(os.path.join(self.out, 'FlowFM_DFMWAQ_source.ext'), 'r') as cpfile:
                lines = cpfile.readlines()
                for line in lines:
                    nmff.write(line)
            with open(os.path.join(self.out, 'FlowFM_DFMWAQ_open.ext'), 'r') as cpfile:
                lines = cpfile.readlines()
                for line in lines:
                    nmff.write(line)
        # clean files
        os.remove(os.path.join(self.out, 'FlowFM_DFMWAQ_extra.ext'))
        os.remove(os.path.join(self.out, 'FlowFM_DFMWAQ_source.ext'))
        os.remove(os.path.join(self.out, 'FlowFM_DFMWAQ_open.ext'))


    def run_set_file(self):
        '''
        write the shell script to run the model
        '''
        if self.nodes is not None and self.version is not None:
            with open(self.mdu[:find_last(self.mdu, '\\')] + 'run_set.sh', 'w') as shfile:
                shfile.write('#!/bin/bash\n')
                shfile.write('module load dflowfm\n')
                shfile.write('run_dflowfm.sh -v %s --partition:ndomains=%i:icgsolver=6 %s\n' %
                             (self.version, self.nodes * self.threads, self.mdu[find_last(self.mdu, '\\'):]))
                shfile.write('#Start qsub \n')
                shfile.write('submit_dflowfm.sh -v %s -m %s -n %i -c %i --processlibrary proc_def.dat -j %s\n' %
                             (self.version, self.mdu[find_last(self.mdu, '\\'):], self.nodes, self.threads,
                              self.mdu[find_last(self.mdu, '\\'):]))


    def process_loads_data(self, files):
        '''
        reads in a deltashell format data file for creation of tim files
        needs to be processed as mass per time, so it can be interpolated, 
        and converted to a concentration based on the flow in the tim file

        Arguments:
            file {str} -- path to deltashel .csv file
        '''
        for file in files:
            df = pd.read_csv(file)
            ret = {}    
            for lind, loc in enumerate(df['location'].unique()):
                print('processing ' + loc)
                # will be unique tim
                dfl = df[df['location'] == loc]
                subs = dfl['substance'].unique()

                if 'FLOW' not in subs and 'flow' not in subs:
                    print('ERROR: flow rate needed for each bound')
                    raise
                else:
                    if 'FLOW' in subs:
                        flowname = 'FLOW'
                    else:
                        flowname = 'flow'

                    subs = [ii for ii in subs if ii != 'FLOW' and ii != 'flow']
                    
                    piv = dfl.pivot_table(index = 'timeBlock', columns = 'substance', values = 'value')
                    piv['time'] = piv.index.values
                    piv['dtime'] = pd.to_datetime(piv['time'])   
                    
                    piv.set_index('dtime', inplace = True, drop = False)
                    piv.sort_index(inplace = True)
                    ret_df = pd.DataFrame(columns = subs, index = piv['dtime'])

                    ret_df['time'] = piv['dtime']
                    for sub in subs:
                    # g/s
                        if lind == 0:
                            piv['OXY']
                        ret_df[sub] = piv[sub] * piv[flowname]

                    ret[loc] = ret_df

        return ret


    def process_bounds_data(self, files):
        '''
        reads in a deltashell format data file for creation of tim files
        needs to be processed as loads per time, so it can be interpolated, 
        and converted to a concentration based on the flow in the tim file

        Arguments:
            file {str} -- path to deltashel .csv file
        '''
        for file in files:
            df = pd.read_csv(file)
            ret = {}
            for loc in df['location'].unique():
                print('processing ' + loc)
                # will be unique tim
                dfl = df[df['location'] == loc]
                subs = dfl['substance'].unique()

                piv = dfl.pivot_table(index='timeBlock', columns='substance', values='value')
                piv['time'] = piv.index.values
                piv['dtime'] = pd.to_datetime(piv['time'])

                piv.set_index('time', inplace=True, drop=False)
                piv.sort_index(inplace=True)
                ret_df = pd.DataFrame(columns=subs, index=piv['time'])
                ret_df['time'] = piv['dtime']
                for sub in subs:
                    # g/m3
                    ret_df[sub] = piv[sub]

                ret[loc] = ret_df

        return ret


class Grevelingen(DFMWAQModel):
    # example subclass
    def __init__(self):
        # mdu to copy
        mdu = r'p:\11203715-006-d-hydro-grevelingen\communicatie\201908XX_verzonden_aan_RWS\model\2008\computations\run01\Grevelingen-FM_save.mdu'
        # boundaries for that mdu
        ext = [r'p:\11203715-006-d-hydro-grevelingen\WAQ\DFMWAQ\model_2008\computations\run01\Grevelingen-FM_bnd.ext']
        # location of new model
        new_bnd_dir = 'p:\\11203715-006-d-hydro-grevelingen\\WAQ\\DFMWAQ\\model_2008\\computations\\run01\\'
        # sub file to use
        subfile = SubFile(
            r'p:\11203715-006-d-hydro-grevelingen\WAQ\DFMWAQ\model_2008\computations\run01\Marine_Algae_parsGrev.sub')
        # initial conditions
        ini = {'OXY': 7}
        # kernel version
        v = '1.2.59.64457'
        # cores to run on h6
        cores = [2, 4]
        super().__init__(mdu, ext, subfile, new_bnd_dir, ini, v, cores)


class Guayaquil(DFMWAQModel):
    # example subclass
    def __init__(self):
        mdu = r'p:\11201302-guayaquil\02_flow\03_baseCase\R26\guayas.mdu'
        ext = [r'p:\11201302-guayaquil\02_flow\03_baseCase\R26\plant_loads_current_local.ext',
               r'p:\11201302-guayaquil\02_flow\03_baseCase\R26\sea_riv_boundary_local_bc.ext']
        new_bnd_dir = 'p:\\11201302-guayaquil\\02_flow\\03_baseCase\\DFMWAQ\\'
        subfile = SubFile(r'p:\11201302-guayaquil\03_waterquality\03_baseCase\01_substances\guayas_V11.sub')
        ini = {'OXY': 7}
        v = '1.2.59.64457'
        cores = [2, 4]
        super().__init__(mdu, ext, subfile, new_bnd_dir, ini, v, cores)
