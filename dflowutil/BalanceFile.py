from .utils import row2array
import numpy as np
import pandas as pd

class BalanceFile():
    '''
    DFMWAQ bal.txt file
    usage:
    file = r'd:\HK-FMWAQ_0000_wq_proc_bal.txt'
    prn = dflowutil.PrnFile(file)
    prn.extract_balances(['BalArea1'],['365-730'])
    df = prn.areas['BalArea1']['OXY']['Inflows']
    '''
    def __init__(self, path):
        self.header_length = 10000

        self.substances = []
        self.areas = {}
        self.fluxes = []
        self.time_ind = []
        self.unique_balance_ind = []
        self.rep_sub_ind = {}
        self.balance_sources = {}

        self.path = path
        print('loading file...')
        with open(self.path, 'r') as prn:
            self.lines = prn.readlines()
        print('file loaded')

        # metadata
        print('parsing file...')

        lines = self.lines[0:self.header_length]
        for ind, line in enumerate(lines):
            if 'Overview of mass balance areas' in line:
                self.area_ind = ind
            if 'Overview of substances' in line:
                self.substances_ind = ind
            if 'Overview of fluxes' in line:
                self.flux_ind = ind
            if 'Mass balances for' in line:
                self.balance_ind = ind

        # number of times and indices of unique balances (for a given substance, time, and area)
        for ind, line in enumerate(self.lines):
            if 'Mass balances for' in line:
                self.time_ind.append(ind)
            if 'Mass balance area             : ' in line:
                self.unique_balance_ind.append(ind)

        self.time_ind = np.array(self.time_ind, dtype=np.int)

        self.get_substances()
        self.get_areas()
        self.get_fluxes()
        self.get_representative_sub_ind()
        self.get_balance_sources()
        self.initialize_df()
        print('file parsed')

    def get_areas(self):
        '''
        discovers all unique areas in file
        '''
        lines = self.lines[0:self.header_length]
        for ind in range(self.area_ind + 3, self.substances_ind):
            line = lines[ind].strip().replace('\n', '').split(' ')
            if len(line) > 1:
                self.areas[line[-1]] = {}

    def get_substances(self):
        '''
        discovers all unique substances in file
        '''
        lines = self.lines[0:self.header_length]
        for ind in range(self.substances_ind + 4, self.flux_ind):
            line = lines[ind].strip().replace('\n', '').split(' ')
            if len(line) > 1 and ':' not in line:
                self.substances.append(line[-1])

    def get_fluxes(self):
        '''
        discovers all unique fluxes in file
        '''
        lines = self.lines[0:self.header_length]
        for ind in range(self.flux_ind + 4, self.balance_ind):
            line = lines[ind].strip().replace('\n', '').split(' ')
            if len(line) > 1:
                self.fluxes.append(line[-2])

    def get_representative_sub_ind(self):
        '''
        gets a representative index of a balance that can be used to determine a substance's balance sources
        A.K.A the first balance output for each substance
        '''
        for sub in self.substances:
            for ind in self.unique_balance_ind:
                curr_sub = self.lines[ind + 2].split(' ')[2].strip()
                if curr_sub not in self.rep_sub_ind.keys():
                    self.rep_sub_ind[curr_sub] = ind
                    self.balance_sources[curr_sub] = []

    def get_balance_sources(self):
        '''
        discovers all unique mass sources and sinks in file, serves as column names in dataframes
        unique to each substance
        '''
        go_subs = list(self.rep_sub_ind.keys())
        for sind, sub in enumerate(go_subs):
            start_ind = self.rep_sub_ind[go_subs[sind]]
            try:
                end_ind = self.rep_sub_ind[go_subs[sind + 1]]
            except:
                # is last, add length between last two
                end_ind = self.rep_sub_ind[go_subs[sind]] + self.rep_sub_ind[go_subs[sind]] - self.rep_sub_ind[
                    go_subs[sind - 1]]

            work_array = self.lines[start_ind + 12: end_ind - 11]
            for line in work_array:
                vals = line.split('  ')
                self.balance_sources[sub].append(vals[0])

    def initialize_df(self):
        '''
        initializes the dataframes containing the balance information
        each substance needs a different column structure
        '''
        for area in self.areas:
            for sub in self.substances:
                self.areas[area][sub] = {}
                self.areas[area][sub]['Inflows'] = pd.DataFrame(columns=self.balance_sources[sub])
                self.areas[area][sub]['Outflows'] = pd.DataFrame(columns=self.balance_sources[sub])

    def extract_balances(self, areas, times):
        """
        appends the data for this period, substance, and area to dataframe within the areas dictionary

        Arguments:
            areas {list} -- list of areas to extract data for
            times {list} -- list of times to extract data for, based on 1, days from start date
                            for example ['365-366','365-730'] represents two balances for the period
                            day 365 to day 366 from reference date, and day 365 to day 730 from reference date
        """

        # balance index is the row number in the file that every new balance (unique time, sub, and are) begins
        for bal_no, bal_ind in enumerate(self.unique_balance_ind):
            if 'Whole model' not in self.lines[bal_ind]:
                area = self.lines[bal_ind].strip().replace('\n', '').split(' ')[-1]
                t1 = self.lines[bal_ind - 3].split(':')[1].strip().replace('d', '').replace(' 0', '').strip()
                t2 = self.lines[bal_ind - 2].split(':')[1].strip().replace('d', '').replace(' 0', '').strip()
                if area in areas and (t1 + '-' + t2) in times:
                    sub = self.lines[bal_ind + 2].replace('Mass substance', '').replace('Begin', '').replace('End',
                                                                                                             '').replace(
                        '\n', '').strip()
                    print('%s balance for area %s and time range %s to %s' % (sub, area, t1, t2))

                    df_in = pd.DataFrame(columns=self.balance_sources[sub], index=[t1 + '-' + t2])
                    df_out = pd.DataFrame(columns=self.balance_sources[sub], index=[t1 + '-' + t2])

                    # for the chunk of this file where the current balance resides
                    params = []
                    for line_no in range(self.unique_balance_ind[bal_no] + 12,
                                         self.unique_balance_ind[bal_no + 1] - 11):
                        param = self.lines[line_no].split('  ')[0].strip()
                        params.append(param)
                        val = row2array(self.lines[line_no].replace('\n', ''))
                        try:
                            val_in = float(val[0])
                            val_out = float(val[1])
                        except:
                            # NOTE: Change after balance file bug is fixed!
                            val_in = 0.0
                            val_out = 0.0

                        df_in[param] = val_in
                        df_out[param] = val_out

                    self.areas[area][sub]['Inflows'] = pd.concat([self.areas[area][sub]['Inflows'], df_in], axis=0,
                                                                 ignore_index=False)
                    self.areas[area][sub]['Outflows'] = pd.concat([self.areas[area][sub]['Outflows'], df_out], axis=0,
                                                                  ignore_index=False)
