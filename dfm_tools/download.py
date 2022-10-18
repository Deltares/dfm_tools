# -*- coding: utf-8 -*-
"""
Created on Tue Oct 18 15:09:26 2022

@author: veenstra
"""

import subprocess
import pandas as pd
import datetime as dt
from pathlib import Path

def download_CMEMS(username, password, #register at: https://resources.marine.copernicus.eu/registration-form' 
                   dir_output='.', #default to pwd
                   longitude_min=-180, longitude_max=180, latitude_min=-90, latitude_max=90,
                   date_min='2010-01-01', date_max='2010-01-03', #'%Y-%m-%d'
                   varlist=['bottomT'], #['thetao','so','zos','bottomT','uo','vo'], ['o2','no3','po4','si','nppv','chl'],
                   reanalysis_forecast=None, quantitygroup=None,
                   motu_url=None, service=None, product=None, #optionally provided with motu_url_dict and source_dict
                   timeout=30, #in seconds #TODO: set timeout back to 300?
                   max_tries=2):
    
    import motuclient #used in motu_commands, so has to be importable. conda install -c conda-forge motuclient #TODO: move to top of script (then make dependency of dfm_tools)
    
    motu_url_dict = {'reanalysis':'http://my.cmems-du.eu/motu-web/Motu',# multiyear reanalysis data (01-01-1993 12:00 till 31-05-2020 12:00)
                     'forecast':'http://nrt.cmems-du.eu/motu-web/Motu'} # operational forecast data (01-01-2019 12:00 till now + several days)
                    #http://my.cmems-du.eu2/motu-web/Motu # invalid url, TimeoutExpired + othererror OUT: [WARNING] Warning: CAS connection failed
                    #http://my.cmems-du.eu/motu-web/Motu2 # invalid url, CalledProcessError + othererror OUT: [ERROR] Execution failed: HTTP Error 404
                    #http://nrt.cmems-du.eu/motu-web/Motu # nrt instead of my, othererror OUT: [ERROR] 010-30 : The requested service is unknown: 'GLOBAL_MULTIYEAR_PHY_001_030-TDS'
    source_dict =   {'reanalysis':{'physchem':{'service': 'GLOBAL_MULTIYEAR_PHY_001_030-TDS',
                                               'product': 'cmems_mod_glo_phy_my_0.083_P1D-m'},
                                   'bio':     {'service': 'GLOBAL_MULTIYEAR_BGC_001_029-TDS',
                                               'product': 'cmems_mod_glo_bgc_my_0.25_P1D-m'}},
                     'forecast':  {'physchem':{'service': 'GLOBAL_ANALYSIS_FORECAST_PHY_001_024-TDS',
                                               'product': 'global-analysis-forecast-phy-001-024'},
                                   'bio':     {'service': 'GLOBAL_ANALYSIS_FORECAST_BIO_001_028-TDS',
                                               'product': 'global-analysis-forecast-bio-001-028-daily'}}
                     }
    
    if motu_url is None:
        motu_url = motu_url_dict[reanalysis_forecast]
    if service is None:
        service = source_dict[reanalysis_forecast][quantitygroup]['service']
    if product is None:
        product = source_dict[reanalysis_forecast][quantitygroup]['product']
    
    date_range = pd.date_range(dt.datetime.strptime(date_min, '%Y-%m-%d'),dt.datetime.strptime(date_max, '%Y-%m-%d'), freq='D')
    
    for var in varlist:
        for date in date_range: #retrieve data per day
            date_str = date.strftime('%Y-%m-%d')
            name_output = f'cmems_{var}_{date_str}.nc' #TODO: add 12h to filename
            check_file = Path(dir_output,name_output)
            tryno = 0
            while not check_file.is_file():
                tryno += 1
                print(f'retrieving variable {var} for {date_str}: try {tryno}')
                
                motu_command = ' '.join(['motuclient', '--motu', motu_url, '--service-id', service, '--product-id', product,
                                         '--longitude-min', str(longitude_min), '--longitude-max', str(longitude_max),
                                         '--latitude-min', str(latitude_min), '--latitude-max', str(latitude_max),
                                         '--date-min', date_str, '--date-max', date_str,#+' 12:00:00',
                                         '--depth-min', '0', '--depth-max', '2e31',
                                         '--variable', str(var),
                                         '--out-dir', dir_output, '--out-name', name_output,
                                         '--user', username, '--pwd', password])
                try:
                    out = subprocess.run(f'python -m {motu_command}', capture_output=True, check=True, universal_newlines=True, timeout=timeout)
                    #stdin=None, input=None, stdout=None, stderr=None, 
                    #capture_output=False, shell=False, cwd=None, 
                    #encoding=None, errors=None, text=None, env=None, universal_newlines=None
                except subprocess.TimeoutExpired as e:
                    out = e
                    print(f'TimeoutExpired: {e} Check above logging.')
                    if tryno >= max_tries:
                        raise Exception(f'TimeoutExpired {max_tries} times: maybe increase max_tries ({max_tries}) or timeout ({timeout} seconds)')
                except subprocess.CalledProcessError as e:
                    out = e
                    raise Exception(f'CalledProcessError: {e} Check above logging.')
                finally:
                    if ('ERROR' in out.stdout) or ('WARNING' in out.stdout): #catch all other errors, and the relevant information in TimeoutExpired and CalledProcessError
                        raise Exception(f'othererror:\nOUT:{out.stdout}\nERR:{out.stderr}')
                    #else:
                    #    print(f'OUT:{out.stdout}\nERR:{out.stderr}')
                
                if tryno >= max_tries:
                    raise Exception(f'max tries ({max_tries}) reached, this should not happen since it is already catched at timeout and errors are raised at others')
            if check_file.is_file():
                print(f'file available after {tryno} tries, continuing to next var/time')
    print('done')
    return