# -*- coding: utf-8 -*-
"""
Created on Mon Oct 17 16:55:46 2022

@author: veenstra
"""

#TODO: source: https://github.com/c-scale-community/use-case-hisea/blob/main/scripts/download/download_cmems_physics.py

import os
import click
from datetime import timedelta, datetime
from pathlib import Path
import subprocess
import xarray as xr # note dependencies: dask, netCDF4

@click.command()
@click.option(    '--username',
        default='',
        help='register at: https://resources.marine.copernicus.eu/registration-form')
@click.option(    '--password',
        default='',
        help='register at: https://resources.marine.copernicus.eu/registration-form')
@click.option(    '--longitude_min',
        type=(float),
        help='Minimum longitude for region of interest',
        default=-180,
        show_default=True)
@click.option(    '--longitude_max',
        type=(float),
        help='Maximum longitude for region of interest',
        default=180,
        show_default=True)
@click.option(    '--latitude_min',
        type=(float),
        help='Minimum latitude for region of interest',
        default=-90,
        show_default=True)
@click.option(    '--latitude_max',
        type=(float),
        help='Maximum latitude for region of interest',
        default=90,
        show_default=True)
@click.option(    '--date_min',
        type=(str),
        help='Start date for data download. Format: YYYY-MM-DD',
        default=(datetime.now()).strftime('%Y-%m-%d'),
        show_default=True)
@click.option(    '--date_max',
        type=(str),
        help='End date for data download. Format: YYYY-MM-DD',
        default=(datetime.now()).strftime('%Y-%m-%d'),
        show_default=True)
@click.option(    '--vars',
        multiple=True,
        help='available vars: https://catalogue.marine.copernicus.eu/documents/PUM/CMEMS-GLO-PUM-001-024.pdf',
        default=(    'thetao',
                'bottomT',
                'so',
                'zos',
                'uo',
                'vo'),
        show_default=True)

def runcommand(username, password, longitude_min, longitude_max, latitude_min, latitude_max, date_min, date_max, vars):
    #make the /data/tmp directory if it does not exist
    dir_output = './cmems_temp'
    Path(dir_output).mkdir(parents=True, exist_ok=True)
    delta = datetime.strptime(date_max, '%Y-%m-%d') - datetime.strptime(date_min, '%Y-%m-%d')
    for var in vars:
        for i in range(delta.days+1):
            max_runs = 2
            run = 0
            day = datetime.strptime(date_min, '%Y-%m-%d').date() + timedelta(days=i)
            check_file = Path(dir_output,str(var)+'_'+str(day)+'.nc')
            while not check_file.is_file():
                while run < max_runs:
                    try:
                        subprocess.run(['python', '-m', 'motuclient',
                            '--motu', 'https://nrt.cmems-du.eu/motu-web/Motu',
                            '--service-id', 'GLOBAL_ANALYSIS_FORECAST_PHY_001_024-TDS',
                            '--product-id', 'global-analysis-forecast-phy-001-024',
                            '--longitude-min',str(longitude_min),
                            '--longitude-max', str(longitude_max),
                            '--latitude-min', str(latitude_min),
                            '--latitude-max', str(latitude_max),
                            '--date-min', str(day)+' 12:00:00',
                            '--date-max', str(day)+' 12:00:00',
                            '--depth-min', '0.493',
                            '--depth-max', '5727.918000000001',
                            '--variable', str(var),
                            '--out-dir', dir_output,
                            '--out-name', 'cmems_'+str(var)+'_'+str(day)+'.nc',
                            '--user', username,
                            '--pwd', password],
                            check=True,
                            timeout=300)
                    except subprocess.TimeoutExpired as e:
                        print(var)
                        print(e.stdout)
                        print(e.stderr)
                        continue
                    else:
                        break
                    finally:
                        run += 1
        ds = xr.open_mfdataset(os.path.join(dir_output,f'{var}_*.nc'), combine='by_coords', decode_times=False)
        ds.to_netcdf(os.path.join(dir_output,f'{var}.nc'))


if __name__ == '__main__':
    runcommand()