# -*- coding: utf-8 -*-
"""
Created on Thu Jul  3 14:23:06 2025

@author: veenstra
"""

import os
from pathlib import Path

if os.name == "nt":
    # windows drive letter should include trailing slash
    # https://github.com/Deltares/dfm_tools/issues/1084
    PDRIVE = "p:/"
else:
    PDRIVE = "/p"


# the paths below can be overwritten with the following code in the user script
# import dfm_tools as dfmt; dfmt.settings.PATH_FES2014 = "updated/path/to/source"
# however, there is no validation whatsoever, so beware of typos

# tide models
PATH_FES2012 = Path(PDRIVE, 'metocean-data', 'open', 'FES2012', 'data','*_FES2012_SLEV.nc') # is probably also licensed
PATH_FES2014 = Path(PDRIVE, 'metocean-data', 'licensed', 'FES2014', '*.nc') # ocean_tide_extrapolated
PATH_EOT20 = Path(PDRIVE, 'metocean-data', 'open', 'EOT20', 'ocean_tides','*_ocean_eot20.nc')
PATH_GTSMv41 = Path(PDRIVE, '1230882-emodnet_hrsm', 'GTSMv3.0EMODnet', 'EMOD_MichaelTUM_yearcomponents', 'GTSMv4.1_yeartide_2014_2.20.06', 'compare_fouhis_fouxyz_v4', 'GTSMv4.1_tide_2014_*_rasterized.nc')
PATH_GTSMv41_opendap = 'https://opendap.deltares.nl/thredds/dodsC/opendap/deltares/GTSM/GTSMv4.1_tide/GTSMv4.1_tide_2014_*_rasterized.nc'
PATH_tpxo80_opendap = 'https://opendap.deltares.nl/thredds/dodsC/opendap/deltares/delftdashboard/tidemodels/tpxo80/tpxo80.nc'

# observation data
PATH_GESLA3 = Path(PDRIVE, "metocean-data", "licensed", "GESLA3")
