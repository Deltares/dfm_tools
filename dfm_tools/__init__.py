"""
.. include:: ../README.md
"""

__author__ = """Jelmer Veenstra"""
__email__ = 'jelmer.veenstra@deltares.nl'
__version__ = '0.8.13'

from dfm_tools.download import *
from dfm_tools.get_nc import *
from dfm_tools.get_nc_helpers import *
from dfm_tools.hydrolib_helpers import *
from dfm_tools.interpolate_grid2bnd import *
from dfm_tools.linebuilder import *
from dfm_tools.modplot import *
from dfm_tools.regulargrid import *
from dfm_tools.ugrid import *
from dfm_tools.xarray_helpers import *
from dfm_tools.energy_dissipation import *

import warnings
warnings.filterwarnings('always',category=DeprecationWarning)

#add plotmethod to xugrid grid object
import xugrid as xr
def plot(self,**kwargs): #TODO: maybe add to xugrid directly
    xr.plot.line(self,**kwargs) #uds.ugrid.grid
    pass
ug2d = xr.ugrid.ugrid2d.Ugrid2d
ug2d.plot = plot
