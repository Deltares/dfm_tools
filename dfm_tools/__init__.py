"""
.. include:: ../README.md
"""

__author__ = """Jelmer Veenstra"""
__email__ = "jelmer.veenstra@deltares.nl"
__version__ = "0.21.0"

from dfm_tools.deprecated import *
from dfm_tools.download import *
from dfm_tools.get_nc import *
from dfm_tools.get_nc_helpers import *
from dfm_tools.hydrolib_helpers import *
from dfm_tools.meshkernel_helpers import *
from dfm_tools.interpolate_grid2bnd import *
from dfm_tools.linebuilder import *
from dfm_tools.modplot import *
from dfm_tools.xarray_helpers import *
from dfm_tools.xugrid_helpers import *
from dfm_tools.energy_dissipation import *
from dfm_tools.bathymetry import *
from dfm_tools.coastlines import *
from dfm_tools import data
from dfm_tools.modelbuilder import *
from dfm_tools.observations import *

import warnings
warnings.filterwarnings(action="always", category=DeprecationWarning)

# suppress chunks warning: https://github.com/Deltares/dfm_tools/issues/679
warnings.filterwarnings("ignore", category=UserWarning, message='The specified chunks separate the stored chunks along dimension "time" starting at index 1. This could degrade performance. Instead, consider rechunking after loading.')
