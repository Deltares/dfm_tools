"""
.. include:: ../README.md
"""

__author__ = """Jelmer Veenstra"""
__email__ = "jelmer.veenstra@deltares.nl"
__version__ = "0.16.2"

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

import warnings
warnings.filterwarnings(action="always", category=DeprecationWarning)

# python 3.8 DeprecationWarning
import sys
version_info = sys.version_info
if (version_info.major == 3) & (version_info.minor <= 8):
    warnings.warn(
        DeprecationWarning("You are using Python 3.8 or lower, please note that "
                           "dfm_tools will soon drop support for Python 3.8: "
                           "https://github.com/Deltares/dfm_tools/issues/267")
        )
