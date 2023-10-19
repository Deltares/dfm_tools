"""
.. include:: ../README.md
"""

__author__ = """Jelmer Veenstra"""
__email__ = "jelmer.veenstra@deltares.nl"
__version__ = "0.15.0"

# TODO workaround for hcdfm calling of deprecated enum: https://github.com/Deltares/HYDROLIB-core/blob/bf75f024e48b515d4cd7c93cd3f997de4fc5b1ef/hydrolib/core/dflowfm/net/models.py#L1185
import meshkernel as mk
mk.DeleteMeshOption.ALL_FACE_CIRCUMCENTERS = 999

from dfm_tools.deprecated import *
from dfm_tools.errors import *
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
warnings.filterwarnings("always",category=DeprecationWarning)
