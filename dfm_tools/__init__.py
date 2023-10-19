"""
.. include:: ../README.md
"""

__author__ = """Jelmer Veenstra"""
__email__ = "jelmer.veenstra@deltares.nl"
__version__ = "0.15.0"

# TODO: workaround for hcdfm calling of deprecated enum: https://github.com/Deltares/HYDROLIB-core/blob/bf75f024e48b515d4cd7c93cd3f997de4fc5b1ef/hydrolib/core/dflowfm/net/models.py#L1185
import meshkernel as mk
mk.DeleteMeshOption.ALL_FACE_CIRCUMCENTERS = 999

# TODO: workaround for hcdfm.FMModel() calling Network/Networkmodel with mk.MeshKernel(is_geographic=is_geographic)
import hydrolib.core.dflowfm as hcdfm
from hydrolib.core.basemodel import ParsableFileModel
from pydantic import Field
def Network_init(self, is_geographic: bool = False) -> None:
    self.meshkernel = mk.MeshKernel(projection=is_geographic)
    # Monkeypatch the meshkernel object, because the "is_geographic" is not saved
    # otherwise, and needed for reinitializing the meshkernel
    self.meshkernel.is_geographic = is_geographic
    self._mesh1d = hcdfm.Mesh1d(meshkernel=self.meshkernel)
    self._mesh2d = hcdfm.Mesh2d(meshkernel=self.meshkernel)
    self._link1d2d = hcdfm.Link1d2d(meshkernel=self.meshkernel)
hcdfm.net.models.Network.__init__ = Network_init
class NetworkModel(ParsableFileModel):
    """Network model representation."""
    network: hcdfm.Network = Field(default_factory=hcdfm.Network)
hcdfm.net.models.NetworkModel = NetworkModel

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
