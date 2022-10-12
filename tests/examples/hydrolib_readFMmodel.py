# -*- coding: utf-8 -*-
"""
Created on Mon Oct  3 12:07:18 2022

@author: veenstra
"""

from pathlib import Path
from hydrolib.core.io.mdu.models import FMModel, NetworkModel, ExtModel, StructureModel
from hydrolib.core.io.bc.models import ForcingModel
from dfm_tools.hydrolib_helpers import forcingobject_to_dataframe
import datetime as dt
import matplotlib.pyplot as plt
plt.close('all')

dtstart = dt.datetime.now()
"""
file_mdu = Path(r'p:\11206813-006-kpp2021_rmm-2d\C_Work\31_RMM_FMmodel\computations\model_setup\run_206_HYDROLIB\RMM_dflowfm.mdu') #model with all but one structure and all but one lateral commented, reduces validation errors from >200 to 5. TODO: resolve validation errors
#file_mdu = Path(r'c:\DATA\dfm_tools_testdata\DFM_3D_z_Grevelingen\computations\run01\Grevelingen-FM.mdu')
#fm = FMModel(file_mdu) #TODO: currently crashes on issues below, and is quite slow since all files are being read


file_struct = Path(r'p:\11206813-006-kpp2021_rmm-2d\C_Work\31_RMM_FMmodel\computations\model_setup\run_206_HYDROLIB\RMM_structures.ini')
#structs = StructureModel(file_struct)
#TODO: pli in structures.ini is currently not supported: https://github.com/Deltares/HYDROLIB-core/issues/353 (use *_original file to test after fix)
#TODO: single structure in structures.ini currently crashes because of missing make_list_validator: https://github.com/Deltares/HYDROLIB-core/pull/352 (use *_original file to test after fix)


file_network = Path(r'p:\11206813-006-kpp2021_rmm-2d\C_Work\31_RMM_FMmodel\computations\model_setup\run_206_HYDROLIB\rmm_v1p7_net.nc')
#network = NetworkModel(file_network) #TODO: what is this used for?
"""

#file_extnew = Path(r'p:\11206813-006-kpp2021_rmm-2d\C_Work\31_RMM_FMmodel\computations\model_setup\run_206_HYDROLIB\RMM_bnd_5bnds.ext')
file_extnew = Path(r'p:\11206813-006-kpp2021_rmm-2d\C_Work\31_RMM_FMmodel\computations\model_setup\run_206_HYDROLIB\RMM_bnd_course.ext')
#ext = ExtModel(fm.external_forcing.extforcefilenew) #TODO: also possible to read from FMmodel?
ext = ExtModel(file_extnew) #TODO: laterals xycoordinates float is not yet supported (int is prescribed in ext model): https://github.com/Deltares/HYDROLIB-core/pull/351 (use *_original file to test after fix)

time_passed = (dt.datetime.now()-dtstart).total_seconds()
print(f'>>time passed: {time_passed:.2f} sec')

max_extforcings = 6 #None for all?

ext_boundaries = ext.boundary
ext_laterals = ext.lateral
for iEB, extbnd in enumerate(ext_boundaries+ext_laterals): #TODO: waterlevelbnd for rivers are present three times: https://github.com/Deltares/HYDROLIB-core/issues/354
    if hasattr(extbnd,'forcingfile'):
        extbnd_filepath = extbnd.forcingfile.filepath
        extbnd_forcings = extbnd.forcingfile.forcing
    elif hasattr(extbnd,'discharge'):
        extbnd_filepath = extbnd.discharge.filepath
        extbnd_forcings = extbnd.discharge.forcing
    else:
        raise Exception('ERROR: not forcingfile/discharge present (boundary/lateral') #TODO: this is not intuitive, make issue?
    print(f'boundary {iEB+1} of {len(ext_boundaries)}: {extbnd_filepath}')
    fig,ax = plt.subplots(figsize=(12,6))
    leglabels_new = []
    for iEBF, forcing in enumerate(extbnd_forcings[:max_extforcings]):
        print(f'forcing {iEBF+1} of {len(extbnd_forcings)}: {forcing.name} ({forcing.function}) ({forcing.quantityunitpair[1].quantity})')
        forcing_pd = forcingobject_to_dataframe(forcing)
        ax.set_title(f'{extbnd_filepath}')
        pc = forcing_pd.plot(ax=ax) #TODO: see CMEMS_interpolate_example.py for pcolormesh in case of verticalpositions
        leglabels = pc.get_legend_handles_labels()[1]
        leglabels_new.append(f'{forcing.name} ({forcing.function}) {leglabels[-1]}')
    ax.legend(leglabels_new,fontsize=8)
    fig.tight_layout()


#dimr = DIMR("dimr_config.xml")

