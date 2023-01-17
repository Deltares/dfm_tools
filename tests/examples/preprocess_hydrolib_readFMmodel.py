# -*- coding: utf-8 -*-
"""
Created on Mon Oct  3 12:07:18 2022

@author: veenstra
"""

from pathlib import Path
try: #0.3.1 release
    from hydrolib.core.io.mdu.models import FMModel, NetworkModel, ExtModel, StructureModel
    from hydrolib.core.io.bc.models import ForcingModel
except: #main branch and next release #TODO: move to easy imports after https://github.com/Deltares/HYDROLIB-core/issues/410
    from hydrolib.core.dflowfm.mdu.models import FMModel, NetworkModel, ExtModel, StructureModel
    from hydrolib.core.dflowfm.bc.models import ForcingModel
#TODO: #import hydrolib.core.dflowfm as hcdfm # (and add hydrolib-core>=0.4.0 to dependencies)
import dfm_tools as dfmt
import datetime as dt
import matplotlib.pyplot as plt
plt.close('all')

dtstart = dt.datetime.now()


file_mdu = Path(r'p:\11206813-006-kpp2021_rmm-2d\C_Work\31_RMM_FMmodel\computations\model_setup\run_206_HYDROLIB_structbc\RMM_dflowfm_general.mdu')
mdu_contents = open(str(file_mdu),'r').readlines()
if '[model]' in mdu_contents[0]:
    print('WARNING: [model] found in mdufile, this should be [general]')
fm = FMModel(file_mdu) #TODO: works, but many mdu-lines are comented, so uncomment things one by one


file_struct = Path(r'p:\11206813-006-kpp2021_rmm-2d\C_Work\31_RMM_FMmodel\computations\model_setup\run_206_HYDROLIB\RMM_structures.ini')
structs = StructureModel(file_struct) #TODO SOLVED with : pli in structures.ini is currently not supported: https://github.com/Deltares/HYDROLIB-core/issues/353
for struct in structs.structure:
    print(struct.id)
#structs.save('tst.ini') #TODO: how to get xydata from plifile in structuremodel?
#structs.structure[0].__dict__


file_network = Path(r'p:\11206813-006-kpp2021_rmm-2d\C_Work\31_RMM_FMmodel\computations\model_setup\run_206_HYDROLIB\rmm_v1p7_net.nc')
#file_network = Path(r'p:\1230882-emodnet_hrsm\GTSMv5.0\runs\reference_GTSMv4.1_wiCA\step11_global_1p25eu_net.nc')
#file_network = Path(r'p:\1230882-emodnet_hrsm\GTSMv5.0\runs\GM50_2000m_eu0900m_ITfac5p5_wx\gtsm_200s_2000m_eu0900m_ca2000m_v4_net.nc')
#network = NetworkModel(file_network) #TODO: what is this used for? plotting network/map is easier with xugrid


#file_extnew = Path(r'p:\11206813-006-kpp2021_rmm-2d\C_Work\31_RMM_FMmodel\computations\model_setup\run_206_HYDROLIB\RMM_bnd.ext') #TODO: waterlevelbnd for rivers are present three times: https://github.com/Deltares/HYDROLIB-core/issues/354
file_extnew = Path(r'p:\11206813-006-kpp2021_rmm-2d\C_Work\31_RMM_FMmodel\computations\model_setup\run_206_HYDROLIB\RMM_bnd_course.ext')
#file_extnew = Path(r'p:\1230882-emodnet_hrsm\GTSMv5.0\SO_NHrivGTSM\computations\BD013_4par_mildslope_wflowdis_JV\gtsm_forcing_bc.ext')
ext = ExtModel(file_extnew)


time_passed = (dt.datetime.now()-dtstart).total_seconds()
print(f'>>time passed: {time_passed:.2f} sec')

max_extforcings = 6 #None for all

ext_boundaries = ext.boundary
ext_laterals = ext.lateral
for iEB, extbnd in enumerate(ext_boundaries+ext_laterals): 
    if hasattr(extbnd,'forcingfile'):
        extbnd_filepath = extbnd.forcingfile.filepath
        extbnd_forcings = extbnd.forcingfile.forcing
    elif hasattr(extbnd,'discharge'):
        extbnd_filepath = extbnd.discharge.filepath
        extbnd_forcings = extbnd.discharge.forcing
    else:
        raise Exception('ERROR: not forcingfile/discharge present (boundary/lateral')
    print(f'boundary {iEB+1} of {len(ext_boundaries)}: {extbnd_filepath}')
    fig,ax = plt.subplots(figsize=(12,6))
    for iEBF, forcing in enumerate(extbnd_forcings[:max_extforcings]):
        print(f'forcing {iEBF+1} of {len(extbnd_forcings)}: {forcing.name} ({forcing.function}) ({forcing.quantityunitpair[1].quantity})')
        forcing_xr = dfmt.forcinglike_to_Dataset(forcing)
        data_vars = list(forcing_xr.data_vars.keys()) #mostly one variable, except for astronomic/uxuy bnd
        ax.set_title(f'{extbnd_filepath}')
        pc = forcing_xr[data_vars[0]].plot(ax=ax, label=f'{forcing.name} ({forcing.function}) ({forcing.quantityunitpair[1].quantity})') # see CMEMS_interpolate_example.py for pcolormesh in case of verticalpositions
    ax.legend(fontsize=8)
    fig.tight_layout()


#dimr = DIMR("dimr_config.xml")

