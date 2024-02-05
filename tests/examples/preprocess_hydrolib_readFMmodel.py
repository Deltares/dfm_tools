# -*- coding: utf-8 -*-
"""
Created on Mon Oct  3 12:07:18 2022

@author: veenstra
"""

from pathlib import Path
import hydrolib.core.dflowfm as hcdfm
# from hydrolib.core.dimr import DIMR
import dfm_tools as dfmt
import datetime as dt
import matplotlib.pyplot as plt
plt.close('all')

#file_mdu = Path(r'p:\archivedprojects\11208053-004-kpp2022-rmm1d2d\C_Work\09_Validatie2018_2020\dflowfm2d-rmm_vzm-j19_6-v2d\computations\validation\2018_HYDROLIB_JV\RMM_VZM.mdu') #TODO: contains multiline obs/crs files, maybe raise more clear validationerror than "File: `\\directory.intra\Project\` not found, skipped parsing. (type=value_error)"
file_mdu = Path(r'p:\archivedprojects\11208053-004-kpp2022-rmm1d2d\C_Work\09_Validatie2018_2020\dflowfm2d-rmm_vzm-j19_6-v2d\computations\validation\2018_HYDROLIB_JV\RMM_VZM_MLnoslash.mdu') #TODO: contains multiline obs/crs files without slash, only reads first file but should raise an error: https://github.com/Deltares/HYDROLIB-core/issues/485#issuecomment-1494053471
file_mdu = Path(r'p:\archivedprojects\11208053-004-kpp2022-rmm1d2d\C_Work\09_Validatie2018_2020\dflowfm2d-rmm_vzm-j19_6-v2d\computations\validation\2018_HYDROLIB_JV\RMM_VZM_noML.mdu')
print('>> opening FMModel: ',end='')
dtstart = dt.datetime.now()
fm = hcdfm.FMModel(file_mdu, recurse=False) #TODO: works, but many mdu-lines are comented, so uncomment things one by one
print(f'{(dt.datetime.now()-dtstart).total_seconds():.2f} sec')

#TODO: enable caching possible to save time on second load?
#TODO: some things are not loaded, what happens when saving at new location?
#TODO: is model validated? (e.g. presence of diskonlyfilemodels or other non-loaded things)
# file_mdu_rewrite = str(file_mdu).replace('.mdu','_rewritepy.mdu')
# fm.save(file_mdu_rewrite)
#TODO: when writing mdu again, are all keywords written: https://github.com/Deltares/HYDROLIB-core/issues/495
structs = fm.geometry.structurefile
if structs is not None: #is empty when commented
    for struct in structs[0].structure: #TODO: why is structs a list of StructureModel instead of just StructureModel?
        print(struct.id)
#TODO: add newext/oldext printer
obs_pd_list = [dfmt.pointlike_to_DataFrame(x) for x in fm.output.obsfile] #TODO: "AttributeError: 'ObservationPointModel' object has no attribute 'points'" >> might be because of double quotes in names in combination with union type mess-up >> results in empty ObservationPointModel instead of XYNModel: https://github.com/Deltares/HYDROLIB-core/issues/519
crs_pd_list = [dfmt.PolyFile_to_geodataframe_linestrings(x) for x in fm.output.crsfile]


file_struct = Path(r'p:\archivedprojects\11206813-006-kpp2021_rmm-2d\C_Work\31_RMM_FMmodel\computations\model_setup\run_206_HYDROLIB\RMM_structures.ini')
file_struct = r'p:\archivedprojects\11205258-006-kpp2020_rmm-g6\C_Work\08_RMM_FMmodel\structures_toRTC\RMM_structures_open_cl10_coeff10.ini'
file_struct = r'p:\archivedprojects\11205258-006-kpp2020_rmm-g6\C_Work\08_RMM_FMmodel\structures_toRTC\RMM_structures_ts_cl10_coeff10.ini' #TODO: incorrect timfiles (e.g. with duplicate times) are now read as ForcingModel instead: https://github.com/Deltares/HYDROLIB-core/issues/519 #TODO: some of the timfiles contain duplicate times, but these are not validated (probably due due to union type)
print('>> opening structurefile: ',end='')
dtstart = dt.datetime.now()
structs = hcdfm.StructureModel(file_struct, recurse=False)
print(f'{(dt.datetime.now()-dtstart).total_seconds():.2f} sec')
for struct in structs.structure:
    print("  ", struct.id)
    print("    ", type(struct.gateloweredgelevel))
#structs.save('tst.ini') #TODO: how to get xydata from plifile in structuremodel? >> dfmt.pointlike_to_DataFrame() per object?
#structs.structure[0].__dict__


#file_extnew = Path(r'p:\archivedprojects\11206813-006-kpp2021_rmm-2d\C_Work\31_RMM_FMmodel\computations\model_setup\run_206_HYDROLIB\RMM_bnd.ext') #TODO: waterlevelbnd for rivers are present three times: https://github.com/Deltares/HYDROLIB-core/issues/354
file_extnew = Path(r'p:\archivedprojects\11206813-006-kpp2021_rmm-2d\C_Work\31_RMM_FMmodel\computations\model_setup\run_206_HYDROLIB\RMM_bnd_course.ext')
#file_extnew = Path(r'p:\1230882-emodnet_hrsm\GTSMv5.0\SO_NHrivGTSM\computations\BD013_4par_mildslope_wflowdis_JV\gtsm_forcing_bc.ext')
#file_extnew = Path(r'p:\archivedprojects\11208053-004-kpp2022-rmm1d2d\C_Work\09_Validatie2018_2020\dflowfm2d-rmm_vzm-j19_6-v2d\computations\validation\2018_HYDROLIB_JV\RMM_bnd.ext') #slow because of waterlevelbnd seaside
print('>> opening ExtModel: ',end='')
dtstart = dt.datetime.now()
ext = hcdfm.ExtModel(file_extnew, recurse=True)
print(f'{(dt.datetime.now()-dtstart).total_seconds():.2f} sec')

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
    print(f'  boundary {iEB+1} of {len(ext_boundaries)}: {extbnd_filepath}')
    fig,ax = plt.subplots(figsize=(12,6))
    for iEBF, forcing in enumerate(extbnd_forcings[:max_extforcings]):
        print(f'    forcing {iEBF+1} of {len(extbnd_forcings)}: {forcing.name} ({forcing.function}) ({forcing.quantityunitpair[1].quantity})')
        forcing_xr = dfmt.forcinglike_to_Dataset(forcing)
        data_vars = list(forcing_xr.data_vars.keys()) #mostly one variable, except for astronomic/uxuy bnd
        ax.set_title(f'{extbnd_filepath}')
        pc = forcing_xr[data_vars[0]].plot(ax=ax, label=f'{forcing.name} ({forcing.function}) ({forcing.quantityunitpair[1].quantity})') # see CMEMS_interpolate_example.py for pcolormesh in case of verticalpositions
    ax.legend(fontsize=8)
    fig.tight_layout()


# file_dimr = r"p:\archivedprojects\11206813-006-kpp2021_rmm-2d\C_Work\31_RMM_FMmodel\computations\model_setup\run_206_HYDROLIB\dimr_config.xml"
# dimr = DIMR(file_dimr, recurse=False)

