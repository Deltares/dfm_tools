# -*- coding: utf-8 -*-
"""
Created on Thu Oct  6 09:18:27 2022

@author: veenstra
"""
import os
import hydrolib.core.dflowfm as hcdfm
import numpy as np
from pathlib import Path

#TODO: bc files are probably not yet supported for structures. If they are, this is probably the way to couple them.

dir_out = r'p:\archivedprojects\11208053-004-kpp2022-rmm1d2d\C_Work\09_Validatie2018_2020\dflowfm2d-rmm_vzm-j19_6-v2d\computations\validation\2018_HYDROLIB_JV_structbc'

file_struct = os.path.join(dir_out,'RMM_structures_j19.ini')
structs = hcdfm.StructureModel(file_struct)
list_structures = [struct.id for struct in structs.structure]

ext_new = hcdfm.ExtModel()
ForcingModel_object = hcdfm.ForcingModel()

arr_time = np.array([0,10,20,60])
arr_GOW = np.array([0,360,360,0])
arr_GLEL = np.array([0,-5.5,-5.5,0])

for structname in list_structures:
    for quantity in ['general_structure_gateOpeningWidth','general_structure_gateLowerEdgeLevel']:
        if quantity=='general_structure_gateOpeningWidth':
            if 'Maeslan' not in structname:
                continue
        if quantity=='general_structure_gateOpeningWidth':
            datablock_incltime = np.concatenate([arr_time[:,np.newaxis],arr_GOW[:,np.newaxis]],axis=1)
        else:
            datablock_incltime = np.concatenate([arr_time[:,np.newaxis],arr_GLEL[:,np.newaxis]],axis=1)
            
        ts_one = hcdfm.TimeSeries(name=structname,
                                  quantityunitpair=[hcdfm.QuantityUnitPair(quantity="time", unit='minutes since 2018-01-20 00:00:00'),
                                                    hcdfm.QuantityUnitPair(quantity=quantity, unit='m')],
                                  timeinterpolation='linear',
                                  datablock=datablock_incltime.tolist(), 
                                  )
        
        ForcingModel_object.forcing.append(ts_one)
        # ForcingModel_object.save(Path(dir_out,'test_struct.bc')) #TODO: has to be saved all the time, otherwise forcingFile in ext is empty
        
        file_pli = f'../../../geometry_j19_6-v2/structures/{structname}.pli'
        boundary_object = hcdfm.Boundary(quantity=quantity,
                                         locationfile=Path(file_pli),
                                         forcingfile=ForcingModel_object,
                                         )
        ext_new.boundary.append(boundary_object)


# ext_new.save(Path(dir_out,'test_struct.ext')) #seems not necesary


#strmod = StructureModel(Path(r'p:\11206813-006-kpp2021_rmm-2d\C_Work\31_RMM_FMmodel\computations\model_setup\run_206_HYDROLIB_structbc\RMM_structures.ini'))

