'''
A tutorial on how to make a DFMWAQ model from a D-FLOW model

python setup.py sdist bdist_wheel
conda env create -f dflowutil.yml
conda activate dflowutil
pip install dist/dflowutil-0.1.3-py3-none-any.whl

'''

from dflowutil.DFMWAQModel import DFMWAQModel
from dflowutil.SubFile import SubFile
import datetime


'''
model with boundary data
'''

process_path= r'd:\projects\dflowutil\tests\DSD\01_substances\proc_def.dat'
system = 'windows'
mdu = r'd:\projects\dflowutil\tests\DSD\00_src\current_situation.mdu'
sub_file = r'd:\projects\dflowutil\tests\DSD\01_substances\guayas_V11.sub'

# boundaries for the mdu
ext = [r'd:\projects\dflowutil\tests\DSD\00_src\plant_loads_current_local.ext',
       r'd:\projects\dflowutil\tests\DSD\00_src\sea_riv_boundary_local_bc.ext']

# water quality data
loads_data  = [r'd:\projects\dflowutil\tests\DSD\03_loads\full_model_ds.csv']
bounds_data = [r'd:\projects\dflowutil\tests\DSD\02_boundaries\river_ds.csv',
               r'd:\projects\dflowutil\tests\DSD\02_boundaries\sea_ds.csv']

# location of new model
new_model_dir = 'd:\\projects\\dflowutil\\tests\\DSD\\R01\\'

# sub file to use
substances = SubFile(sub_file)

# initial conditions
ini = {'OXY' : 7}
tref = datetime.datetime(2000, 1, 1)
ecuador_model = DFMWAQModel(mdu, ext, substances, new_model_dir, tref, 
                            bounds_data = bounds_data, loads_data = loads_data, 
                            run_sys = system, process_path = process_path)
ecuador_model.build()

'''
make a second model to test restart functionality
'''

new_model_dir = 'd:\\projects\\dflowutil\\tests\\DSD\\R02\\'
ecuador_model = DFMWAQModel(mdu, ext, substances, new_model_dir, tref, 
                            bounds_data = bounds_data, loads_data = loads_data, 
                            run_sys = system, process_path = process_path)
ecuador_model.build()


