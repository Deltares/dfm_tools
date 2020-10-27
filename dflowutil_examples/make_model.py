'''
A tutorial on how to make a DFMWAQ model from a D-FLOW model

python setup.py sdist bdist_wheel
conda env create -f dflowutil.yml
conda activate dflowutil
pip install dist/dflowutil-0.1.3-py3-none-any.whl

'''

#import pytest
#import inspect
import os

#dir_tests = os.path.join(os.path.realpath(__file__), os.pardir)
#dir_testoutput = os.path.join(dir_tests,'test_output')
#if not os.path.exists(dir_testoutput):
#    os.mkdir(dir_testoutput)
dir_testinput = os.path.join(r'c:/DATA/werkmap','dfm_tools_testdata')

from dflowutil.DFMWAQModel import DFMWAQModel
from dflowutil.SubFile import SubFile
import datetime as dt


'''
model with boundary data
'''

process_path= os.path.join(dir_testinput,'DSD\\01_substances\\proc_def.dat')
system = 'windows'
mdu = os.path.join(dir_testinput,'DSD\\00_src\\current_situation.mdu')
sub_file = os.path.join(dir_testinput,'DSD\\01_substances\\guayas_V11.sub')

# boundaries for the mdu
ext = [os.path.join(dir_testinput,'DSD\\00_src\\plant_loads_current_local.ext'),
       os.path.join(dir_testinput,'DSD\\00_src\\sea_riv_boundary_local_bc.ext')]

# water quality data
loads_data  = [os.path.join(dir_testinput,r'DSD\03_loads\full_model_ds.csv')]
bounds_data = [os.path.join(dir_testinput,r'DSD\02_boundaries\river_ds.csv'),
               os.path.join(dir_testinput,r'DSD\02_boundaries\sea_ds.csv')]

# location of new model
new_model_dir = os.path.join(dir_testinput,r'DSD\\R01')

# sub file to use
substances = SubFile(sub_file)

# initial conditions
ini = {'OXY' : 7}
tref = dt.datetime(2000, 1, 1)
ecuador_model = DFMWAQModel(mdu, ext, substances, new_model_dir, tref, 
                            bounds_data = bounds_data, loads_data = loads_data, 
                            run_sys = system, process_path = process_path)
ecuador_model.build()

'''
make a second model to test restart functionality
'''

new_model_dir = os.path.join(dir_testinput,r'DSD\\R02')
ecuador_model = DFMWAQModel(mdu, ext, substances, new_model_dir, tref, 
                            bounds_data = bounds_data, loads_data = loads_data, 
                            run_sys = system, process_path = process_path)
ecuador_model.build()


