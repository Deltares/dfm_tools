'''
A tutorial on how to make initial conditions for a DFMWAQ model

'''
#import pytest
#import inspect
import os

#dir_tests = os.path.join(os.path.realpath(__file__), os.pardir)
#dir_testoutput = os.path.join(dir_tests,'test_output')
#if not os.path.exists(dir_testoutput):
#    os.mkdir(dir_testoutput)
dir_testinput = os.path.join(r'c:/DATA/werkmap','dfm_tools_testdata')


from dflowutil.utils import rst_to_xyz


'''
turn the last model into a restart
'''

sub_file = os.path.join(dir_testinput,'DSD\\01_substances\\guayas_V11.sub')
map_dir = os.path.join(dir_testinput,'DSD\\R01\\DFM_OUTPUT_current_situation')
new_model_dir = os.path.join(dir_testinput,'DSD\\R02')
rst_to_xyz(map_dir, sub_file, -1, new_model_dir)

