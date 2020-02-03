'''
A tutorial on how to make initial conditions for a DFMWAQ model

'''

from dflowutil.utils import rst_to_xyz


'''
turn the last model into a restart
'''

sub_file = r'd:\projects\dflowutil\tests\DSD\01_substances\guayas_V11.sub'
map_dir = 'd:\\projects\\dflowutil\\tests\\DSD\\R01\\DFM_OUTPUT_current_situation\\'
new_model_dir = 'd:\\projects\\dflowutil\\tests\\DSD\\R02\\'
rst_to_xyz(map_dir, sub_file, -1, new_model_dir)

