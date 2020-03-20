# -*- coding: utf-8 -*-
"""
Created on Thu Mar 19 15:27:12 2020

@author: veenstra
"""


def getmakeoutputdir(script_dir, function_name):
    import os
    dir_tests = os.path.join(os.path.realpath(script_dir), os.pardir)
    dir_testoutput = os.path.join(dir_tests,'test_output')
    if not os.path.exists(dir_testoutput):
        os.mkdir(dir_testoutput)
    scriptname = os.path.basename(script_dir).replace('.','')
    dir_output = os.path.join(dir_testoutput,'%s - %s'%(scriptname,function_name))
    if not os.path.exists(dir_output):
        os.mkdir(dir_output)
    return dir_output
