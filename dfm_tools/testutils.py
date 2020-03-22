# -*- coding: utf-8 -*-
"""
GNU GENERAL PUBLIC LICENSE
	      Version 3, 29 June 2007

dfm_tools are post-processing tools for Delft3D FM
Copyright (C) 2020 Deltares

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.


INFORMATION
This script is part of dfm_tools: https://github.com/openearth/dfm_tools
Check the README.rst on github for other available functions
Check the tests folder on github for example scripts (this is the pytest testbank of dfm_tools)
Check the pptx and example figures in (created by the testbank): N:/Deltabox/Bulletin/veenstra/info dfm_tools

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
