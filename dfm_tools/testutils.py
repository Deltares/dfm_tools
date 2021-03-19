# -*- coding: utf-8 -*-
"""
dfm_tools are post-processing tools for Delft3D FM
Copyright (C) 2020 Deltares. All rights reserved.

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.

All names, logos, and references to "Deltares" are registered trademarks of
Stichting Deltares and remain full property of Stichting Deltares at all times.
All rights reserved.


INFORMATION
This script is part of dfm_tools: https://github.com/openearth/dfm_tools
Check the README.rst on github for other available functions
Check the tests folder on github for example scripts (this is the dfm_tools pytest testbank)
Check the pptx and example figures in (created by the testbank): N:/Deltabox/Bulletin/veenstra/info dfm_tools

Created on Thu Mar 19 15:27:12 2020

@author: veenstra
"""

def gettestinputdir():
    import os
    
    if 'TEAMCITY_VERSION' in os.environ.keys(): #teamcity path
        dir_testinput = r'\\dfs-trusted.directory.intra\dfs\Teamcity\Testdata\dfm_tools'
    else: #default to this path
        dir_testinput = os.path.join(r'c:\DATA','dfm_tools_testdata')
    
    return dir_testinput



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


def try_importmodule(modulename=None):
    command = '\t- open command window (or anaconda prompt)\n\t- conda activate dfm_tools_env\n\t- conda install -c conda-forge %s'%(modulename)

    try:
        exec('import %s'%(modulename))
    except:
        raise Exception('ERROR: module %s not found, do the following:\n%s'%(modulename, command))
    if modulename == 'shapely':
        try:
            import shapely.geometry
        except:
            raise Exception('ERROR: cannot execute "import shapely.geometry", do the following:\n%s'%(command))
        shpvers = [int(x) for x in shapely.__version__.split('.')]
        correctversion = False
        if shpvers[0] == 1:
            if shpvers[1] <=7:
                correctversion = True
        elif shpvers[0] > 2:
            correctversion = True
        if correctversion is False:
            raise Exception('ERROR: incorrect shapely version (%s), should be 1.7.0 or higher, do the following:\n%s'%(shapely.__version__, command))



            