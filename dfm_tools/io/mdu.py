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

Original author: Koen D. Berends (Deltares / University of Twente)
k.d.berends@utwente.nl / koen.berends@deltares.nl
"""


def read_deltares_ini(filename):
    """
    MDU files are like ini files, but can contain multiple copies of the same 
    headers. Therefore, we adopt a different structure where each 
    section is one entry in a list. 

    read/write Delft3D-FM *.mdu input files to/from dictionary
    
    Example use:
    
    data = mdu.read('ref.mdu')
    
    # Print all sections
    print (data['section'].unique())
    
    # Print all keys in section 'numerics'
    print (data[data['section'] == 'numerics'].key.values)
    
    
    # Print value of key 'TimeStepType' in section 'Numerics'
    print data[(data['section'] == 'numerics') & 
    	   (data['key'] == 'TimeStepType')]
    
    # Change value to 5
    data.loc[(data['section'] == 'numerics') & 
    	 (data['key'] == 'TimeStepType'), 
    	 'value'] = 5
    
    # Write to file
    mdu.write(data, 'ref_mod.mdu')
    """
    from collections import OrderedDict
    import re
    import pandas as pd

    data = pd.DataFrame(columns=('section', 'key', 'value', 'comment'))
    current_section = ''
    comment_re = re.compile(r'((\s+)?\#.+)')
    section_re = re.compile(r'^\s*\[\s*(.+?)\s*\]')
    keyval_re = re.compile(r'^\s*(.+?)\s*=\s*(.*?)\s*$')

    with open(filename, 'r') as f:
        for i, line in enumerate(f):
            data.loc[i] = ['', '', '', '']
            comment = comment_re.search(line)
            section = section_re.search(line)
            keyvalue = keyval_re.search(line)
            
            if comment:
                data.loc[i]['comment'] = comment.group(1)
            else:
                data.loc[i]['comment'] = ''
            if section:
                current_section = section.group(1)
                data.loc[i]['section'] = current_section  
            if keyvalue:
                data.loc[i]['key'] = keyvalue.group(1).split('#')[0]
                data.loc[i]['value'] = keyvalue.group(2).split('#')[0]
                data.loc[i]['section'] = current_section  
    
    return data

def write_deltares_ini(filedata, filename):
    keylen = filedata.key.str.len().max()+1
    f = open(filename,'w')
    for i, row in filedata.iterrows():
        if row['section'] and not row['key'] and not row['comment']:
            f.write('[{}]'.format(row['section']))
        elif row['key']:
            f.write('{: <{}}= {}'.format(row['key'], keylen, row['value']))
        if row['comment']:
            f.write('{}'.format(row['comment']))
        f.write('\n')

