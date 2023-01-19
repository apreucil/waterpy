# -*- coding: utf-8 -*-
"""
Created on Wed Jan 18 14:13:43 2023

@author: AnthonyPreucil
"""

# Run Water Case Studies

import os
import sys
sys.path.append(r'D:\waterpy\waterpy')
from main import waterpy


basin = '100B'
lu = 'agricultural'


input_path = os.path.join(r'D:\WATER_FILES\inputs\BaseV4_2011aULC_inputs',basin)
configfile = os.path.join(input_path,"modelconfigfile_"+lu+".ini")
waterpy(configfile,None)