# -*- coding: utf-8 -*-
"""
Created on Mon Feb 27 09:08:05 2023

@author: AnthonyPreucil
"""

import os
import shutil


# Fix timeseries files
path_to_basins = r"D:\WATER_FILES\outputs\Climbasin_outputs\basin_characteristics"
subbasins = os.listdir(path_to_basins)
for b in subbasins:
    
    if os.path.exists(os.path.join(path_to_basins,b,'f')):
        os.rmdir(os.path.join(path_to_basins,b,'f'))
    else:
        pass

    