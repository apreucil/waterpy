# -*- coding: utf-8 -*-
"""
Created on Wed Jan 18 15:33:08 2023

@author: AnthonyPreucil
"""

# Post Processing for WATER-DRB
import pandas as pd
import os
import re

path_to_output = r"D:\WATER_FILES\outputs\BaseV4_2011aULC_outputs"
basin = '110'

basin_df = pd.DataFrame()
for lu in ['forest','agricultural','developed']:
    print (lu)
    lu_df = pd.read_csv(os.path.join(path_to_output,basin,lu,'output.csv'),index_col='date')
    basin_df = pd.concat([basin_df,lu_df['discharge_predicted average (cfs)']],axis=1)
    
DRB_df = pd.DataFrame()
DRB_df['110'] = basin_df.sum(axis=1)


# OLD WATER
filename = r"D:\WATER_FILES\BaseV4_2011aLULC\110\WATER.txt"
file = open(filename)
lines = file.readlines()
file.close()

data = [re.split("\t|\n",line) for line in lines[6:]]
df = pd.DataFrame(data[1:],columns=data[0])
df['Date'] = pd.to_datetime(df.Date)
df = df.set_index('Date').apply(pd.to_numeric,errors='coerce')