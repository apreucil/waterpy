# -*- coding: utf-8 -*-
"""
Created on Wed Jan 18 15:33:08 2023

@author: AnthonyPreucil
"""

# Post Processing for WATER-DRB
import pandas as pd
import os
import re
from drb_geoprocessing import files

path_to_output = r"D:\WATER_FILES\outputs\climbasin_outputs\basin_characteristics"
basin = '145'

filename = os.path.join(r"D:\WATER_FILES\BaseV4_2011aLULC",str(basin),"WATER.txt")
climbasins = False

if climbasins:
    basins,watersheds = files()
    key = basins.overlay(watersheds[['PSTSubNode','geometry']])
    
    climbasins = list(key[key.PSTSubNode == basin].HydroID.values)
    
    basin_df = pd.DataFrame()
    for b in climbasins:
        for lu in ['forest','agricultural','developed']:
            print (b)
            try:
                lu_df = pd.read_csv(os.path.join(path_to_output,str(b),lu,'output.csv'),index_col='date')
                basin_df = pd.concat([basin_df,lu_df['discharge_predicted average (cfs)']],axis=1)
            except:
                pass
            
        
    DRB_df = pd.DataFrame()
    DRB_df[basin] = basin_df.sum(axis=1)

else:
    for lu in ['forest','agricultural','developed']:
        print (b)
        try:
            lu_df = pd.read_csv(os.path.join(path_to_output,str(b),lu,'output.csv'),index_col='date')
            basin_df = pd.concat([basin_df,lu_df['discharge_predicted average (cfs)']],axis=1)
        except:
            pass
        
    DRB_df = pd.DataFrame()
    DRB_df[basin] = basin_df.sum(axis=1)
        
# OLD WATER
file = open(filename)
lines = file.readlines()
file.close()

data = [re.split("\t|\n",line) for line in lines[6:]]
df = pd.DataFrame(data[1:],columns=data[0])
df['Date'] = pd.to_datetime(df.Date)
df = df.set_index('Date').apply(pd.to_numeric,errors='coerce')