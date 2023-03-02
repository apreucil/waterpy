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
import matplotlib.pyplot as plt
import numpy as np
#%% User set up
basin = '145'
#%% waterpy
# path_to_output = r"D:\WATER_FILES\outputs\BaseV4_2011aULC_outputs"
path_to_output = r"D:\WATER_FILES\outputs\climbasin_outputs\basin_characteristics"

# filename = os.path.join(r"D:\WATER_FILES\CliimBasinParamsFromWATERXML",str(basin),"WATER.txt")
filename = os.path.join(r"D:\WATER_FILES\BaseV4_2011aLULC",str(basin),"WATER.txt")
climbasins = True

if climbasins:
    basins,watersheds = files()
    basins['center'] = basins.centroid
    key = watersheds.sjoin(basins)
    node = key[key.PSTSubNode == basin]
    climbasins = list(node[node.geometry.contains(node.center)].HydroID.values)
    
    fig,ax = plt.subplots()
    watersheds[watersheds.PSTSubNode==basin].plot(facecolor='grey',edgecolor='k',ax=ax)
    node.center.plot(color='r',ax=ax)
    node[node.geometry.contains(node.center)].center.plot(color='green',ax=ax,marker='+')
    plt.show()
    # Manual List
    # climbasins = [613,642,656,667,672,676,680,681,683,685,705,712,718,721,722]
    climbasins = [str(c) for c in climbasins]
    basin_df = pd.DataFrame()
    for b in climbasins:
        for lu in ['forest','agricultural','developed']:
            print (b)
            try:
                lu_df = pd.read_csv(os.path.join(path_to_output,str(b),lu,'output.csv'),index_col='date')
                # basin_df = pd.concat([basin_df,lu_df['discharge_predicted average (cfs)']],axis=1)
                basin_df = pd.concat([basin_df,lu_df['discharge_predicted']],axis=1)

                
            except:
                pass
            
        
    DRB_df = pd.DataFrame()
    DRB_df[basin] = basin_df.sum(axis=1)

else:
    basin_df = pd.DataFrame()
    for lu in ['forest','agricultural','developed']:
        print (basin, lu)
        try:
            lu_df = pd.read_csv(os.path.join(path_to_output,str(basin),lu,'output.csv'),index_col='date')
            # basin_df = pd.concat([basin_df,lu_df['discharge_predicted average (cfs)']],axis=1)
            basin_df = pd.concat([basin_df,lu_df['discharge_predicted']],axis=1)
        except:
            pass
        
    DRB_df = pd.DataFrame()
    DRB_df[basin] = basin_df.sum(axis=1)
DRB_df.index = pd.to_datetime(DRB_df.index)

#%% OLD WATER
file = open(filename)
lines = file.readlines()
file.close()

data = [re.split("\t|\n",line) for line in lines[6:]]
df = pd.DataFrame(data[1:],columns=data[0])
df['Date'] = pd.to_datetime(df.Date)
water_df = df.set_index('Date').apply(pd.to_numeric,errors='coerce')

#%% PST

basedata_path = r"C:\OASIS_DRBC\basedata\basedata_2017_v3.csv"
pst_df = pd.read_csv(basedata_path,index_col=0)
pst_df.index = pd.to_datetime(pst_df.index)

pst_basin = pst_df[basin]


#%% Compare

combine_df = pd.concat([DRB_df,water_df['Discharge (cfs)'],pst_basin],axis=1).dropna()
combine_df.columns = ['waterpy','WATER','PST']

sums = pd.DataFrame(combine_df.sum())
sums = sums.transpose()

#%% Display
fig,ax = plt.subplots(figsize=(10,7))
combine_df.cumsum().plot(ax=ax)
plt.title('PST basin: '+basin)
plt.legend()
ax.table(cellText=sums.round().values,rowLabels=['sum'],colLabels = sums.columns,
         bbox=[0.2, -0.2, 0.6, 0.1])
plt.ylabel('Cumulative Flow (CFS)')
plt.show()
#%%
fig,ax = plt.subplots(figsize=(8,8))
combine_df.describe([p for p in np.linspace(0,1,21)])[4:-2].plot(ax=ax)
plt.title('Cumulative Distribtion, node: '+basin)
plt.grid()
plt.legend(fontsize=12)
plt.ylabel('CFS')
