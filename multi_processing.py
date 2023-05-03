# -*- coding: utf-8 -*-
"""
Created on Wed Jan 18 08:26:50 2023

@author: AnthonyPreucil
"""


import multiprocessing as mp
from GenerateInput import input_main
from drb_geoprocessing import geo_main, files
from worker_file import main
import pandas as pd
import os
#%% For input generation
def gen_input():
    if __name__ == '__main__':
        num_workers = mp.cpu_count()  
        
        pool = mp.Pool(num_workers)
        subbasins = [folder for folder in os.listdir("D:\WATER_FILES\BaseV4_2011aLULC") if len(folder)<=4]
        for basin in subbasins:
            pool.apply_async(input_main, args = (basin,))
        
        pool.close()
        pool.join()
#%% For geoprocessing
def drb_geoprocess(pst_node):
    if __name__ == '__main__':
        
        basins,watersheds = files()
        basins['center'] = basins.centroid
        key = watersheds.sjoin(basins)
        node = key[key.PSTSubNode == b]
        climbasins = list(node[node.geometry.contains(node.center)].HydroID.values)
        
        num_workers = mp.cpu_count()  
        
        pool = mp.Pool(num_workers)
        # climbasins = list(key[key.PSTSubNode == pst_node].HydroID.values)
        for climbasin in climbasins:
            if not os.path.exists(os.path.join(r"D:\WATER_FILES\inputs\Climbasin_inputs\basin_characteristics",str(climbasin))):
                print ('Processing '+str(climbasin)+'...')
                # geo_main(climbasin,basins,watersheds,pst_node) # For Debug
                pool.apply_async(geo_main, args = (climbasin,basins,watersheds,pst_node))
            else:
                pass
        
        pool.close()
        pool.join()
#%% For running the model one basin at a time
# THIS WORKS! NOTE, RUN WHOLE FILE, NOT JUST CELL OR IT WON'T WORK.
def run_water(basin):
    if __name__=='__main__':
        num_workers = mp.cpu_count()  
        
        pool = mp.Pool(num_workers)
        land_use = ['forest','agricultural','developed']
        for lu in land_use:
            
            path = os.path.join(r"D:\WATER_FILES\inputs\BaseV4_2011aULC_inputs",basin)
            pool.apply_async(main, args = (path,lu,))
        
        pool.close()
        pool.join()
#%% For running the entire basin!       
def run_waterDRB(subbasins):
    if __name__=='__main__':
        num_workers = mp.cpu_count()  
        
        pool = mp.Pool(num_workers)
        land_use = ['forest','agricultural','developed']
        # subbasins = [folder for folder in os.listdir("D:\WATER_FILES\BaseV4_2011aLULC") if len(folder)<=4]
        for basin in subbasins:
            for lu in land_use:
                
                path = os.path.join(r"D:\WATER_FILES\inputs\BaseV4_2011aULC_inputs",basin)
                pool.apply_async(main, args = (path,lu,))
        
        pool.close()
        pool.join()
#%% For running smaller climbasins 
 
def run_one_climbasin(climbasin):
    if __name__ == '__main__':
        num_workers = mp.cpu_count()  
        
        pool = mp.Pool(num_workers)
        land_use = ['forest','agricultural','developed']
        for lu in land_use:
            
            path = os.path.join(r"D:\WATER_FILES\inputs\Climbasin_inputs\basin_characteristics",climbasin)
            # main(path,lu) ## For debug
            pool.apply_async(main, args = (path,lu,))
        
        pool.close()
        pool.join()
    
def run_list_climbasin(climbasins):
    if __name__ == '__main__':
        num_workers = mp.cpu_count()  
        
        
        pool = mp.Pool(num_workers)
        land_use = ['forest','agricultural','developed']
        for climbasin in climbasins:
            print (climbasin)
            path = os.path.join(r"D:\WATER_FILES\inputs\Climbasin_inputs\basin_characteristics",str(climbasin))
            for lu in land_use:
                
                # main(path,lu) ## For debug
                pool.apply_async(main, args = (path,lu,))
            
        pool.close()
        pool.join()
    
    
def run_climbasins(pst_basin):
    if __name__ == '__main__':
        basins,watersheds = files()
        # Output is generated based on original watersheds to ensure 
        # correct masks are used
        basins['center'] = basins.centroid
        key = watersheds.sjoin(basins)
        node = key[key.PSTSubNode == pst_basin]
        climbasins = list(node[node.geometry.contains(node.center)].HydroID.values)
        
        num_workers = mp.cpu_count()
        pool = mp.Pool(num_workers)
        
        land_use = ['forest','agricultural','developed']
        for climbasin in climbasins:
            path = os.path.join(r"D:\WATER_FILES\inputs\Climbasin_inputs\basin_characteristics",str(climbasin))
            for lu in land_use:
                try:
                    chars = pd.read_csv(os.path.join(path,'basin_characteristics_'+lu+'.csv'),index_col='name')
                    if chars.loc['basin_area_total'].values != 0:
                        pool.apply_async(main, args = (path,lu,))
                    else:
                        print ('ALERT: Skipping '+str(climbasin)+' - '+lu+': No landuse detected')
                except Exception as e:
                    print ('Error with Climate Basin: '+str(climbasin)+' - '+lu)
                    print (e)
                    
                
        pool.close()
        pool.join()
        
#%% Run
# run_water('1931')
# run_waterDRB(['115A','115B'])
# subbasins = [folder for folder in os.listdir("D:\WATER_FILES\BaseV4_2011aLULC") if len(folder)<=4]
# subbasins = [b for b in subbasins if int(b[:3]) <=365] # To trenton only!
# subbasins = ['100A',
#               '100B']
# for b in subbasins:
#     drb_geoprocess(b)
# gen_input()

# landuse=['forest','agricultural','developed']
# subbasins = []
# for cb in os.listdir(r"D:\WATER_FILES\outputs\climbasin_outputs\basin_characteristics"):
#     ex = 0
#     for lu in landuse:
#         exist = os.path.exists(os.path.join(r"D:\WATER_FILES\outputs\climbasin_outputs\basin_characteristics",str(cb),lu,"output.csv"))
#         if exist:
#             pass
#         else:
#             ex+=1
#     if ex == 3:
#         subbasins.append(cb)
#     else:
#         pass

import geopandas as gpd

gdf = gpd.read_file(r"D:\WATER_FILES\USGS_Gage_Water_PY_test\Watersheds.shp")
# Climate Basins File
gdf_clim = gpd.read_file(r"D:\WATER_FILES\climbasins.shp")

clip_clim_gdf = gpd.sjoin(gdf_clim,gdf,predicate='within')

subbasins = clip_clim_gdf.HydroID.values
print ('Running for USGS gage basins only')
subbasins = subbasins[230:]

run_list_climbasin(subbasins)

# for sb in subbasins:
#     print (sb)
#     run_one_climbasin(str(sb))

# subbasins = [folder for folder in os.listdir("D:\WATER_FILES\BaseV4_2011aLULC") if len(folder)<=4]
# for b in subbasins:
#     print (b)
#     run_climbasins(b)