# -*- coding: utf-8 -*-
"""
Created on Wed Jan 18 08:26:50 2023

@author: AnthonyPreucil
"""

import multiprocessing as mp
from GenerateInput import input_main
from drb_geoprocessing import geo_main, files
from worker_file import main
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
        key = basins.overlay(watersheds[['PSTSubNode','geometry']])
        
        num_workers = mp.cpu_count()  
        
        pool = mp.Pool(num_workers)
        climbasins = list(key[key.PSTSubNode == pst_node].HydroID.values)
        for climbasin in climbasins:
            pool.apply_async(geo_main, args = (climbasin,basins,watersheds))
        
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
    
#%% Run
# run_water('110')
# run_waterDRB(['115A','115B'])
drb_geoprocess('165')
# gen_input()
