# -*- coding: utf-8 -*-
"""
Created on Wed Jan 25 08:06:19 2023

@author: AnthonyPreucil
"""

# Import Libraries

import numpy as np
import pandas as pd
import geopandas as gpd
import os
import sys
sys.path.append('D:/waterpy/waterpy')
# from geospatial import *
import rasterio
from rasterio.features import geometry_mask
# from rasterio.plot import show
from rasterio.features import shapes
from shapely.geometry import shape
import skimage.transform as st
# import matplotlib.pyplot as plt

#%% DRB Specific Functions for area and twi
# Borrowed calculations from geospatial
# Modified to worked with masked arrays given a transform

def calc_area(masked_arr,transform):
    """
    

    Parameters
    ----------
    masked_arr : TYPE
        DESCRIPTION.
    transform : TYPE
        DESCRIPTION.

    Returns
    -------
    tot_area : TYPE
        DESCRIPTION.

    """
    raster_shape = shapes(masked_arr,transform = transform)
    tot_area = 0
    for myshape,value in raster_shape:
        polygon = shape(myshape)
        area = polygon.area
        tot_area += area
    return tot_area

def calc_twi(masked,nbins=30):
    """
    Function copied from geospatial.py
    and modified to work with the masked
    land use twi dataset.

    Parameters
    ----------
    masked : TYPE
        DESCRIPTION.
        
    nbins

    Returns
    -------
    df : TYPE
        DESCRIPTION.

    """
    mx = np.nanmax(masked.filled(np.nan))
    mn = np.nanmin(masked.filled(np.nan))
    intvl = (mx - mn) / (nbins + 1)
    edges = np.arange(mn, mx, intvl)
    histo = np.histogram(masked.filled(np.nan), bins=edges)


    # need mean of each bin.  Get the rest of the stats while there.
    # TWI Mean is the value we need for TopModel Input.

    bins = []

    for i in range(nbins):
        line = []
        bin = i + 1
        if i == 0:
            twi_val = histo[1][i] / 2
        else:
            twi_val = (histo[1][i] + histo[1][i-1]) / 2
        proportion = histo[0][i]/np.sum(histo[0])

        line.append(bin)
        line.append(twi_val)
        line.append(proportion)
        bins.append(line)

    df = pd.DataFrame(bins, columns=['bin', 'twi', 'proportion'])
    df.set_index('bin', inplace=True)

    return df
#%% Calculate characteristics
def drb_lu_params():
    
    # Pre-defined parameters (page 21 https://pubs.usgs.gov/sir/2015/5143/sir20155143.pdf)
    # These are land-use dependent, but not location or seasonally dependent at this time
    params = {}
    params['f'] = {'spatial_coef':0.4,
                   'rooting_depth_factor':0.75,
                   'macropore_fraction':0.15,
                   'impervious_curve_number':90,
                   'snowmelt_temperature_cutoff':0,
                   'snowmelt_rate_coeff':2,
                   'snowmelt_rate_coeff_with_rain':3,
                   'pet_calib_coeff':1.2,
                   'lake_delay':15,
                   'eff_imp':0.7,
                   'imp_delay':0.1,
                   'twi_adj':1,
                   'et_exp_dorm':5,
                   'et_exp_grow':0.5,
                   'grow_trigger':15}
    
    params['a'] = {'spatial_coef':0.3,
                   'rooting_depth_factor':0.25,
                   'macropore_fraction':0.15,
                   'impervious_curve_number':90,
                   'snowmelt_temperature_cutoff':0,
                   'snowmelt_rate_coeff':2,
                   'snowmelt_rate_coeff_with_rain':4,
                   'pet_calib_coeff':1.2,
                   'lake_delay':1.5,
                   'eff_imp':1,
                   'imp_delay':0.5,
                   'twi_adj':1,
                   'et_exp_dorm':5,
                   'et_exp_grow':0.5,
                   'grow_trigger':15}
    
    params['r'] = {'spatial_coef':0.25,
                   'rooting_depth_factor':0.75,
                   'macropore_fraction':0.2,
                   'impervious_curve_number':100,
                   'snowmelt_temperature_cutoff':0,
                   'snowmelt_rate_coeff':4,
                   'snowmelt_rate_coeff_with_rain':6,
                   'pet_calib_coeff':1.2,
                   'lake_delay':1.5,
                   'eff_imp':1,
                   'imp_delay':0.1,
                   'twi_adj':0.5,
                   'et_exp_dorm':5,
                   'et_exp_grow':0.5,
                   'grow_trigger':15}
    return params
    
def drb_characteristics(db_rasters,basin,path_to_masks,basins,watersheds):    
    chars = {}
    twis = {}
    # basin_geometry = basins.loc[[basin],'geometry'].geometry.iloc[0]
    gdf = gpd.GeoDataFrame(basins.loc[[basin],'geometry']).overlay(watersheds[['PSTSubNode','geometry']])
    gdf['area'] = gdf.area.astype(int) # added to filter out broken basins
    gdf = gdf[gdf.area >= 10] # added to filter out broken basins
    for landuse in ['f','r','a']:
        rasters = {}
        if landuse=='a':
            root_depth_factor = .25
        else:
            root_depth_factor = .75
        # Since landuse params do not change, output once
        # params_lu = params['f']
        # gdf = gpd.GeoDataFrame(basins.loc[[basin],'geometry']).overlay(watersheds[['PSTSubNode','geometry']])
        for pst in gdf.PSTSubNode.values:
            basin_geometry = gdf[gdf.PSTSubNode==pst].geometry.iloc[0]
            # fig,ax = plt.subplots()
            # gdf[gdf.PSTSubNode==pst].plot(ax=ax)
            # plt.show()
            mask_path = os.path.join(path_to_masks,str(pst),landuse+"mask\w001001.adf")
            
            with rasterio.open(mask_path) as landuse_src:
                col_off, row_off, lwidth, lheight = landuse_src.window(*basin_geometry.bounds).flatten()
                col_start = int(np.round(col_off))
                row_start = int(np.round(row_off))
                col_stop = col_start + int(np.round(lwidth))
                row_stop = row_start + int(np.round(lheight))
                
                basewindow = ((row_start, row_stop), (col_start, col_stop))
                window = tuple(tuple(max(0, num) for num in inner_tuple) for inner_tuple in basewindow)
                landuse_mask = landuse_src.read(1, masked=True, window=window)
                # show(landuse_mask,title='The mask')
                
                mask = geometry_mask([basin_geometry], landuse_mask.shape, transform=landuse_src.window_transform(window), invert=True)
                lu_basin = np.ma.masked_array(landuse_mask,mask=~mask)
                
            # show(lu_basin,title='Basin Masked')
            
            lu_area = calc_area(lu_basin,landuse_src.transform)

            db_rasters_clip = {}
            db_rasters_clip['area'] = np.ma.array([lu_area],mask=[False])
            for k, v in db_rasters.items():
                # print (k)
                with rasterio.open(v) as src:
                    col_off, row_off, width, height = src.window(*basin_geometry.bounds).flatten()
                    col_start = int(np.round(col_off))
                    row_start = int(np.round(row_off))
                    col_stop = col_start + int(np.round(width))
                    row_stop = row_start + int(np.round(height))                
                    
                    window = ((row_start, row_stop), (col_start, col_stop))
                    data = src.read(1, masked=True, window=window)
                    # show(data)
                    mask = geometry_mask([basin_geometry], data.shape, transform=src.window_transform(window), invert=True)
                    data = np.ma.masked_array(data, mask=~mask)
                    # show(data,title='Before Mask '+k)
                    
                    if data.shape == landuse_mask.shape:
                        lu_data = np.ma.masked_where(landuse_mask!=1,data)
                        # show(lu_data,title='After Mask '+k)
                    else:
                        # imp_transform=d_transform
                        # print (data.shape)
                        # rs_landuse_mask = st.resize(landuse_mask.mask,data.shape)
                        data = np.ma.masked_array(data, mask=~mask)
                        # correct_res_data = st.resize_local_mean(data,lu_basin.shape,grid_mode=False,preserve_range=True)
                        correct_res_data = st.resize(data,lu_basin.shape,mode='constant',cval=0,preserve_range=True)
                        # scale_factor = np.sum(data)/np.sum(correct_res_data)
                        # correct_res_data = resample_imp(v,mask_path)
                        # print (k+ ' was corrected for resolution, corrected to: ')
                        # print (correct_res_data.shape)
                        # show(correct_res_data,title='Resolution Corrected for '+k)
                        # blu_data = np.ma.masked_array(correct_res_data, mask=~mask)
                        lu_data = np.ma.masked_where(lu_basin!=1,correct_res_data)
                        # print(lu_data.sum()/230400*100)
                        # lu_data = mlu_data*scale_factor
                        # show(mlu_data,title='After Mask '+k)
                        # zlu_data = mlu_data.filled(-99999)
                        # lu_data = st.resize(zlu_data,data.shape,mode='constant')
                        # show(lu_data,title='After mask, reduced to original size, filled with 0s '+k)
                        # rs_lu_mask = np.ma.masked_where(rs_landuse_mask==False,rs_landuse_mask)
                        # lu_data = np.ma.masked_where((landuse_mask!=True)|(mlu_data==0),mlu_data)
                        # show(lu_data,title='After Mask '+k)
                
                db_rasters_clip[k] = lu_data
            
            rasters[pst] = db_rasters_clip
                
    
        joined_rasters = {key: np.ma.concatenate([d[key].compressed() for d in rasters.values()]) 
                          for key in rasters[list(rasters.keys())[0]].keys()}

        characteristics_out = {
            "scaling_parameter": {joined_rasters['scaling_parameter'].mean() / 100},
            "saturated_hydraulic_conductivity": {joined_rasters['k_sat'].mean() / 100 * 86.4},
            "saturated_hydraulic_conductivity_multiplier": {joined_rasters['con_mult'].mean() / 100},
            "soil_depth_total": {joined_rasters["soil_thickness"].mean() / 10},
            "soil_depth_roots": {(joined_rasters["soil_thickness"].mean() / 10)*root_depth_factor},
            "field_capacity_fraction": {joined_rasters["field_cap"].mean() / 100},
            "porosity_fraction": {joined_rasters["porosity"].mean() / 100},
            "wilting_point_fraction": {(joined_rasters["field_cap"].mean()/100)  -
                                       (joined_rasters['awc'].mean()/100)},
            "water_holding_capacity": {(joined_rasters["field_cap"].mean() / 100) -
                                       ((joined_rasters["field_cap"].mean()/100)  -
                                        (joined_rasters['awc'].mean()/100))
                                       },
            "latitude": {basins.loc[[basin],'geometry'].to_crs('EPSG:4326').geometry.centroid.y.values[0]},
            "basin_area_total": {joined_rasters['area'].sum() / 10e5},
            "impervious_area_fraction": {joined_rasters['imp'].sum()/joined_rasters['area'].sum()*100},
            "channel_length_max": {2},
            "channel_velocity_avg": {10},
            "flow_initial": {0.1},
            "stream area": {(joined_rasters['snet_10m'].filled(0)>0).sum()*100 / 10e5},
            "lake_area": {0}, #still working
            "up_lake_area" : {0}, #still working
            "rip_area": {(joined_rasters['snet_10m'].filled(0)>0).sum()*100 / 10e5},  # stream_area + lake_area,
        }
    
        df = pd.DataFrame.from_dict(characteristics_out, orient="index")
        df.index.name = "name"
        df.columns = ['value']
    
        twis[landuse] = calc_twi(joined_rasters["twi"])
        chars[landuse] = df

    return chars,twis

#%% Read Files
def files():
    path_to_wbdCatch = r"D:\WATER_FILES\final_2015oct26\USGS\water\water_db\drb\phys\wbdCatch.shp"
    path_to_watersheds = r"D:\WATER_FILES\BaseV4_2011aLULC\Watersheds.shp"
    
    # Read in the clim basins shapefile and convert it 
    basins = gpd.read_file(path_to_wbdCatch)
    basins = basins.to_crs(epsg=5070)
    basins.index = basins.HydroID

    # Larger Watersheds
    watersheds = gpd.read_file(path_to_watersheds)
    watersheds = watersheds.to_crs(epsg=5070)
    
    return basins,watersheds

#%% Main Program

def geo_main(climbasin, basins, watersheds,pst_node):
    
    # Path to basin files
    path_to_db = r'D:\WATER_FILES\final_2015oct26\USGS\water\water_db\drb\phys'
    path_to_masks = r"D:\WATER_FILES\BaseV4_2011aLULC"
    path_to_save = r"D:\WATER_FILES\inputs\Climbasin_inputs\basin_characteristics"
    landuse_path = r"D:\WATER_FILES\inputs\Climbasin_inputs\landuse_params"
    timeseries_path = r"D:\WATER_FILES\inputs\Climbasin_inputs\time_series"

    # Create dictionary of raster paths for use in the calculation of characteristics
    db_rasters = {'awc': os.path.join(path_to_db,'haawc100\w001001.adf'),
                  'con_mult': os.path.join(path_to_db,'haconmult100\w001001.adf'),
                  'field_cap': os.path.join(path_to_db,'hafc100\w001001.adf'),
                  'k_sat': os.path.join(path_to_db,'haksat100\w001001.adf'),
                  'scaling_parameter': os.path.join(path_to_db,'hammilli100\w001001.adf'),
                  'soil_thickness': os.path.join(path_to_db,'hadepth100\w001001.adf'),
                  'porosity': os.path.join(path_to_db,'hapor100\w001001.adf'),
                  'imp': os.path.join(path_to_db,'imp_2011\w001001.adf'),
                  'snet_10m': r"D:\WATER_FILES\final_2015oct26\USGS\water\water_db\drb\topo\snet\w001001.adf",
                  'twi': os.path.join(path_to_db,"twi\w001001.adf"),
                  'stream': r"D:\WATER_FILES\final_2015oct26\USGS\water\water_db\drb\topo\snet\w001001.adf"
                  }
    
    # Solve for the characteristics
    chars, twis = drb_characteristics(db_rasters, climbasin, path_to_masks, basins, watersheds)
    # params = drb_lu_params()
    
    # Dictionary of landuse
    lu_to_landuse = {'f':'forest','r':'developed','a':'agricultural'}
    
    # Create path to save the input files
    if os.path.exists(os.path.join(path_to_save,str(climbasin))):
        pass
    else:
        os.mkdir(os.path.join(path_to_save,str(climbasin)))
              
    # Loop through the landuse params and save files
    for lu in chars.keys():
        # Save the characteristics and the twi
        chars[lu].to_csv(os.path.join(path_to_save,str(climbasin),'basin_characteristics_'+lu_to_landuse[lu]+'.csv'))
        twis[lu].to_csv(os.path.join(path_to_save,str(climbasin),'twi_'+lu_to_landuse[lu]+'.csv'))
        
        # Copy the template for the modelconfig file
        template_ini = open(os.path.join(r'D:\WATER_FILES\inputs','modelconfig_climbasins_template.ini'))
        content = template_ini.read()
        template_ini.close()
        
        # Create the output path, insert it into the model config file
        output_path = os.path.join(path_to_save,str(climbasin))
        content = content.replace('path_to_input_directory',output_path)
        
        # Link to the location of landuse params
        content = content.replace('path_to_landuse_directory',landuse_path)
        
        # Create the output directories for the climbasin
        output_dir = output_path.replace('inputs','outputs')
        if not os.path.exists(os.path.join(output_dir,lu_to_landuse[lu])):
            os.makedirs(os.path.join(output_dir,lu_to_landuse[lu]))
        else:
            pass
        
        content = content.replace('lu',lu_to_landuse[lu])
        content = content.replace('path_to_output_directory',output_dir)
        
        # Link to the location of the climate timeseries data
        content = content.replace('path_to_timeseries',timeseries_path)
        content = content.replace('pstbasin',str(pst_node))
        
        # Save the new modelconfig file
        lu_ini = open(os.path.join(output_path,'modelconfigfile_'+lu_to_landuse[lu]+'.ini'),'w')
        lu_ini.write(content)
        lu_ini.close()
        
        
        
#%% Uncomment to test one basin
# basins,watersheds = files()
# geo_main(712,basins,watersheds,145)