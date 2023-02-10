# -*- coding: utf-8 -*-
"""
Created on Wed Jan 25 08:06:19 2023

@author: AnthonyPreucil
"""

# Set up shapefiles

import numpy as np
import geopandas as gpd
import os
import sys
sys.path.append('D:/waterpy/waterpy')
from geospatial import *
import rioxarray as riox
import rasterio
from rasterio.features import geometry_mask
from rasterio.plot import show
from rasterio.features import shapes
from shapely.geometry import shape

#%% Path to climate basins
path_to_wbdCatch= r"D:\WATER_FILES\final_2015oct26\USGS\water\water_db\drb\phys\wbdCatch.shp"
#%% Read in the clim basins shapefile and convert it 
basins = gpd.read_file(path_to_wbdCatch)
basins = basins.to_crs(epsg=5070)
#%% This demostrates the clipping process for a land use
mask_path = r"D:\WATER_FILES\BaseV4_2011aLULC\170\rmask\w001001.adf"
data_path = r'D:\WATER_FILES\final_2015oct26\USGS\water\water_db\drb\phys\haawc100\w001001.adf'

# Read rasters using rioxarray
raster = riox.open_rasterio(mask_path,masked=True)
raster_data = riox.open_rasterio(data_path,masked=True)

# Use polygon in clip method of rioxarray object to clip rasters to climate basin
clipped_mask = raster.rio.clip(basins.loc[[845],'geometry'].geometry)
clipped_data = raster_data.rio.clip(basins.loc[[845],'geometry'].geometry)

mean_lu = np.nanmean(clipped_data.data[clipped_mask.data == 1])


#%%
def calc_area(masked_arr,transform):
    raster_shape = shapes(masked_arr,transform = transform)
    tot_area = 0
    for myshape,value in raster_shape:
        polygon = shape(myshape)
        area = polygon.area
        tot_area += area
    return tot_area
#%%
def drb_characteristics(db_rasters, shp, basin):
    for landuse in ['f','r','a']:
        mask_path = os.path.join(r"D:\WATER_FILES\BaseV4_2011aLULC\170",landuse+"mask\w001001.adf")

        db_rasters_clip = {}
        for k, v in db_rasters.items():
            with rasterio.open(v) as src:
                basin_geometry = basins.loc[[basin],'geometry'].geometry.iloc[0]
                col_off, row_off, width, height = src.window(*basin_geometry.bounds).flatten()
                col_start = int(np.round(col_off))
                row_start = int(np.round(row_off))
                col_stop = col_start + int(np.round(width))
                row_stop = row_start + int(np.round(height))                
                
                window = ((row_start, row_stop), (col_start, col_stop))
                data = src.read(1, masked=True, window=window)
                show(data)
                mask = geometry_mask([basin_geometry], data.shape, transform=src.window_transform(window), invert=True)
                data = np.ma.masked_array(data, mask=~mask)
                show(data,title='Before Mask '+k)
                
                with rasterio.open(mask_path) as landuse_src:
                    col_off, row_off, lwidth, lheight = landuse_src.window(*basin_geometry.bounds).flatten()
                    col_start = int(np.round(col_off))
                    row_start = int(np.round(row_off))
                    col_stop = col_start + int(np.round(lwidth))
                    row_stop = row_start + int(np.round(lheight))
                    
                    window = ((row_start, row_stop), (col_start, col_stop))
                    landuse_mask = landuse_src.read(1, masked=True, window=window)
                    show(landuse_mask,title='The mask')

                    if data.shape == landuse_mask.shape:
                        lu_data = np.ma.masked_where(landuse_mask!=1,data)
                        show(lu_data,title='After Mask '+k)
                    else:
                        correct_res_data = st.resize(data,(590,590),mode='constant') 
                        show(correct_res_data,title='Resolution Corrected for '+k)
                        lu_data = np.ma.masked_where(landuse_mask!=1,correct_res_data)
                        show(lu_data,title='After Mask '+k)
                    
                    lu_basin = np.ma.masked_array(landuse_mask,mask=~mask)
                    lu_area = calc_area(lu_basin,landuse_src.transform)
                    
                db_rasters_clip[k] = lu_data

        characteristics_out = {
            "scaling_parameter": {db_rasters_clip['scaling_parameter'].mean() / 100},
            "saturated_hydraulic_conductivity": {db_rasters_clip['k_sat'].mean() / 100 * 86.4},
            "saturated_hydraulic_conductivity_multiplier": {db_rasters_clip['con_mult'].mean() / 100},
            "soil_depth_total": {db_rasters_clip["soil_thickness"].mean() / 10},
            "field_capacity_fraction": {db_rasters_clip["field_cap"].mean() / 10000},
            "porosity_fraction": {db_rasters_clip["porosity"].mean() / 10000},
            "wilting_point_fraction": {(db_rasters_clip["field_cap"].mean() / 10000) -
                                       (db_rasters_clip['awc'].mean() / 100)},
            "latitude": {basins.loc[[basin],'geometry'].to_crs('EPSG:4326').geometry.centroid.y.values[0]},
            "basin_area_total": {lu_area / 10e6},
            "impervious_area_fraction": {(db_rasters_clip['imp'].filled(np.nan)>0).sum()*100/lu_area*100},
            "channel_length_max": {2},
            "channel_velocity_avg": {10},
            "flow_initial": {0.1},
            "stream area": {(db_rasters_clip['snet_10m'].filled(0)>0).sum()*100 / 10e5},
            "lake_area": {0}, #still working
            "up_lake_area" : {0}, #still working
            "rip_area": {(db_rasters_clip['snet_10m'].filled(0)>0).sum()*100 / 10e5},  # stream_area + lake_area,
            "lake_delay": {0},
            "eff_imp": {0.7},
            "imp_delay": {0.1},
            "twi_adj": {1},
            "et_exp_dorm": {0.5},
            "et_exp_grow": {0.5},
            "grow_trigger": {15},
        }
        units = ["mm", "mm/day","unitless","mm","fraction","fraction","fraction","degrees", "sq km", "percentage",
                 "km", "km/day", "mm/day", "sq km", "sq km", "sq km", "sq km", "days", "fraction", "days", "unitless",
                 "unitless", "unitless", "temp C"]
    
        description = ['controls the rate of decline of transmissivity in the soil profile',
                       'saturated hydraulic conductivity of the C horizon of the soil',
                       'multiplier to apply to saturated hydraulic conductivity ', 'soil depth',
                       'fraction of soil moisture or water content in the soil after excess water has drained away',
                       'fraction of soil that is porous and is always larger than field_capacity_fraction',
                       'fraction amount of the minimal amount of water in the soil that plants require not to wilt',
                       'centroid latitude of basin', 'total basin area', 'fraction of impervious area of basin',
                       'maximum channel length', 'average channel velocity', 'initial river flow',
                       'total stream surface area', 'total waterbody area', 'total waterbody area upstream',
                       'total riparian area', 'estimated time for water to move through lake',
                       'percentage of impervious area connection to stream network',
                       'estimated delay for impervious runoff to reach stream network',
                       'Adjustment for magnitude of TWI - must be >= 1.',
                       'evapotranspiration Exponent for non-growing season.',
                       'evapotranspiration Exponent for growing season.',
                       'Temperature (C) transition to/from growing season for ET Exp and AMC.']
    
        df = pd.DataFrame.from_dict(characteristics_out, orient="index")
        df.index.name = "name"
        df.columns = ['value']
        df['units'] = units
        df['description'] = description

    return df




#%% Code from geospatial.py
path_to_db = r'D:\WATER_FILES\final_2015oct26\USGS\water\water_db\drb\phys'
# karst_raster = Raster(path="database//sinks.tif")
# karst_shp = ''
shp = Shp(path=r"D:\WATER_FILES\final_2015oct26\USGS\water\water_db\drb\shape_files\basin_near_hawley.shp")


# db_rasters = {'awc': riox.open_rasterio(os.path.join(path_to_db,'haawc100\w001001.adf'),masked=True),
#               'con_mult': riox.open_rasterio(os.path.join(path_to_db,'haconmult100\w001001.adf'),masked=True),
#               'field_cap': riox.open_rasterio(os.path.join(path_to_db,'hafc100\w001001.adf'),masked=True),
#               'k_sat': riox.open_rasterio(os.path.join(path_to_db,'haksat100\w001001.adf'),masked=True),
#               'scaling_parameter': riox.open_rasterio(os.path.join(path_to_db,'hammilli100\w001001.adf'),masked=True), # ??? hammilli100?
#               'soil_thickness': riox.open_rasterio(os.path.join(path_to_db,'hadepth100\w001001.adf'),masked=True),
#               'porosity': riox.open_rasterio(os.path.join(path_to_db,'hapor100\w001001.adf'),masked=True),
#               'imp': riox.open_rasterio(os.path.join(path_to_db,'imp\w001001.adf'),masked=True),
#               'snet_10m': riox.open_rasterio(r"D:\WATER_FILES\final_2015oct26\USGS\water\water_db\drb\topo\snet\w001001.adf",masked=True),
#               'twi': riox.open_rasterio(os.path.join(path_to_db,"twi\w001001.adf"),masked=True),
#               'stream': riox.open_rasterio(r"D:\WATER_FILES\final_2015oct26\USGS\water\water_db\drb\topo\snet\w001001.adf",masked=True)
#               }
db_rasters = {'awc': os.path.join(path_to_db,'haawc100\w001001.adf'),
              'con_mult': os.path.join(path_to_db,'haconmult100\w001001.adf'),
              'field_cap': os.path.join(path_to_db,'hafc100\w001001.adf'),
              'k_sat': os.path.join(path_to_db,'haksat100\w001001.adf'),
              'scaling_parameter': os.path.join(path_to_db,'hammilli100\w001001.adf'), # ??? hammilli100?
              'soil_thickness': os.path.join(path_to_db,'hadepth100\w001001.adf'),
              'porosity': os.path.join(path_to_db,'hapor100\w001001.adf'),
              'imp': os.path.join(path_to_db,'imp_2011\w001001.adf'),
              'snet_10m': r"D:\WATER_FILES\final_2015oct26\USGS\water\water_db\drb\topo\snet\w001001.adf",
              'twi': os.path.join(path_to_db,"twi\w001001.adf"),
              'stream': r"D:\WATER_FILES\final_2015oct26\USGS\water\water_db\drb\topo\snet\w001001.adf"
              }
# shp.karst_flag = karst_detection(karst_raster, shp)

# out_df = characteristics(db_rasters, shp)
# out_twi = twi_bins(db_rasters["twi"], shp)