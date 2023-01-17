# -*- coding: utf-8 -*-
"""
Created on Fri Jan  6 11:50:04 2023

@author: AnthonyPreucil
"""

# Run Water
from waterpy.main import waterpy
import pandas as pd
import numpy as np
import os

#%% Read/Parse XML

xmlfile = r"D:\WATER_FILES\BaseV4_2011aLULC\100A\WATERSimulation.xml"

# Note - may need to set the second param using a keywork, 'xpath' for future version of pandas
node_df = pd.read_xml(xmlfile,'Project/Study/StudySimulation/SimulationFeatures')

twi_df = pd.read_xml(xmlfile,'Project/Study/StudySimulation/SimulationTopographicWetnessIndex')
df_list = [df.reset_index(drop=True).add_suffix(f'_{i+1}') for i, (name, df) in enumerate(twi_df.groupby('SimulID'))]
node_twi_df = pd.concat(df_list, axis=1, join='outer')

temp_data = pd.read_xml(xmlfile,'Project/Study/StudySimulation/ClimaticTemperatureSeries')

precip_data = pd.read_xml(xmlfile,'Project/Study/StudySimulation/ClimaticPrecipitationSeries')

#%% Prepare Land Use and Basin Characteristics Params

# Notes: 1 = Forest, 2 = Agriculture, 3 = Developed

xml_to_py_path = r"D:\WATER_FILES\inputs\xml_waterpy_key.xlsx"
output_path = "D:\WATER_FILES\inputs\BaseV4_2011aULC_inputs"
xml_to_py_key = pd.read_excel(xml_to_py_path)
# drop missing values from the waterpy sample to help with mapping later
xml_to_py_key.dropna(subset=['Name From waterpy sample'], inplace=True)

path_to_inputs = "D:\WATER_FILES\inputs\inputs_from_waterpy"
# initialze empty dataframes for input params
basin_df = pd.DataFrame()
lu_df = pd.DataFrame()
ag_df = pd.DataFrame()
dev_df = pd.DataFrame()
for f in os.listdir(path_to_inputs):
    print (f)
    df = pd.read_csv(os.path.join(path_to_inputs,f))
    if 'forest' in f:
        print ('Forest Params')
        lu_df.index = df.name
        
        # use map function to map values of node_df to lu_df
        lu_df['xml'] = lu_df.index.map(xml_to_py_key.set_index('Name From waterpy sample')['Name From WATER Original XML'])

        lu_df['value'] = lu_df['xml'].map(node_df[node_df.SimulID==1].set_index('AttName')['AttMeanVal'])
        lu_df.drop('xml',axis=1).to_csv(os.path.join(output_path,f))
        
    elif 'agricultural' in f:
        print ('Agriculture Params')
        ag_df.index = df.name
        
        # use map function to map values of node_df to ag_df
        ag_df['xml'] = ag_df.index.map(xml_to_py_key.set_index('Name From waterpy sample')['Name From WATER Original XML'])

        ag_df['value'] = ag_df['xml'].map(node_df[node_df.SimulID==2].set_index('AttName')['AttMeanVal'])
        ag_df.drop('xml',axis=1).to_csv(os.path.join(output_path,f))
        
    elif 'developed' in f:
        print ('Developed Params')
        dev_df.index = df.name
        
        # use map function to map values of node_df to dev_df
        dev_df['xml'] = dev_df.index.map(xml_to_py_key.set_index('Name From waterpy sample')['Name From WATER Original XML'])

        dev_df['value'] = dev_df['xml'].map(node_df[node_df.SimulID==3].set_index('AttName')['AttMeanVal'])
        dev_df.drop('xml',axis=1).to_csv(os.path.join(output_path,f))
    elif 'characteristics' in f:
        print ('Basin Characteristics Params')
        basin_df.index = df.name
        
        # use map function to map values of node_df to basin_df
        basin_df['xml'] = basin_df.index.map(xml_to_py_key.set_index('Name From waterpy sample')['Name From WATER Original XML'])

        basin_df['value'] = basin_df['xml'].map(node_df[node_df.SimulID==1].set_index('AttName')['AttMeanVal'])
        basin_df.drop('xml',axis=1).to_csv(os.path.join(output_path,f[:-4]+'_forest.csv'))
        
        basin_df['xml'] = basin_df.index.map(xml_to_py_key.set_index('Name From waterpy sample')['Name From WATER Original XML'])

        basin_df['value'] = basin_df['xml'].map(node_df[node_df.SimulID==2].set_index('AttName')['AttMeanVal'])
        basin_df.drop('xml',axis=1).to_csv(os.path.join(output_path,f[:-4]+'_agricultural.csv'))
        
        basin_df['xml'] = basin_df.index.map(xml_to_py_key.set_index('Name From waterpy sample')['Name From WATER Original XML'])

        basin_df['value'] = basin_df['xml'].map(node_df[node_df.SimulID==3].set_index('AttName')['AttMeanVal'])
        basin_df.drop('xml',axis=1).to_csv(os.path.join(output_path,f[:-4]+'_developed.csv'))
    
    else:
        print ('file '+f+' does not have an input equivalent for waterpy')
    
#%% Set up twi 

twi_forest = node_twi_df.filter(regex='1$',axis=1)
twi_forest.columns = ['bin','simul','twi','proportion']
twi_forest = (twi_forest.set_index('bin',drop=True)
              .drop('simul',axis=1)
              .assign(cells = lambda x: x.proportion*100))
twi_forest.to_csv(os.path.join(output_path,'twi_forest.csv'))

twi_agricultural = node_twi_df.filter(regex='2$',axis=1)
twi_agricultural.columns = ['bin','simul','twi','proportion']
twi_agricultural = (twi_agricultural.set_index('bin',drop=True)
              .drop('simul',axis=1)
              .assign(cells = lambda x: x.proportion*100))
twi_agricultural.to_csv(os.path.join(output_path,'twi_agricultural.csv'))

twi_developed = node_twi_df.filter(regex='3$',axis=1)
twi_developed.columns = ['bin','simul','twi','proportion']
twi_developed = (twi_developed.set_index('bin',drop=True)
              .drop('simul',axis=1)
              .assign(cells = lambda x: x.proportion*100))
twi_developed.to_csv(os.path.join(output_path,'twi_developed.csv'))

#%% Set up the time series from climate data using xml data
# Note: will not need to pull from the xml for climate runs, just build the time series csvs from the data

# Forest
temp = temp_data[temp_data['SimulID']==1][['SeriesDate','SeriesValue']]
temp = temp.set_index(pd.to_datetime(temp['SeriesDate'], format='%Y-%m-%dT%H:%M:%S%z',utc=True).dt.date).drop(['SeriesDate'], axis=1)
temp.columns = ['temperature (celsius)']

precip = precip_data[precip_data['SimulID']==1][['SeriesDate','SeriesValue']]
precip = precip.set_index(pd.to_datetime(precip['SeriesDate'], format='%Y-%m-%dT%H:%M:%S%z',utc=True).dt.date).drop(['SeriesDate'], axis=1)
precip.columns = ['precipitation (mm/day)']

time_series = pd.concat([temp,precip],axis=1)
time_series.index = time_series.index.rename('date')

time_series['flow_observed (mm/day)'] = [np.nan]*len(time_series)
time_series.to_csv(os.path.join(output_path,'timeseries_forest.csv'))

# Agricultural
temp = temp_data[temp_data['SimulID']==2][['SeriesDate','SeriesValue']]
temp = temp.set_index(pd.to_datetime(temp['SeriesDate'], format='%Y-%m-%dT%H:%M:%S%z',utc=True).dt.date).drop(['SeriesDate'], axis=1)
temp.columns = ['temperature (celsius)']

precip = precip_data[precip_data['SimulID']==2][['SeriesDate','SeriesValue']]
precip = precip.set_index(pd.to_datetime(precip['SeriesDate'], format='%Y-%m-%dT%H:%M:%S%z',utc=True).dt.date).drop(['SeriesDate'], axis=1)
precip.columns = ['precipitation (mm/day)']

time_series = pd.concat([temp,precip],axis=1)
time_series.index = time_series.index.rename('date')

time_series['flow_observed (mm/day)'] = [np.nan]*len(time_series)
time_series.to_csv(os.path.join(output_path,'timeseries_agricultural.csv'))

# Developed
temp = temp_data[temp_data['SimulID']==3][['SeriesDate','SeriesValue']]
temp = temp.set_index(pd.to_datetime(temp['SeriesDate'], format='%Y-%m-%dT%H:%M:%S%z',utc=True).dt.date).drop(['SeriesDate'], axis=1)
temp.columns = ['temperature (celsius)']

precip = precip_data[precip_data['SimulID']==3][['SeriesDate','SeriesValue']]
precip = precip.set_index(pd.to_datetime(precip['SeriesDate'], format='%Y-%m-%dT%H:%M:%S%z',utc=True).dt.date).drop(['SeriesDate'], axis=1)
precip.columns = ['precipitation (mm/day)']

time_series = pd.concat([temp,precip],axis=1)
time_series.index = time_series.index.rename('date')

time_series['flow_observed (mm/day)'] = [np.nan]*len(time_series)
time_series.to_csv(os.path.join(output_path,'timeseries_developed.csv'))

#%% Set up model config files
# Note - set up one for each land use (01/17/2023)
output_dir = output_path.replace('inputs','outputs')
if not os.path.exists(output_dir):
    os.mkdir(output_dir)
else:
    pass

# Plan - run with files the way they are now, come back and fix var problems as they happen?

for lu in ['forest','agricultural','developed']:
    template_ini = open(os.path.join(r'D:\WATER_FILES\inputs','modelconfig_template.ini'))
    content = template_ini.read()
    template_ini.close()
    
    content = content.replace('path_to_input_directory',output_path)
    content = content.replace('lu',lu)
    content = content.replace('path_to_output_directory',output_dir)
    
    lu_ini = open(os.path.join(output_path,'modelconfigfile_'+lu+'.ini'),'w')
    lu_ini.write(content)
    lu_ini.close()
#%% Run WATER with current set up
'''
configfile = r"D:\waterpy\example\input\modelconfig.ini"

waterpy(configfile,None)
'''