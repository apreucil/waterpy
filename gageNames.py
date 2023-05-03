# -*- coding: utf-8 -*-
"""
Created on Mon Mar 27 10:53:26 2023

@author: AnthonyPreucil
"""

# USGS ID to Name

import requests
import padnas as pd
import re

def GetNameFromIDs(gage_ids):

    url = f"https://waterservices.usgs.gov/nwis/site/?format=rdb&sites={gage_ids}"
    response = requests.get(url)
    
    return response
    
    
        
#%% main program

url = r'https://waterservices.usgs.gov/nwis/dv/?format=rdb&bBox=-76.400000,38.800000,-74.400000,42.300000&siteStatus=all'

response = requests.get(url)
content = response.text

meta_re = r'#\s{4}USGS\s(\d*)\s(.*)'
results = re.findall(meta_re,content)
meta_df = pd.DataFrame(results,columns = ['ID','GageName'])
meta_df['GageName'] = [re.sub(',','',gn.upper()) for gn in meta_df.GageName]

att = pd.read_excel(r"D:\WATER_FILES\GIS_Working_Attributes.xlsx")
att['Gage Name'] = [re.sub(',','',gn.upper()) for gn in att['Gage Name']]
join_df = pd.merge(att,meta_df,left_on='Gage Name',right_on='GageName')

# gage_list = ['0'+str(g) for g in att['Gage Number']]

gage_list = ','.join(join_df.ID)

response = GetNameFromIDs(gage_list)
content = response.text

info_re = r'USGS\t(\d*)\t([^\t]*)\tST\t([\d\.]*)\t([-\d\.]*)'
results = re.findall(info_re,content)

df = pd.DataFrame(results)
df.columns = ['ID','Name','Lat','Lon']
df['Name'] = [re.sub(',','',gn.upper()) for gn in df.Name]

final_df = pd.merge(join_df,df,left_on='GageName',right_on='Name')

final_df.to_clipboard()




