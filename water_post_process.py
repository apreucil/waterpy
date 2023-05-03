# -*- coding: utf-8 -*-
"""
Created on Wed Jan 18 15:33:08 2023

@author: AnthonyPreucil
"""

# Post Processing for WATER-DRB
import pandas as pd
import geopandas as gpd
import os
import re
from drb_geoprocessing import files
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np
#%% User set up
# basin = '225'

runall = True
if runall:
    subbasins = [folder for folder in os.listdir("D:\WATER_FILES\BaseV4_2011aLULC") if len(folder)<=4]
    subbasins.remove('255')
    subbasins = [b for b in subbasins if int(b[:3]) <=365] # To trenton only!

    # subbasins = ['100A','100B']
    grouped = {}
    for s in subbasins:
        key = s[:3]
        if key not in grouped:
            grouped[key] = []
        grouped[key].append(s)
    
    result = list(grouped.values())
pdf = PdfPages('output_results.pdf')
for sb in result:
    DRB_df = pd.DataFrame()
    water_df = pd.DataFrame()
    for b in sb:
        print ('Status: '+b)
        #%% waterpy
        # path_to_output = r"D:\WATER_FILES\outputs\BaseV4_2011aULC_outputs"
        path_to_output = r"D:\WATER_FILES\outputs\climbasin_outputs\basin_characteristics"
        
        # filename = os.path.join(r"D:\WATER_FILES\CliimBasinParamsFromWATERXML",str(basin),"WATER.txt")
        filename = os.path.join(r"D:\WATER_FILES\BaseV4_2011aLULC",str(b),"WATER.txt")
        climbasins = True
        
        if climbasins:
            basins,watersheds = files()
            basins['center'] = basins.centroid
            key = watersheds.sjoin(basins)
            node = key[key.PSTSubNode == b]
            climbasins = list(node[node.geometry.contains(node.center)].HydroID.values)
            pst_clim = gpd.GeoDataFrame(pd.DataFrame(node[node.geometry.contains(node.center)].HydroID,columns=['HydroID']).join(basins,on='HydroID',lsuffix='_'))

            fig1,ax = plt.subplots(figsize=(10,10))
            watersheds[watersheds.PSTSubNode==b].plot(facecolor='grey',edgecolor='k',linewidth=2,ax=ax,)
            pst_clim.plot(facecolor='b',edgecolor='k',alpha=0.4,ax=ax)
            node.center.plot(color='r',ax=ax)
            node[node.geometry.contains(node.center)].center.plot(color='w',ax=ax,marker='+',markersize=80)
            plt.title(b,fontsize=20)
            ax.get_xaxis().set_visible(False)
            ax.get_yaxis().set_visible(False)
            legend = plt.legend(['Climate Basin Center','Climate Basin Included in Sum'],frameon = 1,bbox_to_anchor=(0.5,0))
            frame = legend.get_frame()
            frame.set_facecolor('grey')
            
            pdf.savefig(fig1)
            # Manual List
            # climbasins = [613,642,656,667,672,676,680,681,683,685,705,712,718,721,722]
            climbasins = [str(c) for c in climbasins]
            basin_df = pd.DataFrame()
            for cb in climbasins:
                for lu in ['forest','agricultural','developed']:
                    # print (cb)
                    try:
                        lu_df = pd.read_csv(os.path.join(path_to_output,str(cb),lu,'output.csv'),index_col='date')
                        lu_df.index = pd.to_datetime(lu_df.index)
                        # basin_df = pd.concat([basin_df,lu_df['discharge_predicted average (cfs)']],axis=1)
                        basin_df = pd.concat([basin_df,lu_df['discharge_predicted']],axis=1)
        
                        
                    except:
                        print ('issue with '+str(cb)+ ': '+lu)
                    
                
            sb_df = pd.DataFrame()
            sb_df[b] = basin_df.sum(axis=1)
        
        else:
            basin_df = pd.DataFrame()
            for lu in ['forest','agricultural','developed']:
                print (b, lu)
                try:
                    lu_df = pd.read_csv(os.path.join(path_to_output,str(b),lu,'output.csv'),index_col='date')
                    # basin_df = pd.concat([basin_df,lu_df['discharge_predicted average (cfs)']],axis=1)
                    basin_df = pd.concat([basin_df,lu_df['discharge_predicted']],axis=1)
                except:
                    pass
                
            sb_df = pd.DataFrame()
            sb_df[b] = basin_df.sum(axis=1)
        
        sb_df.index = pd.to_datetime(sb_df.index)
        DRB_df = pd.concat([DRB_df,sb_df],axis=1)
        
        #%% OLD WATER
        file = open(filename)
        lines = file.readlines()
        file.close()
        
        data = [re.split("\t|\n",line) for line in lines[6:]]
        df = pd.DataFrame(data[1:],columns=data[0])
        df['Date'] = pd.to_datetime(df.Date)
        waterbasin_df = df.set_index('Date').apply(pd.to_numeric,errors='coerce')
        water_df = pd.concat([water_df,waterbasin_df['Discharge (cfs)']],axis=1)
        
        #%% PST
        
    basedata_path = r"C:\OASIS_DRBC\basedata\basedata_2017_v3.csv"
    pst_df = pd.read_csv(basedata_path,index_col=0)
    pst_df.index = pd.to_datetime(pst_df.index)
    
    pst_basin = pst_df[b[:3]]
        
        
    #%% Compare
    
    combine_df = pd.concat([DRB_df.sum(axis=1),water_df.sum(axis=1),pst_basin],axis=1).dropna()
    combine_df.columns = ['waterpy','WATER','PST']
    combine_df = combine_df[combine_df.index <= pd.to_datetime('20060601')]
    
    sums = pd.DataFrame(combine_df.sum())
    sums = sums.transpose()
    
    annual = pd.DataFrame(combine_df.resample('Y').mean().mean(),columns=['Annual'])
    monthly = combine_df.resample('M').mean()
    out_stats = annual
    for m in range(12):
        m+=1
        out_stats = pd.concat([out_stats,pd.DataFrame(monthly[monthly.index.month==m].resample('Y').mean().mean(),columns = [m])],axis=1)
    out_stats.to_excel('average_cfs_compare_'+b[:3]+'.xlsx')
    
    fig = plt.figure(figsize=(15, 4))
    ax = plt.Axes(fig, [0., 0., 1., 1.])
    ax.set_axis_off()
    fig.add_axes(ax)
    tab = plt.table(cellText=np.asarray(out_stats.round(1)), rowLabels = out_stats.index,colWidths=[0.1]*13, colLabels=out_stats.columns, cellLoc='center', 
                    loc='center',bbox=[0.1,0.1,0.8,0.75])
    tab.auto_set_font_size(False)
    tab.set_fontsize(12)
    tab.scale(0.5, 1.5)
    ax.text(0.25,0.9,'Table of Average CFS - node '+b[:3],fontsize=20,fontweight='bold')
    # plt.tight_layout()
    pdf.savefig(fig)
    
    #%% Display
    fig2,ax = plt.subplots(figsize=(10,7))
    combine_df.cumsum().plot(ax=ax)
    plt.title('PST basin: '+b[:3])
    plt.legend()
    ax.table(cellText=sums.round().values,rowLabels=['sum'],colLabels = sums.columns,
             bbox=[0.2, -0.2, 0.6, 0.1])
    plt.ylabel('Cumulative Flow (CFS)')
    plt.tight_layout()
    # plt.show()
    #%%
    fig3,ax = plt.subplots(figsize=(8,8))
    combine_df.describe([p for p in np.linspace(0,1,21)])[4:-2].plot(ax=ax)
    plt.title('Cumulative Distribtion, node: '+b[:3])
    plt.grid()
    plt.legend(fontsize=12)
    plt.ylabel('CFS')
    pdf.savefig(fig2)
    pdf.savefig(fig3)
pdf.close()
