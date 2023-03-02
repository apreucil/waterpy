# -*- coding: utf-8 -*-
"""
Created on Wed Jan 18 11:36:42 2023

@author: AnthonyPreucil
"""
import os
import sys
sys.path.append(r'D:\waterpy\waterpy')
from main import waterpy

def main(input_path,lu):
    print("Run Water, land use: "+lu)
    # Run WATER with current set up

    configfile = os.path.join(input_path,"modelconfigfile_"+lu+".ini")
    print (configfile)
    waterpy(configfile,None)
    
