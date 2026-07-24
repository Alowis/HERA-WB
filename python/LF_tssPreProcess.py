# -*- coding: utf-8 -*-
"""
Created on Tue Mar 14 18:30:37 2023

@author: tilloal
"""


import numpy as np
import os
import sys
from netCDF4 import Dataset
from matplotlib import pyplot as plt
#from hydroeval import *
import pandas as pd
from pandas import ExcelFile
from pandas.plotting import register_matplotlib_converters
register_matplotlib_converters()
from datetime import date
from datetime import datetime
from datetime import timedelta
import datetime
import json
import csv
import inspect

InitPath='D:/tilloal/Documents/LFRuns_utils/LFPostProcess/Diagnostics/'
os.chdir(InitPath)
os.getcwd()
from lf_readplot_tss import *

# 
# dir_path = os.path.dirname(os.path.realpath(__file__))
# # Change working directory
# os.chdir(dir_path)

# The following five entries can be summarized in a settings_plots.txt file ...

#SubCatchmentPath='/BGFS/COMMON/tilloal/lisflood-fd/out/'
YearPath='D:/tilloal/Documents/LFRuns_utils/data/tss'
os.chdir(YearPath)
#catchments=np.arange(1,6300) # ALL the GloFAS IDs


# Find whether it is leap year or not
def checkYear(year):
    return (((year % 4 == 0) and
             (year % 100 != 0)) or
             (year % 400 == 0))
  
    
Yrstart=1950
Yrend=2020
years=range(Yrstart,Yrend)
#total_num_steps = 14975 # this is the number of computational steps between forcings start and forcings end, be mindfull of the daily and 6hours options
#need to find a way to get this automatically from an input file
total_num_steps = 1460
#efas outlets
#total_num_gauges = 2743
#hybas out
total_num_gauges = 2992
plots_storage_folder = 'D:/tilloal/Documents/LFRuns_utils/Figures/'
suffix_fig_filename = 'tplot_' # optional deatils of the file name

# stationfile = YearPath+'stations_efas_meta.csv'
# stationdata = pd.read_csv(stationfile, sep=";", index_col=0,encoding= 'unicode_escape')
# stationdata.index=stationdata.index.values
# allcat=stationdata.index.values
# catchments=list(allcat)# selection of subcatchments IDs
# with open('matchcatch.csv', newline='') as file:
#     reader = csv.reader(file)
#     catchments = list(reader)[0]
   
   
xml='D:/tilloal/Documents/LFRuns_utils/data/tss/Settings/settings_sixhourlyLiteRuns_warmstart_1955.xml'
tss='D:/tilloal/Documents/LFRuns_utils/data/tss/HERA_SocCF/disWin/disWin_1955.tss'
# timexml=xml_timeinfo(xml)
# tssdata=read_tss(tss,xml)
outputfilenames = ['rainUpsX','snowUpsX','scovUps', 'snowMeltUpsX','frostUpsX', 'tAvgUpsX','theta1totalX','theta2totalX', 'theta3totalX',
                    'theta1totalX','qUzUpsX', 'qLZUpsX','surfaceRunoffUpsX','qLakeIn', 'qLakeOut','resfill','etUpsX','qresin',
                    'qresout','infUpsX','qLakeIn','qLakeOut','ActEvapo','disWin']
#outputfilenames = ['ActEvapo','disWin'] 
# plots: 1'rainUpsX', 2'snowUpsX', 3'snowMeltUpsX', 4'frostUps', 5'actEvapo', 6'theta1totalX', 7'theta2totalX', 8'theta3totalX , 11'qLzUpsX', 12'percUZLZUpsX', 13'dis', 14'frost' , 16'surfaceRunoffUpsX' , 17'gwLossUpsX', 18'lzUpsX']
num_var = len(outputfilenames)
#UpsXtss=np.zeros((total_num_steps*num_var))
ts="6h"
Gaugedat={}
ind=-1
yrlist=range(Yrstart,(Yrend)+1)
sfold="Settings"
sce="HERA_RWCF"
output_folder=os.path.join(YearPath,sce)
#%%
for outfn in outputfilenames:
  UpsXtss=pd.DataFrame()
  print("analysis file: "+ outfn)
  for yr in yrlist:
    print(str(yr))
    yr = pd.to_numeric(yr)  # Convert yr to a numeric value if it's not already
    yr = str(yr)
    path_subcatch = YearPath
    ofold = outfn
    tssfile=os.path.join(path_subcatch,sce,ofold,outfn + '_' + str(yr) + '.tss')
    if yr == '1950':
        xmlfile=os.path.join(path_subcatch,sfold,'settings_sixhourlyFullRuns_coldstart_' + str(yr) + '.xml')
    else:
        xmlfile=os.path.join(path_subcatch,sfold,'settings_sixhourlyLiteRuns_warmstart_' + str(yr) + '.xml')
    if os.path.exists(tssfile):
      print("outputs for year "+ str(yr) + " exists")  
      tssfile_d = read_tss(tssfile,xmlfile)
      UpsXtss =  pd.concat([tssfile_d,UpsXtss])
  UpsXtss.to_csv(os.path.join(output_folder, outfn + '_1951_2020.csv'))
        
