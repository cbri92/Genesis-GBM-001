# -*- coding: utf-8 -*-
"""
Created on Tue Feb 21 01:42:19 2023

@author: Caterina Brighi
"""

#%% Import functions 

import matplotlib.pyplot as plt
import matplotlib as mpl
import SimpleITK as sitk
import numpy as np
import pandas as pd
import datetime
import os
import glob
import gzip
import shutil
import xlsxwriter
import time
from scipy.stats.stats import pearsonr
import dicom2nifti
from multiprocessing.pool import ThreadPool
from functools import partial
from ImageAnalysisFunctions import *

def dice(pred, true, k = 1):
    intersection = np.sum(pred[true==k]) * 2.0
    dice = intersection / (np.sum(pred) + np.sum(true))
    return dice

#%% Set Working directory
        
data_supradir = 'Path to data directory' #Set working directory

subjs_path = [ f.path for f in os.scandir(data_supradir) if f.is_dir() ] #Create a list of the paths to the subjects directories
subjs_name = [ f.name for f in os.scandir(data_supradir) if f.is_dir() ] #Create a lisdt of subjects names
# subjs_name = ['AB001']

#Generate excel file containing voxelwise results
Results = pd.ExcelWriter(data_supradir +'Results_stats.xlsx')

#Create empty dataframes to populate as going through the loop
images =['FET', 'PSMA']
d={}
for image in images:
    d[image] = pd.DataFrame(columns=['Subject_ID', 'BTV Volume [mm3]', 'BTV Mean intensity [SUV]', 'BTV Std of Mean [SUV]', 'BTV Max intensity [SUV]', 'CTRL VOI Volume [mm3]', 'CTRL VOI Mean intensity [SUV]', 'CTRL VOI Std of Mean [SUV]', 'CTRL VOI Max intensity [SUV]', 'BTV TBR Mean', 'BTV TBR Std of Mean', 'BTV TBR Max'])

d['CE Tumour'] = pd.DataFrame(columns=['Subject_ID', 'CE Tumour volume [mm3]'])
d['Metrics'] = pd.DataFrame(columns=['Subject_ID', 'DICE PSMA-FET BTV', 'PSMA_VOI/FET_VOI'])

#%%Create a for loop to perform image analysis on each subject sequentially

for current in subjs_name:
    subj_dir = data_supradir+current
    subj_name = current
    
    #Create empty Dictionary with results of stats to populate as going through the loop
    Results_dict = {}
    
    #Set paths to subfolder    
    reg_imgs = subj_dir +'/Registered'
    ROIs = subj_dir +'/ROIs'

    #Read images
    PET_FET = sitk.ReadImage(reg_imgs +'/PET_FET_inFETCT.nii', sitk.sitkFloat32) #read PET FET image
    PET_TER = sitk.ReadImage(reg_imgs +'/PET_PSMA_inFETCT.nii', sitk.sitkFloat32) #read PET PSMA image   
    TBR_FET = sitk.ReadImage(reg_imgs +'/PET_FET_TBR_map.nii', sitk.sitkFloat32) #read PET FET TBR image
    TBR_TER = sitk.ReadImage(reg_imgs +'/PET_PSMA_TBR_map.nii', sitk.sitkFloat32) #read PET PSMA TBR image
    
    #Read VOIs
    FET_BTV = sitk.ReadImage(ROIs+'/PET_FET_BTV.nii') #Read FET BTV
    TER_BTV = sitk.ReadImage(ROIs+'/PET_PSMA_BTV.nii') #Read PSMA BTV
    FET_CTRL = sitk.ReadImage(ROIs+'/PET_FET_CTRL_VOI.nii') #Read FET CTRL VOI
    TER_CTRL = sitk.ReadImage(ROIs+'/PET_PSMA_CTRL_VOI.nii') #Read PSMA CTRL VOI
    if os.path.exists(ROIs+'/CEtum.nii'):
        CEtum = sitk.ReadImage(ROIs+'/CEtum.nii') #Read CE Tum VOI se esiste

#%% Calculate stats for FET images
    
    #Calculate BTV stats 
    BTV_stats = getStatsRoi(FET_BTV, PET_FET) #Calculate BTV stats
    BTV_MeanSUV = BTV_stats.get('Mean intensity [SUV]') #Calculate BTV mean value of SUV
    BTV_StdSUV = BTV_stats.get('Std of Mean [SUV]') #Calculate BTV Std value of mean SUV
    BTV_MaxSUV = BTV_stats.get('Max intensity [SUV]') #Calculate BTV max value of SUV
    BTV_volume = int(BTV_stats.get('Volume [mm3]')) #Calculate BTV volume in mm3
                
    #Calculate CTRL VOI stats
    CTRL_VOI_stats = getStatsRoi(FET_CTRL, PET_FET) #Calculate new CTRL_VOI stats
    CTRL_VOI_MeanSUV = CTRL_VOI_stats.get('Mean intensity [SUV]') #Calculate CTRL_VOI mean value of SUV
    CTRL_VOI_StdSUV = CTRL_VOI_stats.get('Std of Mean [SUV]') #Calculate CTRL_VOI Std value of mean SUV
    CTRL_VOI_MaxSUV = CTRL_VOI_stats.get('Max intensity [SUV]') #Calculate CTRL_VOI max value of SUV
    CTRL_VOI_volume = int(CTRL_VOI_stats.get('Volume [mm3]')) #Calculate CTRL_VOI volume in mm3
                
    #Calculate BTV TBR stats
    BTV_TBR_stats = getStatsRoi(FET_BTV, TBR_FET)
    BTV_MeanTBR = BTV_TBR_stats.get('Mean intensity [SUV]') #Calculate BTV TBR mean value
    BTV_StdTBR = BTV_TBR_stats.get('Std of Mean [SUV]') #Calculate BTV Std value of mean TBR
    BTV_MaxTBR = BTV_TBR_stats.get('Max intensity [SUV]') #Calculate BTV max value of SUV
    
    #Save statistics of various VOIs            
    sheet_df = {'Subject_ID': subj_name, 'BTV Volume [mm3]': BTV_volume, 'BTV Mean intensity [SUV]': BTV_MeanSUV, 'BTV Std of Mean [SUV]': BTV_StdSUV, 'BTV Max intensity [SUV]':BTV_MaxSUV, 'CTRL VOI Volume [mm3]': CTRL_VOI_volume, 'CTRL VOI Mean intensity [SUV]': CTRL_VOI_MeanSUV, 'CTRL VOI Std of Mean [SUV]': CTRL_VOI_StdSUV, 'CTRL VOI Max intensity [SUV]': CTRL_VOI_MaxSUV,'BTV TBR Mean': BTV_MeanTBR, 'BTV TBR Std of Mean': BTV_StdTBR, 'BTV TBR Max':BTV_MaxTBR}
    Results_dict['FET'] = pd.DataFrame(sheet_df, index=[0], columns=['Subject_ID', 'BTV Volume [mm3]', 'BTV Mean intensity [SUV]', 'BTV Std of Mean [SUV]', 'BTV Max intensity [SUV]', 'CTRL VOI Volume [mm3]', 'CTRL VOI Mean intensity [SUV]', 'CTRL VOI Std of Mean [SUV]', 'CTRL VOI Max intensity [SUV]','BTV TBR Mean', 'BTV TBR Std of Mean', 'BTV TBR Max']) 
    d['FET'] = d['FET'].append(Results_dict['FET'], ignore_index=True)
                
#%% Calculate stats for PSMA images
    
    #Calculate BTV stats 
    BTV_stats = getStatsRoi(TER_BTV, PET_TER) #Calculate BTV stats
    BTV_MeanSUV = BTV_stats.get('Mean intensity [SUV]') #Calculate BTV mean value of SUV
    BTV_StdSUV = BTV_stats.get('Std of Mean [SUV]') #Calculate BTV Std value of mean SUV
    BTV_MaxSUV = BTV_stats.get('Max intensity [SUV]') #Calculate BTV max value of SUV
    t_BTV_volume = int(BTV_stats.get('Volume [mm3]')) #Calculate BTV volume in mm3
                
    #Calculate CTRL VOI stats
    CTRL_VOI_stats = getStatsRoi(TER_CTRL, PET_TER) #Calculate new CTRL_VOI stats
    CTRL_VOI_MeanSUV = CTRL_VOI_stats.get('Mean intensity [SUV]') #Calculate CTRL_VOI mean value of SUV
    CTRL_VOI_StdSUV = CTRL_VOI_stats.get('Std of Mean [SUV]') #Calculate CTRL_VOI Std value of mean SUV
    CTRL_VOI_MaxSUV = CTRL_VOI_stats.get('Max intensity [SUV]') #Calculate CTRL_VOI max value of SUV
    CTRL_VOI_volume = int(CTRL_VOI_stats.get('Volume [mm3]')) #Calculate CTRL_VOI volume in mm3
                
    #Calculate BTV TBR stats
    BTV_TBR_stats = getStatsRoi(TER_BTV, TBR_TER)
    BTV_MeanTBR = BTV_TBR_stats.get('Mean intensity [SUV]') #Calculate BTV TBR mean value
    BTV_StdTBR = BTV_TBR_stats.get('Std of Mean [SUV]') #Calculate BTV Std value of mean TBR
    BTV_MaxTBR = BTV_TBR_stats.get('Max intensity [SUV]') #Calculate BTV max value of SUV
    
    #Save statistics of various VOIs            
    sheet_df = {'Subject_ID': subj_name, 'BTV Volume [mm3]': t_BTV_volume, 'BTV Mean intensity [SUV]': BTV_MeanSUV, 'BTV Std of Mean [SUV]': BTV_StdSUV, 'BTV Max intensity [SUV]':BTV_MaxSUV, 'CTRL VOI Volume [mm3]': CTRL_VOI_volume, 'CTRL VOI Mean intensity [SUV]': CTRL_VOI_MeanSUV, 'CTRL VOI Std of Mean [SUV]': CTRL_VOI_StdSUV, 'CTRL VOI Max intensity [SUV]': CTRL_VOI_MaxSUV,'BTV TBR Mean': BTV_MeanTBR, 'BTV TBR Std of Mean': BTV_StdTBR, 'BTV TBR Max':BTV_MaxTBR}
    Results_dict['PSMA'] = pd.DataFrame(sheet_df, index=[0], columns=['Subject_ID', 'BTV Volume [mm3]', 'BTV Mean intensity [SUV]', 'BTV Std of Mean [SUV]', 'BTV Max intensity [SUV]', 'CTRL VOI Volume [mm3]', 'CTRL VOI Mean intensity [SUV]', 'CTRL VOI Std of Mean [SUV]', 'CTRL VOI Max intensity [SUV]','BTV TBR Mean', 'BTV TBR Std of Mean', 'BTV TBR Max']) 
    d['PSMA'] = d['PSMA'].append(Results_dict['PSMA'], ignore_index=True)   
 
#%% Calculate CE Tumour volume 

    if os.path.exists(ROIs+'/CEtum.nii'):
        CE_tum_stats = getStatsRoi(CEtum, PET_FET)
        CE_volume = int(CE_tum_stats.get('Volume [mm3]')) #Calculate CE tumour volume in mm3
        
        #Save metrics to results dataframe
        sheet_df = {'Subject_ID': subj_name, 'CE Tumour volume [mm3]': CE_volume}
        Results_dict['CE Tumour'] = pd.DataFrame(sheet_df, index=[0], columns=['Subject_ID', 'CE Tumour volume [mm3]']) 
        d['CE Tumour'] = d['CE Tumour'].append(Results_dict['CE Tumour'], ignore_index=True) 
        
#%% Calculate metrics of comparison

    mtx_TER_BTV = sitk.GetArrayFromImage(TER_BTV)
    mtx_FET_BTV = sitk.GetArrayFromImage(FET_BTV)

    dice_coeff = dice(mtx_TER_BTV, mtx_FET_BTV, k=1) #Calculate the dice similarity coefficient between the PSMA BTV and the FET BTV
    Vol_ratio = (t_BTV_volume/BTV_volume)*100 #Calculate the volumetric ratio between the PSMA BTV and the FET BTV
    
    #Save metrics to results dataframe
    sheet_df = {'Subject_ID': subj_name, 'DICE PSMA-FET BTV':dice_coeff, 'PSMA_VOI/FET_VOI':Vol_ratio}
    Results_dict['Metrics'] = pd.DataFrame(sheet_df, index=[0], columns=['Subject_ID', 'DICE PSMA-FET BTV', 'PSMA_VOI/FET_VOI']) 
    d['Metrics'] = d['Metrics'].append(Results_dict['Metrics'], ignore_index=True) 

#%%Sava data to excel spreadsheet    
for df_name, df in d.items():
    df.to_excel(Results, sheet_name=str(df_name), index=False)

Results.save()
