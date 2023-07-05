# -*- coding: utf-8 -*-
"""
Created on Mon Feb  6 21:01:02 2023

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

#%% Set Working directory
        
data_supradir = 'Path to data directory' #Set working directory

subjs_path = [ f.path for f in os.scandir(data_supradir) if f.is_dir() ] #Create a list of the paths to the subjects directories
subjs_name = [ f.name for f in os.scandir(data_supradir) if f.is_dir() ] #Create a lisdt of subjects names
# subjs_name = ['KF006']

#Set path to seeds list excel file
seeds_file_path = data_supradir +'Seeds_list.xlsx'
seeds_dfs = pd.read_excel(seeds_file_path, sheet_name=None, index_col='Patient ID')
n_sheets = len(seeds_dfs)


#%%Create a for loop to perform image analysis on each subject sequentially

for current in subjs_name:
    subj_dir = data_supradir+current
    subj_name = current
      
    #Set paths to subfolder    
    reg_imgs = subj_dir +'/Registered'
    ROIs = subj_dir +'/ROIs'
    

    #Read images
    PET_FET = sitk.ReadImage(reg_imgs +'/PET_FET_inFETCT.nii', sitk.sitkFloat32) #read PET FET image
    PET_TER = sitk.ReadImage(reg_imgs +'/PET_PSMA_inFETCT.nii', sitk.sitkFloat32) #read PET PSMA image
    brain_mask = sitk.ReadImage(reg_imgs +'/brain_mask.nii') #read brain mask
    
    PET_FET = generate_mask(PET_FET, brain_mask)
    PET_TER = generate_mask(PET_TER, brain_mask)
    
    # sitk.WriteImage(PET_FET, reg_imgs+'/FET_PET_bet.nii')
    # sitk.WriteImage(PET_TER, reg_imgs+'/FET_PSMA_bet.nii')
    
    images = [PET_FET, PET_TER]   
    
    print('Start automated segmentation for '+subj_name)
        
    start = time.time()
    
    for df,img in zip(seeds_dfs,images):
        
        if df == 'FET PET':
            PET_image = 'PET_FET'      
            
        elif df == 'PSMA PET':
            PET_image = 'PET_PSMA'
            
        BTV = sitk.ReadImage(ROIs+'/'+PET_image+'_BTV.nii')
                
        #Generate CTRL VOI
        CTRL_VOI = flip(BTV)
        CTRL_VOI = CTRL_VOI-BTV #Remove potential areas of overlay between BTV and CTRL VOI
        CTRL_VOI = sitk.BinaryThreshold(CTRL_VOI, 1, 1, 1, 0)    
        CTRL_VOI_stats = getStatsRoi(CTRL_VOI, img) 
        CTRL_VOI_MeanSUV = CTRL_VOI_stats.get('Mean intensity [SUV]') #Calculate CTRL VOI mean value of SUV
        sitk.WriteImage(CTRL_VOI, ROIs+'/'+PET_image+'_CTRL_VOI.nii')
            
        #Generate PET TBR map
        PET_TBR_map = img/CTRL_VOI_MeanSUV 
        sitk.WriteImage(PET_TBR_map, reg_imgs+'/'+PET_image+'_TBR_map.nii')      
 
    
