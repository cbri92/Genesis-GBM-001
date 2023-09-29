# -*- coding: utf-8 -*-
"""
Created on Tue Feb 28 16:28:17 2023

@author: cbri3325
"""

import SimpleITK as sitk
import pandas as pd
import os
from ImageAnalysisFunctions import getMaxRoi,getMeanRoi

data_supradir = 'Path to data directory'

subjs_path = [ f.path for f in os.scandir(data_supradir) if f.is_dir() ] #Create a list of the paths to the subjects directories
subjs_name = [ f.name for f in os.scandir(data_supradir) if f.is_dir() ] #Create a lisdt of subjects names

n_subj = len(subjs_name) #Total number of subjects

results = pd.ExcelWriter(data_supradir +'TSG_results.xlsx')
res_df = pd.DataFrame(columns=['Subject ID', 'Salivary Glands PSMA mean', 'Salivary Glands PSMA max', 'Tumour PSMA mean', 'Tumour PSMA max', 'TSG PSMA mean', 'TSG PSMA max'])

#%%Create a for loop to perform image analysis on each subject sequentially

for current in subjs_name:
    subj_dir = data_supradir+current
    subj_name = current
    
    PSMA_PET = sitk.ReadImage(subj_dir +'/Registered/PET_PSMA_inFETCT.nii') #read PET PSMA image
    Tumour = sitk.ReadImage(subj_dir +'/ROIs/PET_PSMA_BTV.nii') #read tumour roi
    
    Tumour_Mean = getMeanRoi(Tumour, PSMA_PET)
    Tumour_Max = getMaxRoi(Tumour, PSMA_PET)
    
    if os.path.exists(subj_dir +'/ROIs/SG_VOI.nii'):
        SalivaryGlands = sitk.ReadImage(subj_dir +'/ROIs/SG_VOI.nii') #read salivary glands roi
    
        SG_Mean = getMeanRoi(SalivaryGlands, PSMA_PET)
        SG_Max = getMaxRoi(SalivaryGlands, PSMA_PET)
        
        TSG_Mean = Tumour_Mean/SG_Mean
        TSG_Max = Tumour_Max/SG_Max
    
    else:
        SG_Mean = 0
        SG_Max = 0
        
        TSG_Mean = 0
        TSG_Max = 0
        
    
    df = {'Subject ID':subj_name, 'Salivary Glands PSMA mean':SG_Mean, 'Salivary Glands PSMA max':SG_Max,'Tumour PSMA mean':Tumour_Mean, 'Tumour PSMA max':Tumour_Max, 'TSG PSMA mean':TSG_Mean, 'TSG PSMA max':TSG_Max}
    res_df = res_df.append(df, ignore_index=True)
    
res_df.to_excel(results, sheet_name='TSG', index=False)
results.save()
    
