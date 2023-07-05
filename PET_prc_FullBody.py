# -*- coding: utf-8 -*-
"""
Created on Tue Feb 28 16:28:17 2023

@author: Caterina Brighi
"""

import SimpleITK as sitk
import pandas as pd
import os
from ImageAnalysisFunctions import getMaxRoi,getMeanRoi

data_supradir = 'Path to data directory'

subjs_path = [ f.path for f in os.scandir(data_supradir) if f.is_dir() ] #Create a list of the paths to the subjects directories
subjs_name = [ f.name for f in os.scandir(data_supradir) if f.is_dir() ] #Create a lisdt of subjects names

n_subj = len(subjs_name) #Total number of subjects

results = pd.ExcelWriter(data_supradir +'TLR_results.xlsx')
res_df = pd.DataFrame(columns=['Subject ID', 'Liver PSMA mean', 'Liver PSMA max', 'Tumour PSMA mean', 'Tumour PSMA max', 'TLR PSMA mean', 'TLR PSMA max'])

#%%Create a for loop to perform image analysis on each subject sequentially

for current in subjs_name:
    subj_dir = data_supradir+current
    subj_name = current
    
    # InfoFile=pd.read_excel(subj_dir +'/'+'InfoFile_Subject_'+subj_name+'.xlsx', index_col=0)

    # #Multiply PET raw images by conversion factor
    # GaconvFact = InfoFile.loc["Units conversion factor", "Ga"]
    
    # PET_PSMA_raw = sitk.ReadImage(subj_dir +'/'+ subj_name +'_PSMA_PET_raw.nii') #read PET PSMA raw image
    # PET_PSMA_conv = sitk.ShiftScale(PET_PSMA_raw, shift = 0, scale = GaconvFact) #multiply PET PSMA image by Ga conversion factor
    # sitk.WriteImage(PET_PSMA_conv, subj_dir +'/'+ subj_name +'_PET_PSMA_FullBody.nii') #write converted image as a nifti file PET PSMA
    
    PSMA_PET = sitk.ReadImage(subj_dir +'/'+ subj_name +'_PET_PSMA_FullBody.nii') #read PET PSMA image
    Tumour = sitk.ReadImage(subj_dir +'/'+ subj_name +'_Tumour.nii') #read tumour roi
    Liver = sitk.ReadImage(subj_dir +'/'+ subj_name +'_Liver.nii') #read tumour roi
    
    Liver_Mean = getMeanRoi(Liver, PSMA_PET)
    Liver_Max = getMaxRoi(Liver, PSMA_PET)
    Tumour_Mean = getMeanRoi(Tumour, PSMA_PET)
    Tumour_Max = getMaxRoi(Tumour, PSMA_PET)
    
    TLR_Mean = Tumour_Mean/Liver_Mean
    TLR_Max = Tumour_Max/Liver_Max
    
    df = {'Subject ID':subj_name, 'Liver PSMA mean':Liver_Mean, 'Liver PSMA max':Liver_Max, 'Tumour PSMA mean':Tumour_Mean, 'Tumour PSMA max':Tumour_Max, 'TLR PSMA mean':TLR_Mean, 'TLR PSMA max':TLR_Max}
    res_df = res_df.append(df, ignore_index=True)
    
res_df.to_excel(results, sheet_name='TLR', index=False)
results.save()
    
