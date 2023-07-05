# -*- coding: utf-8 -*-
"""
Created on Mon Feb  6 21:01:02 2023

@author: Caterina Brighi
"""

#%% Import functions 

import SimpleITK as sitk
import numpy as np
import os
from ImageAnalysisFunctions import *

#%% Set Working directory
        
data_supradir = 'Path to data directory' #Set working directory

subjs_path = [ f.path for f in os.scandir(data_supradir) if f.is_dir() ] #Create a list of the paths to the subjects directories
subjs_name = [ f.name for f in os.scandir(data_supradir) if f.is_dir() ] #Create a lisdt of subjects names

#%%Create a for loop to perform image analysis on each subject sequentially

for current in subjs_name:
    subj_dir = data_supradir+current
    subj_name = current
    
    #Set paths to subfolders    
    orig_imgs = subj_dir +'/Original'
    reg_imgs = subj_dir +'/Registered'
    reg_mtrx = subj_dir +'/Transforms'
    ROIs = subj_dir +'/ROIs'

    #Read Original images
    PET_FET = sitk.ReadImage(orig_imgs +'/PET_FET.nii', sitk.sitkFloat32) #read PET FET image
    CT_FET = sitk.ReadImage(orig_imgs +'/CT_FET.nii', sitk.sitkFloat32) #read CT FET image
    PET_TER = sitk.ReadImage(orig_imgs +'/PET_PSMA.nii', sitk.sitkFloat32) #read PET TER image
    if os.path.exists(orig_imgs +'/CT_PSMA.nii'):
        CT_TER = sitk.ReadImage(orig_imgs +'/CT_PSMA.nii', sitk.sitkFloat32) #read CT TER image
    
    
    #%%Registration and Resampling to FET_CT image and resolution
    print('Start registration for '+current)
    
    #Resample PET FET image into CT FET space
    PET_FET_inFETCT = Resample_image(PET_FET, CT_FET)
        
    #Register CT TER to CT FET
    tfm = sitk.ReadTransform(reg_mtrx+'/tfm.txt')

    if os.path.exists(orig_imgs +'/CT_PSMA.nii'):
        CT_TER_inFETCT = sitk.Resample(CT_TER, CT_FET, tfm, sitk.sitkLinear)
    PET_TER_inFETCT = sitk.Resample(PET_TER, CT_FET, tfm, sitk.sitkLinear)
    
    #Align brains with physical axes
    Euler3Dtfm = sitk.ReadTransform(reg_mtrx+'/Euler3Dtfm.txt')
    CT_FET = resample(CT_FET, Euler3Dtfm)
    PET_FET_inFETCT = resample(PET_FET_inFETCT, Euler3Dtfm)
    PET_TER_inFETCT = resample(PET_TER_inFETCT, Euler3Dtfm)
    if os.path.exists(orig_imgs +'/CT_PSMA.nii'):
        CT_TER_inFETCT = resample(CT_TER_inFETCT, Euler3Dtfm)
    
    #Save image
    sitk.WriteImage(CT_FET, reg_imgs +'/CT_FET_aligned.nii')
    sitk.WriteImage(PET_FET_inFETCT, reg_imgs +'/PET_FET_inFETCT.nii')
    if os.path.exists(orig_imgs +'/CT_PSMA.nii'):
        sitk.WriteImage(CT_TER_inFETCT, reg_imgs +'/CT_PSMA_inFETCT.nii')
    sitk.WriteImage(PET_TER_inFETCT, reg_imgs +'/PET_PSMA_inFETCT.nii')
