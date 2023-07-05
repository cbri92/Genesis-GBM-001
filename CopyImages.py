# -*- coding: utf-8 -*-
"""
Created on Mon Feb  6 13:21:43 2023

@author: Caterina Brighi
"""

#%% Import functions 

import os
import shutil

#%% Set Working directory
        
data_supradir = 'Path to data directory' #Set working directory

subjs_path = [ f.path for f in os.scandir(data_supradir) if f.is_dir() ] #Create a list of the paths to the subjects directories
subjs_name = [ f.name for f in os.scandir(data_supradir) if f.is_dir() ] #Create a lisdt of subjects names
subjs_name.remove('Analysis Results Plots')

n_subj = len(subjs_name) #Total number of subjects

dest_dir = 'C:/Users/cbri3325/OneDrive - The University of Sydney (Staff)/Caterina Brighi/Data/test/'

#%%Create a for loop to perform image analysis on each subject sequentially

for current in subjs_name:
    subj_dir = data_supradir+current
    subj_name = current
    
    print('Copying images for '+current)
    
    #Make new patient directory
    if not os.path.exists(dest_dir+current):#if it does not already exist, create a directory of patient data
        os.mkdir(dest_dir+current)
    pt_dir = dest_dir+current
    
    if not os.path.exists(pt_dir+'/Original/'):#if it does not already exist, create a directory of patient data
        os.mkdir(pt_dir+'/Original/')
    orig_dir = pt_dir+'/Original/'
    
    if not os.path.exists(pt_dir+'/ROIs/'):#if it does not already exist, create a directory of patient data
        os.mkdir(pt_dir+'/ROIs/')
    roi_dir = pt_dir+'/ROIs/'
    
    if not os.path.exists(pt_dir+'/Registered/'):#if it does not already exist, create a directory of patient data
        os.mkdir(pt_dir+'/Registered/')
    reg_dir = pt_dir+'/Registered/'
    
    if not os.path.exists(pt_dir+'/Transforms/'):#if it does not already exist, create a directory of patient data
        os.mkdir(pt_dir+'/Transforms/')
    tfm_dir = pt_dir+'/Transforms/'
    
    #Copy images into test folder
    shutil.copy2(subj_dir+'/Original images/'+current+'_CT_FET.nii', orig_dir+'CT_FET.nii')
    shutil.copy2(subj_dir+'/Original images/'+current+'_CT_TER.nii', orig_dir+'CT_PSMA.nii')
    shutil.copy2(subj_dir+'/Original images/'+current+'_PET_FET.nii', orig_dir+'PET_FET.nii')
    shutil.copy2(subj_dir+'/Original images/'+current+'_PET_TER.nii', orig_dir+'PET_PSMA.nii')
    
    # shutil.copy2(subj_dir+'/Registered images/'+current+'_PET_FET_inFETCT.nii', img_dir+'PET_FET.nii')
    # shutil.copy2(subj_dir+'/Registered images/'+current+'_CT_TER_inFETCT.nii', img_dir+'CT_PSMA.nii')
    # shutil.copy2(subj_dir+'/Registered images/'+current+'_PET_TER_inFETCT.nii', img_dir+'PET_PSMA.nii')
   
