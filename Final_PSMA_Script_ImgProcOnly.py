# -*- coding: utf-8 -*-
"""
Created on Wed May 13 08:41:56 2020

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
from scipy.stats.stats import pearsonr
import dicom2nifti
from multiprocessing.pool import ThreadPool
from functools import partial
from ImageAnalysisFunctions import *

#%% Set Working directory
        
data_supradir = 'Path to data directory' #Set working directory

subjs_path = [ f.path for f in os.scandir(data_supradir) if f.is_dir() ] #Create a list of the paths to the subjects directories
subjs_name = [ f.name for f in os.scandir(data_supradir) if f.is_dir() ] #Create a lisdt of subjects names

n_subj = len(subjs_name) #Total number of subjects


#%%Create a for loop to perform image analysis on each subject sequentially

for current in subjs_name:
    subj_dir = data_supradir+current
    subj_name = current
    
    #Unzip all nii.gz files and delete original gzipped files
    for filename in glob.glob(subj_dir +'/' +'*.gz'):
        gunzip_shutil(filename, filename[:-3])
    
    # #Rename all the image files in format "subjectName_scanType.nii"
    rename_img_files(subj_dir, subj_name)
    
    #PET pre-processing
    PET_preprocessing(subj_dir, subj_name)
    
    # Generate subfolders into the patient folder   
    os.mkdir(subj_dir +'/Original images')
    os.mkdir(subj_dir +'/Registered images')
    os.mkdir(subj_dir +'/Registration matrices')
    os.mkdir(subj_dir +'/ROIs')
    os.mkdir(subj_dir +'/InfoTextFiles')
    os.mkdir(subj_dir +'/Results')
    os.mkdir(subj_dir +'/Temp')
    
    #Set paths to subfolders    
    orig_imgs = subj_dir +'/Original images'
    reg_imgs = subj_dir +'/Registered images'
    reg_mtrx = subj_dir +'/Registration matrices'
    ROIs = subj_dir +'/ROIs'
    info = subj_dir +'/InfoTextFiles'
    results = subj_dir +'/Results'
    temp = subj_dir +'/Temp'
    
    # Move files into subfolders   
    for file in glob.glob(subj_dir +'/' +'*.txt'):
        shutil.move(os.path.join(subj_dir, file), info)
    for file in glob.glob(subj_dir +'/' +'*.xlsx'):
        shutil.move(os.path.join(subj_dir, file), info)
    for file in glob.glob(subj_dir +'/' +'*.nii'):
        shutil.move(os.path.join(subj_dir, file), orig_imgs)
    
    #Read Original images
 
    PET_FET = sitk.ReadImage(orig_imgs +'/'+ subj_name +'_PET_FET.nii', sitk.sitkFloat32) #read PET FET image
    CT_FET = sitk.ReadImage(orig_imgs +'/'+ subj_name +'_CT_FET.nii', sitk.sitkFloat32) #read CT FET image
    PET_TER = sitk.ReadImage(orig_imgs +'/'+ subj_name +'_PET_TER.nii', sitk.sitkFloat32) #read PET TER image
    CT_TER = sitk.ReadImage(orig_imgs +'/'+ subj_name +'_CT_TER.nii', sitk.sitkFloat32) #read CT TER image
    
   
    #Registration and Resampling to FET_CT image and resolution
    print('Start registration for '+current)
    
    #Resample PET FET image into CT FET space
    PET_FET_inFETCT = Resample_image(PET_FET, CT_FET)
    sitk.WriteImage(PET_FET_inFETCT, reg_imgs +'/'+ subj_name +'_PET_FET_inFETCT.nii')

    
    #Resample PET TER image into CT TER space
    PET_TER_inTERCT = Resample_image(PET_TER, CT_TER)
    sitk.WriteImage(PET_TER_inTERCT, reg_imgs +'/'+ subj_name +'_PET_TER_inTERCT.nii')

    
    ##Registration of T1CE, PET TER and CT TER to CT FET
    fixed_image = CT_FET

    moving_images = [CT_TER, PET_TER]
        
    for moving_image in moving_images:
        
        #Register and resample image
        if moving_image is CT_TER:
            moving_resampled, final_transform = RegisterResample_image(fixed_image, moving_image, 'same', 'Multires')
        elif moving_image is PET_TER:
            moving_resampled, final_transform = RegisterResample_image(fixed_image, moving_image, 'different', 'ExplorExploit')
    
        #Save image and transform to file

        if moving_image is CT_TER:
            sitk.WriteImage(moving_resampled, reg_imgs +'/'+ subj_name +'_CT_TER_inFETCT.nii')
            sitk.WriteTransform(final_transform, reg_mtrx +'/'+ subj_name +'_CTTER_2_FETCT.tfm')
        elif moving_image is PET_TER:
            sitk.WriteImage(moving_resampled, reg_imgs +'/'+ subj_name +'_PET_TER_inFETCT.nii')
            sitk.WriteTransform(final_transform, reg_mtrx +'/'+ subj_name +'_PETTER_2_FETCT.tfm')
    
    #Read registered/resampled images   
    
    PETFET_inFETCT = sitk.ReadImage(reg_imgs +'/'+ subj_name +'_PET_FET_inFETCT.nii') #read PET FET in FET CT image
    PETTER_inFETCT = sitk.ReadImage(reg_imgs +'/'+ subj_name +'_PET_TER_inFETCT.nii') #read PET TER in FET CT image
    CTTER_inFETCT = sitk.ReadImage(reg_imgs +'/'+ subj_name +'_CT_TER_inFETCT.nii') #read CT TER in FET CT image
    PETTER_inTERCT = sitk.ReadImage(reg_imgs +'/'+ subj_name +'_PET_TER_inTERCT.nii') #read PET TER in TER CT image
    

    
    #Generate Gaussian filter smoothed PET images (using sigma = 2mm)
    print('Generate Gussian smoothed PET images for '+current)
    
    gaussFETPET = gauss_smooth(PETFET_inFETCT, 2.0)
    sitk.WriteImage(gaussFETPET, reg_imgs +'/'+ subj_name +'__PET_FET_inFETCTgauss.nii')
    
    gaussTERPET = gauss_smooth(PETTER_inFETCT, 2.0)    
    sitk.WriteImage(gaussTERPET, reg_imgs +'/'+ subj_name +'__PET_TER_inFETCTgauss.nii')


##%%Generate result images for visualization purposes
#import PSMA_Trial_ImgResultsVisual
