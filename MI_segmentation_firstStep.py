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
# subjs_name = [ f.name for f in os.scandir(data_supradir) if f.is_dir() ] #Create a lisdt of subjects names
subjs_name = ['GS008']

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
    
    for filename in glob.glob(ROIs+'/*.nii'):
        os.remove(filename)

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
            
        seeds_df = seeds_dfs.get(df)
        
        #Get seed from excel file
        subj_seeds = seeds_df.loc[subj_name].values        
        
        #Determine how many seeds are there 
        n_seeds = int(np.count_nonzero(~np.isnan(subj_seeds))/3)
        
        if n_seeds != 1:
            print(str(n_seeds) +' seeds found for '+PET_image)
    
        elif n_seeds == 1:     
            print('Only one seed found for '+PET_image)
            
        for i in range(0, n_seeds):
            print('Analysing seed '+str(i)+' for '+PET_image+' image')
            seed = (int(subj_seeds[3*i]), int(subj_seeds[3*i+1]), int(subj_seeds[3*i+2]))
        
            if df == 'FET PET':
                BTV = sitk.ConnectedThreshold(img, seedList=[seed], lower=2.2, upper=20) #Region growing VOI from the seed with connected component method
            elif df == 'PSMA PET':
                BTV = sitk.ConnectedThreshold(img, seedList=[seed], lower=1.5, upper=15) #Region growing VOI from the seed with connected component method
            
            # BTV = sitk.BinaryMorphologicalClosing(BTV, (1,1,1), sitk.sitkBall, 1.0) #fill holes in initial VOI
            # BTV = sitk.BinaryMorphologicalOpening(BTV, (1,1,1), sitk.sitkBall, 0.0, 1.0) #remove small structures in initial VOI
            sitk.WriteImage(BTV, ROIs +'/'+PET_image+'_BTV'+str(i)+'.nii')
            
        
        n_lesions = len(glob.glob(ROIs+'/'+PET_image+'_BTV*'))
                
        if n_lesions == 0:
            print('No lesion found for '+PET_image)
            break
        
        elif n_lesions == 1:
            os.rename(ROIs+'/'+PET_image+'_BTV0.nii', ROIs+'/'+PET_image+'_BTVinit.nii')
            
        else:             
            BTV0 = sitk.ReadImage(ROIs +'/'+PET_image+'_BTV0.nii')
            BTV0_dim = BTV0.GetSize()
            BTV = sitk.Image(BTV0_dim, sitk.sitkUInt8)
            BTV.CopyInformation(BTV0)
            for filename in glob.glob(ROIs+'/'+PET_image+'_BTV*'):
                imag = sitk.ReadImage(filename)
                BTV = BTV+imag 
            
            BTV = sitk.BinaryThreshold(BTV, 1, n_seeds*2, 1, 0)
            for filename in glob.glob(ROIs+'/'+PET_image+'_BTV*'):
                os.remove(filename)
            sitk.WriteImage(BTV, ROIs +'/'+PET_image+'_BTVinit.nii')
            
        BTV = sitk.ReadImage(ROIs+'/'+PET_image+'_BTVinit.nii')
        
        #Generate CTRL VOI
        CTRL_VOI = flip(BTV)
        CTRL_VOI = CTRL_VOI-BTV #Remove potential areas of overlay between BTV and CTRL VOI
        CTRL_VOI = sitk.BinaryThreshold(CTRL_VOI, 1, 1, 1, 0)    
        CTRL_VOI_stats = getStatsRoi(CTRL_VOI, img) 
        CTRL_VOI_MeanSUV = CTRL_VOI_stats.get('Mean intensity [SUV]') #Calculate CTRL VOI mean value of SUV
        
        #Generate PET TBR map
        PET_TBR_map = img/CTRL_VOI_MeanSUV 
        
        #Generate BTV based on TBR map   
        for i in range(0, n_seeds):
            print('Refining BTV by analysing seed '+str(i)+' for '+PET_image+' TBR map')
            seed = (int(subj_seeds[3*i]), int(subj_seeds[3*i+1]), int(subj_seeds[3*i+2]))
               
            if df == 'FET PET':
                BTV = sitk.ConnectedThreshold(PET_TBR_map, seedList=[seed], lower=1.7, upper=50) #Region growing VOI from the seed with connected component method
            elif df == 'PSMA PET':
                BTV = sitk.ConnectedThreshold(PET_TBR_map, seedList=[seed], lower=4, upper=150) #Region growing VOI from the seed with connected component method
              
            # BTV = sitk.BinaryMorphologicalClosing(BTV, (1,1,1), sitk.sitkBall, 1.0) #fill holes in initial VOI    
            # BTV = sitk.BinaryMorphologicalOpening(BTV, (1,1,1), sitk.sitkBall, 0.0, 1.0) #remove small structures in initial VOI
            sitk.WriteImage(BTV, ROIs +'/'+PET_image+'_BTV_TBR'+str(i)+'.nii')
        
        n_lesions = len(glob.glob(ROIs+'/'+PET_image+'_BTV_TBR*'))
                
        if n_lesions == 0:
            print('No lesion found for '+PET_image+' TBR map')
            break
        
        elif n_lesions == 1:
            os.rename(ROIs+'/'+PET_image+'_BTV_TBR0.nii', ROIs+'/'+PET_image+'_BTVsec.nii')
            
        else:             
            BTV0 = sitk.ReadImage(ROIs +'/'+PET_image+'_BTV_TBR0.nii')
            BTV0_dim = BTV0.GetSize()
            BTV = sitk.Image(BTV0_dim, sitk.sitkUInt8)
            BTV.CopyInformation(BTV0)
            for filename in glob.glob(ROIs+'/'+PET_image+'_BTV_TBR*'):
                imag = sitk.ReadImage(filename)
                BTV = BTV+imag 
            
            BTV = sitk.BinaryThreshold(BTV, 1, n_seeds*2, 1, 0)
            for filename in glob.glob(ROIs+'/'+PET_image+'_BTV_TBR*'):
                os.remove(filename)
            sitk.WriteImage(BTV, ROIs +'/'+PET_image+'_BTVsec.nii')
            
        BTV = sitk.ReadImage(ROIs+'/'+PET_image+'_BTVsec.nii')
            
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
        
        #Generate BTV based on TBR map   
        for i in range(0, n_seeds):
            print('Finalising BTV by analysing seed '+str(i)+' for '+PET_image+' TBR map')
            seed = (int(subj_seeds[3*i]), int(subj_seeds[3*i+1]), int(subj_seeds[3*i+2]))
               
            if df == 'FET PET':
                BTV = sitk.ConnectedThreshold(PET_TBR_map, seedList=[seed], lower=1.7, upper=50) #Region growing VOI from the seed with connected component method
            elif df == 'PSMA PET':
                BTV = sitk.ConnectedThreshold(PET_TBR_map, seedList=[seed], lower=5, upper=150) #Region growing VOI from the seed with connected component method
              
            # BTV = sitk.BinaryMorphologicalClosing(BTV, (1,1,1), sitk.sitkBall, 1.0) #fill holes in initial VOI    
            # BTV = sitk.BinaryMorphologicalOpening(BTV, (2,2,2), sitk.sitkBall, 0.0, 1.0) #remove small structures in initial VOI
            sitk.WriteImage(BTV, ROIs +'/'+PET_image+'_BTV_TBRfin'+str(i)+'.nii')
        
        n_lesions = len(glob.glob(ROIs+'/'+PET_image+'_BTV_TBRfin*'))
                
        if n_lesions == 0:
            print('No lesion found for '+PET_image+' TBR map')
            break
        
        elif n_lesions == 1:
            # BTV = sitk.ReadImage(ROIs +'/'+PET_image+'_BTV_TBRfin0.nii')
            # BTV = sitk.BinaryMorphologicalOpening(BTV, (1,1,1), sitk.sitkBall, 0.0, 1.0) #remove small structures in initial VOI
            os.rename(ROIs+'/'+PET_image+'_BTV_TBRfin0.nii', ROIs+'/'+PET_image+'_BTV.nii')            
            
        else:             
            BTV0 = sitk.ReadImage(ROIs +'/'+PET_image+'_BTV_TBRfin0.nii')
            BTV0_dim = BTV0.GetSize()
            BTV = sitk.Image(BTV0_dim, sitk.sitkUInt8)
            BTV.CopyInformation(BTV0)
            for filename in glob.glob(ROIs+'/'+PET_image+'_BTV_TBRfin*'):
                imag = sitk.ReadImage(filename)
                BTV = BTV+imag 
            
            BTV = sitk.BinaryThreshold(BTV, 1, n_seeds*2, 1, 0)
            # BTV = sitk.BinaryMorphologicalOpening(BTV, (1,1,1), sitk.sitkBall, 0.0, 1.0) #remove small structures in initial VOI
            for filename in glob.glob(ROIs+'/'+PET_image+'_BTV_TBRfin*'):
                os.remove(filename)
            sitk.WriteImage(BTV, ROIs +'/'+PET_image+'_BTV.nii')
               
        #Delete temporary files
        if os.path.isfile(ROIs +'/'+PET_image+'_BTVinit.nii'):
            os.remove(ROIs +'/'+PET_image+'_BTVinit.nii')
        if os.path.isfile(ROIs +'/'+PET_image+'_BTVsec.nii'):
            os.remove(ROIs +'/'+PET_image+'_BTVsec.nii')
            
