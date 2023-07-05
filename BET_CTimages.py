# -*- coding: utf-8 -*-
"""
Created on Tue Nov 30 06:01:05 2021

@author: cbri3325
"""
import SimpleITK as sitk
import os
from ImageAnalysisFunctions import *

data_supradir = 'Path to data directory' #Set working directory

subjs_path = [ f.path for f in os.scandir(data_supradir) if f.is_dir() ] #Create a list of the paths to the subjects directories
subjs_name = [ f.name for f in os.scandir(data_supradir) if f.is_dir() ] #Create a lisdt of subjects names

n_subj = len(subjs_name) #Total number of subjects



for current in subjs_name:
    subj_dir = data_supradir+current
    subj_name = current
    
    # PETFET= sitk.ReadImage(subj_dir +'/'+ subj_name +'_PET_FET_inFETCT.nii')
    # CTFET= sitk.ReadImage(subj_dir +'/'+ subj_name +'_CT_FET.nii')
    
    # brain=sitk.BinaryThreshold(CTFET, 0.0, 100.0)
    # brain2=sitk.BinaryFillhole(brain)
    # brain2=sitk.BinaryMorphologicalOpening(brain2, (8,8,8), sitk.sitkBall, 0.0, 1.0)
    # sitk.WriteImage(brain2, subj_dir+'/brain_mask.nii')

      #Perform brain extraction on PET image
    print('Perform brain extraction on PET image for '+subj_name)
    PET_FET = sitk.ReadImage(subj_dir +'/'+ subj_name +'_PET_FET_inFETCT.nii', sitk.sitkFloat32) #read PET FET image
    brain_mask = sitk.ReadImage(subj_dir +'/brain_mask.nii') #read brain mask
    
    brain_mask_filled = sitk.BinaryMorphologicalClosing(brain_mask, (1,1,1), sitk.sitkBall, 1.0) #fill holes in brain mask
    # sitk.WriteImage(brain_mask_filled, subj_dir +'/brain_mask_filled.nii') #save filled mask
    
    brain_mask_cleaned = sitk.BinaryMorphologicalOpening(brain_mask_filled, (1,1,1), sitk.sitkBall, 0.0, 1.0) #remove small structures in brain mask filled
    # sitk.WriteImage(brain_mask_cleaned, subj_dir +'/brain_mask_cleaned.nii') #save cleaned mask
    
    # brain_mask_cleaned = sitk.ReadImage(subj_dir +'/brain_mask_cleaned.nii')
    
    PET_FET_bet = generate_mask(PET_FET, brain_mask_cleaned)
    sitk.WriteImage(PET_FET_bet, subj_dir +'/' + subj_name + '_PET_FET_bet.nii')
