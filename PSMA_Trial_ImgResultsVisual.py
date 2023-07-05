# -*- coding: utf-8 -*-
"""
Created on Wed May 13 10:54:41 2020

@author: Caterina Brighi
"""

import numpy as np
import pandas as pd
import SimpleITK as sitk
import os
from ImageResultsVisualizationFunctions import *

data_supradir = 'Path to data directory'

subjs_path = [ f.path for f in os.scandir(data_supradir) if f.is_dir() ] #Create a list of the paths to the subjects directories
subjs_name = [ f.name for f in os.scandir(data_supradir) if f.is_dir() ] #Create a lisdt of subjects names

#Make sure that the data_supradir contains an excel file with the number of z_slice with tumour in CT_FET space for each patient

#Make a directory for the images
os.mkdir(data_supradir +'/Visualization Results Images')
imgs_dir = data_supradir +'/Visualization Results Images' #Set paths to plot dir

#Import z_slices values into a dataframe
z_slices = pd.read_excel(data_supradir + 'z_slices.xlsx', sheet_name='Sheet1', index_col='Subject_ID')

for current in subjs_name:
    if current != 'Analysis Results Plots':
        
        subj_dir = data_supradir+current
        subj_name = current
        
        #Set paths to subfolders 
        ROIs = subj_dir +'/ROIs'
        reg_imgs = subj_dir +'/Registered images' 
        reg_ROIs = ROIs +'/Registered ROIs'
        MCTRL_ROIs = ROIs +'/Mirror image CTRL ROIs'
        CCTRL_ROIs = ROIs +'/Crescent CTRL ROIs'
        subCE_ROIs = ROIs +'/subCEtum ROIs'
        
        #Read registered/resampled images   
        if os.path.exists(reg_imgs +'/'+ subj_name +'_T1CE_inFETCT.nii'):
            T1CE_inFETCT = sitk.ReadImage(reg_imgs +'/'+ subj_name +'_T1CE_inFETCT.nii') #read T1CE in FET CT image
        PETFET_inFETCT = sitk.ReadImage(reg_imgs +'/'+ subj_name +'_PET_FET_inFETCT.nii') #read PET FET in FET CT image
        PETTER_inFETCT = sitk.ReadImage(reg_imgs +'/'+ subj_name +'_PET_TER_inFETCT.nii') #read PET TER in FET CT image
        CTTER_inFETCT = sitk.ReadImage(reg_imgs +'/'+ subj_name +'_CT_TER_inFETCT.nii') #read CT TER in FET CT image
        PETTER_inTERCT = sitk.ReadImage(reg_imgs +'/'+ subj_name +'_PET_TER_inTERCT.nii') #read PET TER in TER CT image
        PETFET_inFETCTgauss = sitk.ReadImage(reg_imgs +'/'+ subj_name +'__PET_FET_inFETCTgauss.nii') #read Gaussian PET FET in FET CT image
        PETTER_inFETCTgauss = sitk.ReadImage(reg_imgs +'/'+ subj_name +'__PET_TER_inFETCTgauss.nii') #read Gaussian TER FET in FET CT image
        
        #Read registered and resampled ROIs   
        if os.path.exists(reg_ROIs +'/'+ subj_name +'_CEtum_inFETCT.nii'):
            CEtum_inFETCT = sitk.ReadImage(reg_ROIs +'/'+ subj_name +'_CEtum_inFETCT.nii') #read CEtum roi in FET CT image
        FETroi_inFETCT = sitk.ReadImage(reg_ROIs +'/'+ subj_name +'_FETroi_inFETCT.nii') #read FET roi in FET CT image
        TERroi_inFETCT = sitk.ReadImage(reg_ROIs +'/'+ subj_name +'_TERroi_inFETCT.nii') #read TER roi in FET CT image
        TERroi_inTERCT = sitk.ReadImage(reg_ROIs +'/'+ subj_name +'_TERroi_inTERCT.nii') #read TER roi in TER CT image
        
        #Read mirror image CTRL rois    
        FETroi_inFETCT_CTRL = sitk.ReadImage(MCTRL_ROIs +'/'+ subj_name +'_FETroi_inFETCT_CTRL.nii') #read FET CTRL roi in FET CT image
        TERroi_inFETCT_CTRL = sitk.ReadImage(MCTRL_ROIs +'/'+ subj_name +'_TERroi_inFETCT_CTRL.nii') #read TER CTRL roi in FET CT image
        TERroi_inTERCT_CTRL = sitk.ReadImage(MCTRL_ROIs +'/'+ subj_name +'_TERroi_inTERCT_CTRL.nii') #read TER CTRL roi in TER CT image
        
        #Read crescent-shaped CTRL rois    
        crescentFET = sitk.ReadImage(CCTRL_ROIs +'/'+ subj_name +'_crescentFET.nii')
        crescentTER = sitk.ReadImage(CCTRL_ROIs +'/'+ subj_name +'_crescentTER.nii')
        
        #Read FET thresholded values sub-rois of CEtum roi
        if os.path.exists(reg_ROIs +'/'+ subj_name +'_CEtum_inFETCT.nii'):
            FETminCEtum = sitk.ReadImage(subCE_ROIs +'/'+ subj_name +'_FETminCEtum.nii')
            FETmaxCEtum = sitk.ReadImage(subCE_ROIs +'/'+ subj_name +'_FETmaxCEtum.nii')
        
        #Make a folder of current in visualisation image folder
        os.mkdir(imgs_dir +'/'+ current)
        img_current = imgs_dir +'/'+ current
        
        #Generate images of PET + respective rois
        slice_n = int(z_slices.loc[current,'z_slice'])
        
        FETPET_TERPET_ovl = colormap_imgs_overlay(PETFET_inFETCT, PETTER_inFETCT, slice_n, 'Hot', alpha_value=0.85)    
        sitk.WriteImage(FETPET_TERPET_ovl, img_current +'/'+ subj_name +'_FETTERPEToverlay.png')
        
        FETroi_onFETPET = roi_to_img_overlay(PETFET_inFETCT, FETroi_inFETCT, slice_n, 'red', opacity=0.2)
        sitk.WriteImage(FETroi_onFETPET, img_current +'/'+ subj_name +'_FETroi_onFETPET.png')
        
        TERroi_onTERPET = roi_to_img_overlay(PETTER_inFETCT, TERroi_inFETCT, slice_n, 'red', opacity=0.2)
        sitk.WriteImage(TERroi_onTERPET, img_current +'/'+ subj_name +'_TERroi_onTERPET.png')
        
        FETroi_onTERPET = roi_to_img_overlay(PETTER_inFETCT, FETroi_inFETCT, slice_n, 'red', opacity=0.2)
        sitk.WriteImage(FETroi_onTERPET, img_current +'/'+ subj_name +'_FETroi_onTERPET.png')
        
        if os.path.exists(reg_imgs +'/'+ subj_name +'_T1CE_inFETCT.nii'):
            CEtum_onT1CE = roi_to_img_overlay(T1CE_inFETCT, CEtum_inFETCT, slice_n, 'yellow', opacity=0.2)
            sitk.WriteImage(CEtum_onT1CE, img_current +'/'+ subj_name +'_CEtum_onT1CE.png')
        
        FETroi_onFETPETgauss = roi_to_img_overlay(PETFET_inFETCTgauss, FETroi_inFETCT, slice_n, 'red', opacity=0.2)
        sitk.WriteImage(FETroi_onFETPETgauss, img_current +'/'+ subj_name +'_FETroi_onFETPETgauss.png')
        
        FETroi_onTERPETgauss = roi_to_img_overlay(PETTER_inFETCTgauss, FETroi_inFETCT, slice_n, 'red', opacity=0.2)
        sitk.WriteImage(FETroi_onTERPETgauss, img_current +'/'+ subj_name +'_FETroi_onTERPETgauss.png')
        
        FETroi_plusCTRLMI = FETroi_inFETCT+FETroi_inFETCT_CTRL
        FETPETroisMI = labels_to_img_overlay(PETFET_inFETCT, FETroi_plusCTRLMI, slice_n, opacity=0.2)
        sitk.WriteImage(FETPETroisMI, img_current +'/'+ subj_name +'_FETPETroisMI.png')
        
        FETroi_plusCTRLCS = FETroi_inFETCT+crescentFET
        FETPETroisCS = labels_to_img_overlay(PETFET_inFETCT, FETroi_plusCTRLCS, slice_n, opacity=0.2)
        sitk.WriteImage(FETPETroisCS, img_current +'/'+ subj_name +'_FETPETroisCS.png')
        
        TERroi_plusCTRLMI = TERroi_inFETCT+TERroi_inFETCT_CTRL
        TERPETroisMI = labels_to_img_overlay(PETTER_inFETCT, TERroi_plusCTRLMI, slice_n, opacity=0.2)
        sitk.WriteImage(TERPETroisMI, img_current +'/'+ subj_name +'_TERPETroisMI.png')
        
        TERroi_plusCTRLCS = TERroi_inFETCT+crescentFET
        TERPETroisCS = labels_to_img_overlay(PETTER_inFETCT, TERroi_plusCTRLCS, slice_n, opacity=0.2)
        sitk.WriteImage(TERPETroisCS, img_current +'/'+ subj_name +'_TERPETroisCS.png')
        
        if os.path.exists(reg_imgs +'/'+ subj_name +'_T1CE_inFETCT.nii'):
            FETmin_plusFETmax = FETminCEtum+FETmaxCEtum
            FETminmaxrois_onT1CE = labels_to_img_overlay(T1CE_inFETCT, FETmin_plusFETmax, slice_n, opacity=0.2)
            sitk.WriteImage(FETminmaxrois_onT1CE, img_current +'/'+ subj_name +'_FETminmaxrois_onT1CE.png')
            FETminmaxrois_onTERPET = labels_to_img_overlay(PETTER_inFETCT, FETmin_plusFETmax, slice_n, opacity=0.2)
            sitk.WriteImage(FETminmaxrois_onTERPET, img_current +'/'+ subj_name +'_FETminmaxrois_onTERPET.png')
            
