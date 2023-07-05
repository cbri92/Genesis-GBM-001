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

#Generate excel file containing voxelwise results

voxelwise_results = pd.ExcelWriter(data_supradir +'voxelwise_results.xlsx')
Group_stats_results = pd.ExcelWriter(data_supradir +'Group_stats_results.xlsx')

#Create empty dataframes to populate as going through the loop

CEtum_onT1CEdf = pd.DataFrame(columns=['Subject_ID', 'Volume [mm3]', 'Mean intensity [SUV]', 'Std of Mean [SUV]', 'Median intensity [SUV]', 'Max intensity [SUV]', 'Min intensity [SUV]'])
CEtumFETCT_onT1CEdf = pd.DataFrame(columns=['Subject_ID', 'Volume [mm3]', 'Mean intensity [SUV]', 'Std of Mean [SUV]', 'Median intensity [SUV]', 'Max intensity [SUV]', 'Min intensity [SUV]'])

FETroi_onFETPETdf = pd.DataFrame(columns=['Subject_ID', 'Volume [mm3]', 'Mean intensity [SUV]', 'Std of Mean [SUV]', 'Median intensity [SUV]', 'Max intensity [SUV]', 'Min intensity [SUV]'])
FETroi_onTERPETdf = pd.DataFrame(columns=['Subject_ID', 'Volume [mm3]', 'Mean intensity [SUV]', 'Std of Mean [SUV]', 'Median intensity [SUV]', 'Max intensity [SUV]', 'Min intensity [SUV]'])
CTRLFETroi_onFETPETdf = pd.DataFrame(columns=['Subject_ID', 'Volume [mm3]', 'Mean intensity [SUV]', 'Std of Mean [SUV]', 'Median intensity [SUV]', 'Max intensity [SUV]', 'Min intensity [SUV]'])
CTRLFETroi_onTERPETdf = pd.DataFrame(columns=['Subject_ID', 'Volume [mm3]', 'Mean intensity [SUV]', 'Std of Mean [SUV]', 'Median intensity [SUV]', 'Max intensity [SUV]', 'Min intensity [SUV]'])

TERroi_onFETPETdf = pd.DataFrame(columns=['Subject_ID', 'Volume [mm3]', 'Mean intensity [SUV]', 'Std of Mean [SUV]', 'Median intensity [SUV]', 'Max intensity [SUV]', 'Min intensity [SUV]'])
TERroi_onTERPETdf = pd.DataFrame(columns=['Subject_ID', 'Volume [mm3]', 'Mean intensity [SUV]', 'Std of Mean [SUV]', 'Median intensity [SUV]', 'Max intensity [SUV]', 'Min intensity [SUV]'])
CTRLTERroi_onFETPETdf = pd.DataFrame(columns=['Subject_ID', 'Volume [mm3]', 'Mean intensity [SUV]', 'Std of Mean [SUV]', 'Median intensity [SUV]', 'Max intensity [SUV]', 'Min intensity [SUV]'])
CTRLTERroi_onTERPETdf = pd.DataFrame(columns=['Subject_ID', 'Volume [mm3]', 'Mean intensity [SUV]', 'Std of Mean [SUV]', 'Median intensity [SUV]', 'Max intensity [SUV]', 'Min intensity [SUV]'])

TERroi_onTERPETorigdf = pd.DataFrame(columns=['Subject_ID', 'Volume [mm3]', 'Mean intensity [SUV]', 'Std of Mean [SUV]', 'Median intensity [SUV]', 'Max intensity [SUV]', 'Min intensity [SUV]'])
CTRLTERroi_onTERPETorigdf = pd.DataFrame(columns=['Subject_ID', 'Volume [mm3]', 'Mean intensity [SUV]', 'Std of Mean [SUV]', 'Median intensity [SUV]', 'Max intensity [SUV]', 'Min intensity [SUV]'])

CCRLFETroi_onFETPETdf = pd.DataFrame(columns=['Subject_ID', 'Volume [mm3]', 'Mean intensity [SUV]', 'Std of Mean [SUV]', 'Median intensity [SUV]', 'Max intensity [SUV]', 'Min intensity [SUV]'])
CCRLFETroi_onTERPETdf = pd.DataFrame(columns=['Subject_ID', 'Volume [mm3]', 'Mean intensity [SUV]', 'Std of Mean [SUV]', 'Median intensity [SUV]', 'Max intensity [SUV]', 'Min intensity [SUV]'])
CCRLTERroi_onTERPETdf = pd.DataFrame(columns=['Subject_ID', 'Volume [mm3]', 'Mean intensity [SUV]', 'Std of Mean [SUV]', 'Median intensity [SUV]', 'Max intensity [SUV]', 'Min intensity [SUV]'])

FETminCE_onTERPETdf = pd.DataFrame(columns=['Subject_ID', 'Volume [mm3]', 'Mean intensity [SUV]', 'Std of Mean [SUV]', 'Median intensity [SUV]', 'Max intensity [SUV]', 'Min intensity [SUV]'])
FETmaxCE_onTERPETdf = pd.DataFrame(columns=['Subject_ID', 'Volume [mm3]', 'Mean intensity [SUV]', 'Std of Mean [SUV]', 'Median intensity [SUV]', 'Max intensity [SUV]', 'Min intensity [SUV]'])

TBRdf = pd.DataFrame(columns=['Subject_ID', 'FET TBRmean MI', 'FET stdTBRmean MI', 'FET TBRmax MI', 'FET TBRmean CS', 'FET stdTBRmean CS', 'FET TBRmax CS', 'TER TBRmean MI', 'TER stdTBRmean MI', 'TER TBRmax MI', 'TER TBRmean CS', 'TER stdTBRmean CS', 'TER TBRmax CS'])

VolRatiodf = pd.DataFrame(columns=['Subject_ID', 'TER tum volume/FET tum volume'])

prsdf = pd.DataFrame(columns=['Subject_ID', 'W/O Gauss Pearson corr coeff', 'W/O Gauss p-value', 'W+ Gauss Pearson corr coeff', 'W+ Gauss p-value']) #generate an empty dataframe to populate with Pearson correl coeff


#%%Create a for loop to perform image analysis on each subject sequentially

for current in subjs_name:
    subj_dir = data_supradir+current
    subj_name = current
       
    #Set paths to subfolders    
    orig_imgs = subj_dir +'/Original images'
    reg_imgs = subj_dir +'/Registered images'
    reg_mtrx = subj_dir +'/Registration matrices'
    ROIs = subj_dir +'/ROIs'
    info = subj_dir +'/InfoTextFiles'
    results = subj_dir +'/Results'
    temp = subj_dir +'/Temp'

    #Read Original images
    T1CE = 0
    if os.path.exists(orig_imgs +'/'+ subj_name +'_T1CE.nii'):
        T1CE = sitk.ReadImage(orig_imgs +'/'+ subj_name +'_T1CE.nii', sitk.sitkFloat32) #read T1CE image
    PET_FET = sitk.ReadImage(orig_imgs +'/'+ subj_name +'_PET_FET.nii', sitk.sitkFloat32) #read PET FET image
    CT_FET = sitk.ReadImage(orig_imgs +'/'+ subj_name +'_CT_FET.nii', sitk.sitkFloat32) #read CT FET image
    PET_TER = sitk.ReadImage(orig_imgs +'/'+ subj_name +'_PET_TER.nii', sitk.sitkFloat32) #read PET TER image
    CT_TER = sitk.ReadImage(orig_imgs +'/'+ subj_name +'_CT_TER.nii', sitk.sitkFloat32) #read CT TER image
       
    #Read registered/resampled images   
    if os.path.exists(reg_imgs +'/'+ subj_name +'_T1CE_inFETCT.nii'):
        T1CE_inFETCT = sitk.ReadImage(reg_imgs +'/'+ subj_name +'_T1CE_inFETCT.nii') #read T1CE in FET CT image
    PETFET_inFETCT = sitk.ReadImage(reg_imgs +'/'+ subj_name +'_PET_FET_inFETCT.nii') #read PET FET in FET CT image
    PETTER_inFETCT = sitk.ReadImage(reg_imgs +'/'+ subj_name +'_PET_TER_inFETCT.nii') #read PET TER in FET CT image
    CTTER_inFETCT = sitk.ReadImage(reg_imgs +'/'+ subj_name +'_CT_TER_inFETCT.nii') #read CT TER in FET CT image
    PETTER_inTERCT = sitk.ReadImage(reg_imgs +'/'+ subj_name +'_PET_TER_inTERCT.nii') #read PET TER in TER CT image
    
    #Read transformations    
    if os.path.exists(reg_mtrx +'/'+ subj_name +'_T1CE_2_FETCT.tfm'):
        T1CE2FETCT_tfm = sitk.ReadTransform(reg_mtrx +'/'+ subj_name +'_T1CE_2_FETCT.tfm')
    PETTER2FETCT_tfm = sitk.ReadTransform(reg_mtrx +'/'+ subj_name +'_PETTER_2_FETCT.tfm')
    
    #Set paths to subfolders   
    orig_ROIs = ROIs +'/Original ROIs'
    reg_ROIs = ROIs +'/Registered ROIs'
    MCTRL_ROIs = ROIs +'/Mirror image CTRL ROIs'
    CCTRL_ROIs = ROIs +'/Crescent CTRL ROIs'
    subCE_ROIs = ROIs +'/subCEtum ROIs'
    
    #Read original ROIs
    if os.path.exists(orig_ROIs +'/CEtum.nii'):
        CEtum = sitk.ReadImage(orig_ROIs +'/CEtum.nii') #read CEtum roi
        
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
    
    #Read Gaussian filter smoothed PET images   
    PETFET_inFETCTgauss = sitk.ReadImage(reg_imgs +'/'+ subj_name +'__PET_FET_inFETCTgauss.nii')
    PETTER_inFETCTgauss = sitk.ReadImage(reg_imgs +'/'+ subj_name +'__PET_TER_inFETCTgauss.nii')
    
    ##Calculate general statistics from each ROI
    print('Calculate general stats for '+current)
    
    def update_stats_dict(roi, image):
        roi_onImage_stats = getStatsRoi(roi, image)
        roi_onImage_stats.update({'Subject_ID': subj_name})
        return roi_onImage_stats
    
    if os.path.exists(reg_ROIs +'/'+ subj_name +'_CEtum_inFETCT.nii'):
        CEtum_onT1CE_stats = update_stats_dict(CEtum, T1CE)
        CEtumFETCT_onT1CE_stats = update_stats_dict(CEtum_inFETCT, T1CE_inFETCT)
    
    FETroi_onFETPET_stats = update_stats_dict(FETroi_inFETCT, PETFET_inFETCT)
    FETroi_onTERPET_stats = update_stats_dict(FETroi_inFETCT, PETTER_inFETCT)
    CTRLFETroi_onFETPET_stats = update_stats_dict(FETroi_inFETCT_CTRL, PETFET_inFETCT)
    CTRLFETroi_onTERPET_stats = update_stats_dict(FETroi_inFETCT_CTRL, PETTER_inFETCT)
    
    TERroi_onFETPET_stats = update_stats_dict(TERroi_inFETCT, PETFET_inFETCT)
    TERroi_onTERPET_stats = update_stats_dict(TERroi_inFETCT, PETTER_inFETCT)
    CTRLTERroi_onFETPET_stats = update_stats_dict(TERroi_inFETCT_CTRL, PETFET_inFETCT)
    CTRLTERroi_onTERPET_stats = update_stats_dict(TERroi_inFETCT_CTRL, PETTER_inFETCT)
    
    TERroi_onTERPETorig_stats = update_stats_dict(TERroi_inTERCT, PETTER_inTERCT)
    CTRLTERroi_onTERPETorig_stats = update_stats_dict(TERroi_inTERCT_CTRL, PETTER_inTERCT)
    
    CCRLFETroi_onFETPET_stats = update_stats_dict(crescentFET, PETFET_inFETCT)
    CCRLFETroi_onTERPET_stats = update_stats_dict(crescentFET, PETTER_inFETCT)
    CCRLTERroi_onTERPET_stats = update_stats_dict(crescentTER, PETTER_inTERCT)
    
    if os.path.exists(reg_ROIs +'/'+ subj_name +'_CEtum_inFETCT.nii'):
        FETminCE_onTERPET_stats = update_stats_dict(FETminCEtum, PETTER_inFETCT)
        FETmaxCE_onTERPET_stats = update_stats_dict(FETmaxCEtum, PETTER_inFETCT)
    
    #Calculating TBRmean and TBRmax        
    FET_TBRmeanMI, FET_stdTBRmeanMI, FET_TBRmaxMI = TBR(FETroi_onFETPET_stats, CTRLFETroi_onFETPET_stats)    
    FET_TBRmeanCS, FET_stdTBRmeanCS, FET_TBRmaxCS = TBR(FETroi_onFETPET_stats, CCRLFETroi_onFETPET_stats)     
    TER_TBRmeanMI, TER_stdTBRmeanMI, TER_TBRmaxMI = TBR(TERroi_onTERPET_stats, CTRLTERroi_onTERPET_stats)
    TER_TBRmeanCS, TER_stdTBRmeanCS, TER_TBRmaxCS = TBR(TERroi_onTERPET_stats, CCRLTERroi_onTERPET_stats)    
    TBR_dict = {'Subject_ID': subj_name, 'FET TBRmean MI': FET_TBRmeanMI, 'FET stdTBRmean MI': FET_stdTBRmeanMI, 'FET TBRmax MI': FET_TBRmaxMI, 'FET TBRmean CS': FET_TBRmeanCS, 'FET stdTBRmean CS': FET_stdTBRmeanCS, 'FET TBRmax CS': FET_TBRmaxCS, 'TER TBRmean MI': TER_TBRmeanMI, 'TER stdTBRmean MI': TER_stdTBRmeanMI, 'TER TBRmax MI': TER_TBRmaxMI, 'TER TBRmean CS': TER_TBRmeanCS, 'TER stdTBRmean CS': TER_stdTBRmeanCS, 'TER TBRmax CS': TER_TBRmaxCS}
    
    #Calculate ratio between TER roi volume and FET roi volume
    
    RoisVolRatio = rois_volume_ratio(TERroi_onTERPET_stats,FETroi_onFETPET_stats)
    VolRatio_dict = {'Subject_ID': subj_name, 'TER tum volume/FET tum volume': RoisVolRatio}
    
    #Export avg stats values to dataframe then excel spreadsheet
    print('Export general stats for '+current)
    
    if os.path.exists(reg_ROIs +'/'+ subj_name +'_CEtum_inFETCT.nii'):
        CEtum_onT1CEdf = CEtum_onT1CEdf.append(CEtum_onT1CE_stats, ignore_index=True)
        CEtumFETCT_onT1CEdf = CEtumFETCT_onT1CEdf.append(CEtumFETCT_onT1CE_stats, ignore_index=True)
    
    FETroi_onFETPETdf = FETroi_onFETPETdf.append(FETroi_onFETPET_stats, ignore_index=True)
    FETroi_onTERPETdf = FETroi_onTERPETdf.append(FETroi_onTERPET_stats, ignore_index=True)
    CTRLFETroi_onFETPETdf = CTRLFETroi_onFETPETdf.append(CTRLFETroi_onFETPET_stats, ignore_index=True)
    CTRLFETroi_onTERPETdf = CTRLFETroi_onTERPETdf.append(CTRLFETroi_onTERPET_stats, ignore_index=True)
    
    TERroi_onFETPETdf = TERroi_onFETPETdf.append(TERroi_onFETPET_stats, ignore_index=True)
    TERroi_onTERPETdf = TERroi_onTERPETdf.append(TERroi_onTERPET_stats, ignore_index=True)
    CTRLTERroi_onFETPETdf = CTRLTERroi_onFETPETdf.append(CTRLTERroi_onFETPET_stats, ignore_index=True)
    CTRLTERroi_onTERPETdf = CTRLTERroi_onTERPETdf.append(CTRLTERroi_onTERPET_stats, ignore_index=True)
    
    TERroi_onTERPETorigdf = TERroi_onTERPETorigdf.append(TERroi_onTERPETorig_stats, ignore_index=True) 
    CTRLTERroi_onTERPETorigdf = CTRLTERroi_onTERPETdf.append(CTRLTERroi_onTERPETorig_stats, ignore_index=True)
    
    CCRLFETroi_onFETPETdf = CCRLFETroi_onFETPETdf.append(CCRLFETroi_onFETPET_stats, ignore_index=True)
    CCRLFETroi_onTERPETdf = CCRLFETroi_onTERPETdf.append(CCRLFETroi_onTERPET_stats, ignore_index=True)
    CCRLTERroi_onTERPETdf = CCRLTERroi_onTERPETdf.append(CCRLTERroi_onTERPET_stats, ignore_index=True)
    
    if os.path.exists(reg_ROIs +'/'+ subj_name +'_CEtum_inFETCT.nii'):
        FETminCE_onTERPETdf = FETminCE_onTERPETdf.append(FETminCE_onTERPET_stats, ignore_index=True)
        FETmaxCE_onTERPETdf = FETmaxCE_onTERPETdf.append(FETmaxCE_onTERPET_stats, ignore_index=True)
    
    TBRdf = TBRdf.append(TBR_dict, ignore_index=True)
    
    VolRatiodf = VolRatiodf.append(VolRatio_dict, ignore_index=True)
    
    ##Calculate pixel intensity values for all pixels in roi
    print('Calculate voxelwise stats for '+current)
    
    #PET images without Gaussian filter
    FETroi_onPETFET_vxls = []
    FETroi_onPETTER_vxls = []
    FETroi_onPETFET_vxls = allPixInt(PETFET_inFETCT, FETroi_inFETCT, 'FETroionFETPETmask', temp)
    FETroi_onPETTER_vxls = allPixInt(PETTER_inFETCT, FETroi_inFETCT, 'FETroionTERPETmask', temp)
    r, p = pearsonr(FETroi_onPETFET_vxls, FETroi_onPETTER_vxls)
    
    #PET images with Gaussian filter sigma = 2mm
    FETroi_onPETFETgauss_vxls = []
    FETroi_onPETTERgauss_vxls = []
    FETroi_onPETFETgauss_vxls = allPixInt(PETFET_inFETCTgauss, FETroi_inFETCT, 'FETroionFETPETgaussmask', temp)
    FETroi_onPETTERgauss_vxls = allPixInt(PETTER_inFETCTgauss, FETroi_inFETCT, 'FETroionTERPETgaussmask', temp)
    r1, p1 = pearsonr(FETroi_onPETFETgauss_vxls, FETroi_onPETTERgauss_vxls)
    
    ##Export pixel intensity values and Pearson correl coeff to dataframe and then excel spreadsheet
    print('Export voxelwise stats for '+current)
    
    data = {'FET roi on FETPET': FETroi_onPETFET_vxls, 'FET roi on TERPET': FETroi_onPETTER_vxls, 'FET roi on FETPET gauss': FETroi_onPETFETgauss_vxls, 'FET roi on TERPET gauss': FETroi_onPETTERgauss_vxls}
    df = pd.DataFrame(data, columns=['FET roi on FETPET', 'FET roi on TERPET', 'FET roi on FETPET gauss', 'FET roi on TERPET gauss']) #generate dataframe of pixel intensity values
    df.to_excel(voxelwise_results, sheet_name=subj_name, index=False)
    
    pearsc = {'Subject_ID': subj_name, 'W/O Gauss Pearson corr coeff': r, 'W/O Gauss p-value': p, 'W+ Gauss Pearson corr coeff': r1, 'W+ Gauss p-value': p1}
    prsdf = prsdf.append(pearsc, ignore_index=True)
    
    #Generate correlation plots and save in Results folder
    print('Generate voxelwise correl plots for '+current)
    
    plt.scatter(FETroi_onPETFET_vxls, FETroi_onPETTER_vxls, s=2, c='m', marker='o')
    plt.xlim(0, 8)
    plt.ylim(0, 10)
    plt.axes().xaxis.set_major_locator(plt.MultipleLocator(1))
    plt.axes().yaxis.set_major_locator(plt.MultipleLocator(1))
    plt.axes().xaxis.set_minor_locator(plt.MultipleLocator(0.5))
    plt.axes().yaxis.set_minor_locator(plt.MultipleLocator(0.5))
    plt.xlabel('FET SUV [g mL\u207b\u00b9]')
    plt.ylabel('PSMA SUV [g mL\u207b\u00b9]')
    plt.savefig(results + '/'+ subj_name +'_FETroiCorrPlt.png')     
    plt.close()
    
    plt.scatter(FETroi_onPETFETgauss_vxls, FETroi_onPETTERgauss_vxls, s=2, c='c', marker='o')
    plt.xlim(0, 7)
    plt.ylim(0, 7)
    plt.axes().xaxis.set_major_locator(plt.MultipleLocator(1))
    plt.axes().yaxis.set_major_locator(plt.MultipleLocator(1))
    plt.axes().xaxis.set_minor_locator(plt.MultipleLocator(0.5))
    plt.axes().yaxis.set_minor_locator(plt.MultipleLocator(0.5))
    plt.xlabel('FET SUV [g mL\u207b\u00b9]')
    plt.ylabel('PSMA SUV [g mL\u207b\u00b9]')
    plt.savefig(results + '/'+ subj_name +'_FETroiCorrPlt_Gauss.png')
    plt.close()
    
    
#Save all dataframes to excel files here
print('Save all results to excel files')

CEtum_onT1CEdf.to_excel(Group_stats_results, sheet_name='CEtum_onT1CE', index=False)
CEtumFETCT_onT1CEdf.to_excel(Group_stats_results, sheet_name='CEtumFETCT_onT1CE', index=False)

FETroi_onFETPETdf.to_excel(Group_stats_results, sheet_name='FETroi_onFETPET', index=False)
FETroi_onTERPETdf.to_excel(Group_stats_results, sheet_name='FETroi_onTERPET', index=False)
CTRLFETroi_onFETPETdf.to_excel(Group_stats_results, sheet_name='CTRLFETroi_onFETPET', index=False)
CTRLFETroi_onTERPETdf.to_excel(Group_stats_results, sheet_name='CTRLFETroi_onTERPET', index=False)

TERroi_onFETPETdf.to_excel(Group_stats_results, sheet_name='TERroi_onFETPET', index=False)
TERroi_onTERPETdf.to_excel(Group_stats_results, sheet_name='TERroi_onTERPET', index=False)
CTRLTERroi_onFETPETdf.to_excel(Group_stats_results, sheet_name='CTRLTERroi_onFETPET', index=False)
CTRLTERroi_onTERPETdf.to_excel(Group_stats_results, sheet_name='CTRLTERroi_onTERPET', index=False)

TERroi_onTERPETorigdf.to_excel(Group_stats_results, sheet_name='TERroi_onTERPETorig', index=False)
CTRLTERroi_onTERPETorigdf.to_excel(Group_stats_results, sheet_name='CTRLTERroi_onTERPETorig', index=False)

CCRLFETroi_onFETPETdf.to_excel(Group_stats_results, sheet_name='CCRLFETroi_onFETPET', index=False)
CCRLFETroi_onTERPETdf.to_excel(Group_stats_results, sheet_name='CCRLFETroi_onTERPET', index=False)
CCRLTERroi_onTERPETdf.to_excel(Group_stats_results, sheet_name='CCRLTERroi_onTERPET', index=False)

FETminCE_onTERPETdf.to_excel(Group_stats_results, sheet_name='FETminCE_onTERPET', index=False)
FETmaxCE_onTERPETdf.to_excel(Group_stats_results, sheet_name='FETmaxCE_onTERPET', index=False)

TBRdf.to_excel(Group_stats_results, sheet_name='TBR', index=False)
VolRatiodf.to_excel(Group_stats_results, sheet_name='Volumes ratio', index=False)
prsdf.to_excel(voxelwise_results, sheet_name='Pearson correl coeff', index=False)

Group_stats_results.save()    
voxelwise_results.save()

#%%Execute statistical analysis
print('Perform statistical analysis')

import StatisticalAnalysisScript

#%%Generate result images for visualization purposes
#Remember to include z_slices.xcls file in data_supradir before running the following script
#import PSMA_Trial_ImgResultsVisual
#
##Put all results into a results folder
#os.mkdir(data_supradir +'/Trial Results') 
#trial_results = data_supradir +'/Trial Results'
#
#for file in glob.glob(data_supradir +'/' +'*.xlsx'):
#    shutil.move(os.path.join(data_supradir, file), trial_results)
#shutil.move(data_supradir + '/Analysis Results Plots', trial_results)
#shutil.move(data_supradir +'/Visualization Results Images', trial_results)
