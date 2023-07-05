# -*- coding: utf-8 -*-
"""
Created on Mon Feb 27 08:41:56 2023

@author: Caterina Brighi
"""

#%% Import functions 

import matplotlib.pyplot as plt
import SimpleITK as sitk
import numpy as np
import pandas as pd
import os
from scipy.stats.stats import pearsonr

def allVoxInt(image, roi):
    
    '''Returns a flatten array (in 2D) of the intensity values from all pixels in roi applied to an image.'''
    
    image = sitk.GetArrayFromImage(image)
    roi = sitk.GetArrayFromImage(roi)
    mask = np.array(roi,dtype='bool')
    mask = np.invert(mask)
    masked = image[~np.array(mask)]
    return masked

#%% Set Working directory
        
data_supradir = 'Path to data directory' #Set working directory

subjs_path = [ f.path for f in os.scandir(data_supradir) if f.is_dir() ] #Create a list of the paths to the subjects directories
subjs_name = [ f.name for f in os.scandir(data_supradir) if f.is_dir() ] #Create a lisdt of subjects names

n_subj = len(subjs_name) #Total number of subjects

#Generate excel file containing voxelwise results
voxelwise_results = pd.ExcelWriter(data_supradir +'PSMA_FET_correlation.xlsx')

#Create empty dataframes to populate as going through the loop
prsdf = pd.DataFrame(columns=['Subject_ID', ' SUV Pearson corr coeff', ' SUV p-value', ' TBR Pearson corr coeff', ' TBR p-value']) #generate an empty dataframe to populate with Pearson correl coeff

#%%Create a for loop to perform image analysis on each subject sequentially

for current in subjs_name:
    
    subj_dir = data_supradir+current
    subj_name = current

    #Set paths to subfolder    
    reg_imgs = subj_dir +'/Registered'
    ROIs = subj_dir +'/ROIs'
    
    #Make a directory for new inverse dose images
    if not os.path.exists(data_supradir+current+'/Results'):#if it does not already exist, create a directory where theresult plots will be saved
        os.mkdir(data_supradir+current+'/Results')
    results = data_supradir+current+'/Results'

    #Read images
    PET_FET = sitk.ReadImage(reg_imgs +'/PET_FET_inFETCT.nii', sitk.sitkFloat32) #read PET FET image
    PET_TER = sitk.ReadImage(reg_imgs +'/PET_PSMA_inFETCT.nii', sitk.sitkFloat32) #read PET PSMA image   
    TBR_FET = sitk.ReadImage(reg_imgs +'/PET_FET_TBR_map.nii', sitk.sitkFloat32) #read PET FET TBR image
    TBR_TER = sitk.ReadImage(reg_imgs +'/PET_PSMA_TBR_map.nii', sitk.sitkFloat32) #read PET PSMA TBR image
    
    #Read VOIs
    FET_BTV = sitk.ReadImage(ROIs+'/PET_FET_BTV.nii', sitk.sitkUInt8) #Read FET BTV
    TER_BTV = sitk.ReadImage(ROIs+'/PET_PSMA_BTV.nii', sitk.sitkUInt8) #Read PSMA BTV

    ROI = (FET_BTV+TER_BTV)>1
    sitk.WriteImage(ROI, ROIs+'/Overlapping_BTV.nii')    
    

#%%Calculate pixel intensity values for all pixels in roi
    print('Calculate voxelwise stats for '+current)
    
    #Calculate Pearson correlation between FET and PSMA SUV
    FET_SUV = allVoxInt(PET_FET, ROI)
    PSMA_SUV = allVoxInt(PET_TER, ROI)    
    r_SUV, p_SUV = pearsonr(FET_SUV, PSMA_SUV)
    
    #Calculate Pearson correlation between FET and PSMA SUV
    FET_TBR = allVoxInt(TBR_FET, ROI)
    PSMA_TBR = allVoxInt(TBR_TER, ROI)    
    r_TBR, p_TBR = pearsonr(FET_TBR, PSMA_TBR)
    
#%%Export pixel intensity values and Pearson correl coeff to dataframe and then excel spreadsheet
    print('Export voxelwise stats for '+current)
    
    data = {'FET SUV in overlapping ROI': FET_SUV, 'PSMA SUV in overlapping ROI': PSMA_SUV, 'FET TBR in overlapping ROI': FET_TBR, 'PSMA TBR in overlapping ROI': PSMA_TBR}
    df = pd.DataFrame(data, columns=['FET SUV in overlapping ROI', 'PSMA SUV in overlapping ROI', 'FET TBR in overlapping ROI', 'PSMA TBR in overlapping ROI']) #generate dataframe of pixel intensity values
    df.to_excel(voxelwise_results, sheet_name=subj_name, index=False)
    
    pearsc = {'Subject_ID': subj_name, 'SUV Pearson corr coeff': r_SUV, 'SUV p-value': p_SUV, 'TBR Pearson corr coeff': r_TBR, 'TBR p-value': p_TBR}
    prsdf = prsdf.append(pearsc, ignore_index=True)
    
#%%Generate correlation plots and save in Results folder
    print('Generate voxelwise correl plots for '+current)

    # Plot regression line
    plt.scatter(FET_SUV, PSMA_SUV, s=1, c='g', marker='o')
    # plt.xlim(0, 6)
    # plt.ylim(0, 12)
    # plt.axes().xaxis.set_major_locator(plt.MultipleLocator(1))
    # plt.axes().yaxis.set_major_locator(plt.MultipleLocator(1))
    # plt.axes().xaxis.set_minor_locator(plt.MultipleLocator(0.5))
    # plt.axes().yaxis.set_minor_locator(plt.MultipleLocator(0.5))
    plt.xlim([min(FET_SUV)-0.2, max(FET_SUV)+0.5])
    plt.ylim([min(PSMA_SUV)-0.5, max(PSMA_SUV)+0.5])
    plt.xlabel('FET SUV [g mL\u207b\u00b9]',fontsize=20)
    plt.ylabel('PSMA SUV [g mL\u207b\u00b9]',fontsize=20)
    plt.xticks(fontsize=20)
    plt.yticks(fontsize=20)
    plt.tight_layout()
    plt.savefig(results + '/'+ subj_name +'_SUV_correl.png')     
    plt.close()
    # plt.show()
    
    plt.scatter(FET_TBR, PSMA_TBR, s=1, c='b', marker='o')
    # plt.xlim(0, 7)
    # plt.ylim(0, 80)
    # plt.axes().xaxis.set_major_locator(plt.MultipleLocator(1))
    # plt.axes().yaxis.set_major_locator(plt.MultipleLocator(1))
    # plt.axes().xaxis.set_minor_locator(plt.MultipleLocator(0.5))
    # plt.axes().yaxis.set_minor_locator(plt.MultipleLocator(0.5))
    plt.xlim([min(FET_TBR)-0.2, max(FET_TBR)+0.5])
    plt.ylim([min(PSMA_TBR)-0.5, max(PSMA_TBR)+0.5])
    plt.xlabel('FET TBR',fontsize=20)
    plt.ylabel('PSMA TBR',fontsize=20)
    plt.xticks(fontsize=20)
    plt.yticks(fontsize=20)
    plt.tight_layout()
    plt.savefig(results + '/'+ subj_name +'_TBR_correl.png')    
    plt.close()
    # plt.show()
        
#Save all dataframes to excel files here
print('Save all results to excel files')

prsdf.to_excel(voxelwise_results, sheet_name='Pearson correl coeff', index=False)
  
voxelwise_results.save()
