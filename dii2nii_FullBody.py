# -*- coding: utf-8 -*-
"""
Created on Mon Feb 27 19:03:04 2023

@author: Caterina Brighi
"""

import os
import SimpleITK as sitk
import json


def DICOMseries_toNII(dcm_dir_path, nii_filepath, headers_filepath):
    '''This function converts DICOM series to nii images and saves also a json file with the dicom headers.'''
    reader = sitk.ImageSeriesReader()
    dicom_names = reader.GetGDCMSeriesFileNames(dcm_dir_path)
    reader.SetFileNames(dicom_names)
    image = reader.Execute()
    sitk.WriteImage(image, nii_filepath+'.nii')
    
    reader1 = sitk.ImageFileReader()
    reader1.SetFileName(dicom_names[0])
    reader1.LoadPrivateTagsOn()
    reader1.ReadImageInformation()    
    headers ={}
    for k in reader1.GetMetaDataKeys():
        v = reader1.GetMetaData(k)
        headers[k]=v
        # print(f'({k}) = "{v}"')
    
    with open(headers_filepath+'.json', 'w') as fp:
        json.dump(headers, fp)

    
data_dir = 'Path to data directory'

subjs_path = [ f.path for f in os.scandir(data_dir) if f.is_dir() ] #Create a list of the paths to the subjects directories
subjs_name = [ f.name for f in os.scandir(data_dir) if f.is_dir() ] #Create a list of subjects names

for current in subjs_name:
    
    subj_dir=data_dir+current
    DICOMseries_toNII(subj_dir, data_dir+'/'+current, data_dir+'/'+current)
