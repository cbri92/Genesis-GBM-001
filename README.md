# Code for analysis of Genesis GBM 001 phase I/II study

**Author:** *Caterina Brighi*

This repository contains python scripts used for the analysis of the Genesis GBM 001 phase I/II study, which aimed at evaluating the role of [68Ga]Ga-PSMA-617 as an imaging biomarker in adult recurrent glioblastoma patients.

## Setup/Build/Install

Beside standard python packages for data analysis (numpy, pandas, glob, matplotlib, scipy, etc...) running the code will require the following image analysis packages, which will have to be downloaded in python:
- SimpleITK
- dicom2nifti
- json

## Usage

To run most of the script, make sure to have in the same repository from which the script is run the following two scripts, containing general image analysis and visualisation functions:
- *ImageAnalysisFunctions.py*
- *ImageResultsVisualizationFunctions.py*

## Citation

If you use this repository in your work, please cite 10.5281/zenodo.8124206 
