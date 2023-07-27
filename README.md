# USC-Radiomics
Code repository for the USC Radiomics Lab
This is a MATLAB based calculation of IBSI compliant image features. A binary ROI segmentation method is used in the image overlay function. 
parRadCalc.m is the main script that calls the individual functions, which are separated by feature groups Intensity, Intensity-Volume Histogram,
Gray Level Co-occurence Matrix (GLCM), Gray Level Run Length Matrix (GLRLM), Gray Level Distance Zone Matrix (GLDZM), Gray Level Size Zone Matrix (GLSZM), 
Neighborhood Grey Tone Difference Matrix (NGTDM), and Neighborhood Grey Level Dependence Matrix (NGTDM). Note, the subfolder GLRL was authored by Xunkai Wei
of the Beijing Aeronautical Technology Research Center and is used in the GLRLM function. 

All function and feature calculations are based on the following directory: https://ibsi.readthedocs.io/en/latest/03_Image_features.html#. 
