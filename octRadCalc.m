%% testing OCTAVE code
% download the image package
% download the statistics package
% include the GrayscaleQuantification folder and subfolders in your directory
%clc;
%clear;

digitalPhantom = zeros(4,5,4);
digitalPhantom(:,:,1) = [1,4,4,1,1; 1, 4, 6, 1, 1; 4,1,6,4,1; 4,4,6,4,1];
digitalPhantom(:,:,2) =[1,4,4,1,1; 1,1,6,1,1; NaN,1,3,1,1; 4,4,6,1,1];
digitalPhantom(:,:,3) =[1,4,4,NaN,NaN; 1,1,1,1,1;1,1,NaN,1,1; 1,1,6,1,1];
digitalPhantom(:,:,4) =[1,4,4,NaN,NaN;1,1,1,1,1;1,1,1,1,1;1,1,6,1,1];




%glrlmResults = octGLRLM(digitalPhantom);
%gldzmResults = octGLDZM(digitalPhantom);
%glcmResults = octGLCM(digitalPhantom);
intensityResults = octIntensity(digitalPhantom);
histogramResults = octHisto(digitalPhantom,6,0);
glszmResults = octGLSZM(digitalPhantom);
ngtdmResults = octNGTDM(digitalPhantom);
ngldmResults = octNGLDM(digitalPhantom);
