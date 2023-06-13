%% IBSI Texture Analysis
% main compilation of the texture analysis functions. Can change the number
% of discrete gray levels, originally set at 6. Included GLDM from USC
% Radiomics calculations, not delineated in IBSI core function
% https://ibsi.readthedocs.io/en/latest/03_Image_features.html#grey-level-run-length-based-features
clc;
clear;
%% Obtain User Input

% filename grab - preexisting function, check grayscale quant
filename = 'roiImageTitleTest.xlsx';
opts = detectImportOptions(filename);
imageVar = opts.SelectedVariableNames{1};
roiVar = opts.SelectedVariableNames{2};

% ImageList = readmatrix(filename,opts,'Range','A:A');
% ROIList = readmatrix(filename,opts,'Range','B:B');

% checkbox GUI for final data set, initialize "req", or request, structure
% tie struct format to check box conditions
req.intens = 1;
req.histo = 1; 
req.GLCM = 1;
req.GLRLM = 1;
req.GLDZM = 1;
req.GLSZM = 1;
req.GLDM = 1;
req.NGTDM = 1;
req.NGLDM = 1; 

%% Feature Calculation GUI

ImageTechnique = ImageProcessingTechniqueGUI('IBSI');

%% Image Upload

image = niftiread('ROI_602_Tumor_a.nii'); % change to string var ver of input image
ROI = niftiread('r0501_ARTERIAL.nii'); % change to string var ver of input ROI

imageInfoN = niftiinfo('ROI_602_Tumor_a.nii'); % change to string var ver of input image
ROIinfoN = niftiinfo('r0501_ARTERIAL.nii'); % change to string var ver of input ROI

% check orientation
var1 = imageInfoN.Transform.T;
var2 = ROIinfoN.Transform.T;

if var1 ~= var2
    print('Error, orientation of ROI incorrect')
else
    windowedImage = imageOverlay(image,ROI);
end


%% Discretize

% check if any val from text mat checkbox is 1, anything but
% intensity/histo stats (histo might involve discretization, double check

minimum = min(min(min(image)));
maximum = max(max(max(image)));
range = maximum - minimum;
% DLG
 textMatStructs = [req.GLCM, req.GLRLM, req.NGTDM, req.NGLDM,req.GLSZM, req.GLDZM, req.GLDM];
 textTest = sum(textMatStructs);
  if (textTest > 0)
      disAlg = questdlg( ...
          'Choose discretization algorithm: ','Discretization','Fixed Bin Number','Fixed Bin Size','Fixed Bin Size');
% % 
      switch disAlg
          case 'Fixed Bin Number'
              numLevels = str2double(inputdlg('Number of bins (number of grey levels): '));
              binWidth = ceil(range/numLevels);
              winImg = FBN(windowedImage,numLevels);
              radiomicsNFBN.disAlg = 'FBN';
          case 'Fixed Bin Size'
              binWidth = str2double(inputdlg('Bin width: '));
              numLevels = ceil(range/binWidth);
              winImg = FBS(windowedImage,binWidth);
              radiomicsNFBS.disAlg = 'FBS';
      end
 else
      numLevels = max(max(max(image))); % for regular runs
      binWidth = 0; % so the histogram stats will register as needing to
% %     discretize, not used in other conditionals
       % defaults to FBN for canceled window
end

numLevels = 10;
binWidth = ceil(range/numLevels);
disAlg = 'Fixed Bin Size';
% based on number of levels
winImgFBN = FBN(windowedImage,numLevels); % change uploaded image here

% based on a fixed bin width
winImgFBS = FBS(windowedImage,binWidth);

% store type of discretization algorithm and number bins/levels
radiomicsNFBS.numLevels = ceil(range/binWidth);
radiomicsNFBS.binWidth = binWidth;

radiomicsNFBN.numLevels = numLevels;
radiomicsNFBN.binWidth = binWidth;

radiomics.numLevels = 6;
radiomics.binWidth = binWidth;

%% image phantom used in IBSI database

ROIlayer1 = [1 4 4 1 1; 1 4 6 1 1; 4 1 6 4 1; 4 4 6 4 1];
ROIlayer2 = [1 4 4 1 1; 1 1 6 1 1; NaN 1 3 1 1; 4 4 6 1 1];
ROIlayer3 = [1 4 4 NaN NaN; 1 1 1 1 1; 1 1 NaN 1 1; 1 1 6 1 1];
ROIlayer4 = [1 4 4 NaN NaN; 1 1 1 1 1; 1 1 1 1 1; 1 1 6 1 1];
% % 
phantom(:,:,1) = ROIlayer1;
phantom(:,:,2) = ROIlayer2; 
phantom(:,:,3) = ROIlayer3; 
phantom(:,:,4) = ROIlayer4; 

% create nifti file for phantom, only needs to be done once 
%niftiwrite(phantom,'IBSIPhantom'); 

%% Calculations for IBSI Phantom

% Aggregating Methods Summary for GLCM and GLRLM
% Rows 1-6 correspond to:
% 1 -- 2D, by slice, w/o merging
% 2 -- 2D, by slice, w/ merging by slice
% 3 -- 2.5D, by slice, w/ merging by direction
% 4 -- 2.5D, by slice, w/ full merging
% 5 -- 3D, as volume, w/o merging
% 6 -- 3D, as volume, w/ full merging
featuresGLCM = GLCM(phantom,6); %gray level co-occurance matrix, data in form of struct
featuresGLRLM = GLRLM(phantom,6); %gray level run length matrix

% Aggregating Methods Summary for GLSZM, GLDZM, NGTDM, NGLDM
% Rows 1-3 correspond to:
% 1 -- 2D, by slice, without merging
% 2 -- 2.5D, by slice, with merging
% 3 -- 3D: as volume

featuresGLSZM = GLSZM(phantom,6); % gray level size zone matrix
featuresGLDZM = GLDZM(phantom,6); % gray level distance zone matrix
featuresNGTDM = NGTDM(phantom,6); % neighborhood grey tone difference matrix
featuresNGLDM = NGLDM(phantom,6); % neighborhood grey level difference matrix
% 
% USC specific, row 1 corresponds to 2D features and row 2 to 3D features
featuresGLDM = GLDM(phantom,6); % gray level dependence matrix, USC only
% 

%% Output/Data Export

% calculate intensity statistics
if req.intens == 1
    radiomics.intensity = IntensityStats(phantom);
    radiomicsNFBN.intensity = IntensityStats(winImgFBN);
    radiomicsNFBS.intensity = IntensityStats(winImgFBS);
end

% calculate histogram statistics
if req.histo == 1
    radiomics.histogram = HistoStats(phantom,6,0); % set one of the two parameters to zero, depending on whether or not it is user defined

    if strcmpi(radiomicsNFBN.disAlg,'FBS') == 1
    radiomicsNFBS.histogram = HistoStats(windowedImage,0,binWidth);
    elseif strcmpi(radiomicsNFBS.disAlg,'FBN') == 1
    radiomicsNFBN.histogram = HistoStats(windowedImage,numLevels,0);

    end
end

% calculate GLCM statistics
if req.GLCM == 1
    radiomics.GLCM =featuresGLCM;
    radiomicsNFBN.GLCM = GLCM(winImgFBN,numLevels);
    radiomicsNFBS.GLCM = GLCM(winImgFBS,numLevels);
end

% calculate GLRLM statistics
if req.GLRLM == 1
    radiomics.GLRLM = featuresGLRLM;
    radiomicsNFBN.GLRLM = GLRLM(winImgFBN,numLevels);
    radiomicsNFBS.GLRLM = GLRLM(winImgFBS,numLevels);
end

% calculate GLSZM statistics
if req.GLSZM == 1
    radiomics.GLSZM = featuresGLSZM;
    radiomicsNFBN.GLSZM = GLSZM(winImgFBN,numLevels);
    radiomicsNFBS.GLSZM = GLSZM(winImgFBS,numLevels);
end

% calculate GLDZM statistics
if req.GLDZM == 1
    radiomics.GLDZM = featuresGLDZM;
    radiomicsNFBN.GLDZM = GLDZM(winImgFBN,numLevels);
    radiomicsNFBS.GLDZM = GLDZM(winImgFBS,numLevels);
end

% calculate NGTDM statistics
if req.NGTDM == 1
    radiomics.NGTDM = featuresNGTDM;
    radiomicsNFBN.NGTDM = NGTDM(winImgFBN,numLevels);
    radiomicsNFBS.NGTDM = NGTDM(winImgFBS,numLevels);
end

% calculate NGLDM statistics
if req.NGLDM == 1
    radiomics.NGLDM = featuresNGLDM;
    radiomicsNFBN.NGLDM = NGLDM(winImgFBN,numLevels);
    radiomicsNFBS.NGLDM = NGLDM(winImgFBS,numLevels);
end

% calculate GLDM statistics
if req.GLDM == 1
    radiomics.GLDM = GLDM(phantom,6);
    radiomicsNFBN.GLDM = GLDM(winImgFBN,numLevels);
    radiomicsNFBS.GLDM = GLDM(winImgFBS,numLevels);
end

% save .mat file, struct form
save("radiomics.mat","radiomics");
save("radiomicsNFBN.mat","radiomicsNFBN");
save("radiomicsNFBS.mat","radiomicsNFBS");

% 

%% Discretization Algorithms

% discretization based on fixed number of bins
function disImg = FBN(image,numLevels)
minimum = min(min(min(image)));
maximum = max(max(max(image)));
range = maximum - minimum;
interval = ceil(range/numLevels);
image(image > (minimum+(interval*(numLevels-1)))) = numLevels;
image(image <= minimum) = 1;

edges = [minimum:interval:maximum];
disImg = discretize(image,edges);

end


function disImg = FBS(image,binWidth)
minimum = min(min(min(image)));
maximum = max(max(max(image)));
range = maximum - minimum;
numLevels = ceil(range/binWidth);
interval = binWidth;
% 
edges = [minimum:interval:maximum];
disImg = discretize(image,edges);

end

%% ROI Segmentation

function windowedImage = imageOverlay(actual,ROI)
% initial image overlay, isolates ROI, NaN fills gaps
roiInfo=struct();
roiInfo.img = ROI;
roiInfo.roiSize=size(roiInfo.img);

% Converts ROI file to binary mask
if min(ROI(:))<0
    roiInd=find(ROI~=min(ROI(:)));
else
    roiInd=find(ROI>0);
end
roiInfo.roiInd=roiInd;

valueFile = actual;
tempImg=ROIImageCombiner(roiInfo,valueFile,'3D');


function tempImg =ROIImageCombiner(roiInfo,value_file,type)
valueImg=actual;
% Tests to see if there is a size mismatch between the ROI and Vaule files.
if sum(roiInfo.roiSize-size(valueImg))~=0
    errordlg('ROI and Value Images are mismatched in size.');
    imageStack=[];
    return;
end
roiInd=roiInfo.roiInd;
tempImg=nan(roiInfo.roiSize);
tempImg(roiInd)=valueImg(roiInd);
end

%
 % assuming ROI read in as 0's for excluded voxels

% limit window 
for k=1:3
    tempROI=shiftdim(tempImg,k);
    roiInd=find(tempROI>0);
    tempVal=shiftdim(tempImg,k);
    
    tempImg=nan(size(tempROI));
    tempImg(roiInd)=tempVal(roiInd);

    [indX,indY,indZ]=ind2sub(size(tempROI),roiInd);
end
        

windowedImage(:,:,:)=tempImg(min(indX):max(indX),min(indY):max(indY),...
       min(indZ):max(indZ));


end

%% GUI Functions


function ImageTechnique=ImageProcessingTechniqueGUI(type)
 
ImageTechnique=[];
 %if nargin==0
 PosID=[];
 %end
 techList=GUIList(type);
 listWidth=DescriptSize(techList);
 if listWidth<20
     listWidth=20;
 end
 fieldNum=length(techList);
 set(0,'Units','characters');
 scr_size=get(0,'ScreenSize');
 x_size=(listWidth+10);
 y_size=fieldNum*2+4;
% 
 main_fig = figure('Visible','on','Name','Select Technique',...
     'Units','characters','Position',[scr_size(3)/2-(listWidth+10)/2,...
     scr_size(4)/2-(fieldNum*2+4)/2,x_size,...
     y_size],'MenuBar',...
     'none','BusyAction','cancel','NumberTitle','off');
 if nargin~=0
     set(main_fig,'Name',['Select ',PosID]);
 end
% 
 counter=fieldNum;
% 
 for k=1:fieldNum
     c(k)=uicontrol('Style','checkbox','Units','characters',...
         'position',[1 (k-1)*2+3.5 5 2],'parent',main_fig,'Min',0,...
         'Max',1,'Value',techList(counter).def);
     u(k)= uicontrol('Style','text','String',...
         sprintf('%s',techList(counter).descript),'Units','characters',...
         'pos',[6 (k-1)*2+3 listWidth+4 2],'parent',main_fig,'HandleVisibility',...
         'off','Tag',num2str(counter),'HorizontalAlignment','left');
     counter=counter-1;
     uicontrol('Style','pushbutton','String','OK','Units','characters',...
         'pos',[1 .5 x_size/2-1.5 2],'parent',main_fig,'Callback',@btn_Callback);
     uicontrol('Style','pushbutton','String','Cancel','Units','characters',...
         'pos',[x_size/2-0.5 .5 x_size/2-1.5 2],'parent',main_fig,'Callback',@btn_Cancel);
 end
% 
     uiwait(gcf);
% 
     function btn_Callback(~,~)
         cnt=1;
         for q=1:length(c)
             if get(c(q),'Value')==1
                 ImageTechnique=[ImageTechnique; techList(str2double(get(u(q),'Tag')))];
                 cnt=cnt+1;
             end
         end
         if ~isempty(ImageTechnique)
             ImageTechnique=struct2table(ImageTechnique);
         end
         close(main_fig);
     end
     function btn_Cancel(~,~)
         ImageTechnique=[];
         close(main_fig);
     end
 
 end
% 
 function techList=GUIList(type)
 counter=1;
 switch type
     case 'All'
         techList(counter).descript='Intensity';
         techList(counter).tech='Intensity';
         techList(counter).type='All';
         techList(counter).def=1;
         counter=counter+1;
         techList(counter).descript='Histogram';
         techList(counter).tech='GLCM';
         techList(counter).type='3D';
         techList(counter).def=1;
         counter=counter+1;
        techList(counter).descript='GLCM';
        techList(counter).tech='GLCM';
         techList(counter).type='6';
         techList(counter).def=1;
         counter=counter+1;
         techList(counter).descript='GLRLM';
         techList(counter).tech='GLRLM';
         techList(counter).type='6';
         techList(counter).def=1;
         counter=counter+1;
         techList(counter).descript='GLSZM';
         techList(counter).tech='GLSZM';
         techList(counter).type='3';
         techList(counter).def=1;
         counter=counter+1;
         techList(counter).descript='GLDZM';
         techList(counter).tech='GLDZM';
         techList(counter).type='3';
         techList(counter).def=1;
         counter=counter+1;
         techList(counter).descript='NGTDM';
         techList(counter).tech='NGTDM';
         techList(counter).type='3';
         techList(counter).def=1;
         counter=counter+1;
         techList(counter).descript='NGLDM';
         techList(counter).tech='NGLDM';
         techList(counter).type='3';
         techList(counter).def=1;
         counter=counter+1;
         techList(counter).descript='GLDM';
         techList(counter).tech='GLDM';
         techList(counter).type='2';
         techList(counter).def=1;

% 
     case 'IBSI'
%         techList(counter).descript='Intensity';
         techList(counter).tech='Enhance';
         techList(counter).type='All';
         techList(counter).def=1;
         counter=counter+1;
         techList(counter).descript='Histogram';
         techList(counter).tech='GLCM';
         techList(counter).type='3D';
         techList(counter).def=1;
         counter=counter+1;
        techList(counter).descript='GLCM';
        techList(counter).tech='GLCM';
         techList(counter).type='6';
         techList(counter).def=1;
         counter=counter+1;
         techList(counter).descript='GLRLM';
         techList(counter).tech='GLRLM';
         techList(counter).type='6';
         techList(counter).def=1;
         counter=counter+1;
         techList(counter).descript='GLSZM';
         techList(counter).tech='GLSZM';
         techList(counter).type='3';
         techList(counter).def=1;
         counter=counter+1;
         techList(counter).descript='GLDZM';
         techList(counter).tech='GLDZM';
         techList(counter).type='3';
         techList(counter).def=1;
         counter=counter+1;
         techList(counter).descript='NGTDM';
         techList(counter).tech='NGTDM';
         techList(counter).type='3';
         techList(counter).def=1;
         counter=counter+1;
         techList(counter).descript='NGLDM';
         techList(counter).tech='NGLDM';
         techList(counter).type='3';
         techList(counter).def=1;
         counter=counter+1;
% 
 end
 end

function MaxSize=DescriptSize(techList)
MaxSize=0;
for k=1:length(techList)
    temp=length(techList(k).descript);
    if MaxSize<temp
        MaxSize=temp;
    end
end
end
    
