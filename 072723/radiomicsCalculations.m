%% IBSI Image Feature Calculation
% main compilation of the texture analysis functions. Can change the number
% of discrete gray levels, originally set at 6 for phantom. Included GLDM from USC
% Radiomics calculations, not delineated in IBSI core function
% https://ibsi.readthedocs.io/en/latest/03_Image_features.html#grey-level-run-length-based-features
version = '1.1'; % version of code package

%% Feature Calculation GUI
textCase = 'IBSI'; % either 'IBSI' or 'All'
ImageTechnique = ImageProcessingTechniqueIBSIGUI('IBSI');

reqLabels = ImageTechnique.descript;
reqVals = num2cell(ImageTechnique.val);
for name = 1:size(reqLabels)
    indivLabel = reqLabels{name};
    req.(indivLabel) = reqVals{name};
end
% checkbox GUI for final data set, initialize "req", or request, structure tie struct format to check box conditions

tic
%% File Upload
% CSV
filename = inputdlg('Enter filename (CSV): ','Image Name and ROI List');
fileList = CSVImporter(string(filename));
radiomics = struct();
for numFiles = 1:size(fileList,1) % loop through each of the files
    imageSTRT = (fileList(numFiles,1));
    ROISTRT = (fileList(numFiles,2));
    imageSTR = string(imageSTRT.(1){1});
    ROISTR = string(ROISTRT.(1){1});

    nifti = true;
    dicom = false;

    radiomics.info.Creator = getenv('USER');
    radiomics.info.date = datetime;

    % PHANTOM CHECK, CSV
    if strcmpi(imageSTR,'phantom') == 1 || strcmpi(ROISTR,'phantom') == 1
        windowedImage = niftiread('IBSIPhantom.nii');
        image = windowedImage;
        ROI = windowedImage;
        nifti = false;
        winImg = windowedImage;

        %
    end
    %% image upload and preprocessing

    if nifti == true
        % nifti image
        image = niftiread(imageSTR);
        imageInfoN = niftiinfo(imageSTR);
        cifar = 1;
        if strcmpi(ROISTR,'null') == true
            ROI = image;
            windowedImage = image;
        elseif cifar == 1 % if direct img/gif upload
            ROI = image;
            windowedImage = image;
        else
            ROI = niftiread(ROISTR);
            ROIinfoN = niftiinfo(ROISTR);
            % Check Orientation
            var1 = imageInfoN.Transform.T;
            var2 = ROIinfoN.Transform.T;

            if var1 ~= var2
                print('Error, orientation of ROI incorrect')
            else
                windowedImage = imageOverlay(image,ROI);
            end

        end
          ROI = image;
          windowedImage = image;
        radiomics.metrics{numFiles,1}.Image = imageSTR;
        radiomics.metrics{numFiles,1}.ROI = ROISTR;
    % 
    elseif dicom == true % not used extensively @ USC, may need to edit

        newInfo = dicominfo(imageSTR);
        ROIdcm = double(dicomread(ROISTR));
        valueFile = double(dicomread(imageSTR));
        ROIdcm(ROIdcm~=0) = 1; % convert to a binary mask

        for k = 1: size(valueFile,3)
            windowedImage(:,:,k) = ROIdcm.*valueFile(:,:,k);
        end

        image = windowedImage;
        % will want to discretize w/ FBS & 25 HU
        radiomics.metrics{numFiles,1}.Image = imageSTR;
        radiomics.metrics{numFiles,1}.ROI = ROISTR;

    else
        radiomics.metrics{numFiles,1}.Image = imageSTR;
        radiomics.metrics{numFiles,1}.ROI = ROISTR;

    end
    radiomics.version = version;


    %% Discretize

    % check if any val from text mat checkbox is 1, anything but intensity/histo stats

    minimum = min(min(min(image)));
    maximum = max(max(max(image)));
    range = maximum - minimum;
    % % % DLG
    textMatStructs = [req.GLCM, req.GLRLM, req.NGTDM, req.NGLDM,req.GLSZM, req.GLDZM];
    textTest = sum(textMatStructs);
    if (textTest > 0)
        disAlg = questdlg( ...
            'Choose discretization algorithm: ','Discretization','Fixed Bin Number','Fixed Bin Size','Fixed Bin Size');
        % %
        switch disAlg
            case 'Fixed Bin Number'
                numLevels = str2double(inputdlg('Number of bins (number of grey levels): '));
                winImg = FBN(windowedImage,numLevels);
                binWidth = ceil((max(max(max(winImg)))-min(min(min(winImg))))/numLevels);
                radiomics.metrics{numFiles,1}.disAlg.type = 'FBN';
                radiomics.metrics{numFiles,1}.disAlg.numLevels = numLevels;
                radiomics.metrics{numFiles,1}.disAlg.binWidth = binWidth;
            case 'Fixed Bin Size'
                binWidth = str2double(inputdlg('Bin width: '));
                winImg = FBS(windowedImage,binWidth);
                numLevels = max(max(max(winImg))); % max of the discretized image = number of discretized bins
                radiomics.metrics{numFiles,1}.disAlg.type = 'FBS';
                radiomics.metrics{numFiles,1}.disAlg.numLevels = numLevels;
                radiomics.metrics{numFiles,1}.disAlg.binWidth = binWidth;
        end
    else
        numLevels = max(max(max(image))); % for regular runs
        binWidth = 0; % so the histogram stats will register as needing to
        % %     discretize, not used in other functions
        radiomics.metrics{numFiles,1}.disAlg.type = 'Not discretized';
        radiomics.metrics{numFiles,1}.disAlg.numLevels = numLevels;
        radiomics.metrics{numFiles,1}.disAlg.binWidth = ceil(range/numLevels);
        winImg = windowedImage;
        % defaults to FBN for canceled window
    end


    %% BRODATZ
    % used to bypass discretize section, if want all the same discretization algorithm and
    % number of bins

    % windowedImage = double(windowedImage);
    % image = double(image);
    % binWidth = 0.2;
    % winImg = FBS(windowedImage,binWidth);
    % numLevels = max(max(max(winImg)));
    % 
    % 
    % radiomics.metrics{numFiles,1}.disAlg.type = 'FBS';
    % radiomics.metrics{numFiles,1}.disAlg.numLevels = numLevels;
    % radiomics.metrics{numFiles,1}.disAlg.binWidth = binWidth;
    % 

    %% CIFAR
    % used to bypass discretize section, if want all the same discretization algorithm and
    % number of bins

    % windowedImage = double(windowedImage);
    % image = double(image);
    % numLevels = 20;
    % winImg = FBN(windowedImage,numLevels);
    % binWidth = ceil(range/numLevels);
    % 
    % radiomics.metrics{numFiles,1}.disAlg.type = 'FBN';
    % radiomics.metrics{numFiles,1}.disAlg.numLevels = numLevels;
    % radiomics.metrics{numFiles,1}.disAlg.binWidth = binWidth;

    %% image phantom used in IBSI database

    % ROIlayer1 = [1 4 4 1 1; 1 4 6 1 1; 4 1 6 4 1; 4 4 6 4 1];
    % ROIlayer2 = [1 4 4 1 1; 1 1 6 1 1; NaN 1 3 1 1; 4 4 6 1 1];
    % ROIlayer3 = [1 4 4 NaN NaN; 1 1 1 1 1; 1 1 NaN 1 1; 1 1 6 1 1];
    % ROIlayer4 = [1 4 4 NaN NaN; 1 1 1 1 1; 1 1 1 1 1; 1 1 6 1 1];
    % % %
    % phantom(:,:,1) = ROIlayer1;
    % phantom(:,:,2) = ROIlayer2;
    % phantom(:,:,3) = ROIlayer3;
    % phantom(:,:,4) = ROIlayer4;

    % create nifti file for phantom, only needs to be done once
    % niftiwrite(phantom,'IBSIPhantom');

    %% Feature Calculation Based on Requested Data
    % tied to feature GUI

    % calculate intensity statistics
    if req.Intensity == 1
        radiomics.metrics{numFiles,1}.intensity = IntensityStats(windowedImage);
    end

    % calculate histogram statistics
    if req.Histogram == 1

        if strcmpi(radiomics.metrics{numFiles,1}.disAlg,'FBS') == 1
            radiomics.metrics{numFiles,1}.histogram = HistoStats(windowedImage,0,binWidth);
        elseif strcmpi(radiomics.metrics{numFiles,1}.disAlg,'FBN') == 1
            radiomics.metrics{numFiles,1}.histogram = HistoStats(winImg,numLevels,0);
        else
            radiomics.metrics{numFiles,1}.histogram = HistoStats(windowedImage,6,0); % set one of the two parameters to zero, depending on whether or not it is user defined
        end
    end

    % calculate GLCM statistics
    if req.GLCM == 1
        radiomics.metrics{numFiles,1}.GLCM =GLCM(winImg,numLevels,binWidth);
    end

    % calculate GLRLM statistics
    if req.GLRLM == 1
        radiomics.metrics{numFiles,1}.GLRLM = GLRLM(winImg,numLevels,binWidth);
    end

    % calculate GLSZM statistics
    if req.GLSZM == 1
        radiomics.metrics{numFiles,1}.GLSZM = GLSZM(winImg,numLevels,binWidth);
    end

    % calculate GLDZM statistics
    if req.GLDZM == 1
        radiomics.metrics{numFiles,1}.GLDZM = GLDZM(winImg,numLevels,binWidth);
    end

    % calculate NGTDM statistics
    if req.NGTDM == 1
        radiomics.metrics{numFiles,1}.NGTDM = NGTDM1(winImg,numLevels,binWidth);
    end

    % calculate NGLDM statistics
    if req.NGLDM == 1
        radiomics.metrics{numFiles,1}.NGLDM = NGLDM(winImg,numLevels,binWidth);
    end

    % calculate GLDM statistics
    if strcmpi(textCase,'All') == 1 && req.GLDM == 1
        radiomics.metrics{numFiles,1}.GLDM = GLDM(winImg,numLevels);
    end


end % end line for file loop, use if CSV import
toc
%% Save/Output Features

% save .mat file, struct form
save("radiomics.mat","radiomics");

outputFileName = string(inputdlg('Name out output file (.csv)?','Output Title'));
csvExporter(fileList,req,radiomics,outputFileName);