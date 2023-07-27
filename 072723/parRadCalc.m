%% Parallel Radiomics Calculations

%% IBSI Image Feature Calculation
% main compilation of the texture analysis functions. Can change the number
% of discrete gray levels, originally set at 6 for phantom. Included GLDM from USC
% Radiomics calculations, not delineated in IBSI core function
% https://ibsi.readthedocs.io/en/latest/03_Image_features.html#grey-level-run-length-based-features
clc;
clear;
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
radiomics = struct();

%% File Upload
% CSV
filename = inputdlg('Enter filename (CSV): ','Image Name and ROI List');
fileList = CSVImporter(string(filename));

tic
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
        windowedImage = double(niftiread('IBSIPhantom.nii'));
        image = windowedImage;
        ROI = windowedImage;
        nifti = false;
        winImg = windowedImage;
        numLevels = 6;
        binWidth = 0;
    end


    %% image upload and preprocessing

    if nifti == true
        % nifti image
        image = double(niftiread(imageSTR));
        imageInfoN = niftiinfo(imageSTR);
        cifar = 1;
        if strcmpi(ROISTR,'null') == true
            ROI = image;
            windowedImage = image;

        elseif cifar == 1
            ROI = image;
            windowedImage = image;
        else
            ROI = double(niftiread(ROISTR));
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

        radiomics.metrics{numFiles,1}.Image = imageSTR;
        radiomics.metrics{numFiles,1}.ROI = ROISTR;

    else

        radiomics.metrics{numFiles,1}.Image = imageSTR;
        radiomics.metrics{numFiles,1}.ROI = ROISTR;

    end

    radiomics.version = version;
    reqTab = struct2table(req);

%% Discretize

    % check if any val from text mat checkbox is 1, anything but intensity/histo stats

    minimum = min(min(min(image)));
    maximum = max(max(max(image)));
    range = maximum - minimum;
    % % DLG
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
                binWidth = 0;
                radiomics.metrics{numFiles,1}.disAlg.type = 'FBN';
                radiomics.metrics{numFiles,1}.disAlg.numLevels = numLevels;
                radiomics.metrics{numFiles,1}.disAlg.binWidth = ceil((max(max(max(winImg)))-min(min(min(winImg))))/numLevels);
            case 'Fixed Bin Size'
                binWidth = str2double(inputdlg('Bin width: '));
                winImg = FBS(windowedImage,binWidth);
                numLevels = 0; % max of the discretized image = number of discretized bins
                radiomics.metrics{numFiles,1}.disAlg.type = 'FBS';
                radiomics.metrics{numFiles,1}.disAlg.numLevels = max(max(max(winImg)));w
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

    %% Feature Calculation Based on Requested Data
    % tied to feature GUI
    
    funcList = {@NGLDM,@NGTDM1,@GLDZM,@GLSZM,@GLRLM,@GLCM,@HistoStats,@IntensityStats};
    parfor reqTabInd = 1:size(funcList,2)
        if reqTab{1,reqTabInd} == 1
            result{reqTabInd} = funcList{reqTabInd}(winImg,numLevels,binWidth);
        end
    end


end % end line for file loop, use if CSV import

toc
%% Save/Output Features

% save .mat file, struct form
save("radiomics.mat","radiomics");

outputFileName = string(inputdlg('Name out output file (.csv)?','Output Title'));
csvExporter(fileList,req,radiomics,outputFileName);
nonTextOutput(fileList,radiomics,req,outputFileName);
