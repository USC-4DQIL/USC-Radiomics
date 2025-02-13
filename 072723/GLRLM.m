%% Gray Level Run Length Matrix (GLRLM)
function GLRLMFeatures = GLRLM(winImg,numLevels,binWidth)

numLevels = max(max(max(winImg)));
featuresA = GLRLMa(winImg,numLevels); % 2D, by slice, w/o merging
featuresB = GLRLMb(winImg,numLevels); % 2D, by slice, w/ merging by slice
featuresC = GLRLMc(winImg,numLevels); % 2.5D, by slice, w/ merging by direction
featuresD = GLRLMd(winImg,numLevels); % 2.5D, by slice, w/ full merging
featuresE = GLRLMe(winImg,numLevels);  % 3D, as volume, w/o merging
featuresF = GLRLMf(winImg,numLevels);  % 3D, as volume, w/ full merging

featuresAF = [featuresA,featuresB,featuresC,featuresD,featuresE,featuresF];

agMethods = {'A','B','C','D','E','F'};
conCat = [agMethods' struct2table(featuresAF)]; % labeling step
GLRLMFeatures = renamevars(conCat,'Var1','Aggregating_Methods');

%USCGLRLM2D =  GLRLM_auto2D(winImg(:,:,1),numLevels);
%USCGLRLM3D = GLRLM_auto3D(winImg,numLevels);
%featuresG = GLRLMStats(USCGLRLM3D,Nv); % based on 3D averaging ignoring the corner voxels

%% IBSI Aggregating Methods

% a -- 2D, by slice, w/o merging
% b -- 2D, by slice, w/ merging by slice
% c -- 2.5D, by slice, w/ merging by direction
% d -- 2.5D, by slice, w/ full merging
% e -- 3D, as volume, w/o merging
% f -- 3D, as volume, w/ full merging

    function featuresA = GLRLMa(winImg,numLevels)
        offset=[3; 1; 2; 4];
        zDim = size(winImg,3);
        numDir = size(offset,1);
        matrixA = zeros(16,numDir);

        for k = 1:zDim
            slice2D = winImg(:,:,k);
            for q=1:size(offset,1)
                x = numel(slice2D(:));
                voxelstoRemove = sum(isnan(slice2D(:)));
                Nv = x - voxelstoRemove;
                temp=grayrlmatrix(slice2D,'GrayLimits',[min(winImg(:)),...
                    max(winImg(:))],'NumLevels',numLevels,'Offset',offset(q,:));
                tempGLRLM = temp{1};
                fresult = GLRLMStats(tempGLRLM,Nv);
                biomarkers = fieldnames(fresult);
                cellResults = cell2mat(struct2cell(fresult));
                matrixA(:,q) = cellResults; % stores directional matrices, column 1 == direction 1, etc.
            end
            averageofSlice = sum(matrixA,2)./numDir;
            matrixB(:,k) = averageofSlice;
        end
        average = sum(matrixB,2)./zDim;
        dim = size(average,1);
        intermediate = num2cell(average,dim);
        featuresA = cell2struct(intermediate,biomarkers,1);
    end

    function featuresB = GLRLMb(winImg,numLevels)
        offset=[3; 1; 2; 4];
        zDim = size(winImg,3);
        dimens = [size(winImg,1);size(winImg,2);size(winImg,3)];
        maxRunLength = ceil(sqrt(2)*max(dimens));
        numDir = size(offset,1);
        for k = 1:zDim
            mergeGLRLM = zeros(numLevels,maxRunLength); % intentionally too large, zero pad for ability to add.
            slice2D = winImg(:,:,k);
            for q=1:numDir
                temp=grayrlmatrix(slice2D,'GrayLimits',[min(winImg(:)),...
                    max(winImg(:))],'NumLevels',numLevels,'Offset',offset(q,:));
                tempGLRLM = temp{1};

                % change to make it possible to add directly
                % parse through padded GLRLM
                rowSum = sum(tempGLRLM,1); % for runs
                counter = 0;
                % remove excess
                for runIndex = size(rowSum,2):-1:1
                    if (rowSum(1,runIndex) == 0) && counter == 0
                        tempGLRLM(:,runIndex) = [];
                    else
                        counter = 1;
                    end
                end
                % add extra
                tempMaxRunLength = size(tempGLRLM,2);
                numColumnsToAdd = maxRunLength - tempMaxRunLength;
                newColumns = zeros(numLevels,numColumnsToAdd);
                padding = [tempGLRLM newColumns];
                % merge
                mergeGLRLM = mergeGLRLM + padding;
            end
            totNv = numel(slice2D(:));
            voxelstoRemove = sum(isnan(slice2D(:)));
            Nv = numDir*(totNv - voxelstoRemove);
            fresult = GLRLMStats(mergeGLRLM,Nv);
            biomarkers = fieldnames(fresult);
            cellResults = cell2mat(struct2cell(fresult));
            matrixB(:,k) = cellResults;
        end
        average = sum(matrixB,2)./zDim;
        dim = size(average,1);
        intermediate = num2cell(average,dim);
        featuresB = cell2struct(intermediate,biomarkers,1);
    end

    function featuresC = GLRLMc(winImg,numLevels)
        offset=[3; 1; 2; 4];
        zDim = size(winImg,3);
        numDir = size(offset,1);
        dimens = [size(winImg,1);size(winImg,2);size(winImg,3)];
        maxRunLength = ceil(sqrt(2)*max(dimens));
        Nv = numel(winImg(:,:,:));
        voxelsToSubtract = sum(sum(sum(isnan(winImg(:,:,:)))));
        totNv = Nv - voxelsToSubtract;

        for q=1:numDir
            dir = offset(q);
            directionalGLRLM = zeros(numLevels,maxRunLength);
            for k = 1:zDim
                slice2D = winImg(:,:,k);
                temp=grayrlmatrix(slice2D,'GrayLimits',[min(winImg(:)),...
                    max(winImg(:))],'NumLevels',numLevels,'Offset',dir);
                tempGLRLM = temp{1};

                % change to make it possible to add directly
                % parse through padded GLRLM
                rowSum = sum(tempGLRLM,1); % for runs
                counter = 0;
                % remove excess
                for runIndex = size(rowSum,2):-1:1
                    if (rowSum(1,runIndex) == 0) && counter == 0
                        tempGLRLM(:,runIndex) = [];
                    else
                        counter = 1;
                    end
                end
                % add extra
                tempMaxRunLength = size(tempGLRLM,2);
                numColumnsToAdd = maxRunLength - tempMaxRunLength;
                newColumns = zeros(numLevels,numColumnsToAdd);
                padding = [tempGLRLM newColumns];
                % merge
                directionalGLRLM = directionalGLRLM + padding;
            end
            indivDir = GLRLMStats(directionalGLRLM,totNv);
            biomarkers = fieldnames(indivDir);
            cellResults = cell2mat(struct2cell(indivDir));
            matrixB(:,q) = cellResults;
        end
        average = sum(matrixB,2)./numDir; % average over directions
        dim = size(average,1);
        intermediate = num2cell(average,dim);
        featuresC = cell2struct(intermediate,biomarkers,1);
    end

    function featuresD = GLRLMd(winImg,numLevels)
        offset=[3; 1; 2; 4];
        zDim = size(winImg,3);
        dimens = [size(winImg,1);size(winImg,2);size(winImg,3)];
        maxRunLength = ceil(sqrt(2)*max(dimens));
        numDir = size(offset,1);
        mergeGLRLMB = zeros(numLevels,maxRunLength);
        Nv = numel(winImg);
        voxelsToSubtract = sum(sum(sum(isnan(winImg))));
        totNv = zDim*(Nv - voxelsToSubtract);
        for q=1:numDir
            dir = offset(q);
            mergeGLRLMA = zeros(numLevels,maxRunLength);
            for k = 1:zDim
                slice2D = winImg(:,:,k);
                temp=grayrlmatrix(slice2D,'GrayLimits',[min(winImg(:)),...
                    max(winImg(:))],'NumLevels',numLevels,'Offset',dir);
                tempGLRLM = temp{1};

                % change to make it possible to add directly
                % parse through padded GLRLM
                rowSum = sum(tempGLRLM,1); % for runs
                counter = 0;
                % remove excess
                for runIndex = size(rowSum,2):-1:1
                    if (rowSum(1,runIndex) == 0) && counter == 0
                        tempGLRLM(:,runIndex) = [];
                    else
                        counter = 1;
                    end
                end
                % add extra
                tempMaxRunLength = size(tempGLRLM,2);
                numColumnsToAdd = maxRunLength - tempMaxRunLength;
                newColumns = zeros(numLevels,numColumnsToAdd);
                padding = [tempGLRLM newColumns];
                % merge
                mergeGLRLMA = mergeGLRLMA + padding;
            end
            mergeGLRLMB = mergeGLRLMB + mergeGLRLMA;
        end
        featuresD = GLRLMStats(mergeGLRLMB,totNv);
    end

    function featuresE = GLRLMe(winImg,numLevels)
        FScounter = 1;
        xDim = size(winImg,1);
        yDim = size(winImg,2);
        zDim = size(winImg,3);
        featureStorage = zeros(16,13); % 13 direction matrices, 16 features calculated
        counter = 1; % used as a placeholder, determines which direction, which column to add to not overwrite data

        %% Nondiagonal GLRLM Features
        % vector only changes in max 2 dimensions at a time, not all three
        dimens = [size(winImg,1);size(winImg,2);size(winImg,3)];
        maxRunLength = ceil(sqrt(2)*max(dimens));

        for v = 1:3
            if v == 1
                offsetNonDiag  = [1;2;3;4]; % 0, 45, 90
            elseif v == 2
                offsetNonDiag = [1;2;4]; % 0, 45, 135
            elseif v == 3
                offsetNonDiag = [2;4]; % 45, 135
            end

            shiftedImage = shiftdim(winImg,v); %shift image to slice along all three directions

            for vect = 1:size(offsetNonDiag)
                perSlice = zeros(numLevels,maxRunLength);
                sliceVox = numel(shiftedImage) - sum(sum(sum(isnan(shiftedImage))));
                for h = 1:size(shiftedImage,3)
                    slice(:,:) = shiftedImage(:,:,h); % 2D input for grayrlmatrix function
                    padding = zeros(numLevels,maxRunLength); % allows for easy merging w/ different max run lengths
                    temp = grayrlmatrix(slice,'GrayLimits',[min(winImg(:)),...
                        max(winImg(:))],'NumLevels',numLevels,'Offset',offsetNonDiag(vect));
                    tempGLRLM = temp{1}; % change array type
                    % change to make it possible to add directly
                    % parse through padded GLRLM
                    rowSum = sum(tempGLRLM,1); % for runs

                    % remove excess
                    for runIndex = size(rowSum,2):-1:1
                        if (rowSum(1,runIndex) == 0) && counter == 0
                            tempGLRLM(:,runIndex) = [];
                        else
                            counter = 1;
                        end
                    end
                    % add extra
                    tempMaxRunLength = size(tempGLRLM,2);
                    numColumnsToAdd = maxRunLength - tempMaxRunLength;
                    newColumns = zeros(numLevels,numColumnsToAdd);
                    padding = [tempGLRLM newColumns];
                    % merge
                    perSlice = perSlice + padding;
                end
                features = GLRLMStats(perSlice,sliceVox);
                biomarkers = fieldnames(features);
                featureStorage(:,FScounter) = cell2mat(struct2cell(features));
                FScounter = FScounter + 1;

                clear slice % will be overwritten, additional step
                clear perSlice
            end

            clear offsetNonDiag % additional step for sanity check
        end

        %% Diagonal GLRLM Features

        %blank new images, NaN pad to not affect run length
        niA = nan(2*xDim-1, 2*yDim-1, zDim);
        niB = nan(2*xDim-1, 2*yDim-1, zDim);
        niC = nan(2*xDim-1, 2*yDim-1, zDim);
        niD = nan(2*xDim-1, 2*yDim-1, zDim);

        for k = 1:zDim
            %set different possible dimensions
            sxda = k:(k+xDim-1);
            sxdb = (xDim + 1 - k):(2*xDim - k);
            syda= k:1:(k+yDim-1);
            sydb = (yDim - k + 1):(2*yDim - k);

            % stagger 4 images for 4 3D diag direction matrices
            niA(sxda,syda,k) = winImg(:,:,k); % staggers upper L down
            niB(sxdb,syda,k) = winImg(:,:,k); % staggers bottom L up
            niC(sxda,sydb,k) = winImg(:,:,k); % staggers upper R down
            niD(sxdb,sydb,k) = winImg(:,:,k); % staggers bottom R up
        end

        dim = ceil(sqrt(2)*zDim);

        stagVox = numel(niA) - sum(sum(sum(isnan(niA))));

        % slightly repetitive, calculate and store features for each diagonal GLRLM
        featuresI = GLRLMStats(staggeredImageCount(niA,numLevels,dim),stagVox);
        featureStorage(:,FScounter) = cell2mat(struct2cell(featuresI));
        FScounter = FScounter +1;

        featuresII = GLRLMStats(staggeredImageCount(niB,numLevels,dim),stagVox);
        featureStorage(:,FScounter) = cell2mat(struct2cell(featuresII));
        FScounter = FScounter +1;

        featuresIII = GLRLMStats(staggeredImageCount(niC,numLevels,dim),stagVox);
        featureStorage(:,FScounter) = cell2mat(struct2cell(featuresIII));
        FScounter = FScounter +1;

        featuresIV = GLRLMStats(staggeredImageCount(niD,numLevels,dim),stagVox);
        featureStorage(:,FScounter) = cell2mat(struct2cell(featuresIV));

        average = sum(featureStorage,2)./size(featureStorage,2); % average features across the 13 directions
        dim = size(average,1);
        intermediate = num2cell(average,dim);
        featuresE = cell2struct(intermediate,biomarkers,1); % change array data type
    end

    function featuresF = GLRLMf(winImg,numLevels)
        xDim = size(winImg,1);
        yDim = size(winImg,2);
        zDim = size(winImg,3);

        %% diagonal GLRLM
        GLRLMdiag = diag3D(winImg,numLevels); % for the four 3D direction vectors

        %% Nondiagonal GLRLM
        % vector only changes in max 2 dimensions at a time, not all three
        GLRLMnonDiag = zeros(numLevels,ceil(sqrt(2)*xDim)); % initial padding, note there's an excess of columns. accounted for in the GLRLM stats function
        for v = 1:3
            if v == 1
                offsetNonDiag  = [1;2;3;4]; % 0, 45, 90
            elseif v == 2
                offsetNonDiag = [1;2;4]; % 0, 45, 135
            elseif v == 3
                offsetNonDiag = [2;4]; % 45, 135
            end

            shiftedImage = shiftdim(winImg,v); %shift image to slice along all three directions
            for h = 1:size(shiftedImage,3)
                slice(:,:) = shiftedImage(:,:,h); % 2D input for grayrlmatrix function
                for vect = 1:size(offsetNonDiag)
                    padding = zeros(numLevels,ceil(sqrt(2)*xDim)); % allows for easy merging w/ different max run lengths
                    temp = grayrlmatrix(slice,'GrayLimits',[min(winImg(:)),...
                        max(winImg(:))],'NumLevels',numLevels,'Offset',offsetNonDiag(vect));
                    padding(1:numel(temp{1})) = temp{1}; % change array type
                    GLRLMnonDiag = GLRLMnonDiag + padding; % add individual slice's GLRLM to total sum
                end
                clear slice % will be overwritten, additional step
            end
            clear offsetNonDiag % additional step for sanity check
        end

        rowSum = sum(GLRLMnonDiag,1); % for runs
        counter = 0;
        % remove excess
        for runIndex = size(rowSum,2):-1:1
            if (rowSum(1,runIndex) == 0) && counter == 0
                GLRLMnonDiag(:,runIndex) = [];
            else
                counter = 1;
            end
        end


        rowSum = sum(GLRLMdiag,1); % for runs
        counter = 0;
        % remove excess
        for runIndex = size(rowSum,2):-1:1
            if (rowSum(1,runIndex) == 0) && counter == 0
                GLRLMdiag(:,runIndex) = [];
            else
                counter = 1;
            end
        end

        difference = size(GLRLMnonDiag,2) - size(GLRLMdiag,2);
        columnsToAdd = zeros(numLevels,difference);
        if difference > 0
            GLRLMdiag = [GLRLMdiag columnsToAdd];
        elseif difference < 0
            GLRLMnonDiag = [GLRLMnonDiag columnnsToAdd];

        end

        total = GLRLMdiag + GLRLMnonDiag; % merge the 13 directions' texture matrices

        Nv = numel(winImg); % take total number of voxels
        voxelsToSubtract = sum(sum(sum(isnan(winImg)))); % sum NaN values
        totNv = Nv - voxelsToSubtract; % subtract NaN to get number of voxels in ROI

        featuresF = GLRLMStats(total,totNv); % calculate image features from merged texture matrix
    end

% additional functions for featuresE and featuresF functions

    function GLRLM = diag3D(winImg,numLevels)
        xDim = size(winImg,1);
        yDim = size(winImg,2);
        zDim = size(winImg,3);

        %blank new images, NaN pad to not affect run length
        niA = nan(2*xDim-1, 2*yDim-1, zDim);
        niB = nan(2*xDim-1, 2*yDim-1, zDim);
        niC = nan(2*xDim-1, 2*yDim-1, zDim);
        niD = nan(2*xDim-1, 2*yDim-1, zDim);

        for k = 1:zDim
            %set different possible dimensions
            sxda = k:(k+xDim-1);
            sxdb = (xDim + 1 - k):(2*xDim - k);
            syda= k:1:(k+yDim-1);
            sydb = (yDim - k + 1):(2*yDim - k);

            % stagger 4 images for 4 3D diag direction matrices
            niA(sxda,syda,k) = winImg(:,:,k); % staggers upper L down
            niB(sxdb,syda,k) = winImg(:,:,k); % staggers bottom L up
            niC(sxda,sydb,k) = winImg(:,:,k); % staggers upper R down
            niD(sxdb,sydb,k) = winImg(:,:,k); % staggers bottom R up
        end
        dimensions = [xDim yDim zDim];
        maxDim = max(dimensions);
        dim = ceil(sqrt(2)*maxDim);
        GLRLMA = staggeredImageCount(niA,numLevels,dim);
        GLRLMB = staggeredImageCount(niB,numLevels,dim);
        GLRLMC = staggeredImageCount(niC,numLevels,dim);
        GLRLMD = staggeredImageCount(niD,numLevels,dim);

        GLRLM = GLRLMA + GLRLMB + GLRLMC + GLRLMD;
    end

    function indivGLRLM = staggeredImageCount(newImage,numLevels,dim)
        indivGLRLM = zeros(numLevels,dim);
        
        for j = 1:size(newImage,2)
            for i = 1:size(newImage,1)
                column = zeros((size(newImage,3)),1); % column will stretch along z
                column(:) = newImage(i,j,:);
                indivGLRLM = indivGLRLM + CRL(column,numLevels,dim); % add each column GLRLM
            end
        end

    end

    function temp = CRL(c,numLevels,dim) % column run length, feed in column of pixels in correct run direction
        temp = zeros(numLevels,dim);
        oldGL = NaN;
        runlength = 1;
        for t = 1:length(c) % iterate thru z column
            GL = c(t); % assign new gray level
            if isnan(GL)
                oldGL = NaN; % resets
                runlength = 1;
            elseif (oldGL == GL)&& (isnan(GL) == false) % if run is incomplete
                temp(GL,runlength) =  temp(GL,runlength) - 1; % avoids double counting
                runlength = runlength + 1;
                temp(GL,runlength) = temp(GL,runlength) + 1;
            else % if run is complete && GL is not NaN
                runlength = 1;
                temp(GL,runlength) = temp(GL,runlength) + 1; % add the new GL to the matrix
                oldGL = GL;
            end
        end
    end



%% GLRLM Stats, Image Feature Calculation from Texture Matrices
    function A = GLRLMStats(GLRLM,Nv)

        Ns = sum(GLRLM(:)); % sum over all elements
        Ng = size(GLRLM,1); % number of gray levels

        % parse through padded GLRLM
        rowSum = sum(GLRLM,1); % for runs
        Nr = size(rowSum,2);
        counter = 0;
        for runIndex = size(rowSum,2):-1:1
            if (rowSum(1,runIndex) == 0) && counter == 0
                Nr = Nr - 1;
            else
                counter = 1;
            end
        end

        %Calculate num pixels based on run lengths
        runLength = [1:size(GLRLM,2)];
        Nv = runLength*rowSum';

        rj = sum(GLRLM,1);
        ri = sum(GLRLM,2);
        p = GLRLM./Ns;

        mui = 0;
        for i = 1:Ng
            for j = 1:Nr
                mui = mui + (i*p(i,j));
            end
        end

        muj = 0;
        for i = 1:Ng
            for j = 1:Nr
                muj = muj + (j*p(i,j));
            end
        end

        A = struct();

        % term = []; good for multiple data types
        term = 0;
        for j = 1:Nr
            for i = 1:Ng
                term = term + (GLRLM(i,j)./(j.^2));
            end
        end
        A.SRE = (1/Ns)*term;

        term = 0;
        for j = 1:Nr
            for i = 1:Ng
                term = term + (GLRLM(i,j)*(j.^2));
            end
        end
        A.LRE = (1/Ns)*term;

        term = 0;
        for j = 1:Nr
            for i = 1:Ng
                term = term + (GLRLM(i,j)./(i.^2));
            end
        end
        A.LGLRE = (1/Ns)*term;

        term = 0;
        for j = 1:Nr
            for i = 1:Ng
                term = term + (GLRLM(i,j)*(i.^2));
            end
        end
        A.HGLRE = (1/Ns)*term;

        term = 0;
        for i = 1:Ng
            for j = 1:Nr
                term = term + (GLRLM(i,j)./(i^2*j^2));
            end
        end
        A.SRLGLE = (1/Ns)* term;

        term = 0;
        for i = 1:Ng
            for j = 1:Nr
                term = term + ((GLRLM(i,j)*i^2)/(j^2));
            end
        end
        A.SRHGLE = (1/Ns)* term;

        term = 0;
        for i = 1:Ng
            for j = 1:Nr
                term = term + ((GLRLM(i,j)*j^2)/(i^2));
            end
        end
        A.LRLGLE = (1/Ns)* term;


        term = 0;
        for i = 1:Ng
            for j = 1:Nr
                term = term + (GLRLM(i,j)*j^2*i^2);
            end
        end
        A.LRHGLE = (1/Ns)* term;

        term = 0;
        for i = 1:Ng
            term = term + ri(i).^2;
        end
        A.GLN = (1/Ns)*term;

        term = 0;
        for i = 1:Ng
            term = term + ri(i).^2;
        end
        A.GLNN = (1/(Ns^2))*term;

        term = 0;
        for j = 1:Nr
            term = term + rj(j).^2;
        end
        A.RLN = (1/Ns)*term;

        term = 0;
        for j = 1:Nr
            term = term + rj(j).^2;
        end
        A.RLNN = (1/(Ns^2))*term;

        A.RP = Ns/Nv;

        term = 0;
        for i =1:Ng
            for j = 1:Nr
                term = term + ((i - mui).^2*p(i,j));
            end
        end
        A.GLV = term;

        term = 0;
        for i =1:Ng
            for j = 1:Nr
                term = term + ((j - muj).^2*p(i,j));
            end
        end
        A.RLV = term;

        term = 0;
        for i =1:Ng
            for j = 1:Nr
                pixel = p(i,j);
                if pixel ~= 0
                    term = term + (p(i,j)*(-log2(p(i,j))));
                else
                    continue;
                end
            end
        end
        A.ENT = term;

    end

%% USC Radiomics Aggregating Methods

    function fRLMatrix=GLRLM_auto2D(winImg,numLevels)
        % USC in house aggregating method, used when ignoring corner voxels
        % offset=[-1 0;0 1; -1 1; -1 -1];
        offset=[3; 1; 2; 4];
        fRLMatrix=zeros([numLevels,numLevels]);

        for q=1:size(offset,1)
            temp=grayrlmatrix(winImg,'GrayLimits',[min(winImg(:)),...
                max(winImg(:))],'NumLevels',numLevels,'Offset',offset(q,:));
            fRLMatrix=MatAdd(fRLMatrix,temp{1});
        end
        % fRLMatrix=fRLMatrix/size(offset,1);
    end

    function fRLMatrix=GLRLM_auto3D(winImg,NumLevels)
        % USC in-house method, ignores corner voxels
        %offset=[-1 0;0 1; -1 1; -1 -1];
        offset=[3; 1; 2; 4];
        fRLMatrix=zeros([NumLevels,NumLevels]);
        for z=1:3
            RotwinImg=shiftdim(winImg,z);
            fMatrix=zeros([NumLevels,NumLevels]);
            for q=z:size(offset,1)
                fresult=zeros([NumLevels,NumLevels]);
                for k=1:size(RotwinImg,3)
                    tempImg(:,:)=RotwinImg(:,:,k);
                    if sum(~isnan(tempImg(:)))~=0
                        temp=grayrlmatrix(tempImg,'GrayLimits',[min(winImg(:)),...
                            max(winImg(:))],'NumLevels',NumLevels,'Offset',offset(q,:));
                        fresult=MatAdd(fresult,temp{1});
                        clear tempImg;

                    end
                end
                fMatrix=MatAdd(fresult,fMatrix);
            end
            fRLMatrix=MatAdd(fRLMatrix,fMatrix);
        end
        % fRLMatrix=fRLMatrix/9;
    end

    function fresult=MatAdd(A,B)
        % matrix addition
        s1=size(A);
        s2=size(B);
        s3=[s1;s2];
        s3=max(s3);
        tempA=zeros(s3);
        tempA(1:size(A,1),1:size(A,2))=A;
        tempB=zeros(s3);
        tempB(1:size(B,1),1:size(B,2))=B;
        fresult=tempA+tempB;
    end

end
