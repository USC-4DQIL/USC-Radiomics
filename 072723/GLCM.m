%% Gray Level Cooccurrence Matrix (GLCM)
function GLCMFeatures = GLCM(winImg,numLevels,binWidth)

numLevels = max(max(max(winImg)));
%% Compute Features
featuresA = GLCMa(winImg,numLevels); % 2D, by slice, w/o merging
featuresB = GLCMb(winImg,numLevels); % 2D, by slice, w/ merging by slice
featuresC = GLCMc(winImg,numLevels); % 2.5D, by slice, w/ merging by direction
featuresD = GLCMd(winImg,numLevels); % 2.5D, by slice, w/ full merging
featuresE = GLCMe(winImg,numLevels); % 3D, as volume, w/o merging
featuresF = GLCMf(winImg,numLevels); % 3D, as volume, w/ full merging

featuresAF = [featuresA,featuresB,featuresC,featuresD,featuresE,featuresF];
agMethods = {'A','B','C','D','E','F'};
conCat = [agMethods' struct2table(featuresAF)]; % labeling step
GLCMFeatures = renamevars(conCat,'Var1','Aggregating_Methods');

% GLCM3D = GLCM_auto3D(winImg,numLevels); % USC Radiomics aggregating
% method, ignores corner voxels
% featuresGLCM3D = GLCMStats(GLCM3D); % compute statistics
% GLCM2D = GLCM_auto2D(winImg(:,:,1),numLevels); % USC radiomics GLCM method, only used for 2D inputs
% featuresGLCM2D = GLCMStats(GLCM2D);

%% Secondary Statistics Generator
    function fresult=GLCMStats(GLCM)
        normGLCM=GLCM./sum(GLCM(:));
        Pxplusy=GLCM3DPxplusy(normGLCM);
        Pxminusy=GLCM3DPxminusy(normGLCM);
        %     Contrast (CON)
        fresult.CON=GLCM3DContrast(normGLCM,'CON');
        %     Disssimilarity (DIS)
        fresult.DIS=GLCM3DContrast(normGLCM,'DIS');
        %     Homogeneity (HOM)
        fresult.HOM=GLCM3DContrast(normGLCM,'HOM');
        %     Angular Second Moment (ASM)
        fresult.ASM=GLCM3DASM(normGLCM);
        %     GLCM Mean
        fresult.Mean=GLCM3DMean(normGLCM);
        %     Texture Variance (SQV)
        fresult.SQV=GLCM3DSQV(normGLCM);
        %     GLCM Standard Deviation
        fresult.STD=sqrt(fresult.SQV);
        %     Energy (also called Uniformity)
        fresult.Uniformity = sqrt(fresult.ASM);
        % Joint Maximum
        fresult.jointMax = GLCM3DjointMax(normGLCM);
        % Joint Average
        fresult.jointAvg = GLCM3DjointAvg(normGLCM);
        %     GLCM Variance
        fresult.jointVar=GLCM3DVar(normGLCM,fresult.Mean);
        % Joint Entropy
        fresult.jointENT = GLCM3DjointEnt(normGLCM);
        % Difference Average
        fresult.difAve=GLCM3DdifAve(Pxminusy); %% Edit this function
        % Difference variance
        fresult.difVar = GLCM3DdifVar(Pxminusy,fresult.DIS);
        %     Difference entropy (DIF_ENT)
        fresult.difENT=GLCM3DsumENT(Pxminusy);
        % Cluster Tendency
        fresult.clusterTend = clusterTend(normGLCM);
        % Cluster Shade
        fresult.clusterShade = clusterShade(normGLCM);
        % Cluster Prominence
        fresult.clusterProminence = clusterProminence(normGLCM);
        %    image average gray texture within each region (SUM_AVE)
        fresult.sumAve=GLCM3DsumAve(fresult.jointAvg);
        %     Sum variance (SUM_VAR)
        fresult.sumVar=GLCM3DsumVar(fresult.clusterTend);
        %     Sum entropy (SUM_ENT)
        fresult.sumENT=GLCM3DsumENT(Pxplusy);
        %     Entropy (ENT)
        fresult.ENT=GLCM3DEntropy(normGLCM);
        %     GLCM Correlation
        fresult.Corr=GLCM3DCorrelation(normGLCM,fresult.Mean,...
            fresult.jointVar);
        % Inverse Difference
        fresult.invDif = GLCM3DinvDif(normGLCM);
        % Normalized Inverse Difference
        fresult.normInvDif = GLCM3DnormInvDif(normGLCM);
        %     inverse difference moment(IDM)
        fresult.IDM=GLCM3DIDM(normGLCM);
        % Normalized IDM
        fresult.normIDM = GLCM3DnormIDM(normGLCM);
        % Inverse Variance
        fresult.invVar = GLCM3DinvVar(normGLCM);
        % Autocorrelation
        fresult.autocorr = autocorrelation(normGLCM);
        %     Information Measure of Correlation 1 (IMCI) indicates the direction
        %     measurement of image texture energy distribution
        fresult.IMC1=GLCM3DIMC1(normGLCM,fresult.ENT);
        %     Information Measure of Correlation 2 (IMC2) indicates that the energy
        %     intensity distribution measurement of the image texture
        fresult.IMC2=GLCM3DIMC2(normGLCM,fresult.jointENT);
        %     Maximum Correlation Coefficient (MCC)
        fresult.MCC=GLCM3DMCC(normGLCM);

    end

%% Diagonal and Cross-Diagonal Calculations, Mean & ENT
    function fresult=GLCM3DPxplusy(normGLCM)
        fresult=zeros(2*size(normGLCM,1),1);
        for k=1:size(normGLCM,1)
            for q=1:size(normGLCM,2)
                fresult(k+q-1,1)=fresult(k+q-1,1)+normGLCM(k,q);
            end
        end
    end

    function fresult=GLCM3DPxminusy(normGLCM)
        fresult=zeros(2*size(normGLCM,1),1);
        for k=1:size(normGLCM,1)
            for q=1:size(normGLCM,2)
                fresult(abs(k-q)+1,1)=fresult(abs(k-q)+1,1)+normGLCM(k,q);
            end
        end
    end

    function fresult=GLCM3DEntropy(normGLCM)
        temp=normGLCM;
        temp(temp~=0)=-log(temp(temp~=0));
        temp=temp.*normGLCM;
        fresult=sum(temp(:));
    end

    function fresult=GLCM3DMean(normGLCM)
        matsize=size(normGLCM,1);
        conweight=repmat((0:matsize-1)',1,matsize);
        temp=conweight.*normGLCM;
        fresult=sum(temp(:));
    end

%% Joint Measures
    function fresult = GLCM3DjointMax(normGLCM)
        fresult = max(max(normGLCM));
    end

    function fresult = GLCM3DjointAvg(normGLCM)
        matsize=size(normGLCM,1);
        conweight=repmat((1:matsize)',1,matsize);
        temp=conweight.*normGLCM;
        fresult=sum(temp(:));
    end

    function fresult = GLCM3DVar(normGLCM,GLCMMean)
        matsize=size(normGLCM,1);
        conweight=repmat((0:matsize-1)',1,matsize);
        temp=(conweight-GLCMMean).^2.*normGLCM;
        fresult=sum(temp(:));
    end

    function fresult = GLCM3DjointEnt(normGLCM)
        temp=normGLCM;
        temp(temp~=0)=-log2(temp(temp~=0));
        fresult = sum(normGLCM(:).*temp(:));
    end

%% Difference Measures
    function fresult=GLCM3DdifAve(Pxminusy)
        fresult=sum((0:(size(Pxminusy,1)-1))'.*Pxminusy);
    end

    function fresult = GLCM3DdifVar(Pxminusy,dis)
        fresult = sum(((0:(size(Pxminusy,1)-1))'-dis).^2.*Pxminusy);
    end

%% Sum Measures
    function fresult=GLCM3DsumAve(jointAvg)
        fresult = 2.*jointAvg;
    end

    function fresult=GLCM3DsumVar(clusterTend)
        fresult=clusterTend;
    end

    function fresult=GLCM3DsumENT(Pxplusy)
        temp=Pxplusy;
        temp(temp~=0)=-log2(temp(temp~=0));
        fresult=sum(temp.*Pxplusy);
    end

%% Contrast, ASM, GLCM-Specific Measures
    function fresult=GLCM3DContrast(normGLCM,conType)
        matsize=size(normGLCM,1);
        conweight=zeros([matsize,matsize]);
        for k=0:matsize-1
            v=ones(matsize-k,1)*GLCM3DContrastWeight(k,conType);
            conweight=conweight+diag(v,k);
        end
        conweight=conweight+triu(conweight,1)';
        temp=conweight.*normGLCM;
        fresult=sum(temp(:));
    end

    function fresult=GLCM3DContrastWeight(k,conType)
        switch conType
            case 'CON'
                fresult=k^2;
            case 'DIS'
                fresult=k;
            case 'HOM'
                fresult=1/(1+k^2);
        end
    end

    function fresult=GLCM3DASM(normGLCM)
        temp=normGLCM.^2;
        fresult=sum(temp(:));
    end

    function fresult=GLCM3DSQV(normGLCM)
        GLCMMean=sum(normGLCM(normGLCM>0));
        matsize=size(normGLCM,1);
        conweight=repmat((0:matsize-1)',1,matsize);
        temp=(conweight-GLCMMean).^2.*normGLCM;
        fresult=sum(temp(:));
    end

    function fresult=GLCM3DCorrelation(normGLCM,GLCMMean,GLCMVariance)
        matsize=size(normGLCM,1);
        conweighti=repmat((0:matsize-1)',1,matsize);
        conweightj=conweighti';
        temp=((conweighti-GLCMMean).*(conweightj-GLCMMean))./GLCMVariance.*...
            normGLCM;
        fresult=sum(temp(:));
    end

%% Correlation & Inverse Measures
    function fresult = GLCM3DinvDif(normGLCM)
        matsize = size(normGLCM,1);
        iarray = repmat((1:matsize),matsize,1);
        jarray = repmat((1:matsize),matsize,1)';
        weight = abs(iarray-jarray) + 1;
        fresult = sum(sum(normGLCM./weight));
    end

    function fresult = GLCM3DnormInvDif(normGLCM)
        matsize = size(normGLCM,1);
        iarray = repmat((1:matsize),matsize,1);
        jarray = repmat((1:matsize),matsize,1)';
        weight = (abs(iarray-jarray))./matsize + 1;
        fresult = sum(sum(normGLCM./weight));
    end

    function fresult=GLCM3DIDM(normGLCM)
        matsize = size(normGLCM,1);
        iarray = repmat((1:matsize),matsize,1);
        jarray = repmat((1:matsize),matsize,1)';
        weight = ((iarray-jarray).^2) + 1;
        fresult = sum(sum(normGLCM./weight));
    end

    function fresult = GLCM3DnormIDM(normGLCM)
        matsize = size(normGLCM,1);
        iarray = repmat((1:matsize),matsize,1);
        jarray = repmat((1:matsize),matsize,1)';
        weight = (((iarray-jarray).^2)./matsize^2) + 1;
        fresult = sum(sum(normGLCM./weight));
    end

    function fresult = GLCM3DinvVar(normGLCM)
        matsize = size(normGLCM,1);
        iarray = repmat((1:matsize),matsize,1);
        jarray = iarray';
        weight = (iarray-jarray).^2;
        x = normGLCM./(weight);
        x(isnan(x)|isinf(x)) = 0;
        fresult = sum(sum(x));
    end

    function fresult = autocorrelation(normGLCM)
        matsize= size(normGLCM,1);
        iarray = repmat((1:matsize),matsize,1);
        jarray = iarray';
        fresult = sum(sum(iarray.*jarray.*normGLCM));
    end

    function fresult=GLCM3DIMC1(normGLCM,ENT)
        Px=sum(normGLCM,2);
        Py=sum(normGLCM,1);
        temp=Px;
        temp(temp~=0)=-log(temp(temp~=0));
        ENTx=sum(temp.*Px);
        temp=Py;
        temp(temp~=0)=-log(temp(temp~=0));
        ENTy=sum(temp.*Py);
        temp=Px*Py;
        temp2=zeros(size(temp));
        temp2(temp~=0)=-log(temp(temp~=0));
        HXY1=sum(normGLCM(:).*temp2(:));
        HXY2=sum(temp(:).*temp2(:));
        fresult=(ENT-HXY1)/max(ENTx,ENTy);
    end

    function fresult=GLCM3DIMC2(normGLCM,ENT)
        Px=sum(normGLCM,2);
        Py=sum(normGLCM,1);
        temp=Px*Py;
        temp2=zeros(size(temp));
        temp2(temp~=0)=-log2(temp(temp~=0));
        HXY2=sum(temp(:).*temp2(:));
        fresult=sqrt(1-exp(-2.0*(HXY2-ENT)));
    end

    function fresult=GLCM3DMCC(normGLCM)
        Px=sum(normGLCM,2);
        Py=sum(normGLCM,1);
        for k=1:size(normGLCM,1)
            for q=1:size(normGLCM,2)
                if Px(k)==0 || Py(q)==0
                    Q(k,q)=0;
                else
                    Q(k,q)=sum(normGLCM(k,:).*normGLCM(q,:)/(Px(k)*Py(q)));
                end
            end
        end
        lambda=eig(Q);
        lambda=sort(lambda,'descend');
        fresult=sqrt(lambda(2));
    end

%% Cluster Measures
    function fresult = clusterTend(normGLCM)
        nGrayLevels = size(normGLCM,1);
        nglcm = 1;
        sub   = 1:nGrayLevels*nGrayLevels;
        [I,J] = ind2sub([nGrayLevels,nGrayLevels],sub);
        uX = zeros(nglcm,1);
        uY = zeros(nglcm,1);
        currentGLCM = normGLCM;

        for k = 1:nglcm
            uX(k)   = sum(I.*currentGLCM(sub));
            uY(k)   = sum(J.*currentGLCM(sub));
            clusterTend = sum((I+J-uX(k)-uY(k)).^2.*currentGLCM(sub)); %OK
        end
        fresult = clusterTend;
    end

    function fresult = clusterShade(normGLCM)
        nGrayLevels = size(normGLCM,1);
        nglcm = 1;
        sub   = 1:nGrayLevels*nGrayLevels;
        [I,J] = ind2sub([nGrayLevels,nGrayLevels],sub);
        uX = zeros(nglcm,1);
        uY = zeros(nglcm,1);
        currentGLCM = normGLCM;

        for k = 1:nglcm
            uX(k)   = sum(I.*currentGLCM(sub));
            uY(k)   = sum(J.*currentGLCM(sub));
            clusterShade = sum((I+J-uX(k)-uY(k)).^3.*currentGLCM(sub)); %OK
        end

        fresult = clusterShade;
    end

    function fresult = clusterProminence(normGLCM)
        matsize=size(normGLCM,1);
        nGrayLevels = size(normGLCM,1);
        nglcm = 1;
        sub   = 1:nGrayLevels*nGrayLevels;
        [I,J] = ind2sub([nGrayLevels,nGrayLevels],sub);
        uX = zeros(nglcm,1);
        uY = zeros(nglcm,1);
        currentGLCM = normGLCM;
        for k = 1:nglcm
            uX(k)   = sum(I.*currentGLCM(sub));
            uY(k)   = sum(J.*currentGLCM(sub));
            clusterProminence = sum((I+J-uX(k)-uY(k)).^4.*currentGLCM(sub)); %OK
        end
        fresult = clusterProminence;
    end

%% IBSI Based Joint Calculation of GLCM and Features
% separated by aggregating method

% a -- 2D, by slice, w/o merging
% b -- 2D, by slice, w/ merging by slice
% c -- 2.5D, by slice, w/ merging by direction
% d -- 2.5D, by slice, w/ full merging
% e -- 3D, as volume, w/o merging
% f -- 3D, as volume, w/ full merging

    function featuresA = GLCMa(winImg,numLevels)
        offset = [1 0; 1 1; 0 1; -1 1]; % directions 1-4
        numDir = size(offset,1);
        stackNum = size(winImg,3);
        for k = 1:stackNum
            slice2D = winImg(:,:,k);
            for q = 1:numDir % compute direction matrices from 1-4
                tempGLCM = graycomatrix(slice2D,'GrayLimits',[min(winImg(:)),...
                    max(winImg(:))],'NumLevels',numLevels,'Offset',offset(q,:),...
                    'Symmetric',true);
                tempNormGLCM = tempGLCM/sum(tempGLCM(:));
                fresults = GLCMStats(tempNormGLCM);
                biomarkers = fieldnames(fresults);
                cellResults = cell2mat(struct2cell(fresults));
                matrixA(:,q) = cellResults; % stores directional matrices, column 1 == direction 1, etc.
            end
            averageofSlice = sum(matrixA,2)./numDir; % divide by number of directions
            matrixB(:,k) = averageofSlice;  % average GLCM for each direction computed
        end
        average = sum(matrixB,2)./stackNum;
        dim = size(average,1);
        intermediate = num2cell(average,dim);
        featuresA = cell2struct(intermediate,biomarkers,1);  % put back into structure format
    end

    function featuresB = GLCMb(winImg,numLevels)
        offset = [1 0; 1 1; 0 1; -1 1]; % directions 1-4
        numDir = size(offset,1);
        stackNum = size(winImg,3);
        for k = 1:stackNum
            matrixA = zeros(numLevels,numLevels);
            slice2D = winImg(:,:,k);
            for q = 1:numDir % % compute direction matrices from 1-4
                tempGLCM = graycomatrix(slice2D,'GrayLimits',[min(winImg(:)),...
                    max(winImg(:))],'NumLevels',numLevels,'Offset',offset(q,:),...
                    'Symmetric',true);
                matrixA = matrixA + tempGLCM; % merge the 4 directions
            end
            normMatrixA = matrixA/sum(matrixA(:));
            fresults = GLCMStats(normMatrixA);
            biomarkers = fieldnames(fresults);
            sliceResults = cell2mat(struct2cell(fresults));
            matrixB(:,k) = sliceResults;
        end

        averageOverSlices = sum(matrixB,2)./stackNum; % average over slices
        dim = size(averageOverSlices,1);
        intermediate = num2cell(averageOverSlices,dim);
        featuresB = cell2struct(intermediate,biomarkers,1);
    end

    function featuresC = GLCMc(winImg,numLevels)
        offset = [1 0; 1 1; 0 1; -1 1]; % directions 1-4
        numDir = size(offset,1);
        stackNum = size(winImg,3);
        for q = 1:numDir
            matrixA = zeros(numLevels,numLevels);
            for k = 1:stackNum
                slice2D = winImg(:,:,k);
                tempGLCM = graycomatrix(slice2D,'GrayLimits',[min(winImg(:)),...
                    max(winImg(:))],'NumLevels',numLevels,'Offset',offset(q,:),...
                    'Symmetric',true);
                matrixA = matrixA + tempGLCM;
            end
            normMatrixA = matrixA./sum(matrixA(:));
            fresults = GLCMStats(normMatrixA);
            biomarkers = fieldnames(fresults);
            dirResults = cell2mat(struct2cell(fresults));
            matrixB(:,q) = dirResults;
        end

        averageOverDirections = sum(matrixB,2)./numDir;
        dim = size(averageOverDirections,1);
        intermediate = num2cell(averageOverDirections,dim);
        featuresC = cell2struct(intermediate,biomarkers,1);
    end

    function featuresD = GLCMd(winImg,numLevels)
        offset = [1 0; 1 1; 0 1; -1 1]; % directions 1-4
        numDir = size(offset,1);
        stackNum = size(winImg,3);
        matrixB = zeros(numLevels, numLevels);
        for q = 1:numDir
            matrixA = zeros(numLevels,numLevels);
            for k = 1:stackNum
                slice2D = winImg(:,:,k);
                tempGLCM = graycomatrix(slice2D,'GrayLimits',[min(winImg(:)),...
                    max(winImg(:))],'NumLevels',numLevels,'Offset',offset(q,:),...
                    'Symmetric',true);
                matrixA = matrixA + tempGLCM;
            end
            matrixB = matrixB + matrixA;
        end

        normGLCM = matrixB./sum(matrixB(:));
        featuresD = GLCMStats(normGLCM);
    end

    function featuresE = GLCMe(winImg,numLevels)
        xDim = size(winImg,1);
        yDim = size(winImg,2);
        zDim = size(winImg,3);
        temp = zeros(32,13); % first dim == number of biomarkers, second dim == number of direction vectors
        % 13 3D direction vectors
        offset3D = [0 0 1; 0 1 0; 1 0 0; 0 1 1; 0 1 -1; 1 0 1; 1 0 -1; 1 1 0; 1 -1 0; 1 1 1; 1 1 -1; 1 -1 1; 1 -1 -1];
        for dir = 1:size(offset3D,1)
            tempGLCM = zeros(numLevels, numLevels);
            for i = 1:xDim
                for j= 1:yDim
                    for k = 1:zDim
                        
                        i_max = xDim;
                        j_max = yDim;
                        k_max = zDim;

                        referencePixel = winImg(i,j,k);
                        vect = offset3D(dir,:,:);
                        xShift = vect(1);
                        yShift = vect(2);
                        zShift = vect(3);
                        if (i+xShift < 1) || (i+xShift>i_max) || (j+yShift < 1) || (j+yShift>j_max) || (k+zShift<1) ||  (k+zShift > k_max)
                            continue;
                        else
                            neighborPixel = winImg(i+xShift,j+yShift,k+zShift);

                            for x = 1:numLevels
                                for y = 1:numLevels
                                    if referencePixel == x && neighborPixel == y
                                        tempGLCM(x,y) = tempGLCM(x,y) + 1;
                                    end
                                end
                            end
                        end
                    end
                end
            end
            tempNormGLCM = tempGLCM./sum(tempGLCM(:));
            biomarkers = fieldnames(GLCMStats(tempGLCM));
            temp(:,dir) = cell2mat(struct2cell(GLCMStats(tempNormGLCM)));
        end
        intermediate = mean(temp,2); % average across directions
        featuresE = cell2struct(num2cell(intermediate),biomarkers,1);
    end

    function featuresF = GLCMf(winImg,numLevels)
        xDim = size(winImg,1);
        yDim = size(winImg,2);
        zDim = size(winImg,3);
        M = zeros(numLevels,numLevels);
        offset3D = [0,0,1; 0,1,0; 1,0,0; 0,1,1; 0,1,-1; 1,0,1; 1,0,-1; ...
            1,1,0; 1,-1,0; 1,1,1; 1,1,-1; 1, -1, 1; 1, -1, -1];

        for dir = 1:size(offset3D,1)
            tempGLCM = zeros(numLevels, numLevels);
            vect = offset3D(dir,:,:);
            for i = 1:xDim
                for j= 1:yDim
                    for k = 1:zDim
                        referencePixel = winImg(i,j,k);
                        xShift = vect(1);
                        yShift = vect(2);
                        zShift = vect(3);
                        if (i+xShift < 1) || (i+xShift>xDim) || (j+yShift < 1) || (j+yShift>yDim) || (k+zShift<1) ||  (k+zShift > zDim)
                            continue;
                        else
                            neighborPixel = winImg(i+xShift,j+yShift,k+zShift);

                            for x = 1:numLevels
                                for y = 1:numLevels
                                    if referencePixel == x && neighborPixel == y
                                        tempGLCM(x,y) = tempGLCM(x,y) + 1;
                                    end
                                end
                            end
                        end
                    end
                end
            end
            M = M + tempGLCM; % merge all direction matrices
        end
        normM = M./sum(M(:)); % normalize
        featuresF = GLCMStats(normM);
    end

end


%% USCRadiomics Calculated GLCM
function fCoMatrix=GLCM_auto3D(winImg,NumLevels)
offset=[-1 0;0 1; -1 1; -1 -1];
fCoMatrix=zeros([NumLevels,NumLevels]);
for z=1:3
    RotwinImg=shiftdim(winImg,z);
    fMatrix=zeros([NumLevels,NumLevels]);
    for q=z:size(offset,1)
        fresult=zeros([NumLevels,NumLevels]);
        for k=1:size(RotwinImg,3)
            tempImg(:,:)=RotwinImg(:,:,k);
            fresult=fresult+graycomatrix(tempImg,'GrayLimits',[min(winImg(:)),...
                max(winImg(:))],'NumLevels',NumLevels,'Offset',offset(q,:),...
                'Symmetric',true);
            clear tempImg;
        end
        fMatrix=fresult+fMatrix;
    end
    fCoMatrix=fCoMatrix+fMatrix;
end
fCoMatrix=fCoMatrix/9;
end

function fCoMatrix=GLCM_auto2D(winImg,NumLevels)
offset=[-1 0;0 1; -1 1; -1 -1];
fCoMatrix=zeros([NumLevels,NumLevels]);

for q=1:size(offset,1)
    fCoMatrix=fCoMatrix+graycomatrix(winImg,'GrayLimits',[min(winImg(:)),...
        max(winImg(:))],'NumLevels',NumLevels,'Offset',offset(q,:),...
        'Symmetric',true);
end
fCoMatrix=fCoMatrix/size(offset,1);
end
