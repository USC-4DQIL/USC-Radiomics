%% GLCM 2

function GLCMFeatures = GLCM2(winImg,numLevels,binWidth)

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

end
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

function featuresE = GLCMe(img, n_levels)
level_stack = expand_to_level(img);
temp = unique(img);
vect1 = temp(~isnan(temp));
arrayOI = table2array(combinations(vect1,vect1));
featuresSlice = zeros(32,size(arrayOI,1));
offsetArray = [0,0,1; 0,1,0; 1,0,0; 0,1,1; 0,1,-1; 1,0,1; 1,0,-1; ...
    1,1,0; 1,-1,0; 1,1,1; 1,1,-1; 1, -1, 1; 1, -1, -1];
for j = 1:size(offsetArray,1)
    glcm = zeros(n_levels, n_levels);
    for i = 1:size(arrayOI,1)
        two_level_img = level_stack(:,:,:,[arrayOI(i,1),arrayOI(i,2)]);

        filt = filter1(offsetArray(j,1), offsetArray(j,2),offsetArray(j,3));
        conv_out = convn(two_level_img, filt) == 2;
        glcm(arrayOI(i,1), arrayOI(i,2)) = glcm(arrayOI(i,1), arrayOI(i,2)) + sum(conv_out(:)); % merging
    end
    fresults = GLCMStats(glcm);
    biomarkers = fieldnames(fresults);
    featuresSlice(:,j) = cell2mat(struct2cell(fresults));
end

avgFeatures = sum(featuresSlice,2)./size(offsetArray,1);
dim = size(avgFeatures,1);
intermediate = num2cell(avgFeatures,dim);
featuresE = cell2struct(intermediate,biomarkers,1);
end

% Calc GLCM
function featuresF = GLCMf(img, n_levels)
level_stack = expand_to_level(img);
glcm = zeros(n_levels, n_levels);
temp = unique(img);
vect1 = temp(~isnan(temp));
arrayOI = table2array(combinations(vect1,vect1));
for i = 1:size(arrayOI,1)
    two_level_img = level_stack(:,:,:,[arrayOI(i,1),arrayOI(i,2)]);
    offset_lim = [-1,0,1];
    offsetArray = [0,0,1; 0,1,0; 1,0,0; 0,1,1; 0,1,-1; 1,0,1; 1,0,-1; ...
        1,1,0; 1,-1,0; 1,1,1; 1,1,-1; 1, -1, 1; 1, -1, -1];
    for j = 1:size(offsetArray,1)
        filt = filter1(offsetArray(j,1), offsetArray(j,2),offsetArray(j,3));
        conv_out = convn(two_level_img, filt) == 2;
        glcm(arrayOI(i,1), arrayOI(i,2)) = glcm(arrayOI(i,1), arrayOI(i,2)) + sum(conv_out(:)); % merging
    end
end
normM = glcm./sum(glcm(:));
featuresF = GLCMStats(normM);
end

% set up binary image (i,j,k,level)
% given image with discrete levels return binary image for each level
function out = expand_to_level(img)
out = zeros([size(img), length(unique(img))]);
for i = min(img,[],'all'):max(img,[],'all')  % to do vectorize
    out(:,:,:,i) = img == i;
end
end

function out = offset_to_ind(i)
if i >= 0
    first_i = 1;
    second_i = i+1;
else
    first_i = 1-i;
    second_i = 1;
end
out = [first_i, second_i];
end

function out = filter1(i,j,k) %--> matrix with 1's at corners, size i,j,k
out = zeros(abs(i)+1,abs(j)+1,abs(k)+1,2);
temp = offset_to_ind(i); % first_i, second_i = offset_to_ind(i)
first_i = temp(1);
second_i = temp(2);
temp1 = offset_to_ind(j);
first_j = temp1(1);
second_j  = temp1(2);
temp2 = offset_to_ind(k);
first_k = temp2(1);
second_k = temp2(2);
out(first_i, first_j, first_k, 1) = 1;
out(second_i, second_j, second_k, 2) = 1;
end



