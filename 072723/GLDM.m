function featuresGLDM = GLDM(winImg,numLevels)
dynRange = [0,abs(max(winImg(:))-min(winImg(:)))];
flag = 1; % can be 1, 2, or 3
zDim = size(winImg,3);
slice2D = winImg(:,:,round(zDim/2)); % replace with widest image
GLDM2D = GLDM_auto2D(slice2D,dynRange,numLevels);
GLDM3D = GLDM_auto3D(winImg,dynRange,numLevels);

features2D = GLDMStats(GLDM2D);
features3D = GLDMStats(GLDM3D);
featuresGLDM = [features2D,features3D];

%% Aggregating Methods (same as IBSI GLCM)

% a -- 2D, by slice, w/o merging
% b -- 2D, by slice, w/ merging by slice
% c -- 2.5D, by slice, w/ merging by direction
% d -- 2.5D, by slice, w/ full merging
% e -- 3D, as volume, w/o merging
% f -- 3D, as volume, w/ full merging

function featuresA = GLDMA(winImg,dynRange,numLevels)

end

function featuresB = GLDMB(winImg,dynRange,numLevels)
        
end

function featuresC = GLDMC(winImg,dynRange,numLevels)
        
end

function featuresD = GLDMD(winImg,dynRange,numLevels)
        
end

function featuresE = GLDME(winImg,dynRange,numLevels)
        
end


%% 2D USC GLDM Calculation
function fDiffMatrix=GLDM_auto2D(slice2D,dynRange,NumLevels)
fDifMatrix=[];
for q=1:4
    temp=graydiffmatrix2D(slice2D,q,dynRange,NumLevels);
    fDifMatrix=matrixAdd(fDifMatrix,temp);
end
fDiffMatrix=fDifMatrix/4;
end

function fdMatrix=graydiffmatrix2D(imgVol,flag,dynRange,NumLevels)
volSize=size(imgVol);
fdMatrix=[];

 
    switch flag
        case 1
            temp=matrixDiff(imgVol(:,:)',dynRange,NumLevels);
        case 2
            temp=matrixDiff(imgVol(:,:),dynRange,NumLevels);
        case 3
            counter=1;
            tempVol=nan(length(-(volSize(1)-1):volSize(2)-1),...
                min(volSize(1:2)));
            for q=-(volSize(1)-1):volSize(2)-1
                tempVal=diag(imgVol(:,:),q);
                tempVol(counter,1:length(tempVal))=tempVal';
                counter=counter+1;
            end
            temp=matrixDiff(tempVol,dynRange,NumLevels);
        case 4
            counter=1;
            tempVol=nan(length(-(volSize(1)-1):volSize(2)-1),...
                min(volSize(1:2)));
            for q=-(volSize(1)-1):volSize(2)-1
                tempVal=diag(rot90(imgVol(:,:),1),q);
                tempVol(counter,1:length(tempVal))=tempVal';
                counter=counter+1;
            end
            temp=matrixDiff(tempVol,dynRange,NumLevels);
    end
    
    fdMatrix=matrixAdd(fdMatrix,temp);
end

%% 3D USC GLDM Calculation
function fDiffMatrix=GLDM_auto3D(winImg,dynRange,NumLevels)
fDifMatrix=[];
for z=1:3
     RotwinImg=shiftdim(winImg,z);   
     for q=z:4
         temp=graydiffmatrix3D(RotwinImg,q,dynRange,NumLevels);
         fDifMatrix=matrixAdd(fDifMatrix,temp);
     end
end
fDiffMatrix=fDifMatrix/9;
end

function fdMatrix=graydiffmatrix3D(imgVol,flag,dynRange,NumLevels)
volSize=size(imgVol);
fdMatrix=[];

 for k=2:volSize(3)-1
    switch flag
        case 1
            temp=matrixDiff(imgVol(:,:,k)',dynRange,NumLevels);
        case 2
            temp=matrixDiff(imgVol(:,:,k),dynRange,NumLevels);
        case 3
            counter=1;
            tempVol=nan(length(-(volSize(1)-1):volSize(2)-1),...
                min(volSize(1:2)));
            for q=-(volSize(1)-1):volSize(2)-1
                tempVal=diag(imgVol(:,:,k),q);
                tempVol(counter,1:length(tempVal))=tempVal';
                counter=counter+1;
            end
            temp=matrixDiff(tempVol,dynRange,NumLevels);
        case 4
            counter=1;
            tempVol=nan(length(-(volSize(1)-1):volSize(2)-1),...
                min(volSize(1:2)));
            for q=-(volSize(1)-1):volSize(2)-1
                tempVal=diag(rot90(imgVol(:,:,k),1),q);
                tempVol(counter,1:length(tempVal))=tempVal';
                counter=counter+1;
            end
            temp=matrixDiff(tempVol,dynRange,NumLevels);
    end
    
    fdMatrix=matrixAdd(fdMatrix,temp);
end
    
end

function fmatrix=matrixDiff(tempIMG,dynRange,NumLevels)
% Systematically calculates the difference between each number progressing
% for the left of the matrix to the right of the matrix and sorting the
% resulting histogram by distance between the voxels.
fmatrix=[];
finalBin=1;
for k=1:size(tempIMG,2)-3
    temp=abs(tempIMG(:,2:size(tempIMG,2)-k)-tempIMG(:,2+k:size(tempIMG,2)));
    temp=histc(temp(~isnan(temp))',...
        dynRange(1):dynRange(2)/NumLevels:dynRange(2));
        
    % Controls whether or not the last bin is just the last scalar value or
    % if the value gets added to the second to last bin.
    if finalBin==1
        fmatrix(k,:)=temp(1:NumLevels+1);
    else
        fmatrix(k,:)=temp(1:NumLevels);
        fmatrix(k,NumLevels)=fmatrix(k,NumLevels)+temp(NumLevels+1);
    end
end
end

function fmatrix=matrixAdd(oldMat,newMat)
if isempty(oldMat)
    fmatrix=newMat;
elseif isempty(newMat)
    fmatrix=oldMat;
end
oldMatSize=size(oldMat);
newMatSize=size(newMat);
newSize=max([oldMatSize; newMatSize]);
regOld=zeros(newSize);
regNew=zeros(newSize);
regOld(1:oldMatSize(1),1:oldMatSize(2))=oldMat;
regNew(1:newMatSize(1),1:newMatSize(2))=newMat;
fmatrix=regOld+regNew;
end


%% GLDM Feature Calculation
% delineated in pyradiomics online database
%https://pyradiomics.readthedocs.io/en/latest/features.html#radiomics.gldm.RadiomicsGLDM

function A = GLDMStats(GLDM)
    
    Nz = sum(GLDM(:));
    Ng = size(GLDM,1);
    Nd = size(GLDM,2);
    
    pij = GLDM./Nz;
    
    mui = 0;
    for i = 1:Ng
        for j = 1:Nd
            mui = mui + i*pij(i,j);
        end
    end
    
    muj = 0;
    for i = 1:Ng
        for j = 1:Nd
            muj = muj + j*pij(i,j);
        end
    end
    
    % Small Dependence Emphasis
    % A measure of the distribution of small dependencies, with a greater value 
    % indicative of smaller dependence and less homogeneous textures.
    term = 0;
    for i = 1:Ng
        for j = 1:Nd
            term = term + (GLDM(i,j)/(i^2));
        end
    end
    A.SDE = term/Nz;
    
    % Large Dependence Emphasis
    % A measure of the distribution of large dependencies, with a greater value 
    % indicative of larger dependence and more homogeneous textures.
    term = 0;
    for i = 1:Ng
        for j = 1:Nd
            term = term + (GLDM(i,j) * j^2);
        end
    end
    A.LDE = term/Nz;
    
    % Gray Level Non-Uniformity
    % Measures the similarity of gray-level intensity values in the image, 
    % where a lower GLN value correlates with a greater similarity in intensity values.
    term = 0;
    indivsum = 0;
    for i = 1:Ng
        for j = 1:Nd
            indivsum = indivsum + GLDM(i,j);
        end
        term = term + (indivsum)^2;
        indivsum = 0;
    end
    A.GLN = term/Nz;
    
    % Normalized Gray Level Non-Uniformity
    % mathematically equal to First Order - Uniformity
    term = 0;
    indivsum = 0;
    newsum = 0;
    for i = 1:Ng
        for j = 1:Nd
            indivsum = indivsum+GLDM(i,j);
            newsum = newsum + (GLDM(i,j)^2);
        end
        term = term + indivsum^2;
        indivsum = 0;
    end
    A.GLNN = indivsum/newsum;
    
    %Dependence Non-Uniformity
    % Measures the similarity of dependence throughout the image, with 
    %a lower value indicating more homogeneity among dependencies in the image.
    term = 0;
    indivsum = 0;
    for j = 1:Nd
        for i = 1:Ng
            indivsum = indivsum + GLDM(i,j);
        end
        term = term + (indivsum)^2;
        indivsum = 0;
    end
    A.DN = term/Nz;
    
    % Dependence Non-Uniformity Normalized
    A.DNN = term/(Nz^2);
    
    % Gray Level Variance
    term = 0;
    for i = 1:Ng
        for j = 1:Nd
            term = term + pij(i,j)*((i-mui)^2);
        end
    end
    A.GLV = term;
    
    % Dependence Variance
    term = 0;
    for i = 1:Ng
        for j = 1:Nd
            term = term + pij(i,j)*((j - muj)^2);
        end
    end
    A.DV = term;
    
    % Low Grey Level Emphasis
    term = 0;
    for i = 1:Ng
        for j = 1:Nd
            term = term + (GLDM(i,j)/(i^2));
        end
    end
    A.LGLE = term/Nz;
    
    % High Gray Level Emphasis
    term = 0;
    for i = 1:Ng
        for j = 1:Nd
            term = term + (GLDM(i,j)*(i^2));
        end
    end
    A.HGLE = term/Nz;
    
    % Small Dependence Low Grey Level Emphasis
    term = 0;
    for i = 1:Ng
        for j = 1:Nd
            term = term + (GLDM(i,j)/((i^2)*(j^2)));
        end
    end
    A.SDLGLE = term/Nz;
    
    % Large Dependence Low Gray Level Emphasis
    term = 0;
    for i = 1:Ng
        for j = 1:Nd
            term = term + ((GLDM(i,j)*(j^2))/((i^2)));
        end
    end
    A.LDLGLE = term/Nz;
    
    % Large Dependence High Gray Level Emphasis
    term = 0;
    for i = 1:Ng
        for j = 1:Nd
            term = term + (GLDM(i,j)*(i^2)*(j^2));
        end
    end
    A.LDHGLE = term/Nz;
end

end
