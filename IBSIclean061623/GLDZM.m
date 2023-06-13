%% Created from IBSI Definition
% https://ibsi.readthedocs.io/en/latest/03_Image_features.html#grey-level-distance-zone-based-features

function featuresGLDZM = GLDZM(winImg,numLevels)
GLDZMa = GLDZMA(winImg); % 2D, average across slices
GLDZMb = GLDZMB(winImg); % 2.5D, calculate GLDZM per slice, merge GLDZM, calculate features
GLDZMc = GLDZM3DStats(winImg); % 3D

featuresGLDZM = [GLDZMa,GLDZMb,GLDZMc];
%% Functions

function G = GLDZMA(winImg)
    % first aggregating method
    % calculate features for each individual slice, average together for
    % final result
zDim = size(winImg,3);
for k = 1:zDim 
    slice2D = winImg(:,:,k);
    imagePerimeter = bwperim(slice2D,4);
    distmap =~ imagePerimeter + 1;
    tempG = GLDZM2DStats(slice2D,distmap);
    names = fieldnames(tempG);
    cellResults = cell2mat(struct2cell(tempG));
    tot(:,k) = cellResults;
end
average = mean(tot,2);
G = cell2struct(num2cell(average),names,1);   
end

function G = GLDZMB(winImg)
    % merge the GLDZM for each slice, use this to calculate features
Ng = 6;
maxVal = max(winImg(:));
minVal = min(winImg(:));
xDim = size(winImg,1);
yDim = size(winImg,2);
zDim = size(winImg,3);
dimensions = [xDim,yDim,zDim];
max_dim = max(dimensions);
Mtot = zeros(Ng, floor(max_dim/2));
for w = 1:zDim % loop through, calculate M matrix for each slice then merged
    slice2D = winImg(:,:,w);
    imagePerimeter = bwperim(slice2D,4);
    distmap =~ imagePerimeter + 1;
    Ng = 6; % numLevels
    minVal = min(slice2D(:));
    maxVal = max(slice2D(:));
    I_quant = round((Ng-1)*(slice2D-minVal)/(maxVal-minVal))+1;
    reshapedDistMap = reshape(distmap,(xDim*yDim),1);
% Loop through all grey levels
    CC_list = cell(Ng, 1);
    for i = 1 : Ng
        B = I_quant == i;
        CC = bwconncomp(B,8);
        CC_list{i} = CC; 
    end

    maxDist = max(distmap(:));
    M = zeros(Ng, floor(max_dim/2));
    for i = 1:Ng
        structA = CC_list{i};
        structB = structA.PixelIdxList;
        if length(structB) >= 1
            for a = 1: length(structB)
                pixelIDs = structB{a};
                iterations = size(pixelIDs);
                min_distance = 100; % arbitrary
                for c = 1:iterations
                    index = pixelIDs(c);
                    distance = reshapedDistMap(index);
                    if distance < min_distance
                        min_distance = distance;
                    end
            % of the linked pixels, need to find the minimum distance to ROI edge
                end
            M(i,min_distance) = M(i,min_distance) + 1;
            end
        end
    end 
Mtot = Mtot + M;
% merge all the matrices
end
% remove zero columns, since may have fewer key distances than maxDist
unused_distances = find(sum(Mtot, 1) == 0);
Mtot(:,unused_distances) = [];

Nd = size(Mtot,2);
Nv = numel(winImg)-numel(find(isnan(winImg)));
Ns = 0;
for i = 1:Ng
    for j = 1:Nd
        Ns = Ns + Mtot(i,j);
    end
end

G = GLDZMsecondStats(Mtot,winImg);

end

function G = GLDZM2DStats(slice2D,distmap)
xDim = size(slice2D,1);
yDim = size(slice2D,2);
Ng = 6; % numLevels
minVal = min(slice2D(:));
maxVal = max(slice2D(:));
I_quant = round((Ng-1)*(slice2D-minVal)/(maxVal-minVal))+1;
reshapedDistMap = reshape(distmap,(xDim*yDim),1);
% Loop through all grey levels
CC_list = cell(Ng, 1);
for i = 1 : Ng
    B = I_quant == i;
    CC = bwconncomp(B,8);
    CC_list{i} = CC; 
end

maxDist = max(distmap(:));
M = zeros(Ng, maxDist);
for i = 1:Ng
    structA = CC_list{i};
    structB = structA.PixelIdxList;
    if length(structB) >= 1
        for a = 1: length(structB)
            pixelIDs = structB{a};
            iterations = size(pixelIDs);
            min_distance = 100; % arbitrary
            for c = 1:iterations
                index = pixelIDs(c);
                distance = reshapedDistMap(index);
                if distance < min_distance
                    min_distance = distance;
                end
               
            % of the linked pixels, need to find the minimum distance to ROI edge
            end
         M(i,min_distance) = M(i,min_distance) + 1;
        end
    end
end 

% remove zero columns, since may have fewer key distances than maxDist
unused_distances = find(sum(M, 1) == 0);
M(:,unused_distances) = [];

G = GLDZMsecondStats(M,winImg);
end
      
function G = GLDZM3DStats(original)
Ng = 6; % number of gray levels
minVal = min(original(:));
maxVal = max(original(:));

%% Distance Map
zDim = size(original, 3);
distmap = zeros(size(original));
temp = original;
matMin = min(original(:));
for i = 1:size(original,1)
    for j = 1:size(original,2)
       for k = 1:size(original,3)
           if temp(i,j,k) < matMin || (isnan(temp(i,j,k)) == 1)
                temp(i,j,k) = matMin - 1;
           end
       end
   end
end

imagePerimeter = bwperim(temp);
distmap = imerode(temp,imagePerimeter)+1;

% Binning of dynamic range
I_quant = round((Ng-1)*(original-minVal)/(maxVal-minVal))+1;
% G is master struct to store all results
% Build GLSZM in G.GLSZM
% Initialize fields for secondary statistics based on the GLSZM
% Collect all cell arrays of connected components                          
CC_list = cell(Ng, 1);
% Update with biggest connected component size
% Will be used along with number of gray levels being considered to
% preallocate the GLSZM matrix
CC_max = 0;             
% Loop through all grey levels
for i = 1 : Ng

    % Select voxels in image of this gray level
    B = I_quant == i;

    % Get cell array of indices of connected components of this
    % gray level
    if size(original,3)==1
        CC = bwconncomp(B,8);
    else
        CC = bwconncomp(B,26); % was 18, no observed change when shifting to 26-bit connectedness.
        % IBSI recommends 26 bit connectedness for 3D analysis, so implementing here
    end

    CC_list{i} = CC.PixelIdxList;
    

    % Find size of largest connected component found
    largest_cc = max(cellfun(@length, CC_list{i}));

    totalPixels = size(original,1) * size(original,2) * size(original,3);
    newDistMap = reshape(distmap,totalPixels,1);
    % Find mininum distance of connected component
    if numel(CC.PixelIdxList) == 0
        distances(i,:) = 0;
        continue;
    else
        list = CC.PixelIdxList{1};
        for index = 1:size(list,1)
            element = list(index);
            temp_dist = newDistMap(element);
            if index == 1
                min_dist = temp_dist;
            elseif temp_dist < min_dist
                min_dist = temp_dist;
            end
        end
    
        distances(i,:) = min_dist;
    end
    
    % If larger than previously found largest connected component,
    % replace largest_cc value
    if largest_cc > CC_max
        CC_max = largest_cc;
    end
end
P = zeros(Ng, CC_max);
for i = 1 : length(CC_list)
    for j = 1 : length(CC_list{i})
        cc_size = length(CC_list{i}{j});
        P(i, cc_size) = P(i, cc_size) + 1;
    end
end

% variable definition
Ng = size(P,1);
Nd = size(P,2);
Nv = numel(original)-numel(find(isnan(original)));
Ns = 0;
for i = 1:Ng
    for j = 1:Nd
        Ns = Ns + P(i,j);
    end
end
x = sum(P,2);
Pnew = x.*distances; % matrix size num gray levels x max zone distance. W/ digital phantom, max zone distance = 1. 
D = Pnew;

G = GLDZMsecondStats(D,winImg);
end

function G = GLDZMsecondStats(matrix,winImg)
 G = struct();
%% Secondary statistics
% Based on: https://github.com/cerr/CERR/wiki/GLSZM_global_features

numLevels = size(matrix,1); % number of levels with linked pixels
maxDistance = size(matrix,2); % max number of linked pixels
% i = gray level i
% j = number of linked pixels
irange = 1:numLevels;
jrange = 1:maxDistance;

totalZone = sum(sum(matrix(irange,jrange))); %total number of zones

% num levels with linked pixels
Ng = numLevels;
% Largest zone size in image
Nz = maxDistance;
% Total number of voxels in image
Nv = numel(winImg)-numel(find(isnan(winImg)));
% Total number of zones in image
Ns = totalZone;
% Normalized size zone matrix
m = matrix ./ totalZone;
di = sum(matrix,2);
dj = sum(matrix,1);

% mu subscript i, intermediate variable
term = 0;
for i = 1:numLevels
    for j = 1:maxDistance
        term = term + (i * m(i,j));
    end
end
muI = term;
%G.mu = sum(sum(p(irange,jrange) * i_idx));

% mu subscript j, intermediate variable
term = 0;
for i = 1:numLevels
    for j = 1:maxDistance
        term = term + (j * m(i,j));
    end
end
muJ = term;

% Small distance emphasis
% Higher value indicates fine texture %% EDITED to match IBSI version
G.SDE = (1/totalZone).* sum( dj./(jrange.^2));

% Large distance emphasis
% Higher value indicates course texture
G.LDE = (1/totalZone).* sum( dj.*(jrange.^2));

% Gray level nonuniformity
% Measures distribution of zone counts over the grey values. Value
% is low when zone counts are equally distributed along grey
% values.
G.GLN = sum(sum(matrix, 2).^2) / Ns;

% Gray level nonuniformity normalized
G.GLNN = sum(sum(matrix, 2).^2) / (Ns^2);

% Zone distance nonuniformity
% Measures distribution of zone counts over the zone sizes. Value
% is low when zone counts are equally distributed along zone sizes.
G.ZDN = sum(sum(matrix, 1).^2) / Ns;

% Zone distance nonuniformity normalized
G.ZDNN = sum(sum(matrix, 1).^2) / (Ns^2);

% Zone percentage
% Measures the fraction of the number of realized zones and the
% maximum number of potential zones. Uniformity produces low zone
% percentage.
G.ZP = Ns / Nv;

% Low gray level distance emphasis
% Gray level analogue to small zone emphasis. Low gray levels are
% emphasized instead of small zone sizes.
term = 0;
for i = 1:numLevels
    term = term + (di(i)/i.^2);
end
G.LGLDE = (1/totalZone) * term;

% High gray level distance emphasis
% Gray level analogue to large zone emphasis. High gray levels are
% emphasized.
term = 0;
for i = 1:numLevels
    term = term + (di(i)*i.^2);
end
G.HGLDE = (1/totalZone) * term;
%G.HGLZE = sum(sum(P .* i_idx.^2)) / G.Ns;

% Small distance low gray level emphasis
% Emphasizes zone counts in the upper left quadrant of the GLSZM
% where small zone sizes and low gray levels are located.
term = 0;
for i = 1:numLevels
    for j = 1:maxDistance
        term = term + (matrix(i,j) / (i^2 * j^2));
    end
end
G.SDLGLE = (1/totalZone) * term;

% Small distance high gray level emphasis
% Emphasizes zone counts in the lower left quadrant of the GLSZM
% where small zone sizes and high grey levels are loacated.
term = 0;
for i = 1:numLevels
    for j = 1:maxDistance
        term = term + ((matrix(i,j) *i^2) / (j^2));
    end
end
G.SDHGLE = (1/totalZone) * term;
%G.SAHGLE = sum(sum(P .* (i_idx.^2 ./ j_idx.^2))) / G.Ns;

% Large distance low gray level emphasis
% Emphasizes zone counts in upper right quadrant of the GLSZM where
% large zone sizes and low gray levels are located.
term = 0;
for i = 1:numLevels
    for j = 1:maxDistance
        term = term + ((matrix(i,j) *j^2) / (i^2));
    end
end
G.LDLGLE = (1/totalZone) * term;
%G.LALGLE = sum(sum(P .* (j_idx.^2 ./ i_idx.^2))) / G.Ns;

% Large distance high gray level emphasis
% Emphasizes zone counts in lower right quadrant of the GLSZM where
% large zone sizes and high gray levels are locaed.
term = 0;
for i = 1:numLevels
    for j = 1:maxDistance
        term = term + (matrix(i,j) *j^2 * i^2);
    end
end
G.LDHGLE = (1/totalZone) * term;
%G.LAHGLE = sum(sum(P .* (i_idx.^2 .* j_idx.^2))) / G.Ns;

% Grey level variance
% Variance in zone counts for the gray levels.
term = 0;
for i = 1:numLevels
    for j = 1:maxDistance
        term = term + (i - muI).^2 * m(i,j);
    end
end
G.GLV = term;

% Zone distance variance
% Variance in zone counts for the different zone sizes.
term = 0;
for i = 1:numLevels
    for j = 1:maxDistance
        term = term + (j - muJ).^2 * m(i,j);
    end
end
G.ZDV = term;
%G.SZV = sum(sum(p .* (j_idx - G.mu).^2));
NBin=numLevels;

% Zone Distance Entropy
term = 0;
temp=m;
temp(temp~=0)=-log2(temp(temp~=0));
for i = 1:numLevels
    for j = 1:maxDistance
        
        term = term + (m(i,j)*temp(i,j));
    end
end
G.zoneDistENT = term;       
end

end
