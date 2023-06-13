%% IBSI USC Comparison

function featuresNGTDM = NGTDM(winImg,numLevels)
A = num2str(numLevels); % num gray levels
B = '1'; % neighborhood size
C = 'false'; % debug, true or false
varargin =char(A,B,C);
%%
NGTDMa = a(winImg,varargin); % 2D: compute features per 2D slice, average over slices
NGTDMb = b(winImg,varargin); % 2.5D: merge slice NGTDM's, compute features
NGTDMc = NGTDM(winImg,varargin); % 3D

featuresNGTDM = [NGTDMa,NGTDMb,NGTDMc];

function NGTDMa = a(winImg,varargin)
zDim = size(winImg,3);
for k = 1:zDim 
    slice2D = winImg(:,:,k);
    sliceNGTDM = NGTDM(slice2D,varargin);
    names = fieldnames(sliceNGTDM);
    cellResults = cell2mat(struct2cell(sliceNGTDM));
    tot(:,k) = cellResults;
end
average = mean(tot,2);
NGTDMa = cell2struct(num2cell(average),names,1);
end

function NGTDMb = b(original,varargin)
% Calculate the neighborhood graytone difference matrix of an input
% image
%
% The first argument must be the 2D or 3D input data, and the
% following arguments must be char-array, value pairs. The
% available options are:

% Version tracking for function
appVersion='1.0';

% Returns the version information when called.
if nargin==1 && strcmpi(original,'Version')
    A=appVersion;
    return
end

%% Set defaults

% Number of graytone bins to consider. Will be reassigned after the
% matrix is built to remove consideration of graytone bins not
% present in the data.
Ng = 6;

% Neighborhood size to consider around each voxel. d=1 implies a
% 3x3 region centered about the current voxel of interest, not
% counting the value of the center voxel. Will also be used to
% create a border of zeros around the entire dataset to allow
% consideration of edge voxels.
d = 1;

% Don't save debugging information unless set to true.
debug = false;

%% Handle optional arguments

%assert((mod(length(varargin), 2) == 0), 'Must pass optional parameters in name-value pairs.');
%acceptable_param_names = {'Ng', 'NeighborhoodSize', 'Debug'};

%for param = 1 : 2 : length(varargin)-1
    %parameter = varargin{param};
    %val = varargin{param + 1};
    %assert((any(strcmp(acceptable_param_names, parameter))), 'Invalid parameter name passed.');
    %switch parameter
    %  case 'Ng'
    %    Ng = val;
    %  case 'Neighborhoodsize'
    %    d = val;
    %  case 'Debug'
    %    debug = val;
%    end
%end

%% Main processing
depth = size(original,3);
npsTot = size(Ng,3);
sAdded = zeros(Ng,1);
pAdded = zeros(Ng,1);
nAdded = zeros(Ng,1);
for k = 1:depth
    slice2D = original(:,:,k);
    nps = zeros(Ng, 3);
% Quantize input data into graytone bins according to same process
% used in GLSZM algorithm
% Previously used code kept here for posterity, but Matlab's rescale()
% function works nicely even on inputs that contain data whose values don't
% go past the number of bins being considered. The code previously used
% sometimes broke when testing on data that, for example, only contained
% values between 1 and 5 and attempted to bin that into a space of bins
% from 1 to 20.
%max_val = max(original(:));
%min_val = min(original(:));
%original_quant = round((Ng-1)*(original-min_val)/(max_val-min_val))+1;
    original_quant = round(rescale(slice2D, 1, Ng));



% NGTDM matrix. One row per graytone bin and three columns. First
% column, ni, is the number of occurrances in the dataset of that
% graytone. Second column, pi, is the probability of probability of
% that graytone. Third column, si, is t he sum of absolute
% differences for graytone i for each neighborhood around an
% occurrance of that graytone in the dataset.


% Pad the input dataset with zeros for calculating si
% values. Keeping the original dataset around unmodified to
% calculate probabilities and other metrics without confusion over
% the added presence of zeros.
% original_padded = padarray(original_quant, [d, d], NaN);

% Get si values
% Only considering square window neighborhoods, so middle index will always
% be the same, so calculating and saving this index outside the loop

    original_padded = padarray(original_quant, [d, d], NaN);
    middle = ceil(((2 * d + 1)^2)/2);
    for a = 1 + d : size(original_padded, 1) - d
        for b = 1 + d : size(original_padded, 2) - d
            
            % Pull iteration's window from padded matrix
            w = original_padded(a-d:a+d, b-d:b+d);
            
            % Determine which graytone is at center of window for which
            % row of the nps matrix's si value to increment
            graytone = w(middle);
            
            if ~isnan(graytone)
                % Nan out the value at center of window so it's value doesn't
                % contribute to the calculation of the average. We use the
                % nanmean() later to prevent counting NaNs around the border as
                % well as prevent counting this voxel of interest.
                w(middle) = NaN;
                
                % Go to nps matrix and increment the value already present
                % by the new voxel's difference from neighborhood average
                %nps(graytone, 3) = nps(graytone, 3) + abs(graytone - (total/count));
                nps(graytone, 3) = nps(graytone, 3) + abs(graytone - mean(w(:)));                    
            end
        end
    end

%%
% Get ni values
for a = 1 : Ng
    nps(a, 1) = sum(original_quant(:) == a);
end

% Get pi values
voxel_count = sum(~isnan(original), 'all');
nps(:, 2) = nps(:, 1) / voxel_count;

npsTot = npsTot + nps;

nList = nps(:,1,:);
pList = nps(:,2,:);
sList = nps(:,3,:);

nAdded = nAdded + nList;
pAdded = pAdded + pList;
sAdded = sAdded + sList;
end

nps = [nAdded, pAdded, sAdded];

n = nAdded./depth;
p = pAdded./depth;
s = sAdded./depth;
% Delete rows where all data remain zeros
nps(sum(nps, 2) == 0, :) = [];

if debug
    % Print nps matrix for debugging purposes
    disp(nps);
end

% Number of graytones present in the image
nps_deleted = nps;
nps_deleted(sum(nps, 2) == 0, :) = [];

% Create struct to store results
A = struct();
NBin = size(nps_deleted, 1);
Ng = Ng;

irange = 1:Ng;

Nv = sum(nAdded(:));

%% Secondary statistics

% Coarseness
term = 0;
for i = 1:Ng
    term = term + (pAdded(i).*sAdded(i));
end
A.COARSE = 1 / term;

% Contrast
termA = 0;
for i1 = 1:Ng
    for i2 = 1:Ng
        termA = termA + (pAdded(i1)*pAdded(i2)*((i1-i2).^2));
    end
end
sTot = sum(sAdded(:));
constant1 = 1/((NBin)*(NBin - 1));
constant2 = (1/Nv);
A.CON = (constant1*termA)*(constant2*sTot);
% Busyness
% Other libraries calculate the difference_magnitudes value including the
% absolute value, despite there being no absolute value mentioned in the
% paper's original expression. However, excluding the abs() here appears to
% always result in the value being zero and causing busyness being
% infinite, so we include this strategy here.
numerator = 0;
for i = 1:Ng
    numerator = numerator + ( pAdded(i) * sAdded(i));
end

denominator = 0;
for i1 = 1:Ng
    for i2 = 1:Ng
        if pAdded(i1) ~= 0 && sAdded(i2) ~=0
        denominator = denominator + abs(i1 * pAdded(i1) - i2 * pAdded(i2));
        end
    end
end
A.Busyness = numerator/denominator;
% Complexity
constant = (1/Nv);
term = 0;
for i1 = 1:Ng
    for i2 = 1:Ng
        x = abs(i1-i2);
        num = (pAdded(i1)*sAdded(i1))+ (pAdded(i2)*sAdded(i2));
        denom = (pAdded(i1)+pAdded(i2));
        if denom ~= 0 && pAdded(i1) ~= 0 && pAdded(i2) ~= 0 % added denom condition to formula
            term = term + x*(num/denom);
        end
    end
end
A.Complex = constant * term;

% Texture strength
term = 0;
for i1 = 1:Ng
    for i2 = 1:Ng
        if p(i1) ~= 0 && p(i2) ~= 0
            term = term + (p(i1) + p(i2))*((i1-i2).^2);
        end
    end
end
denom = sum(s(:));
A.TextStrength = term/denom;

NGTDMb = A;
end

function NGTDMc = NGTDM(original,varargin)
% Calculate the neighborhood graytone difference matrix of an input
% image
%
% The first argument must be the 2D or 3D input data, and the
% following arguments must be char-array, value pairs. The
% available options are:

% Version tracking for function
appVersion='1.0';

% Returns the version information when called.
if nargin==1 && strcmpi(original,'Version')
    A=appVersion;
    return
end

%% Set defaults

% Number of graytone bins to consider. Will be reassigned after the
% matrix is built to remove consideration of graytone bins not
% present in the data.
Ng = 6;

% Neighborhood size to consider around each voxel. d=1 implies a
% 3x3 region centered about the current voxel of interest, not
% counting the value of the center voxel. Will also be used to
% create a border of zeros around the entire dataset to allow
% consideration of edge voxels.
d = 1;

% Don't save debugging information unless set to true.
debug = false;

%% Handle optional arguments

%assert((mod(length(varargin), 2) == 0), 'Must pass optional parameters in name-value pairs.');
%acceptable_param_names = {'Ng', 'NeighborhoodSize', 'Debug'};

%for param = 1 : 2 : length(varargin)-1
    %parameter = varargin{param};
    %val = varargin{param + 1};
    %assert((any(strcmp(acceptable_param_names, parameter))), 'Invalid parameter name passed.');
    %switch parameter
    %  case 'Ng'
    %    Ng = val;
    %  case 'Neighborhoodsize'
    %    d = val;
    %  case 'Debug'
    %    debug = val;
%    end
%end

%% Main processing

% Quantize input data into graytone bins according to same process
% used in GLSZM algorithm
% Previously used code kept here for posterity, but Matlab's rescale()
% function works nicely even on inputs that contain data whose values don't
% go past the number of bins being considered. The code previously used
% sometimes broke when testing on data that, for example, only contained
% values between 1 and 5 and attempted to bin that into a space of bins
% from 1 to 20.
%max_val = max(original(:));
%min_val = min(original(:));
%original_quant = round((Ng-1)*(original-min_val)/(max_val-min_val))+1;
original_quant = round(rescale(original, 1, Ng));


% Create struct to store results
A = struct();

% NGTDM matrix. One row per graytone bin and three columns. First
% column, ni, is the number of occurrances in the dataset of that
% graytone. Second column, pi, is the probability of probability of
% that graytone. Third column, si, is t he sum of absolute
% differences for graytone i for each neighborhood around an
% occurrance of that graytone in the dataset.
nps = zeros(Ng, 3);

% Pad the input dataset with zeros for calculating si
% values. Keeping the original dataset around unmodified to
% calculate probabilities and other metrics without confusion over
% the added presence of zeros.
% original_padded = padarray(original_quant, [d, d], NaN);

% Get si values
% Only considering square window neighborhoods, so middle index will always
% be the same, so calculating and saving this index outside the loop
if size(original,3)==1
    original_padded = padarray(original_quant, [d, d], NaN);
    middle = ceil(((2 * d + 1)^2)/2);
    for a = 1 + d : size(original_padded, 1) - d
        for b = 1 + d : size(original_padded, 2) - d
            
            % Pull iteration's window from padded matrix
            w = original_padded(a-d:a+d, b-d:b+d);
            
            % Determine which graytone is at center of window for which
            % row of the nps matrix's si value to increment
            graytone = w(middle);
            
            if ~isnan(graytone)
                % Nan out the value at center of window so it's value doesn't
                % contribute to the calculation of the average. We use the
                % nanmean() later to prevent counting NaNs around the border as
                % well as prevent counting this voxel of interest.
                w(middle) = NaN;
                
                % Go to nps matrix and increment the value already present
                % by the new voxel's difference from neighborhood average
                %nps(graytone, 3) = nps(graytone, 3) + abs(graytone - (total/count));
                nps(graytone, 3) = nps(graytone, 3) + abs(graytone - mean(w(:)));                    
            end
        end
    end
else
    original_padded = padarray(original_quant, [d, d, d], NaN);
    middle = ceil((3*(2 * d + 1)^2)/2);
    for a = 1 + d : size(original_padded, 1) - d
        for b = 1 + d : size(original_padded, 2) - d
            for c = 1 + d : size(original_padded, 3) - d
                % Pull iteration's window from padded matrix
                w = original_padded(a-d:a+d, b-d:b+d, c-d:c+d);
                
                % Determine which graytone is at center of window for which
                % row of the nps matrix's si value to increment
                graytone = w(middle);

                if ~isnan(graytone)
                    % Nan out the value at center of window so it's value doesn't
                    % contribute to the calculation of the average. We use the
                    % nanmean() later to prevent counting NaNs around the border as
                    % well as prevent counting this voxel of interest.
                    w(middle) = NaN;
                    
                    % Go to nps matrix and increment the value already present
                    % by the new voxel's difference from neighborhood average
                    %nps(graytone, 3) = nps(graytone, 3) + abs(graytone - (total/count));
                    nps(graytone, 3) = nps(graytone, 3) + abs(graytone - mean(w, 'all'));
                end
            end
        end
    end
end


% Get ni values
for a = 1 : Ng
    nps(a, 1) = sum(original_quant(:) == a);
end

% Get pi values
voxel_count = sum(~isnan(original), 'all');
nps(:, 2) = nps(:, 1) / voxel_count;

% Delete rows where all data remain zeros
% nps(sum(nps, 2) == 0, :) = [];

if debug
    % Print nps matrix for debugging purposes
    disp(nps);
end

% Number of graytones present in the image
nps_deleted = nps;
nps_deleted(sum(nps, 2) == 0, :) = [];
NBin = size(nps_deleted, 1);
Ng = Ng;

n = nps(:,1,:);
Nv = sum(n(:));

%% Secondary statistics

% Matrix sum preparation
p_i = nps(:, 2);
p_j = nps(:, 2);
s_i = nps(:, 3);
s_j = nps(:, 3);
ii = (1:numel(p_i))';
jj = (1:numel(p_j))';

zero_row = sum(nps, 2) == 0;
p_i(zero_row, :) = [];
p_j(zero_row, :) = [];
s_i(zero_row, :) = [];
s_j(zero_row, :) = [];
ii(zero_row, :) = [];
jj(zero_row, :) = [];

copies = numel(p_i);

p_i = repelem(p_i, copies);
s_i = repelem(s_i, copies);
ii = repelem(ii, copies);

p_j = repmat(p_j, copies, 1);
jj = repmat(jj, copies, 1);
s_j = repmat(s_j, copies, 1);

% Coarseness
% A.coarseness = 1 / sum(nps(:, 2) .* nps(:, 3));
A.COARSE = 1 / sum(nps(:, 2) .* nps(:, 3));

% Contrast
p = nps(:,2,:);
s = nps(:,3,:);
termA = 0;
for i1 = 1:Ng
    for i2 = 1:Ng
        termA = termA + (p(i1)*p(i2)*((i1-i2).^2));
    end
end
sTot = sum(s(:));
constant1 = 1/((NBin)*(NBin - 1));
constant2 = (1/Nv);
A.CON = (constant1*termA)*(constant2*sTot);

% Busyness
% Other libraries calculate the difference_magnitudes value including the
% absolute value, despite there being no absolute value mentioned in the
% paper's original expression. However, excluding the abs() here appears to
% always result in the value being zero and causing busyness being
% infinite, so we include this strategy here.
spacial_derivative = sum(nps(:, 2) .* nps(:, 3));
difference_magnitudes = sum(abs((ii .* p_i) - (jj .* p_j)));
% A.busyness = spacial_derivative / difference_magnitudes;
A.Busyness = spacial_derivative / difference_magnitudes;

% Complexity
% A.complexity = sum((abs(ii-jj)./((voxel_count^2)*(p_i+p_j))).*((p_i.*s_i)+(p_j.*s_j)));
A.Complex = sum((abs(ii-jj)./((voxel_count)*(p_i+p_j))).*((p_i.*s_i)+(p_j.*s_j)));

% Texture strength
% A.texture_strength = sum((p_i+p_j).*((ii-jj).^2)) / sum(nps(:, 3));
A.TextStrength = sum((p_i+p_j).*((ii-jj).^2)) / sum(nps(:, 3));

NGTDMc = A;
end

end