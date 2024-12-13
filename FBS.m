%% FBS
% discretization based on a fixed bin size
% assuming input is 3D and already in array form
function disImg = FBS(image,binWidth)
maximum = max(max(max(image)));
minimum = min(min(min(image)));
rangeImage = maximum-minimum;
numLevels = ceil(rangeImage/binWidth); % the number of bins and size of the bins is related

disImg = discretize(image,ceil(numLevels+binWidth)); % fixed number of levels approach
% MATLAB discretization function < not <=, need to add an extra interval to
% not lop off the maximum

end
