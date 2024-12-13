%% FBN
% discretization based on a fixed bin number
% assuming input is 3D and already in array form
function disImg = FBN(image,numLevels)
minimum = min(min(min(image)));
maximum = max(max(max(image)));
range = maximum - minimum;
interval = ceil(range/numLevels); % the number of bins and size of the bins is related

% top level is a < not <=; need to add an extra bin above it
disImg = discretize(image,numLevels);

end
