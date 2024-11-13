%% FBN
% discretization based on a fixed bin number
% assuming input is 3D and already in array form
function disImg = FBN(image,numLevels)
minimum = min(min(min(image)));
maximum = max(max(max(image)));
range = maximum - minimum;
interval = ceil(range/numLevels); % the number of bins and size of the bins is related

disImg = floor(numLevels*(image-minimum)/(maximum-minimum))+1;
disImg(disImg>numLevels) = numLevels;
end
