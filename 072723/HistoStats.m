function fresult = HistoStats(imgData,numLevels,binWidth)
XGL = imgData; % undiscretized image
range = max(max(max(XGL))) - min(min(min(XGL))); % range BEFORE discretization
NV = sum(sum(sum(~isnan(XGL)))); % number of pixels in whole image

%   Building the Histogram
if numLevels == 0 && binWidth ~= 0
    XDGL = FBS(XGL,binWidth); % array of discretized intensities
    numLevels = max(max(max(XDGL)));
    hist = histogram(imgData,numLevels,'HandleVisibility','off');
    
    nG = hist.NumBins;
    BinCenters = 1 : binWidth : (max(max(max(XDGL)))); %gives the edges
    binFreq =hist.Values;

    binInd = 1:nG;
    pI = binFreq/NV;
    imgDatanew = repelem(binInd,binFreq);
    

elseif numLevels ~= 0 && binWidth == 0
    XDGL = FBN(XGL,numLevels); % array of discretized intensities
    binWidth = ceil(range/numLevels);
    hist = histogram(imgData,numLevels,'HandleVisibility','off');
    
    nG = hist.NumBins;
    BinCenters = 1:binWidth:numLevels;
    binFreq = hist.Values;
    
    binInd = 1:nG;
    pI = binFreq./NV;
    imgDatanew = repelem(binInd,binFreq);
    
end
close all;

%       Building the Histogram Gradient
hGradient = zeros(1,nG);
for k = 1:length(hGradient)
    if k == 1
        hGradient(k) = binFreq(k+1) - binFreq(k);
    elseif k == nG
        hGradient(k) = binFreq(k) - binFreq(k-1);
    else
        hGradient(k) = (binFreq(k+1) - binFreq(k-1))/2;
    end
end

%       Mean discretized intensity
fresult.mean = sum(binInd .* pI);
fresult.mean1 = mean(imgDatanew,"omitnan");

%       Discretized intensity varaince
fresult.var = var(imgDatanew,"omitnan");  %MATLAB
fresult.var1 = sum(((binInd - fresult.mean).^2) .* pI);  %IBSI
%fResult.var2 = sum((imgDatanew-fResult.mean).^2)/nV; %IBSI

if fresult.var == 0
    fresult.skew = 0;
    fresult.kurt = 0;

else

    %       Discretized intensity skewness
    fresult.skew = skewness(imgDatanew);  %MATLAB

    skewNum1 =  sum((binInd - fresult.mean).^3 .* pI);
    skewDen1 = (sum((binInd - fresult.mean).^2 .* pI))^(3/2);
    fresult.skew1 = skewNum1/skewDen1; %IBSI

    %{
        skewNum2 = sum((imgDatanew - fResult.mean).^3)/nV;
        skewDen2 = (sum((imgDatanew - fResult.mean).^2) /nV)^(3/2);
        fResult.skew2 = skewNum2/skewDen2; %IBSI 
    %}

    %       Discretized intesnity kurtosis
    fresult.kurt = kurtosis(imgDatanew); %MATLAB

    kurtNum1 = sum((binInd - fresult.mean).^4 .* pI);
    kurtDen1 = (sum((binInd - fresult.mean).^2 .* pI))^2;
    fresult.kurt1 = kurtNum1/kurtDen1 - 3; %IBSI

    %{
        kurtNum2 = sum((imgDatanew - fResult.mean).^4)/nV;
        kurtDen2 = (sum((imgDatanew - fResult.mean).^2 )/nV)^2;
        fResult.kurt2 = kurtNum2/kurtDen2 - 3; %IBSI
    %}

end

%       Median discretized intensity
fresult.median = median(imgDatanew,"omitnan");

%       Minimum discretized intensity
fresult.min = min(imgDatanew);

%       Intensity Percentile
if NV < 10
    fresult.tenPer = 0;
    fresult.ninePer = 0;
else
    % 10th Intensity Percentile
    tenInd = round(0.1*(NV));
    fresult.tenPer = imgDatanew(tenInd);

    %90th Intensity Percentile
    nineInd = round(0.9*(NV));
    fresult.ninePer = imgDatanew(nineInd);
end

%       Maximum discretized intensity
fresult.max = max(imgDatanew);


%       Intensity histogram mode
imgDataunique = unique(imgDatanew);
varDataFrqmax = find(binFreq==max(binFreq));
if length(varDataFrqmax)>1
    [minValue,closestIndex] = min(abs(imgDataunique(varDataFrqmax)-fresult.mean));
    fresult.mode = imgDataunique(varDataFrqmax(closestIndex));
else
    fresult.mode = max(mode(imgDataunique));
end

%       Discretized intensity interquartile range
quartInd = round(.25*(NV));
quartPer = imgDatanew(quartInd);
sevenInd = round(.75*(NV));
sevenPer = imgDatanew(sevenInd);
fresult.interquart = sevenPer - quartPer;

%       Discretized intensity range
fresult.range = fresult.max - fresult.min;

%       Intensity histogram mean absolute deviation
fresult.ihmad = sum(abs(imgDatanew-fresult.mean))/NV;

%       Intesnity histogram robust mean absolute deviation
if fresult.tenPer == 0 && fresult.ninePer == 0
    fresult.ihrmad = fresult.ihmad;
else
    ind = imgDatanew >= fresult.tenPer & imgDatanew <= fresult.ninePer;
    imgDataihrmad = imgDatanew(ind);
    nVrmad = length(imgDataihrmad);
    fresult.ihrmad = sum(abs(imgDataihrmad - mean(imgDataihrmad)))/nVrmad;
end

%       Intensity histogram median absolute deviation
fresult.ihmedian = sum(abs(imgDatanew - fresult.median))/NV;

%       Intensity histogram coefficient of variation
fresult.std = std(imgDatanew);
fresult.ihcov = fresult.std/fresult.mean;
fresult.ihcov2 = (fresult.var1^(1/2))/fresult.mean;

%       Intensity histogram quartile coefficient of dispersion
fresult.ihqcod = (sevenPer - quartPer)/(sevenPer + quartPer);

%       Discretized intensity of entropy
sumpI = zeros(length(pI));
for i = 1:length(pI)
    if pI(i) > 0
        sumpI(i) = pI(i)*log2(pI(i));
    end
end

fresult.ihentropy = -sum(sumpI,'all');

%       Discretized intensity of uniformity/energy
fresult.ihuniformity = sum(pI.^2);

%       Maximum histogram gradient
fresult.maxgrad = max(max(hGradient));

%       Maximum histogram gradient intensity
maxgradIndex = find(hGradient == max(hGradient));
fresult.maxgradint = maxgradIndex;

%       Minimum histogram gradient
fresult.mingrad = min(hGradient);

%       Minimum histogram gradient intensity
mingradIndex = find(hGradient == min(hGradient));
fresult.mingradint = mingradIndex;


G = 1:numLevels; % range of discretized intensities
minG = min(G);
maxG = max(G);
IVH = zeros(numLevels,3);
IVH(:,1) = G;
for i = 1:numLevels
    newSum = sum(sum(sum(XDGL< i)));
    IVH(i,2) = 1 - (1/NV)*newSum; % fancy looking v, volume fraction
    IVH(i,3) = (i - minG)/(maxG - minG); % gamma, intensity fraction
end


% volume at intensity fraction

check10VI = (IVH(:,3) >= 0.1);

check90VI = (IVH(:,3) >= 0.9);

if sum(check10VI) == 0 
    fresult.ivhV10 = 0;
else
    fresult.ivhV10 = max(IVH(check10VI,2));
end

if sum(check90VI) == 0 
    fresult.ivhV90 = 0;
else
    fresult.ivhV90 = max(IVH(check90VI,2));
end


% intensity at volume fraction

check10IV = (IVH(:,2) < 0.1);

check90IV = (IVH(:,2) < 0.9);

if sum(check10IV) == 0
    fresult.ivhI10 = 0;
else
    fresult.ivhI10 = min(IVH(check10IV,1));
end

if sum(check90IV) == 0
    fresult.ivhI90 = 0;
else
    fresult.ivhI90 = min(IVH(check90IV,1));
end

fresult.V10minusV90 = fresult.ivhV10 - fresult.ivhV90;
fresult.I10minusI90 = fresult.ivhI10 - fresult.ivhI90;

end
