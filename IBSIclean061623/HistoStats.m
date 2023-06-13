

function fresult = HistoStats(imgData,binNum,binWidth)

 

%   Building the Histogram 
    if binNum == 0 && binWidth ~= 0
        
        minGl =  min(imgData(:));
        fresult.minbefore = minGl;
        maxGl = max(imgData(:));
        fresult.maxbefore = maxGl;
        BinCenters = minGl : binWidth : maxGl; %gives the edges 
        binFreq = histcounts(imgData(:),BinCenters);
        fresult.hist = binFreq;
        nG = length(binFreq);
        nV =  nnz(~isnan(imgData));
        binInd = 1:nG;
        pI = binFreq/nV;
        imgDatanew = repelem(binInd,binFreq);

        %{
        minGl =  min(imgData(:));
        maxGl = max(imgData(:));
        imgDatanew = floor((imgData - minGl)/binWidth) + 1;
        fResult.data = imgDatanew;
        nG = max(imgDatanew(:));
        BinCenters = 1:nG
        nV =  nnz(~isnan(imgDatanew))
        binInd = 1:nG;
        binFreq = histcounts(imgDatanew,[BinCenters, BinCenters(end)+1])
        pI = binFreq/nV;
        %}

        %{
        minGl =  min(min(imgData));
        maxGl = max(max(imgData));
        BinCenters = minGl : binWidth : maxGl;
        fResult.hist = hist(imgData(:),36);
        discretInten = 1:length(fResult.hist);
        binFreq = fResult.hist;
        imgDatanew = repelem(discretInten,binFreq);
        nV = sum(binFreq);
        nG = length(fResult.hist);
        pI = fResult.hist/nV;
        binInd = 1:nG;
        %}
        
    elseif binNum ~= 0 && binWidth == 0
      
        fresult.hist = histogram(imgData,binNum);
        minGl = min(min(min(imgData))); maxGl = max(max((max(imgData))));
        BinCenters = minGl:1:maxGl;
        binFreq = fresult.hist.Values;
        discretInten = 1:binNum;
        imgDatanew = repelem(discretInten,binFreq);
        nV = sum(binFreq);
        nG = fresult.hist.NumBins;
        pI = binFreq./nV;
        binInd = 1:nG;
        
    end
   

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
    if nV < 10
        fresult.tenPer = 0;
        fresult.ninePer = 0;
    else
            % 10th Intensity Percentile
        tenInd = round(0.1*(nV));
        fresult.tenPer = imgDatanew(tenInd);
    
            %90th Intensity Percentile
        nineInd = round(0.9*(nV));
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
        fresult.mode = imgDataunique(binFreq==max(binFreq));
    end

%       Discretized intensity interquartile range
    quartInd = round(.25*(nV));
    quartPer = imgDatanew(quartInd);
    sevenInd = round(.75*(nV));
    sevenPer = imgDatanew(sevenInd);
    fresult.interquart = sevenPer - quartPer;

%       Discretized intensity range
    fresult.range = fresult.max - fresult.min;

%       Intensity histogram mean absolute deviation
    fresult.ihmad = sum(abs(imgDatanew-fresult.mean))/nV;

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
    fresult.ihmedian = sum(abs(imgDatanew - fresult.median))/nV;

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
    fresult.maxgrad = max(hGradient);

%       Maximum histogram gradient intensity
    maxgradIndex = find(hGradient == max(hGradient));
    fresult.maxgradint = maxgradIndex;

%       Minimum histogram gradient 
    fresult.mingrad = min(hGradient);

%       Minimum histogram gradient intensity
    mingradIndex = find(hGradient == min(hGradient));
    fresult.mingradint = mingradIndex;

%       Fractional volume
    vSum = zeros(nV,nG);
    for i = 1:nG
        for k = 1:nV
            if imgData(k) < i
            vSum(k,i) = 1;
            end
        end
    end
    fresult.ivhv = 1 - (sum(vSum)/nV);

%       Intensity fraction
    y = 1:nG;
    fresult.ivhy = (y - fresult.min)/ (fresult.max - fresult.min);

%       Volume at 10 % intensity fraction
    vInd = find(fresult.ivhy>=0.10);
    fresult.ivhV10 = max(fresult.ivhv(vInd));

%       Volume at 90 % intensity fraction
    vInd = find(fresult.ivhy>=0.90);
    fresult.ivhV90 = max(fresult.ivhv(vInd));

%       Intensity at 10 % volume fraction
    histraw = histcounts(imgDatanew,[BinCenters, BinCenters(end)+1]);
    HistAc = cumsum(histraw);
    V = 1 - (HistAc ./ HistAc(end));
    V = [1, V(1:end-1)];
    fresult.ivhI10 = BinCenters(min(length(BinCenters) , find(V<=0.1,1,'first')));

%       Intensity at 90 % volume fraction
    fresult.ivhI90 = BinCenters(min(length(BinCenters) , find(V<=0.9,1,'first')));
    
%       Volume fraction difference between intensity fractions
    fresult.V10minuV90 = fresult.ivhV10 - fresult.ivhV90;

%       Intesnity fraction difference between volume fractions
    fresult.I10minusI90 = fresult.ivhI10 - fresult.ivhI90;

end