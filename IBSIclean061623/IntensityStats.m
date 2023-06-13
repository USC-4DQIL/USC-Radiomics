
%
% fresult = IntensityStats(imgData);

function fResult = IntensityStats(imgData)
    imgDatanew = double(sort(imgData(:)));
    nV = length(imgDatanew); %number of voxels
     
%       Mean 
    fResult.mean = mean(imgDatanew,"omitnan");

%       Intensity variance
    fResult.var = var(imgDatanew,"omitnan");
    fResult.var1 = sum((imgDatanew-fResult.mean).^2,"omitnan")/nV;

%       Intensity skewness
    fResult.skew = skewness(imgDatanew);
    
    skewNum1 = sum((imgDatanew - fResult.mean).^3,"omitnan")/nV;
    skewDen1 = (sum((imgDatanew - fResult.mean).^2,"omitnan") /nV)^(3/2);
    fResult.skew1 = skewNum1/skewDen1;

%       Excess intensity kurtosis
    fResult.kurt = kurtosis(imgDatanew);

    kurtNum1 = sum((imgDatanew - fResult.mean).^4,"omitnan")/nV;
    kurtDen1 = (sum((imgDatanew - fResult.mean).^2 ,"omitnan")/nV)^2;
    fResult.kurt1 = kurtNum1/kurtDen1 - 3; %IBSI

%       Median intensity
    fResult.median = median(imgDatanew,"omitnan");

%       Minimum intensity
    fResult.min = min(imgDatanew);

%       Standard Deviation
    fResult.std = std(imgDatanew,"omitnan");

%       Intensity Percentile
    if nV < 10
        fResult.tenPer = 0;
        fResult.ninePer = 0;
    else
        % 10th Intensity Percentile
        tenInd = round(0.1*(nV));
        fResult.tenPer = imgDatanew(tenInd);
    
        %90th Intensity Percentile
        nineInd = round(0.9*(nV));
        fResult.ninePer = imgDatanew(nineInd);
    end
%       Maximum intensity
    fResult.max = max(imgDatanew);

%       Intensity interquartile range
    quartInd = round(.25*(nV));
    quartPer = imgDatanew(quartInd);
    sevenInd = round(.75*(nV));
    sevenPer = imgDatanew(sevenInd);
    fResult.interquart = sevenPer - quartPer;

%       Intensity Range
    maxVal = max(imgDatanew);
    minVal = min(imgDatanew);
    fResult.range = maxVal - minVal;

%       Intensity-based mean absolute deviation
    fResult.mad = sum(abs(imgDatanew - fResult.mean),"omitnan")/nV;

%       Intensity-based robust mean absolute deviation
    if fResult.tenPer == 0 && fResult.ninePer == 0
        fResult.rmad = fResult.mad;
    else
        ind = imgDatanew >= fResult.tenPer & imgDatanew <= fResult.ninePer;
        imgDatarmad = imgDatanew(ind);
        nVrmad = length(imgDatarmad);
        fResult.rmad = sum(abs(imgDatarmad - mean(imgDatarmad)),"omitnan")/nVrmad;
    end

%       Intensity-based median absolute deviation
    fResult.statmedad = sum(abs(imgDatanew - fResult.median),"omitnan")/nV;

%       Intensity-based coeffiecient of variation
    fResult.cov = fResult.std/fResult.mean;
    fResult.cov1 = (fResult.var1^(1/2))/fResult.mean; %IBSI

%       Intensity-based quartile coeffiecient of dispersion
    fResult.qcod = (sevenPer - quartPer)/(sevenPer + quartPer);

%       Intensity-based energy
    fResult.energy = sum(imgDatanew.^2,"omitnan");

%       Root mean squared intensity
    fResult.rms = sqrt(fResult.energy/nV);
end


