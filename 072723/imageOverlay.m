%% ROI Segmentation

function windowedImage = imageOverlay(actual,ROI)
% initial image overlay, isolates ROI, NaN fills gaps
roiInfo=struct();
roiInfo.img = ROI;
roiInfo.roiSize=size(roiInfo.img);

% Converts ROI file to binary mask
if min(ROI(:))<0
    roiInd=find(ROI~=min(ROI(:)));
else
    roiInd=find(ROI>0);
end
roiInfo.roiInd=roiInd;

valueFile = actual;
tempImg=ROIImageCombiner(roiInfo,valueFile,'3D');


function tempImg =ROIImageCombiner(roiInfo,value_file,type)
valueImg=actual;
% Tests to see if there is a size mismatch between the ROI and Vaule files.
if sum(roiInfo.roiSize-size(valueImg))~=0
    errordlg('ROI and Value Images are mismatched in size.');
    imageStack=[];
    return;
end
roiInd=roiInfo.roiInd;
tempImg=nan(roiInfo.roiSize);
tempImg(roiInd)=valueImg(roiInd);
end

%
 % assuming ROI read in as 0's for excluded voxels

% limit window 
for k=1:3
    tempROI=shiftdim(tempImg,k);
    roiInd=find(tempROI>0);
    tempVal=shiftdim(tempImg,k);
    
    tempImg=nan(size(tempROI));
    tempImg(roiInd)=tempVal(roiInd);

    [indX,indY,indZ]=ind2sub(size(tempROI),roiInd);
end
        

windowedImage(:,:,:)=tempImg(min(indX):max(indX),min(indY):max(indY),...
       min(indZ):max(indZ));


end
