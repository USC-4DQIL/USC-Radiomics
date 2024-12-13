function featuresNGLDM = octNGLDM(winImg)
%% Define Key Variables
% winImg is analyzed Image
numLevels = max(max(max(winImg)));
alpha = 0; % typical coarseness parameter, however can be changed
neiRad = floor(sqrt(2)); % neighborhood radius, is typically 1 but can be altered
%ngldm = NGLDM2D(winImg(:,:,1),numLevels,alpha,neiRad);
featuresA = NGLDMa(winImg,numLevels,alpha,neiRad); % 2D, average over slices
featuresB = NGLDMb(winImg,numLevels,alpha,neiRad); % 2.5D, merge each slice texture matrix, then calculate features
Nv = sum(sum(sum(~isnan(winImg))));
featuresC = NGLDMStats(NGLDM3D(winImg,numLevels,alpha,neiRad),Nv);  % 3D

featuresAC = cell(4,18);
featuresAC{1,1} = "AgMethods";
featuresAC{2,1} = "A";
featuresAC{3,1} = "B";
featuresAC{4,1} = "C";
featuresAC(1,2:18) = featuresA(1,1:17);
featuresAC(2,2:18) = featuresA(2,1:17);
featuresAC(3,2:18) = featuresB(2,1:17);
featuresAC(4,2:18) = featuresC(2,1:17);

featuresNGLDM = featuresAC;

    function twoDNGLDM = NGLDM2D(slice2D,numLevels,alpha,neiRad)
        NGLDM = zeros(numLevels,4*numLevels);
        % define neighborhood of the pixel
        for x = 1:size(slice2D,1)
            for y = 1:size(slice2D,2)
                referencePixel = slice2D(x,y);
                if isnan(referencePixel)
                    continue;
                end
                neighborhood = zeros(3);

                m = -neiRad:neiRad;
                for i = 1:length(m)
                    for j = 1:length(m)
                        if x+m(i) > size(slice2D,1) || y+m(j) > size(slice2D,2) ...
                                || x+m(i) < 1 || y+m(j) < 1
                            continue;
                        end
                        if isnan(slice2D(x+m(i),y+m(j)))
                            continue;
                        end
                        neighborhood(i,j) = slice2D(x+m(i),y+m(j));
                    end
                end

                n = numel(find(abs(referencePixel - neighborhood(:,:)) <= alpha));
                if n == 0
                    continue;
                end
                NGLDM(referencePixel,n) = NGLDM(referencePixel,n) + 1;
            end
        end
        twoDNGLDM = NGLDM;
    end

    function threeDNGLDM = NGLDM3D(winImg,numLevels,alpha,neiRad)
        NGLDM = zeros(numLevels,4*numLevels);
        % define neighborhood of the pixel
        for x = 1:size(winImg,1)
            for y = 1:size(winImg,2)
                for z = 1:size(winImg,3)
                    referencePixel = winImg(x,y,z);
                    if isnan(referencePixel)
                        continue;
                    end
                    neighborhood = zeros(3,3,3);

                    m = -neiRad:neiRad;
                    for i = 1:length(m)
                        for j = 1:length(m)
                            for k = 1:length(m)
                                if x+m(i) > size(winImg,1) || y+m(j) > size(winImg,2) ...
                                        || x+m(i) < 1 || y+m(j) < 1 ||  z+m(k) > size(winImg,3) ...
                                        || z+m(k) < 1
                                    continue;
                                end
                                if isnan(winImg(x+m(i),y+m(j),z+m(k)))
                                    continue;
                                end
                                neighborhood(i,j,k) = winImg(x+m(i),y+m(j),z+m(k));
                            end
                        end
                    end


                    n = numel(find(abs(referencePixel - neighborhood(:,:,:)) <= alpha));
                    if n == 0
                        continue;
                    end
                    NGLDM(referencePixel,n) = NGLDM(referencePixel,n) + 1;
                end
            end
        end
        threeDNGLDM = NGLDM;
    end

    function featuresA = NGLDMa(winImg,numLevels,alpha,neiRad)
        for k = 1:size(winImg,3)
            slice2D = winImg(:,:,k);
            tempngldm = NGLDM2D(slice2D,numLevels,alpha,neiRad);
            Nv = numel(slice2D) - sum(sum(isnan(slice2D)));
            tempFeatures = NGLDMStats(tempngldm,Nv);
            matrixA(:,k) = tempFeatures(2,:)';
        end

        numSlices = size(winImg,3);
        avgFeatures = sum(cell2mat(matrixA),2)./numSlices;
        dim = size(avgFeatures,1);
        intermediate = num2cell(avgFeatures,dim);
        featuresA(1,:) = tempFeatures(1,:);
        featuresA(2,:)= intermediate;
    end

    function featuresB = NGLDMb(winImg,numLevels,alpha,neiRad)
        totNv = 0;
        matrixA = zeros(numLevels,4*numLevels);
        for k = 1:size(winImg,3)
            slice2D = winImg(:,:,k);
            tempngldm = NGLDM2D(slice2D,numLevels,alpha,neiRad);
            matrixA = matrixA + tempngldm;
            Nv = numel(slice2D) - sum(sum(isnan(slice2D)));
            totNv = totNv + Nv;
        end

        numSlices = size(winImg,3);
        avgNGLDM = matrixA./numSlices;
        winImgNv = numel(winImg) - sum(sum(sum(isnan(winImg))));
        featuresB = NGLDMStats(matrixA,winImgNv);
       end




    function features = NGLDMStats(ngldm,Nv)
        F = cell(2,17);
        Ns = sum(sum(ngldm));
        si = sum(ngldm,2);
        sj = sum(ngldm,1);
        Nn = size(ngldm,2);
        Ng = size(ngldm,1);

        % low dependence emphasis
        term = 0;
        for j = 1:size(ngldm,2)
            term = term + (sj(j)./j.^2);
        end
        F{1,1} = "LDE";
        F{2,1} = (1/Ns) * term;

        %
        term = 0;
        for j = 1:size(ngldm,2)
            term = term + (sj(j)*j.^2);
        end
        F{1,2} = "HDE";
        F{2,2} = (1/Ns) * term;

        term = 0;
        for i = 1:size(ngldm,1)
            term = term + (si(i)./i.^2);
        end
        F{1,3} = "LDGLCE";
        F{2,3} = (1/Ns) * term;

        term = 0;
        for i = 1:size(ngldm,1)
            term = term + (si(i)*i.^2);
        end
        F{1,4} = "HGLCE";
        F{2,4} = (1/Ns) * term;


        term = 0;
        for i = 1:size(ngldm,1)
            for j = 1:size(ngldm,2)
                term = term + (ngldm(i,j)./(i^2 * j^2));
            end
        end
        F{1,5} = "LDLGLE";
        F{2,5} = (1/Ns) * term;


        term = 0;
        for i = 1:size(ngldm,1)
            for j = 1:size(ngldm,2)
                term = term + ((ngldm(i,j)*(i^2)) ./ j^2);
            end
        end
        F{1,6} = "LDHGLE";
        F{2,6} = (1/Ns) * term;



        term = 0;
        for i = 1:size(ngldm,1)
            for j = 1:size(ngldm,2)
                term = term + ((ngldm(i,j)*(j^2)) ./ i^2);
            end
        end
        F{1,7} = "HDLGLE";
        F{2,7} = (1/Ns) * term;


        term = 0;
        for i = 1:size(ngldm,1)
            for j = 1:size(ngldm,2)
                term = term + (ngldm(i,j)*(j^2)*(i^2));
            end
        end
        F{1,8} = "HDHGLE";
        F{2,8} = (1/Ns) * term;

        % gray level nonuniformity
        term = 0;
        for i = 1:Ng
            term = term + si(i).^2;
        end
        F{1,9} = "GLN";
        F{2,9} = (1/Ns) * term;

        term = 0;
        for i = 1:Ng
            term = term + si(i).^2;
        end
        F{1,10} = "GLNN";
        F{2,10}= (1/(Ns^2)) * term;

        term = 0;
        for j = 1:Nn
            term = term + sj(j).^2;
        end
        F{1,11} = "DCN";
        F{2,11} = (1/Ns) * term;

        term = 0;
        for j = 1:Nn
            term = term + sj(j).^2;
        end
        F{1,12} = "DCNN";
        F{2,12} = (1/(Ns.^2)) * term;

        F{1,13} = "DCP";
        F{2,13} = Ns/Nv;

        term = 0;
        p = ngldm/Ns;
        for i = 1:size(ngldm,1)
            for j = 1:size(ngldm,2)
                term = term + (p(i,j)*i);
            end
        end
        mui = term;

        term = 0;
        for i = 1:size(ngldm,1)
            for j = 1:size(ngldm,2)
                term = term + (p(i,j)*j);
            end
        end
        muj = term;

        term = 0;
        for i = 1:size(ngldm,1)
            for j = 1:size(ngldm,2)
                term = term + ((i - mui).^2*p(i,j));
            end
        end
        F{1,14} = "GLV";
        F{2,14} = term;


        term = 0;
        for i = 1:size(ngldm,1)
            for j = 1:size(ngldm,2)
                term = term + ((j - muj).^2*p(i,j));
            end
        end
        F{1,15} = "DCV";
        F{2,15} = term;

        term = 0;
        temp=p;
        temp(temp~=0)=-log2(temp(temp~=0));
        for i = 1:Ng
            for j = 1:Nn
                term = term + (p(i,j)*temp(i,j));
            end
        end
        F{1,16} = "DCEnt";
        F{2,16} = term;

        term = 0;
        for i = 1:Ng
            for j = 1:Nn
                term = term + p(i,j).^2;
            end
        end
        F{1,17} = "DCEnergy";
        F{2,17} = term;

        features = F;
    end

end
