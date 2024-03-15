function featuresNGTDM = NGTDM1(winImg,numLevels,binWidth)
%% Define Key Variables
% winImg is analyzed Image
numLevels = max(max(max(winImg)));
Nv = numel(winImg) - sum(sum(sum(isnan(winImg))));
alpha = 0; % typical coarseness parameter, however can be changed
neiRad = floor(sqrt(2)); % neighborhood radius, is typically 1 but can be altered
featuresA = NGTDMa(winImg,numLevels,neiRad); % 2D, average over slices
featuresB = NGTDMb(winImg,numLevels,neiRad); % 2.5D, merge each slice texture matrix, then calculate features
featuresC = NGTDMc(winImg,numLevels,neiRad); % 3D

featuresInt = [featuresA,featuresB,featuresC];

agMethods = {'A','B','C'};
featuresNGTDM = renamevars([agMethods' struct2table(featuresInt)],'Var1','Aggregating_Methods'); % labeling step

%% Aggregating Methods
    function featuresA = NGTDMa(winImg,numLevels,neiRad)
        numSlices = size(winImg,3);
        featureStorage = zeros(5,size(winImg,3));

        % si
        % fill in grey levels for initializing NGTDM array
        for z = 1:size(winImg,3)
            siSlice=zeros(numLevels,1);
            slice2D = winImg(:,:,z);
            indivNGTDM = zeros(numLevels,4);
            indivNGTDM(:,1) = 1:numLevels;

            ni = zeros(numLevels,1);
            for index = 1:numLevels
                sliceNumLevels = sum(slice2D(:)==index);
                if sliceNumLevels >0
                    ni(index) = sliceNumLevels;
                else
                    ni(index) =0;
                end
            end
            indivNGTDM(:,2)=ni;
            paddedSlice = padarray(slice2D,[neiRad neiRad], NaN, 'both');


            for x = (1+neiRad):(size(paddedSlice,1)-neiRad)
                for y = (1+neiRad):(size(paddedSlice,2)-neiRad)

                    referencePixel = paddedSlice(x,y);
                    neighborhood = paddedSlice((x-neiRad):(x+neiRad),(y-neiRad):(y+neiRad));
                    editedNei = neighborhood(~isnan(neighborhood));
                    Wk = size(editedNei,1) -1; % careful making it w not wk
                    avgGL = (1/Wk) * (sum(editedNei)-referencePixel);


                    if referencePixel ~= 0 && ~isnan(referencePixel) && ~isnan(avgGL)...
                            && (Wk ~= 0)
                        iverson = abs(referencePixel - avgGL); % reference grey tone minus the average grey level for that neighborhood
                        siSlice(referencePixel) = siSlice(referencePixel) + iverson; %si
                    end

                end
            end


            nvSlice = sum(sum(~isnan(slice2D)));
            indivNGTDM(:,3) = indivNGTDM(:,2)./nvSlice; % pi
            indivNGTDM(:,4) = siSlice;

            sliceFeat = NGTDMStats(indivNGTDM); % calculate features
            biomarkers = fieldnames(sliceFeat);
            featureStorage(:,z) = cell2mat(struct2cell(sliceFeat)); % add together and merge each si column

        end

        arrayFeat = sum(featureStorage,2)./numSlices; % average features from the individual slices
        intermediate = num2cell(arrayFeat);
        featuresA = cell2struct(intermediate,biomarkers,1);

    end


    function featuresB = NGTDMb(winImg,numLevels,neiRad)
        NGTDMtot = zeros(numLevels,4);
        NGTDMtot(:,1) = 1:numLevels;

        % construct pi and ni
        for i = 1:numLevels
            NGTDMtot(i,2) = sum(sum(sum(winImg(:,:) == i)));
        end

        Nvc = sum(sum(sum(~isnan(winImg))));
        NGTDMtot(:,3) = NGTDMtot(:,2)./Nvc; %pi

        % si
        totSI = zeros(numLevels,1);
        % fill in grey levels for initializing NGTDM array
        for z = 1:size(winImg,3)
            siSlice=zeros(numLevels,1);
            slice2D = winImg(:,:,z);
            paddedSlice = padarray(slice2D,[neiRad neiRad], NaN, 'both');
            for x = (1+neiRad):(size(paddedSlice,1)-neiRad)
                for y = (1+neiRad):(size(paddedSlice,2)-neiRad)
                    referencePixel = paddedSlice(x,y);
                    neighborhood = paddedSlice((x-neiRad):(x+neiRad),(y-neiRad):(y+neiRad));
                    editedNei = neighborhood(~isnan(neighborhood));
                    Wk = size(editedNei,1) -1; % careful making it w not wk
                    avgGL = (1/Wk) * (sum(editedNei)-referencePixel);


                    if referencePixel ~= 0 && ~isnan(referencePixel) && ~isnan(avgGL)...
                            && (Wk ~= 0)
                        iverson = abs(referencePixel - avgGL); % reference grey tone minus the average grey level for that neighborhood
                        siSlice(referencePixel) = siSlice(referencePixel) + iverson; %si
                    end

                end
            end
            totSI = totSI + siSlice; % add together and merge each si column
        end
        NGTDMtot(:,4) = totSI;
        featuresB= NGTDMStats(NGTDMtot);
    end


    function featuresC = NGTDMc(winImg,numLevels,neiRad)
        % fill in grey levels for initializing NGTDM array
        NGTDMtot = zeros(numLevels,4);
        NGTDMtot(:,1) = 1:numLevels;

        % construct pi and ni
        for i = 1:numLevels
            NGTDMtot(i,2) = sum(sum(sum(winImg(:,:) == i)));
        end

        Nvc = sum(sum(sum(~isnan(winImg))));
        NGTDMtot(:,3) = NGTDMtot(:,2)./Nvc; %pi

        % si
        paddedWinImg = padarray(winImg,[neiRad neiRad neiRad],NaN,'both');
        % fill in grey levels for initializing NGTDM array
        for x = (1+neiRad):(size(paddedWinImg,1)-neiRad)
            for y = (1+neiRad):(size(paddedWinImg,2)-neiRad)
                for z = (1+neiRad):(size(paddedWinImg,3)-neiRad)
                    referencePixel = paddedWinImg(x,y,z);
                    neighborhood = paddedWinImg((x-neiRad):(x+neiRad),(y-neiRad):(y+neiRad),...
                        (z-neiRad):(z+neiRad));
                    editedNei = neighborhood(~isnan(neighborhood));
                    Wk = size(editedNei,1) -1; % careful making it w not wk
                    avgGL = (1/Wk) * (sum(editedNei)-referencePixel);

                    if referencePixel ~= 0 && ~isnan(referencePixel) && ~isnan(avgGL)...
                            && (Wk ~= 0)
                        iverson = abs(referencePixel - avgGL); % reference grey tone minus the average grey level for that neighborhood
                        NGTDMtot(referencePixel,4) = NGTDMtot(referencePixel,4) + iverson; %si
                    end

                end
            end
        end
        featuresC = NGTDMStats(NGTDMtot);
    end

%% Secondary Statistics Calculations
    function A = NGTDMStats(NGTDM)

        % separate pi,si,ni
        numLevels = size(NGTDM,1);
        s= NGTDM(:,4);
        p = NGTDM(:,3);
        n = NGTDM(:,2);
        A = struct();

        NBin = sum(p(:)>0);
        irange = 1:numLevels;
        Ng =  numLevels;

        Nv = sum(n(:));

        %% Secondary statistics

        % Coarseness
        term = 0;
        for i = 1:Ng
            term = term + (p(i).*s(i));
        end
        A.COARSE = 1 / term;

        % Contrast
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
        numerator = 0;
        for i = 1:Ng
            numerator = numerator + ( p(i) * s(i));
        end

        denominator = 0;
        for i1 = 1:Ng
            for i2 = 1:Ng
                if p(i1) ~= 0 && s(i2) ~=0
                    denominator = denominator + abs(i1 * p(i1) - i2 * p(i2));
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
                num = (p(i1)*s(i1))+ (p(i2)*s(i2));
                denom = (p(i1)+p(i2));
                if denom ~= 0 && p(i1) ~= 0 && p(i2) ~= 0 % added denom condition to formula
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

        % % Coarseness
        % term = 0;
        % for i = 1:numLevels
        %     term = term + (pi(i)*si(i));
        % end
        %
        % if term == 0
        %     A.COARSE = 10^6;
        % else
        %     A.COARSE = 1/term;
        % end
        %
        % % Contrast
        % term = 0;
        % nvn = sum(ni(:));
        % ngp = sum(pi(:)>0);
        % for i1 = 1:numLevels
        %     for i2 = 1:numLevels
        %         term = term + pi(i1)*pi(i2)*((i1-i2)^2);
        %     end
        % end
        %
        % intermediate = siSum*(1/nvn);
        %
        % A.CON = (1/((ngp)*(ngp-1)))*intermediate*term;
        %
        % % Busyness
        %
        % term1 = 0;
        % for i = 1:numLevels
        %     term1 = term1 + pi(i)*si(i);
        % end
        % term2 = 0;
        % for i1 = 1:numLevels
        %     for i2 = 1:numLevels
        %         aa = i1*pi(i1);
        %         bb = i2*pi(i2);
        %         term2 = term2 + abs(aa-bb);
        %     end
        % end
        %
        % if ngp == 1
        %     A.Busyness = 0;
        % else
        %     A.Busyness = term1/term2;
        % end
        %
        % % Complexity
        % term = 0;
        % for i1 = 1:numLevels
        %     for i2 = 1:numLevels
        %         if pi(i1) ~= 0 && pi(i2) ~= 0
        %             numerator = (pi(i1)*si(i1)) + (pi(i2)*si(i2));
        %             denom = pi(i1)+ pi(i2);
        %             if denom ~= 0
        %                 term = term + (abs(i1-i2)*(numerator/denom));
        %             end
        %         end
        %     end
        % end
        %
        % A.Complex = (1/nvn)*term;
        %
        % % Strength
        % term = 0;
        % for i1 = 1:numLevels
        %     for i2 = 1:numLevels
        %         term = term + (pi(i1) + pi(i2))*((i1-i2)^2);
        %     end
        % end
        %
        % if siSum == 0
        %     A.Strength = 0;
        % else
        %     A.Strength = term/siSum;
        % end
        % %
        % end
    end
end