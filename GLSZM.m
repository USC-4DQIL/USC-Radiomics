%% Gray Level Size Zone Matrix
function GLSZMfeatures = GLSZM(winImg,numLevels,binWidth)

numLevels = max(max(max(winImg)));
featuresA = GLSZMA(winImg); % 2D, by slice, w/o merging
featuresB = GLSZMB(winImg); % 2.5D, by slice, with merging
featuresC = GLSZMc(winImg); % 3D: as volume

featuresAC = [featuresA,featuresB,featuresC];

agMethods = {'A','B','C'};
conCat = [agMethods' struct2table(featuresAC)]; % labeling step
GLSZMfeatures = renamevars(conCat,'Var1','Aggregating_Methods');

    function featuresA = GLSZMA(winImg)
        zDim = size(winImg,3);
        for k = 1:zDim
            slice2D = winImg(:,:,k);
            tempG = GLSZMc(slice2D);
            names = fieldnames(tempG);
            cellResults = cell2mat(struct2cell(tempG));
            tot(:,k) = cellResults;
        end
        average = mean(tot,2);
        featuresA = cell2struct(num2cell(average),names,1);
    end

    function G = GLSZMB(original)

        Ng = max(max(max(original)));


        maxVal = max(max(max(original)));
        minVal = min(min(min(original)));


        % G is master struct to store all results
        % Build GLSZM in G.GLSZM
        % Initialize fields for secondary statistics based on the GLSZM
        G = struct();
        I_quant = round((Ng-1)*(original-minVal)/(maxVal-minVal))+1;
        % Collect all cell arrays of connected components
        CC_list = cell(Ng, 1);

        % Update with biggest connected component size
        % Will be used along with number of gray levels being considered to
        % preallocate the GLSZM matrix
        CC_max = 0;
        zDim = size(original,3);
        for k = 1:zDim
            slice2D = original(:,:,k);
            % Loop through all grey levels
            for i = 1 : Ng

                % Select voxels in image of this gray level
                B = I_quant == i;

                % Get cell array of indices of connected components of this
                % gray level
                CC = bwconncomp(B,8);

                CC_list{i} = CC.PixelIdxList;

                % Find size of largest connected component found
                largest_cc = max(cellfun(@length, CC_list{i}));

                % If larger than previously found largest connected component,
                % replace largest_cc value
                if largest_cc > CC_max
                    CC_max = largest_cc;
                end
            end


            % Preallocate and populate the GLSZM matrix
            % Matrix now has a row for each gray level considered and a column
            % for every possible size of connected component in the
            % image. CC_list is a cell array, one element for each gray level
            % considered, and each element in it is also a cell array of the
            % indices of connected components of that gray level. We don't care
            % what the indices were, only how big the connected component
            % was. The gray level is the row, and the size of the connected
            % component is the column, so we increment the matrix's value.
            Pnew = zeros(Ng, CC_max);
            for i = 1 : length(CC_list)
                for j = 1 : length(CC_list{i})
                    cc_size = length(CC_list{i}{j});
                    Pnew(i, cc_size) = Pnew(i, cc_size) + 1;
                end
            end

            if k == 1
                newDim = size(Pnew,2);
                Ptot = zeros(Ng,newDim);
            end

            Ptot = Ptot + Pnew; % merge size zone matrix
        end

        P = Ptot./zDim; % average over number of slices
        Nv = numel(original)-numel(find(isnan(original)));
        G = GLSZMStats(P,Nv);
    end

    function G = GLSZMc(original)

        Ng = max(max(max(original)));
        maxVal = max(original(:));
        minVal = min(original(:));
        I_quant = round((Ng-1)*(original-minVal)/(maxVal-minVal))+1;
        % G is master struct to store all results
        % Build GLSZM in G.GLSZM
        % Initialize fields for secondary statistics based on the GLSZM
        G = struct();

        % Collect all cell arrays of connected components
        CC_list = cell(Ng, 1);

        % Update with biggest connected component size
        % Will be used along with number of gray levels being considered to
        % preallocate the GLSZM matrix
        CC_max = 0;

        % Loop through all grey levels
        for i = 1 : Ng

            % Select voxels in image of this gray level
            B = I_quant == i;

            % Get cell array of indices of connected components of this
            % gray level
            if size(original,3) <= 1
                CC = bwconncomp(B,8);
            else
                CC = bwconncomp(B,26); % was 18, no observed change when shifting to 26-bit connectedness.
                % IBSI recommends 26 bit connectedness for 3D analysis, so implementing here
            end

            CC_list{i} = CC.PixelIdxList;

            % Find size of largest connected component found
            largest_cc = max(cellfun(@length, CC_list{i}));

            % If larger than previously found largest connected component,
            % replace largest_cc value
            if largest_cc > CC_max
                CC_max = largest_cc;
            end
        end

        % Preallocate and populate the GLSZM matrix
        % Matrix now has a row for each gray level considered and a column
        % for every possible size of connected component in the
        % image. CC_list is a cell array, one element for each gray level
        % considered, and each element in it is also a cell array of the
        % indices of connected components of that gray level. We don't care
        % what the indices were, only how big the connected component
        % was. The gray level is the row, and the size of the connected
        % component is the column, so we increment the matrix's value.
        P = zeros(Ng, CC_max);
        for i = 1 : length(CC_list)
            for j = 1 : length(CC_list{i})
                cc_size = length(CC_list{i}{j});
                P(i, cc_size) = P(i, cc_size) + 1;
            end
        end

        Nv = numel(original)-numel(find(isnan(original)));
        G = GLSZMStats(P,Nv);
    end

    function G = GLSZMStats(P,Nv)
        % Based on: https://github.com/cerr/CERR/wiki/GLSZM_global_features
        numLevels = size(P,1); % number of levels with linked pixels
        maxZoneSize = size(P,2); % max number of linked pixels
        % i = gray level i
        % j = number of linked pixels
        irange = 1:numLevels;
        jrange = 1:maxZoneSize;
        %
        totalNumZones = sum(sum(P(irange,jrange))); %total number of zones

        % num levels with linked pixels
        Ng = numLevels;
        % Largest zone size in image
        Nz = maxZoneSize;
        
        % Total number of zones in image
        Ns = totalNumZones;
        % Normalized size zone matrix
        p = P ./ totalNumZones;
        si = sum(P,2);
        sj = sum(P,1);

        % Useful index matrices
        [i_idx, j_idx] = meshgrid(1:Ng,1:Nz);

        % mu
        term = 0;
        for i = 1:numLevels
            for j = 1:maxZoneSize
                term = term + (i * p(i,j));
            end
        end
        muI = term;
        %G.mu = sum(sum(p(irange,jrange) * i_idx));

        term = 0;
        for i = 1:numLevels
            for j = 1:maxZoneSize
                term = term + (j * p(i,j));
            end
        end
        muJ = term;

        % Small area emphasis
        % Higher value indicates fine texture %% EDITED to match IBSI version
        G.SAE = (1/totalNumZones).* sum( sj./(jrange.^2));

        % Large area emphasis
        % Higher value indicates course texture
        G.LAE = (1/totalNumZones).* sum( sj.*(jrange.^2));

        % Gray level nonuniformity
        % Measures distribution of zone counts over the grey values. Value
        % is low when zone counts are equally distributed along grey
        % values.
        G.GLN = sum(sum(P, 2).^2) / Ns;

        % Gray level nonuniformity normalized
        G.GLNN = sum(sum(P, 2).^2) / (Ns^2);

        % Size zone nonuniformity
        % Measures distribution of zone counts over the zone sizes. Value
        % is low when zone counts are equally distributed along zone sizes.
        G.SZN = sum(sum(P, 1).^2) / Ns;

        % Size zone nonuniformity normalized
        G.SZNN = sum(sum(P, 1).^2) / (Ns^2);

        % Zone percentage
        % Measures the fraction of the number of realized zones and the
        % maximum number of potential zones. Uniformity produces low zone
        % percentage.
        G.ZP = Ns / Nv;

        % Low gray level zone emphasis
        % Gray level analogue to small zone emphasis. Low gray levels are
        % emphasized instead of small zone sizes.
        term = 0;
        for i = 1:numLevels
            term = term + (si(i)/i.^2);
        end
        G.LGLZE = (1/totalNumZones) * term;

        % High gray level zone emphasis
        % Gray level analogue to large zone emphasis. High gray levels are
        % emphasized.
        term = 0;
        for i = 1:numLevels
            term = term + (si(i)*i.^2);
        end
        G.HGLZE = (1/totalNumZones) * term;
        %G.HGLZE = sum(sum(P .* i_idx.^2)) / G.Ns;

        % Small area low gray level emphasis
        % Emphasizes zone counts in the upper left quadrant of the GLSZM
        % where small zone sizes and low gray levels are located.
        term = 0;
        for i = 1:numLevels
            for j = 1:maxZoneSize
                term = term + (P(i,j) / (i^2 * j^2));
            end
        end
        G.SALGLE = (1/totalNumZones) * term;
        %G.SALGLE = sum(sum(P ./ (i_idx.^2 .* j_idx.^2))) / G.Nz;

        % Small area high gray level emphasis
        % Emphasizes zone counts in the lower left quadrant of the GLSZM
        % where small zone sizes and high grey levels are loacated.
        term = 0;
        for i = 1:numLevels
            for j = 1:maxZoneSize
                term = term + ((P(i,j) *i^2) / (j^2));
            end
        end
        G.SAHGLE = (1/totalNumZones) * term;
        %G.SAHGLE = sum(sum(P .* (i_idx.^2 ./ j_idx.^2))) / G.Ns;

        % Large area low gray level emphasis
        % Emphasizes zone counts in upper right quadrant of the GLSZM where
        % large zone sizes and low gray levels are located.
        term = 0;
        for i = 1:numLevels
            for j = 1:maxZoneSize
                term = term + ((P(i,j) *j^2) / (i^2));
            end
        end
        G.LALGLE = (1/totalNumZones) * term;
        %G.LALGLE = sum(sum(P .* (j_idx.^2 ./ i_idx.^2))) / G.Ns;

        % Large area high gray level emphasis
        % Emphasizes zone counts in lower right quadrant of the GLSZM where
        % large zone sizes and high gray levels are locaed.
        term = 0;
        for i = 1:numLevels
            for j = 1:maxZoneSize
                term = term + (P(i,j) *j^2 * i^2);
            end
        end
        G.LAHGLE = (1/totalNumZones) * term;
        %G.LAHGLE = sum(sum(P .* (i_idx.^2 .* j_idx.^2))) / G.Ns;

        % Grey level variance
        % Variance in zone counts for the gray levels.
        term = 0;
        for i = 1:numLevels
            for j = 1:maxZoneSize
                term = term + (i - muI).^2 * p(i,j);
            end
        end
        G.GLV = term;

        % Size zone variance
        % Variance in zone counts for the different zone sizes.
        term = 0;
        for i = 1:numLevels
            for j = 1:maxZoneSize
                term = term + (j - muJ).^2 * p(i,j);
            end
        end
        G.SZV = term;
        %G.SZV = sum(sum(p .* (j_idx - G.mu).^2));
        NBin=numLevels;

        term = 0;
        temp=p;
        temp(temp~=0)=-log2(temp(temp~=0));
        for i = 1:numLevels
            for j = 1:maxZoneSize

                term = term + (p(i,j)*temp(i,j));
            end
        end
        G.zoneSizeENT = term;
    end

end
