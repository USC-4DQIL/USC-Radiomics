%% CSV Export
% fills to the right of image file name and ROI name
% needs to put each image feature sequentially, same order as data
% dictionary
% will need [] placeholders that are not zero

function csvExporter(fileList,req,radiomics,outputFileName)

numFiles = size(fileList,1);

%% GLCM output file
if req.GLCM == 1
    % fill row by row --- 6x33 mat filed in
    GLCMfilename = strcat(outputFileName,'_GLCM.csv');
    firstGLCMtable = radiomics.metrics{1,1}.GLCM; % calling for dims, should be same dimensions throughout
    numAgMethods = size(firstGLCMtable,1);
    numFeatures = size(firstGLCMtable,2) - 1; % subtract column delineating ag method
    totNumFeatures = numAgMethods*numFeatures;
    % column 1 - image title
    % column 2 - ROI title
    % column 3 - aggregating method
    % column 4-end feature name

    featList = [];
    for ind = 1:totNumFeatures
        featName = strcat('GLCM_',string(ind));
        featList = [featList, featName];
    end

    varNames = ["Image","Classification",featList];
    x = 2+length(featList);
    varTypes = strings(1,x);
    varTypes(:) = "double";
    varTypes(1,1) = "string";


    outputGLCMTable = table('Size',[numFiles, 2+totNumFeatures],'VariableTypes',varTypes,'VariableNames',varNames);

    for filenum = 1:numFiles
        GLCMtable = radiomics.metrics{filenum,1}.GLCM;
        fileName = fileList(filenum,1);
        classificationName = fileList(filenum,2);

        className = classificationName{1,1};
        className = str2double(className);

        outputGLCMTable(filenum,1) = fileName;
        outputGLCMTable{filenum,2} = className;
        totRow = [];
        GLCMarray = table2array(GLCMtable(1:end,2:end));
        for x = 1:numAgMethods
            totRow = [totRow, GLCMarray(1,1:end)];
        end
        
    
        outputGLCMTable(filenum, 3:end) = array2table(totRow(1,:));
        
    end
    writetable(outputGLCMTable,GLCMfilename);
end

%% GLRLM
if req.GLRLM == 1
        % fill row by row --- 6x33 mat filed in
    GLRLMfilename = strcat(outputFileName,'_GLRLM.csv');
    firstGLRLMtable = radiomics.metrics{1,1}.GLRLM; % calling for dims, should be same dimensions throughout
    numAgMethods = size(firstGLRLMtable,1);
    numFeatures = size(firstGLRLMtable,2) - 1; % subtract column delineating ag method
    totNumFeatures = numAgMethods*numFeatures;
    % column 1 - image title
    % column 2 - ROI title
    % column 3 - aggregating method
    % column 4-end feature name

    featList = [];
    for ind = 1:totNumFeatures
        featName = strcat('GLRLM_',string(ind));
        featList = [featList, featName];
    end

    varNames = ["Image","Classification",featList];
    x = 2+length(featList);
    varTypes = strings(1,x);
    varTypes(:) = "double";
    varTypes(1,1) = "string";


    outputGLRLMTable = table('Size',[numFiles, 2+totNumFeatures],'VariableTypes',varTypes,'VariableNames',varNames);

    for filenum = 1:numFiles
        GLRLMtable = radiomics.metrics{filenum,1}.GLRLM;
        fileName = fileList(filenum,1);
        classificationName = fileList(filenum,2);

        className = classificationName{1,1};
        className = str2double(className);

        outputGLRLMTable(filenum,1) = fileName;
        outputGLRLMTable{filenum,2} = className;
        totRow = [];
        GLRLMarray = table2array(GLRLMtable(1:end,2:end));

        for x = 1:numAgMethods
            totRow = [totRow, GLRLMarray(1,1:end)];
        end
        
    
        outputGLRLMTable(filenum, 3:end) = array2table(totRow(1,:));
        
    end
    writetable(outputGLRLMTable,GLRLMfilename);
end

%% GLDZM

if req.GLDZM == 1
    GLDZMfilename = strcat(outputFileName,'_GLDZM.csv');
    firstGLDZMtable = radiomics.metrics{1,1}.GLDZM; % calling for dims, should be same dimensions throughout
    numAgMethods = size(firstGLDZMtable,1);
    numFeatures = size(firstGLDZMtable,2) - 1; % subtract column delineating ag method
    totNumFeatures = numAgMethods*numFeatures;
    % column 1 - image title
    % column 2 - ROI title
    % column 3 - aggregating method
    % column 4-end feature name

    featList = [];
    for ind = 1:totNumFeatures
        featName = strcat('GLDZM_',string(ind));
        featList = [featList, featName];
    end

    varNames = ["Image","Classification",featList];
    x = 2+length(featList);
    varTypes = strings(1,x);
    varTypes(:) = "double";
    varTypes(1,1) = "string";


    outputGLDZMTable = table('Size',[numFiles, 2+totNumFeatures],'VariableTypes',varTypes,'VariableNames',varNames);

    for filenum = 1:numFiles
        GLDZMtable = radiomics.metrics{filenum,1}.GLDZM;
        fileName = fileList(filenum,1);
        classificationName = fileList(filenum,2);

        className = classificationName{1,1};
        className = str2double(className);

        outputGLDZMTable(filenum,1) = fileName;
        outputGLDZMTable{filenum,2} = className;
        totRow = [];
        GLDZMarray = table2array(GLDZMtable(1:end,2:end));
        for x = 1:numAgMethods
            totRow = [totRow, GLDZMarray(1,1:end)];
        end
    
        outputGLDZMTable(filenum, 3:end) = array2table(totRow(1,:));
        
    end
    writetable(outputGLDZMTable,GLDZMfilename);
end



%% GLSZM

if req.GLSZM == 1

    GLSZMfilename = strcat(outputFileName,'_GLSZM.csv');
    firstGLSZMtable = radiomics.metrics{1,1}.GLSZM; % calling for dims, should be same dimensions throughout
    numAgMethods = size(firstGLSZMtable,1);
    numFeatures = size(firstGLSZMtable,2) - 1; % subtract column delineating ag method
    totNumFeatures = numAgMethods*numFeatures;
    % column 1 - image title
    % column 2 - ROI title
    % column 3 - aggregating method
    % column 4-end feature name

    featList = [];
    for ind = 1:totNumFeatures
        featName = strcat('GLSZM_',string(ind));
        featList = [featList, featName];
    end

    varNames = ["Image","Classification",featList];
    x = 2+length(featList);
    varTypes = strings(1,x);
    varTypes(:) = "double";
    varTypes(1,1) = "string";


    outputGLSZMTable = table('Size',[numFiles, 2+totNumFeatures],'VariableTypes',varTypes,'VariableNames',varNames);

    for filenum = 1:numFiles
        GLSZMtable = radiomics.metrics{filenum,1}.GLSZM;
        fileName = fileList(filenum,1);
        classificationName = fileList(filenum,2);

        className = classificationName{1,1};
        className = str2double(className);

        outputGLSZMTable(filenum,1) = fileName;
        outputGLSZMTable{filenum,2} = className;
        totRow = [];
        GLSZMarray = table2array(GLSZMtable(1:end,2:end));
        for x = 1:numAgMethods
            totRow = [totRow, GLSZMarray(1,1:end)];
        end
        
    
        outputGLSZMTable(filenum, 3:end) = array2table(totRow(1,:));
        
    end
    writetable(outputGLSZMTable,GLSZMfilename);
end


%% NGTDM
if req.NGTDM == 1
    
    NGTDMfilename = strcat(outputFileName,'_NGTDM.csv');
    firstNGTDMtable = radiomics.metrics{1,1}.NGTDM; % calling for dims, should be same dimensions throughout
    numAgMethods = size(firstNGTDMtable,1);
    numFeatures = size(firstNGTDMtable,2) - 1; % subtract column delineating ag method
    totNumFeatures = numAgMethods*numFeatures;
    % column 1 - image title
    % column 2 - ROI title
    % column 3 - aggregating method
    % column 4-end feature name

    featList = [];
    for ind = 1:totNumFeatures
        featName = strcat('NGTDM_',string(ind));
        featList = [featList, featName];
    end

    varNames = ["Image","Classification",featList];
    x = 2+length(featList);
    varTypes = strings(1,x);
    varTypes(:) = "double";
    varTypes(1,1) = "string";
    outputNGTDMTable = table('Size',[numFiles, 2+totNumFeatures],'VariableTypes',varTypes,'VariableNames',varNames);

    for filenum = 1:numFiles
        NGTDMtable = radiomics.metrics{filenum,1}.NGTDM;
        fileName = fileList(filenum,1);
        classificationName = fileList(filenum,2);

        className = classificationName{1,1};
        className = str2double(className);

        outputNGTDMTable(filenum,1) = fileName;
        outputNGTDMTable{filenum,2} = className;
        totRow = [];
        NGTDMarray = table2array(NGTDMtable(1:end,2:end));
        for x = 1:numAgMethods
            totRow = [totRow, NGTDMarray(1,1:end)];
        end
        outputNGTDMTable(filenum, 3:end) = array2table(totRow(1,:)); 
    end
    writetable(outputNGTDMTable,NGTDMfilename);
end


%% NGLDM
if req.NGLDM == 1
       NGLDMfilename = strcat(outputFileName,'_NGLDM.csv');
    firstNGLDMtable = radiomics.metrics{1,1}.NGLDM; % calling for dims, should be same dimensions throughout
    numAgMethods = size(firstNGLDMtable,1);
    numFeatures = size(firstNGLDMtable,2) - 1; % subtract column delineating ag method
    totNumFeatures = numAgMethods*numFeatures;
    % column 1 - image title
    % column 2 - ROI title
    % column 3 - aggregating method
    % column 4-end feature name

    featList = [];
    for ind = 1:totNumFeatures
        featName = strcat('NGLDM_',string(ind));
        featList = [featList, featName];
    end

    varNames = ["Image","Classification",featList];
    x = 2+length(featList);
    varTypes = strings(1,x);
    varTypes(:) = "double";
    varTypes(1,1) = "string";
    outputNGLDMTable = table('Size',[numFiles, 2+totNumFeatures],'VariableTypes',varTypes,'VariableNames',varNames);

    for filenum = 1:numFiles
        NGLDMtable = radiomics.metrics{filenum,1}.NGLDM;
        fileName = fileList(filenum,1);
        classificationName = fileList(filenum,2);

        className = classificationName{1,1};
        className = str2double(className);

        outputNGLDMTable(filenum,1) = fileName;
        outputNGLDMTable{filenum,2} = className;
        totRow = [];
        NGLDMarray = table2array(NGLDMtable(1:end,2:end));
        for x = 1:numAgMethods
            totRow = [totRow, NGLDMarray(1,1:end)];
        end
        outputNGLDMTable(filenum, 3:end) = array2table(totRow(1,:)); 
    end
    writetable(outputNGLDMTable,NGLDMfilename);
end


%% intensity

if req.Intensity == 1
Ifilename = strcat(outputFileName,'_Intensity.csv');
    firstIntenstable = struct2table(radiomics.metrics{1,1}.intensity); % calling for dims, should be same dimensions throughout

    numFeatures = size(firstIntenstable,2); % subtract column delineating ag method

    % column 1 - image title
    % column 2 - ROI title (bin for CIFAR images)
    % column 3 - end feature name
    % column 4 - actual feature value
    featList  =[];
    for ind = 1:numFeatures
        featName = strcat('Inten_',string(ind));
        featList = [featList, featName];
    end


    varNames = ["Image","Classification",featList];
    x = 2+length(featList);
    varTypes = strings(1,x);
    varTypes(:) = "double";
    varTypes(1,1) = "string";
    outputIntenTable = table('Size',[filenum 2+numFeatures],'VariableTypes',varTypes,'VariableNames',varNames);
    
    totRowIndex = 1;
    for filenum = 1:numFiles
        Itable = struct2table(radiomics.metrics{filenum,1}.intensity);
        fileName = fileList(filenum,1);
        classificationName = fileList(filenum,2);

        className = classificationName{1,1};
        className = str2double(className);

        outputIntenTable(filenum,1) = fileName;
        outputIntenTable{filenum,2} = className;


        outputIntenTable(filenum,3:end) = Itable;

    end
    writetable(outputIntenTable,Ifilename);
end


if req.Histogram == 1

    numFiles = size(fileList,1);
    Hfilename = strcat(outputFileName,'_Histogram.csv');
    firstHistotable = struct2table(radiomics.metrics{1,1}.histogram); % calling for dims, should be same dimensions throughout

    numFeatures = size(firstHistotable,2); % subtract column delineating ag method

    % column 1 - image title
    % column 2 - ROI title (bin for CIFAR images)
    % column 3 - end feature name
    % column 4 - actual feature value
    featList  =[];
    for ind = 1:numFeatures
        featName = strcat('Histo_',string(ind));
        featList = [featList, featName];
    end


    varNames = ["Image","Classification",featList];
    x = 2+length(featList);
    varTypes = strings(1,x);
    varTypes(:) = "double";
    varTypes(1,1) = "string";
    outputHistoTable = table('Size',[numFiles 2+numFeatures],'VariableTypes',varTypes,'VariableNames',varNames);
    
    totRowIndex = 1;
    for filenum = 1:numFiles
        Htable = struct2table(radiomics.metrics{filenum,1}.histogram);
        fileName = fileList(filenum,1);
        classificationName = fileList(filenum,2);

        className = classificationName{1,1};
        className = str2double(className);

        outputHistoTable(filenum,1) = fileName;
        outputHistoTable{filenum,2} = className;


        outputHistoTable(filenum,3:end) = Htable;

    end
    writetable(outputHistoTable,Hfilename);
end


end
