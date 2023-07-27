%% CSV Importer

function fileList=CSVImporter(filename)
if isempty(filename)
    [fname,pname,filterindex]=uigetfile(fullfile(pwd,'*.csv'),...
        'Open .csv format file list');
    if filterindex==0
        return;
    end
    filename=fullfile(pname,fname);
elseif isstruct(filename)
    filename=fullfile(filename.pname,filename.fname);
end

opts = detectImportOptions(filename,'Delimiter',',','TextType','char');
for k=1:length(opts.VariableTypes)
    opts.VariableTypes{1,k}='char';
end
% opts.VariableNamesLine=1;
% opts.DataLines=2;

% opts = delimitedTextImportOptions('Delimiter',',');
% opts=setvartype(opts,opts.VariableNames,'char');
fileList=readtable(filename,opts);
existCol=strcmpi('Var1',fileList.Properties.VariableNames);
if ~isempty(existCol(existCol==1))
    for k=1:length(fileList.Properties.VariableNames)
        fileList.Properties.VariableNames{...
            fileList.Properties.VariableNames{k}}=fileList.(...
            fileList.Properties.VariableNames{k}){1,1};
    end
    fileList(1,:)=[];
end
end