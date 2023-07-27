function dirFinal=dirAll(dirName,isuniqueDir)

home_dir=pwd;
if nargin==0
    dirName=home_dir;
    isuniqueDir=0;
elseif nargin==1
    isuniqueDir=0;
end
search_str='';
[pathstr, name, ext] = fileparts(dirName);
if ~isempty(pathstr)
    try
        cd(pathstr);
        if exist([name,ext],'dir')==7
            cd([name,ext]);
            dirFinal=dir;
        else
            dirFinal=dir([name,ext]);
            search_str=[name,ext];
        end
    catch ME
        disp('Error: Directory does not exist');
        dirFinal=0;
        return;
    end
    dirBase=pwd;
    
else
    dirBase=pwd;
    dirFinal=dir([name,ext]);
    search_str=[name,ext];
end

for k=1:length(dirFinal)
    dirFinal(k).dirBase=dirBase;
end


  dirData = dir(dirBase);      %Get the data for the current directory
  dirIndex = [dirData.isdir];  % Find the index for directories
  fileList = {dirData(~dirIndex).name}';  %Get a list of the files
  subDirs = {dirData(dirIndex).name};  %Get a list of the subdirectories
  validIndex = ~ismember(subDirs,{'.','..'});  %Find index of subdirectories
 
  for iDir = find(validIndex)                  %Loop over valid subdirectories
      if isempty(search_str)
          nextDir = fullfile(dirBase,subDirs{iDir});    %Get the subdirectory path
      else
          nextDir = fullfile(dirBase,subDirs{iDir},search_str);    %Get the subdirectory path
      end
      newDir=dirAll(nextDir); %Recursively call dirAll
      if isempty(newDir)
      elseif isnumeric(newDir) && newDir==0
      elseif isempty(dirFinal)
          dirFinal=newDir;
      else
          dirFinal=[dirFinal; newDir;];
      end

  end

  cd(home_dir);
  if isuniqueDir==1
      for k=1:length(dirFinal)
          dirTemp{1,k}=dirFinal(k).folder;
      end
      finalList=unique(dirTemp);
      clear dirFinal;
      for k=1:length(finalList)
          dirFinal(k).folder=finalList{1,k};
      end
  end
end