function SubjectListGen()

dataHome=uigetdir(fullfile(pwd,'*.nii'),...
    'Select Root Data Folder');
if dataHome==0
    return;
end
button=questdlg('Load or Create MRN_List.csv','','Load','Create',...
    'Cancel','Create');
        
switch button
    case 'Load'
        [fname,pname,filterindex]=uigetfile('*.csv',...
            'Select MRN_List.csv');
        if filterindex ==0
            return;
        else
            ProList=CSVImporter(fullfile(pname,fname));
        end
    case 'Create'
        PatientList=dir(dataHome);
        %             Eliminates the first two directories '.' and '..'
        PatientList=PatientList(3:end);
        PatientList=struct2table(PatientList,'AsArray',true);
        PatientList=PatientList(PatientList.isdir==1,:);
        ProList=PatientList(:,'name');
        ProList.Properties.VariableNames{'name'}='PatientID';
        clear PatientList;
        [fname,pname,filterindex]=uiputfile(fullfile(pwd,...
            'MRN_List.csv'),'Save MRN_List.csv');
        if filterindex ~=0
            writetable(ProList,fullfile(pname,fname),'Delimiter',',');
        end
        
    case 'Cancel'
        return;
end
PatientID=ProList{1,'PatientID'};
if iscell(PatientID)
    PatientID=PatientID{1};
end

fileList=dirAll(fullfile(dataHome,PatientID,'*.nii.gz'));

[fname,~,filterindex]=uigetfile(fullfile(fileList(1).folder,...
    '*.nii.gz'),'Select ROI Files','MultiSelect','on');
if filterindex ==0
    return;
end
if iscell(fname)
    for k=1:length(fname)
        ROIFile(k).fname=fname{k};
    end
else
    ROIFile.fname=fname;
end

ROIFile=LabelGUI(ROIFile);

[fname,~,filterindex]=uigetfile(fullfile(fileList(1).folder,...
    '*.nii.gz'),'Select Value Files','MultiSelect','on');
if filterindex ==0
    return;
end
if iscell(fname)
    for k=1:length(fname)
        ValueFile(k).fname=fname{k};
    end
else
    ValueFile.fname=fname;
end

ValueFile=LabelGUI(ValueFile);

counter=1;
counter2=1;
StudyDataList=StudyFolderReturn(ProList,dataHome);
for k=1:length(StudyDataList)
    for r=1:length(ROIFile)
        if exist(fullfile(dataHome,StudyDataList(k).PatientID,...
                StudyDataList(k).StudyID,'NIfTI',ROIFile(r).fname),...
                'file')==0
            MissingFiles(counter).PatientID=StudyDataList(k).PatientID;
            MissingFiles(counter).StudyID=StudyDataList(k).StudyID;
            MissingFiles(counter).NIfTIDir=fullfile(dataHome,...
                StudyDataList(k).PatientID,...
                StudyDataList(k).StudyID,'NIfTI');
            MissingFiles(counter).fname=ROIFile(r).fname;
            counter=counter+1;
        else
            for v=1:length(ValueFile)
                if exist(fullfile(dataHome,StudyDataList(k).PatientID,...
                        StudyDataList(k).StudyID,'NIfTI',...
                        ValueFile(v).fname),'file')==0
                    MissingFiles(counter).PatientID=...
                        StudyDataList(k).PatientID;
                    MissingFiles(counter).StudyID=...
                        StudyDataList(k).StudyID;
                    MissingFiles(counter).NIfTIDir=fullfile(dataHome,...
                        StudyDataList(k).PatientID,...
                        StudyDataList(k).StudyID,'NIfTI');
                    MissingFiles(counter).fname=ValueFile(v).fname;
                    counter=counter+1;
                else
                    finalList(counter2).PatientID=...
                        StudyDataList(k).PatientID;
                    finalList(counter2).StudyID=...
                        StudyDataList(k).StudyID;
                    finalList(counter2).NIfTIDir=...
                        fullfile(dataHome,StudyDataList(k).PatientID,...
                        StudyDataList(k).StudyID,'NIfTI');
                    finalList(counter2).ROIFile=ROIFile(r).fname;
                    finalList(counter2).ValueFile=ValueFile(v).fname;
                    finalList(counter2).ROI=ROIFile(r).Label;
                    finalList(counter2).Phase=ValueFile(v).Label;
                    counter2=counter2+1;
                end
            end
        end
    end
end
if exist('pname','var')==0 || isnumeric(pname)
    pname=pwd;
end
t=datetime('today','format','yyyyMMdd_HHmmss');
if exist('MissingFiles','var')==1
    MissingFiles=struct2table(MissingFiles,'AsArray',true);
    writetable(MissingFiles,fullfile(pname,...
        sprintf('MissingFiles_%s.csv',t)),'Delimiter',',');
end

finalList=struct2table(finalList,'AsArray',true);

writetable(finalList,fullfile(pname,sprintf('InputFiles_%s.csv',t)),...
    'Delimiter',',');
end

function fresult=StudyFolderReturn(MRNList,dataHome)
MRNList=table2struct(MRNList);
counter=1;
for k=1:length(MRNList)
    StudyListTemp=dir(fullfile(dataHome,MRNList(k).PatientID));
    StudyListTemp=StudyListTemp(3:end);
    StudyListTemp=struct2table(StudyListTemp,'AsArray',true);
    StudyListTemp=StudyListTemp(StudyListTemp.isdir==1,:);
    StudyListTemp=table2struct(StudyListTemp);
    for q=1:length(StudyListTemp)
        fresult(counter).PatientID=MRNList(k).PatientID;
        fresult(counter).StudyID=StudyListTemp(q).name;
        counter=counter+1;
    end
end

end

function MaxSize=FileNamesSize(FileNames)
MaxSize=0;
for k=1:length(FileNames)
    temp=length(FileNames(k).fname);
    if MaxSize<temp
        MaxSize=temp;
    end
end
end

function FileNames=LabelGUI(FileNames)
listWidth=FileNamesSize(FileNames);
if listWidth<20
    listWidth=20;
end
fieldNum=length(FileNames);
set(0,'Units','characters');
scr_size=get(0,'ScreenSize');
x_size=(2*listWidth+15);
y_size=(fieldNum+1)*2+4;

main_fig = figure('Visible','on','Name','Select Technique',...
    'Units','characters','Position',[scr_size(3)/2-(listWidth+10)/2,...
    scr_size(4)/2-(fieldNum*2+4)/2,x_size,...
    y_size],'MenuBar',...
    'none','BusyAction','cancel','NumberTitle','off');

counter=1;

for k=fieldNum:-1:1
    FileNames(k).Label='';
    fname=FileNames(counter).fname;
    [~,remain]=strtok(fname,'_');
    [~,remain,~]=fileparts(remain(2:end));
    [~,tempLabel,~]=fileparts(remain);
    
    c(k)=uicontrol('Style','text','Units','characters',...
        'position',[1 (k-1)*2+2.5 listWidth+7 2],'Min',0,...
        'String',fname);
    u(k)= uicontrol('Style','edit','String',...
        tempLabel,'Units','characters',...
        'pos',[listWidth+9 (k-1)*2+3 listWidth+7 2],'parent',main_fig,...
        'HandleVisibility','off','Tag',num2str(counter),'HorizontalAlignment','left');
    counter=counter+1;
    uicontrol('Style','pushbutton','String','OK','Units','characters',...
        'pos',[1 .5 x_size/2-1.5 2],'parent',main_fig,'Callback',@btn_Callback);
    uicontrol('Style','pushbutton','String','Cancel','Units','characters',...
        'pos',[x_size/2-0.5 .5 x_size/2-1.5 2],'parent',main_fig,'Callback',@btn_Cancel);
end
k=fieldNum+1;
uicontrol('Style','text','Units','characters',...
    'position',[1 (k-1)*2+2.5 listWidth+7 2],'Min',0,...
    'String','File Name');
uicontrol('Style','text','String',...
    'Label','Units','characters',...
    'pos',[listWidth+9 (k-1)*2+2.5 listWidth+7 2]);
uicontrol(u(fieldNum));

    uiwait(gcf);
    
    function btn_Callback(~,~)
        cnt=fieldNum;
        for q=1:length(c)
            FileNames(q).Label=get(u(cnt),'String');
            cnt=cnt-1;
        end
        close(main_fig);
    end
    function btn_Cancel(~,~)
        close(main_fig);
    end

end
