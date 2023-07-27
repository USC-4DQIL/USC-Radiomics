%% GUI Functions

function ImageTechnique=ImageProcessingTechniqueIBSIGUI(type)
 
ImageTechnique=[];
 %if nargin==0
 PosID=[];
 %end
 techList=GUIList(type);
 listWidth=DescriptSize(techList);
 if listWidth<20
     listWidth=20;
 end
 fieldNum=length(techList);
 set(0,'Units','characters');
 scr_size=get(0,'ScreenSize');
 x_size=(listWidth+10);
 y_size=fieldNum*2+4;
% 
 main_fig = figure('Visible','on','Name','Select Technique',...
     'Units','characters','Position',[scr_size(3)/2-(listWidth+10)/2,...
     scr_size(4)/2-(fieldNum*2+4)/2,x_size,...
     y_size],'MenuBar',...
     'none','BusyAction','cancel','NumberTitle','off');
 if nargin~=0
     set(main_fig,'Name',['Select ',PosID]);
 end
% 
 counter=fieldNum;
% 
 for k=1:fieldNum
     c(k)=uicontrol('Style','checkbox','Units','characters',...
         'position',[1 (k-1)*2+3.5 5 2],'parent',main_fig,'Min',0,...
         'Max',1,'Value',techList(counter).def);
     u(k)= uicontrol('Style','text','String',...
         sprintf('%s',techList(counter).descript),'Units','characters',...
         'pos',[6 (k-1)*2+3 listWidth+4 2],'parent',main_fig,'HandleVisibility',...
         'off','Tag',num2str(counter),'HorizontalAlignment','left');
     counter=counter-1;
     uicontrol('Style','pushbutton','String','OK','Units','characters',...
         'pos',[1 .5 x_size/2-1.5 2],'parent',main_fig,'Callback',@btn_Callback);
     uicontrol('Style','pushbutton','String','Cancel','Units','characters',...
         'pos',[x_size/2-0.5 .5 x_size/2-1.5 2],'parent',main_fig,'Callback',@btn_Cancel);
 end
% 
     uiwait(gcf);

     function btn_Callback(~,~)
         cnt=1;
         for q=1:length(c)
            ImageTechnique=[ImageTechnique; techList(str2double(get(u(q),'Tag')))];
            ImageTechnique(q).val = get(c(q),"Value");
             cnt=cnt+1;
         end
         if ~isempty(ImageTechnique)
             ImageTechnique=struct2table(ImageTechnique);
         end
         close(main_fig);
     end
     function btn_Cancel(~,~)
         ImageTechnique=[];
         close(main_fig);
     end
 
 end

function techList=GUIList(type)
 counter=1;
 switch type
     case 'All'
         techList(counter).descript='Intensity';
         techList(counter).tech='Intensity';
         techList(counter).type='All';
         techList(counter).def=1;
         techList(counter).val = 1;
         counter=counter+1;
         techList(counter).descript='Histogram';
         techList(counter).tech='GLCM';
         techList(counter).type='3D';
         techList(counter).def=1;
         techList(counter).val = 1;
         counter=counter+1;
        techList(counter).descript='GLCM';
        techList(counter).tech='GLCM';
         techList(counter).type='6';
         techList(counter).def=1; techList(counter).val = 1;
         counter=counter+1;
         techList(counter).descript='GLRLM';
         techList(counter).tech='GLRLM';
         techList(counter).type='6';
         techList(counter).def=1; techList(counter).val = 1;
         counter=counter+1;
         techList(counter).descript='GLSZM';
         techList(counter).tech='GLSZM';
         techList(counter).type='3';
         techList(counter).def=1; techList(counter).val = 1;
         counter=counter+1;
         techList(counter).descript='GLDZM';
         techList(counter).tech='GLDZM';
         techList(counter).type='3';
         techList(counter).def=1; techList(counter).val = 1;
         counter=counter+1;
         techList(counter).descript='NGTDM';
         techList(counter).tech='NGTDM';
         techList(counter).type='3';
         techList(counter).def=1; techList(counter).val = 1;
         counter=counter+1;
         techList(counter).descript='NGLDM';
         techList(counter).tech='NGLDM';
         techList(counter).type='3';
         techList(counter).def=1; techList(counter).val = 1;
         counter=counter+1;
         techList(counter).descript='GLDM';
         techList(counter).tech='GLDM';
         techList(counter).type='2';
         techList(counter).def=1; techList(counter).val = 1;

% 
     case 'IBSI'
         techList(counter).descript='Intensity';
         techList(counter).tech='Intensity';
         techList(counter).type='All';
         techList(counter).def=1; techList(counter).val = 1;
         counter=counter+1;
         techList(counter).descript='Histogram';
         techList(counter).tech='Histogram';
         techList(counter).type='3D';
         techList(counter).def=1; techList(counter).val = 1;
         counter=counter+1;
        techList(counter).descript='GLCM';
        techList(counter).tech='GLCM';
         techList(counter).type='6';
         techList(counter).def=1; techList(counter).val = 1;
         counter=counter+1;
         techList(counter).descript='GLRLM';
         techList(counter).tech='GLRLM';
         techList(counter).type='6';
         techList(counter).def=1; techList(counter).val = 1;
         counter=counter+1;
         techList(counter).descript='GLSZM';
         techList(counter).tech='GLSZM';
         techList(counter).type='3';
         techList(counter).def=1; techList(counter).val = 1;
         counter=counter+1;
         techList(counter).descript='GLDZM';
         techList(counter).tech='GLDZM';
         techList(counter).type='3';
         techList(counter).def=1; techList(counter).val = 1;
         counter=counter+1;
         techList(counter).descript='NGTDM';
         techList(counter).tech='NGTDM';
         techList(counter).type='3';
         techList(counter).def=1; techList(counter).val = 1;
         counter=counter+1;
         techList(counter).descript='NGLDM';
         techList(counter).tech='NGLDM';
         techList(counter).type='3';
         techList(counter).def=1; techList(counter).val = 1;
% 
 end
 end

function MaxSize=DescriptSize(techList)
MaxSize=0;
for k=1:length(techList)
    temp=length(techList(k).descript);
    if MaxSize<temp
        MaxSize=temp;
    end
end
end
