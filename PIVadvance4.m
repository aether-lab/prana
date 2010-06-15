function varargout = PIVadvance4(varargin)
% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @PIVadvance4_OpeningFcn, ...
                   'gui_OutputFcn',  @PIVadvance4_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT

% --- Opening function for figure / variable initialization ---
function PIVadvance4_OpeningFcn(hObject, eventdata, handles, varargin)
handles.syscolor=get(hObject,'color');
if ispc
    handles.loaddirec=[pwd '\'];
else
    handles.loaddirec=[pwd '/'];
end    

handles.data.imdirec=pwd;
handles.data.imbase='Img_';
handles.data.imzeros='6';
handles.data.imext='tif';
handles.data.imcstep='1';
handles.data.imfstep='2';
handles.data.imfstart='1';
handles.data.imfend='1';

handles.data.wrmag='1';
handles.data.wrsamp='1';
handles.data.wrsep='1';
handles.data.batchname='Proc1';
handles.data.datout='1';
handles.data.multiplematout='0';
handles.data.outdirec=pwd;

handles.data.exp_date='';
handles.data.exp_L='';
handles.data.exp_v0='';
handles.data.exp_notes={'Camera Description:' '' 'Lens Description:' '' 'Notes:' ''};
handles.data.exp_density='1000';
handles.data.exp_viscosity='1.308e-3';
handles.data.exp_surfacetension='0.07197';
handles.data.exp_partD='';
handles.data.exp_partdensity='';
handles.data.exp_wavelength='.532';
handles.data.exp_pixelsize='';
handles.data.exp_lensfocal='';
handles.data.exp_lensfnum='';
handles.data.exp_micro='0';
handles.data.exp_NA='';
handles.data.exp_n='';
handles.data.exp_Re='';
handles.data.exp_St='';
handles.data.exp_M='';
handles.data.exp_ROI='';
handles.data.exp_diffractiondiameter='';
handles.data.exp_depthoffocus='';

handles.data.masktype='none';
handles.data.staticmaskname='';
handles.data.maskdirec=pwd;
handles.data.maskbase='maskfor_Img_';
handles.data.maskzeros='6';
handles.data.maskext='tif';
handles.data.maskfstep='1';
handles.data.maskfstart='1';

handles.data.PIV0.winres='32,32';
handles.data.PIV0.winsize='64,64';
handles.data.PIV0.winauto='1';
handles.data.PIV0.gridres='8,8';
handles.data.PIV0.winoverlap='75,75';
handles.data.PIV0.gridtype='1';
handles.data.PIV0.gridbuf='8,8';
handles.data.PIV0.BWO='0,0';
handles.data.PIV0.corr='2';
handles.data.PIV0.RPCd='2.8';
handles.data.PIV0.velsmooth='0';
handles.data.PIV0.velsmoothfilt='2';
handles.data.PIV0.val='0';
handles.data.PIV0.uod='1';
handles.data.PIV0.bootstrap='0';
handles.data.PIV0.thresh='0';
handles.data.PIV0.uod_type='2';
handles.data.PIV0.uod_window='7,7;7,7';
handles.data.PIV0.uod_thresh='3,2';
handles.data.PIV0.bootstrap_percentsampled='15';
handles.data.PIV0.bootstrap_iterations='700';
handles.data.PIV0.bootstrap_passes='12';
handles.data.PIV0.valuthresh='-16,16';
handles.data.PIV0.valvthresh='-16,16';
handles.data.PIV0.valextrapeaks='1';
handles.data.PIV0.savepeakinfo='0';
handles.data.PIV0.corrpeaknum='1';
handles.data.PIV0.savepeakmag='0';
handles.data.PIV0.savepeakvel='0';
handles.data.PIV0.outbase='PIV_';
handles.data.PIV0.write='1';

handles.data.PIV1=handles.data.PIV0;
handles.data.PIV2=handles.data.PIV0;

handles.data.passes='2';
handles.data.cpass=num2str(get(handles.passlist,'Value'));
handles.data.method='1';
handles.data.velinterp='3';
handles.data.iminterp='1';
handles.data.framestep='3';
handles.data.PIVerror='0.1';

handles.data0=handles.data;
handles.Njob=num2str(size(get(handles.joblist,'String'),1));
handles.Cjob=num2str(get(handles.joblist,'String'));
handles=rmfield(handles,'data');
handles=update_data(handles);
handles.output = hObject;
guidata(hObject, handles);

% --- Outputs from this function are returned to the command line ---
function varargout = PIVadvance4_OutputFcn(hObject, eventdata, handles) 
varargout{1} = handles.output;

% --- Job Menu ---
function jobmenu_Callback(hObject, eventdata, handles)
if str2double(handles.Njob)>0
    set(handles.jobmenu_save,'Enable','on')
    set(handles.jobmenu_copy,'Enable','on')
    set(handles.jobmenu_delete,'Enable','on')
else
    set(handles.jobmenu_save,'Enable','off')
    set(handles.jobmenu_copy,'Enable','off')
    set(handles.jobmenu_delete,'Enable','off')
end

% --- Job Menu -> New Job ---
function jobmenu_new_Callback(hObject, eventdata, handles)
newjobbutton_Callback(hObject, eventdata, handles)

% --- Job Menu -> Load Job ---
function jobmenu_load_Callback(hObject, eventdata, handles)
loadjobbutton_Callback(hObject, eventdata, handles)

% --- Job Menu -> Save Job ---
function jobmenu_save_Callback(hObject, eventdata, handles)
savejobbutton_Callback(hObject, eventdata, handles)

% --- Job Menu -> Copy Job ---
function jobmenu_copy_Callback(hObject, eventdata, handles)
copyjobbutton_Callback(hObject, eventdata, handles)

% --- Job Menu -> Delete Job ---
function jobmenu_delete_Callback(hObject, eventdata, handles)
deletejobbutton_Callback(hObject, eventdata, handles)

% --- Execute Menu ---
function executemenu_Callback(hObject, eventdata, handles)
if str2double(handles.Njob)>0
    set(handles.executemenu_currentjob,'Enable','on')
    set(handles.executemenu_alljobs,'Enable','on')
else
    set(handles.executemenu_currentjob,'Enable','off')
    set(handles.executemenu_alljobs,'Enable','off')
end

% --- Execute Menu -> Run Current Job ---
function executemenu_currentjob_Callback(hObject, eventdata, handles)
runcurrent_Callback(hObject, eventdata, handles)

% --- Execute Menu -> Run All Jobs ---
function executemenu_alljobs_Callback(hObject, eventdata, handles)
runall_Callback(hObject, eventdata, handles)

% --- Help Menu ---
function helpmenu_Callback(hObject, eventdata, handles)

% --- Help Menu -> About ---
function helpmenu_about_Callback(hObject, eventdata, handles)
msgbox({'PIVAdvance4 v1.0','','Modified by B.Drew on 5/14/10','','[License / Copyright]'},'About')

% --- Help Menu -> Help Topics ---
function helpmenu_helptopics_Callback(hObject, eventdata, handles)
PIVhelp

% --- Experiment Parameters Tab ---
function outputtoggle_Callback(hObject, eventdata, handles)
set(handles.imagetoggle,'Value',0)
set(handles.processingtoggle,'Value',0)
set(handles.imagepanel,'Visible','off')
set(handles.processingpanel,'Visible','off')
set(handles.outputpanel,'Visible','on')

% --- Image and Data I / O Tab ---
function imagetoggle_Callback(hObject, eventdata, handles)
set(handles.outputtoggle,'Value',0)
set(handles.processingtoggle,'Value',0)
set(handles.outputpanel,'Visible','off')
set(handles.processingpanel,'Visible','off')
set(handles.imagepanel,'Visible','on')

% --- PIV Processing Tab ---
function processingtoggle_Callback(hObject, eventdata, handles)
set(handles.imagetoggle,'Value',0)
set(handles.outputtoggle,'Value',0)
set(handles.imagepanel,'Visible','off')
set(handles.outputpanel,'Visible','off')
set(handles.processingpanel,'Visible','on')

% --- Job List ---
function joblist_Callback(hObject, eventdata, handles)
if str2double(handles.Njob)>0
    Jlist=char(get(handles.joblist,'String'));
    eval(['handles.' Jlist(str2double(handles.Cjob),:) '=handles.data;']);
    handles.Cjob=num2str(get(hObject,'Value'));
    handles=update_data(handles);
    guidata(hObject,handles)
end
function joblist_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Run Current Job Button ---
function runcurrent_Callback(hObject, eventdata, handles)
if str2double(handles.Njob)>0
    Data=handles.data;
    write_expsummary(Data,handles);
    PIVadvance4code(Data);
end

% --- Run All Jobs Button ---
function runall_Callback(hObject, eventdata, handles)
if str2double(handles.Njob)>0
    Jlist=char(get(handles.joblist,'String'));
    eval(['handles.' Jlist(str2double(handles.Cjob),:) '=handles.data;']);
    
%     matlabpool open local 8
    for e=1:size(Jlist,1)
        Data=eval(['handles.' Jlist(e,:)]);
        write_expsummary(Data,handles);
        PIVadvance4code(Data);
    end
%     matlabpool close
end

% --- New Job Button ---
function newjobbutton_Callback(hObject, eventdata, handles)
if str2double(handles.Njob)>0
    Jlist=char(get(handles.joblist,'String'));
    eval(['handles.' Jlist(str2double(handles.Cjob),:) '=handles.data;']);
end
Data=handles.data0;
vn=0;
Data.batchname=char(inputdlg('Job Name?                  ','NEW JOB',1,{['PIV' num2str(str2double(handles.Njob)+1)]}));
if isempty(Data.batchname)
    vn=-1;
end
while vn==0
    if isfield(handles,Data.batchname)
        Data.batchname=char(inputdlg('Job already exists, rename?','NEW JOB',1,{['PIV' num2str(str2double(handles.Njob)+1)]}));
        if isempty(Data.batchname)
            vn=-1;
        end
    else
        handles=setfield(handles,Data.batchname,Data);
        vn=1;
    end
end
if vn~=-1
    if str2double(handles.Njob)>0
        Jlist=char(get(handles.joblist,'String'));
        Jlist={Jlist;Data.batchname};
    else
        Jlist={Data.batchname};
    end
    handles.Njob=num2str(str2double(handles.Njob)+1);
    handles.Cjob=handles.Njob;
    set(handles.joblist,'String',Jlist,'Value',str2double(handles.Cjob));
    handles=update_data(handles);
    guidata(hObject,handles)
end

% --- Load Job Button ---
function loadjobbutton_Callback(hObject, eventdata, handles)
if str2double(handles.Njob)>0
    Jlist=char(get(handles.joblist,'String'));
    eval(['handles.' Jlist(str2double(handles.Cjob),:) '=handles.data;']);
end
handlesP=handles;
[f,handles.loaddirec]=uigetfile('*.mat','LOAD JOB',handles.loaddirec,'multiselect','on');
if ischar(f)==1
    f={f};
else
    f=sort(f);
end
if isnumeric(f)==0
    try
        for pp=1:length(f)
            load([handles.loaddirec char(f(pp))]);
            if exist('Data')~=0
                vn=0;
                while vn==0
                    if isfield(handles,Data.batchname)
                        Data.batchname=char(inputdlg('Job already exists, rename?','LOAD JOB',1,{Data.batchname}));
                        if isempty(Data.batchname)
                            vn=-1;
                        end
                    else
                        handles=setfield(handles,Data.batchname,Data);
                        vn=1;
                    end
                end
                if vn~=-1
                    if str2double(handles.Njob)>0
                        Jlist=char(get(handles.joblist,'String'));
                        Jlist={Jlist;Data.batchname};
                    else
                        Jlist={Data.batchname};
                    end
                    handles.Njob=num2str(str2double(handles.Njob)+1);
                    handles.Cjob=handles.Njob;
                    set(handles.joblist,'String',Jlist,'Value',str2double(handles.Cjob));

                    handles=update_data(handles);
                    guidata(hObject,handles)
                end
            end
        end
    catch
        handles=handlesP;
        handles=update_data(handles);
        guidata(hObject,handles)
        beep,disp('ERROR: Invalid Data Format, Reloading Previous Settings')
    end
end

% --- Save Job Button ---
function savejobbutton_Callback(hObject, eventdata, handles)
if str2double(handles.Njob)>0
    Data=handles.data;
    uisave('Data',[handles.loaddirec handles.data.batchname '.mat']);
end

% --- Copy Job Button ---
function copyjobbutton_Callback(hObject, eventdata, handles)
if str2double(handles.Njob)>0
    Jlist=char(get(handles.joblist,'String'));
    eval(['handles.' Jlist(str2double(handles.Cjob),:) '=handles.data;']);
end
Data=handles.data;
vn=0;
Data.batchname=char(inputdlg('Job Name?                  ','NEW JOB',1,{strtrim(Jlist(str2double(handles.Cjob),:))}));
if isempty(Data.batchname)
    vn=-1;
end
while vn==0
    if isfield(handles,Data.batchname)
        Data.batchname=char(inputdlg('Job already exists, rename?','NEW JOB',1,{strtrim(Jlist(str2double(handles.Cjob),:))}));
        if isempty(Data.batchname)
            vn=-1;
        end
    else
        handles=setfield(handles,Data.batchname,Data);
        vn=1;
    end
end
if vn~=-1
    if str2double(handles.Njob)>0
        Jlist=char(get(handles.joblist,'String'));
        Jlist={Jlist;Data.batchname};
    else
        Jlist={Data.batchname};
    end
    handles.Njob=num2str(str2double(handles.Njob)+1);
    handles.Cjob=handles.Njob;
    set(handles.joblist,'String',Jlist,'Value',str2double(handles.Cjob));
    handles=update_data(handles);
    guidata(hObject,handles)
end

% --- Delete Job Button ---
function deletejobbutton_Callback(hObject, eventdata, handles)
if str2double(handles.Njob)>0
    del=questdlg('Are You Sure?','Delete Job File','yes','no','yes');
    if strcmp(del,'yes')==1
        Jlist=char(get(handles.joblist,'String'));
        handles=rmfield(handles,Jlist(str2double(handles.Cjob),:));
        id=setdiff(1:str2double(handles.Njob),str2double(handles.Cjob));
        Jlist=Jlist(id,:);
        set(handles.joblist,'String',Jlist,'Value',1);
        handles.Cjob='1';
        handles.Njob=num2str(str2double(handles.Njob)-1);
        handles=update_data(handles);
        guidata(hObject,handles)
    end
end

% --- Job Name Text Box ---
function currentjobname_Callback(hObject, eventdata, handles)
if str2double(handles.Njob)>0
    Jlist0=char(get(handles.joblist,'String'));
    handles=rmfield(handles,handles.data.batchname);
    handles.data.batchname = get(hObject,'String');
    Jlist=cell(str2double(handles.Njob),1);
    for e=1:str2double(handles.Njob)
        if e==str2double(handles.Cjob)
            Jlist(e)={handles.data.batchname};
        else
            Jlist(e)={Jlist0(e,:)};
        end
    end    
    handles=setfield(handles,handles.data.batchname,handles.data);
    set(handles.joblist,'String',char(Jlist));
    guidata(hObject,handles)
end
function currentjobname_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Load First Image Button ---
function loadfirstimage_Callback(hObject, eventdata, handles)
if str2double(handles.Njob)>0
    D = handles.data.imdirec;
    B = handles.data.imbase;
    Z = handles.data.imzeros;
    E = handles.data.imext;
    S = handles.data.imfstart;
    F = handles.data.imfend;
    [imname,handles.data.imdirec] = uigetfile('*.*','Select an image file',handles.data.imdirec);
    if handles.data.imdirec==0
        handles.data.imdirec = D;
        handles.data.imbase = B;
        handles.data.imzeros = Z;
        handles.data.imext = E;
        handles.data.imfstart = S;
        handles.data.imfend = F;
    else
        i=strfind(imname,'.');
        handles.data.imext=imname((i+1):end);
        fstart=0;zeros=0;
        while ~isnan(fstart)
            zeros=zeros+1;
            fstart=str2double(imname((i-zeros):(i-1)));
        end
        handles.data.imzeros=num2str(zeros-1);
        fstart=str2double(imname((i-(zeros-1)):(i-1)));
        handles.data.imfstart=num2str(fstart);
        handles.data.imfend=num2str(fstart);
        handles.data.imbase=imname(1:(i-zeros));
    end
    set(handles.imagedirectory,'string',handles.data.imdirec);
    set(handles.imagebasename,'string',handles.data.imbase);
    set(handles.imagezeros,'string',handles.data.imzeros);
    set(handles.imageextension,'string',handles.data.imext);
    set(handles.imageframestart,'string',handles.data.imfstart);
    set(handles.imageframeend,'string',handles.data.imfend);
    if strcmp(handles.data.masktype,'dynamic')
        load_masklist(handles)
    end
    load_imlist(handles);
    guidata(hObject,handles)
end

% --- Image Directory Text Box ---
function imagedirectory_Callback(hObject, eventdata, handles)
if str2double(handles.Njob)>0
    handles.data.imdirec = get(hObject,'String');
    if strcmp(handles.data.masktype,'dynamic')
        load_masklist(handles)
    end
    load_imlist(handles);
    guidata(hObject,handles)
end
function imagedirectory_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Load Image Directory Button ---
function loadimagedirectorybutton_Callback(hObject, eventdata, handles)
if str2double(handles.Njob)>0
    D = handles.data.imdirec;
    handles.data.imdirec = uigetdir(handles.data.imdirec);
    if handles.data.imdirec==0
        handles.data.imdirec = D;
    end
    set(handles.imagedirectory,'string',handles.data.imdirec);
    if strcmp(handles.data.masktype,'dynamic')
        load_masklist(handles)
    end
    load_imlist(handles);
    guidata(hObject,handles)
end

% --- Image Basename Text Box ---
function imagebasename_Callback(hObject, eventdata, handles)
if str2double(handles.Njob)>0
    handles.data.imbase = get(hObject,'String');
    if strcmp(handles.data.masktype,'dynamic')
        load_masklist(handles)
    end
    load_imlist(handles);
    guidata(hObject,handles)
end
function imagebasename_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Image Zeros Text Box ---
function imagezeros_Callback(hObject, eventdata, handles)
if str2double(handles.Njob)>0
    handles.data.imzeros = get(hObject,'String');
    if strcmp(handles.data.masktype,'dynamic')
        load_masklist(handles)
    end
    load_imlist(handles);
    guidata(hObject,handles)
end
function imagezeros_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Image File Extension Text Box ---
function imageextension_Callback(hObject, eventdata, handles)
if str2double(handles.Njob)>0
    handles.data.imext = get(hObject,'String');
    if strcmp(handles.data.imext(1),'.')
        set(hObject,'backgroundcolor','r');
    else
        set(hObject,'backgroundcolor',[1 1 1]);
    end
    if strcmp(handles.data.masktype,'dynamic')
        load_masklist(handles)
    end
    load_imlist(handles);
    guidata(hObject,handles)
end
function imageextension_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Image Correlation Step Text Box ---
function imagecorrelationstep_Callback(hObject, eventdata, handles)
if str2double(handles.Njob)>0
    handles.data.imcstep = get(hObject,'String');
    if strcmp(handles.data.masktype,'dynamic')
        load_masklist(handles)
    end
    load_imlist(handles);
    guidata(hObject,handles)
end
function imagecorrelationstep_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Image Frame Start Text Box ---
function imageframestart_Callback(hObject, eventdata, handles)
if str2double(handles.Njob)>0
    handles.data.imfstart = get(hObject,'String');
    if strcmp(handles.data.masktype,'dynamic')
        load_masklist(handles)
    end
    load_imlist(handles);
    guidata(hObject,handles)
end
function imageframestart_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Image Frame Step Text Box ---
function imageframestep_Callback(hObject, eventdata, handles)
if str2double(handles.Njob)>0
    handles.data.imfstep = get(hObject,'String');
    if strcmp(handles.data.masktype,'dynamic')
        load_masklist(handles)
    end
    load_imlist(handles);
    guidata(hObject,handles)
end
function imageframestep_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Image Frame End Text Box ---
function imageframeend_Callback(hObject, eventdata, handles)
if str2double(handles.Njob)>0
    handles.data.imfend = get(hObject,'String');
    if strcmp(handles.data.masktype,'dynamic')
        load_masklist(handles)
    end
    load_imlist(handles);
    guidata(hObject,handles)
end
function imageframeend_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Dynamic Mask Directory Text Box ---
function maskdirectory_Callback(hObject, eventdata, handles)
if str2double(handles.Njob)>0 && get(handles.dynamicmaskbutton,'Value')==1
    handles.data.maskdirec = get(hObject,'String');
    load_masklist(handles);
    guidata(hObject,handles)
end
function maskdirectory_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Load Dynamic Mask Directory Button ---
function loadmaskdirectorybutton_Callback(hObject, eventdata, handles)
if str2double(handles.Njob)>0 && get(handles.dynamicmaskbutton,'Value')==1
    D = handles.data.maskdirec;
    handles.data.maskdirec = uigetdir(handles.data.maskdirec);
    if handles.data.maskdirec==0
        handles.data.maskdirec = D;
    end
    set(handles.maskdirectory,'string',handles.data.maskdirec);
    load_masklist(handles);
    guidata(hObject,handles)
end
function loadmaskdirectorybutton_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Dynamic Mask Basename Text Box ---
function maskbasename_Callback(hObject, eventdata, handles)
if str2double(handles.Njob)>0 && get(handles.dynamicmaskbutton,'Value')==1
    handles.data.maskbase = get(hObject,'String');
    load_masklist(handles);
    guidata(hObject,handles)
end
function maskbasename_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Dynamic Mask Zeros Text Box ---
function maskzeros_Callback(hObject, eventdata, handles)
if str2double(handles.Njob)>0 && get(handles.dynamicmaskbutton,'Value')==1
    handles.data.maskzeros = get(hObject,'String');
    load_masklist(handles);
    guidata(hObject,handles)
end
function maskzeros_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Dynamic Mask File Extension Text Box ---
function maskextension_Callback(hObject, eventdata, handles)
if str2double(handles.Njob)>0 && get(handles.dynamicmaskbutton,'Value')==1
    handles.data.maskext = get(hObject,'String');
    if strcmp(handles.data.maskext(1),'.')
        set(hObject,'backgroundcolor','r');
    else
        set(hObject,'backgroundcolor',[1 1 1]);
    end
    load_masklist(handles);
    guidata(hObject,handles)
end
function maskextension_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Dynamic Mask Frame Start Text Box ---
function maskframestart_Callback(hObject, eventdata, handles)
if str2double(handles.Njob)>0 && get(handles.dynamicmaskbutton,'Value')==1
    handles.data.maskfstart = get(hObject,'String');
    load_masklist(handles);
    guidata(hObject,handles)
end
function maskframestart_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Dynamic Mask Frame Step Text Box ---
function maskframestep_Callback(hObject, eventdata, handles)
if str2double(handles.Njob)>0 && get(handles.dynamicmaskbutton,'Value')==1
    handles.data.maskfstep = get(hObject,'String');
    load_masklist(handles);
    guidata(hObject,handles)
end
function maskframestep_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Image Processing Order ---
function imagelist_Callback(hObject, eventdata, handles)
function imagelist_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Mask Processing Order ---
function masklist_Callback(hObject, eventdata, handles)
function masklist_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- No Processing Mask Radio Button ---
function nomaskbutton_Callback(hObject, eventdata, handles)
handles.data.masktype='none';
set(hObject,'Value',1)
set(handles.staticmaskbutton,'Value',0)
set(handles.dynamicmaskbutton,'Value',0)
if str2double(handles.Njob)>0
    Jlist=char(get(handles.joblist,'String'));
    eval(['handles.' Jlist(str2double(handles.Cjob),:) '=handles.data;']);
    handles=update_data(handles);
    guidata(hObject,handles)
end

% --- Static Processing Mask Radio Button ---
function staticmaskbutton_Callback(hObject, eventdata, handles)
handles.data.masktype='static';
set(hObject,'Value',1)
set(handles.dynamicmaskbutton,'Value',0)
set(handles.nomaskbutton,'Value',0)
if str2double(handles.Njob)>0
    Jlist=char(get(handles.joblist,'String'));
    eval(['handles.' Jlist(str2double(handles.Cjob),:) '=handles.data;']);
    handles=update_data(handles);
    guidata(hObject,handles)
end

% --- Dynamic Masking Radio Button ---
function dynamicmaskbutton_Callback(hObject, eventdata, handles)
handles.data.masktype='dynamic';
set(hObject,'Value',1)
set(handles.nomaskbutton,'Value',0)
set(handles.staticmaskbutton,'Value',0)
if str2double(handles.Njob)>0
    Jlist=char(get(handles.joblist,'String'));
    eval(['handles.' Jlist(str2double(handles.Cjob),:) '=handles.data;']);
    if get(handles.passtype,'Value')==4 && get(hObject,'Value')==1
        errordlg('Dynamic Masking is not compatible with the Ensemble correlation.','Warning')
    end
    handles=update_data(handles);
    guidata(hObject,handles)
end

% --- Static Mask File Text Box ---
function staticmaskfile_Callback(hObject, eventdata, handles)
if str2double(handles.Njob)>0 && get(handles.staticmaskbutton,'Value')==1
    handles.data.staticmaskname = get(hObject,'String');
    guidata(hObject,handles)
end
function staticmaskfile_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Load Static Mask File Button ---
function loadstaticmaskfile_Callback(hObject, eventdata, handles)
if str2double(handles.Njob)>0 && get(handles.staticmaskbutton,'Value')==1
    D = handles.data.staticmaskname;
    [a,b] = uigetfile('*.*');
    handles.data.staticmaskname = [b a];
    if handles.data.staticmaskname==0
        handles.data.staticmaskname = D;
    end
    guidata(hObject,handles)
    set(handles.staticmaskfile,'string',handles.data.staticmaskname);
end

% --- Preview Image + Mask Button ---
function impreview_Callback(hObject, eventdata, handles)
if str2double(handles.Njob)>0
    try
        im1=double(imread([handles.data.imdirec '\' handles.data.imbase sprintf(['%0.' handles.data.imzeros 'i.' handles.data.imext],str2double(handles.data.imfstart))]));
        im1=flipud(im1)/255;
        try
            if strcmp(handles.data.masktype,'static')
                mask = double(imread(handles.data.staticmaskname));
                mask = flipud(mask);
            elseif strcmp(handles.data.masktype,'dynamic')
                mask = double(imread([handles.data.maskdirec '\' handles.data.maskbase sprintf(['%0.' handles.data.maskzeros 'i.' handles.data.maskext],str2double(handles.data.maskfstart))]));
                mask = flipud(mask);
            else
                mask = ones(size(im1));
            end
            try
                im1(mask==0)=0.5*im1(mask==0);
                e=0;
            catch
                msgbox('Mask / Image Not Compatible');
                e=-1;
            end
        catch
            msgbox('Mask Not Found');
            e=-1;
        end
    catch
        msgbox('Image Frame Not Found');
        e=-1;
    end
    
    if e==0
        figure,hold on

        L=size(im1);
        Ps=get(0,'screensize');
        if (L(1)/Ps(4))>(L(2)/Ps(3))
            dy=0.8*Ps(4);
            dx=dy*(L(2)/L(1));
        else
            dx=0.8*Ps(3);
            dy=dx*(L(1)/L(2));
        end
        if gcf~=gcbf
            set(gcf,'position',[(Ps(3)-dx)/2 (Ps(4)-dy)/2 dx dy])
            set(gcf,'color',0.5*[1 1 1])
            imagesc(im1,[0 1]),axis image,colormap gray,axis off,set(gca,'position',[0 0 1 1])
        end

        A=get(handles.gridbuffer,'String');
        G=[str2double(A(1:(strfind(A,',')-1))) str2double(A((strfind(A,',')+1):end))];
        A=get(handles.gridres,'String');
        S=[str2double(A(1:(strfind(A,',')-1))) str2double(A((strfind(A,',')+1):end))];

        S=[S(2) S(1)];
        G=[G(2) G(1) L(1)-G(2)+1 L(2)-G(1)+1];

        %form grid
        if max(S)==0
            y=(1:L(1))';
            x=1:L(2);
        else
            if G(1)==0 || G(3)<=G(1) || G(3)>L(1)
                y=(ceil((L(1)-(floor(L(1)/S(1))-2)*S(1))/2):S(1):(L(1)-S(1)))';
            else
                y=(G(1):S(1):G(3))';
            end
            if G(2)==0 || G(4)<=G(2) || G(4)>L(2)
                x=ceil((L(2)-(floor(L(2)/S(2))-2)*S(2))/2):S(2):(L(2)-S(2));
            else
                x=(G(2):S(2):G(4));
            end
        end

        %vector2matrix conversion
        X=x(ones(length(y),1),:);
        Y=y(:,ones(1,length(x)));
        
        Y=Y(:);
        X=X(:);

        Eval=zeros(size(X));
        for e=1:length(X)
            if mask(Y(e),X(e))>0
                Eval(e)=1;
            end
        end
        X=X(Eval==1);
        Y=Y(Eval==1);
        if gcf~=gcbf
            plot(X,Y,'r.')
        end
    end
end

% --- Output Directory Text Box ---
function outputdirectory_Callback(hObject, eventdata, handles)
if str2double(handles.Njob)>0
    handles.data.outdirec = get(hObject,'String');
    guidata(hObject,handles)
end
function outputdirectory_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Load Ouput Directory Button ---
function loadoutputdirectorybutton_Callback(hObject, eventdata, handles)
if str2double(handles.Njob)>0
    D = handles.data.outdirec;
    handles.data.outdirec = uigetdir(handles.data.outdirec);
    if handles.data.outdirec==0
        handles.data.outdirec = D;
    end
    set(handles.outputdirectory,'string',handles.data.outdirec);
    guidata(hObject,handles)
end

% --- .dat Filetype Checkbox ---
function datcheckbox_Callback(hObject, eventdata, handles)
if get(hObject,'Value')==1 || get(handles.multiplematcheckbox,'Value')==0
    set(hObject,'Value',1);
    handles.data.datout='1';
    if str2double(handles.Njob)>0
        Jlist=char(get(handles.joblist,'String'));
        eval(['handles.' Jlist(str2double(handles.Cjob),:) '=handles.data;']);
        guidata(hObject,handles)
    end
else    
    handles.data.datout='0';
    if str2double(handles.Njob)>0
        Jlist=char(get(handles.joblist,'String'));
        eval(['handles.' Jlist(str2double(handles.Cjob),:) '=handles.data;']);
        guidata(hObject,handles)
    end
end

% --- .mat Filetype Checkbox ---
function multiplematcheckbox_Callback(hObject, eventdata, handles)
if get(hObject,'Value')==1 || get(handles.datcheckbox,'Value')==0
    set(hObject,'Value',1);
    handles.data.multiplematout='1';
    if str2double(handles.Njob)>0
        Jlist=char(get(handles.joblist,'String'));
        eval(['handles.' Jlist(str2double(handles.Cjob),:) '=handles.data;']);
        guidata(hObject,handles)
    end
else    
    handles.data.multiplematout='0';
    if str2double(handles.Njob)>0
        Jlist=char(get(handles.joblist,'String'));
        eval(['handles.' Jlist(str2double(handles.Cjob),:) '=handles.data;']);
        guidata(hObject,handles)
    end
end

% --- List of Passes ---
function passlist_Callback(hObject, eventdata, handles)
if str2double(handles.Njob)>0
    handles=set_PIVcontrols(handles);
    guidata(hObject,handles)
end
function passlist_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Add Pass Button ---
function addpassbutton_Callback(hObject, eventdata, handles)
if str2double(handles.Njob)>0
    N=str2double(handles.data.passes);
    eval(['handles.data=setfield(handles.data,''PIV' num2str(N+1) ''',handles.data.PIV0);']);
    eval(['handles.data.PIV' num2str(N+1) '.outbase=[''Pass'' num2str(N+1) ''_''];']);
    handles.data.passes=num2str(N+1);
    load_PIVlist(handles);
    handles=set_PIVcontrols(handles);
    guidata(hObject,handles)
end

% --- Delete Pass Button ---
function deletepassbutton_Callback(hObject, eventdata, handles)
if str2double(handles.Njob)>0
    if str2double(handles.data.passes)>1
        p=get(handles.passlist,'Value');
        eval(['handles.data=rmfield(handles.data,''PIV' num2str(p) ''');']);
        for e=(p+1):str2double(handles.data.passes)
            eval(['handles.data=setfield(handles.data,''PIV' num2str(e-1) ''',handles.data.PIV' num2str(e) ');']);
            eval(['handles.data.PIV' num2str(e-1) '.outbase=[''Pass'' num2str(e-1) ''_''];']);
            if e==str2double(handles.data.passes)
                eval(['handles.data=rmfield(handles.data,''PIV' num2str(e) ''');']);
            end
        end
        handles.data.passes=num2str(str2double(handles.data.passes)-1);
        set(handles.passlist,'Value',1);
        load_PIVlist(handles)
        guidata(hObject,handles)
    end
    handles=set_PIVcontrols(handles);
    guidata(hObject,handles)
end

% --- Velocity Interpolation Function Drop-Down Menu ---
function velocityinterptype_Callback(hObject, eventdata, handles)
if str2double(handles.Njob)>0
    handles.data.velinterp=num2str(get(hObject,'Value'));
    guidata(hObject,handles)
end
function velocityinterptype_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Image Interpolation Function Drop-Down Menu ---
function imageinterptype_Callback(hObject, eventdata, handles)
if str2double(handles.Njob)>0
    handles.data.iminterp=num2str(get(hObject,'Value'));
    guidata(hObject,handles)
end
function imageinterptype_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Maximum Framestep Textbox ---
function framestep_Callback(hObject, eventdata, handles)
if str2double(handles.Njob)>0
    handles.data.framestep=get(hObject,'String');
    guidata(hObject,handles)
end
function framestep_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function PIVerror_Callback(hObject, eventdata, handles)
if str2double(handles.Njob)>0
    handles.data.PIVerror=get(hObject,'String');
    guidata(hObject,handles)
end
function PIVerror_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Pass Processing Type Drop-Down Menu ---
function passtype_Callback(hObject, eventdata, handles)
if str2double(handles.Njob)>0
    handles.data.method=num2str(get(hObject,'Value'));
    if get(hObject,'Value')==1
        set(handles.velocityinterptype,'backgroundcolor',0.5*[1 1 1]);
        set(handles.imageinterptype,'backgroundcolor',0.5*[1 1 1]);
        set(handles.smoothingsize,'backgroundcolor',0.5*[1 1 1]);
        set(handles.framestep,'backgroundcolor',0.5*[1 1 1]);
        set(handles.PIVerror,'backgroundcolor',0.5*[1 1 1]);
        for e=1:str2double(handles.data.passes)
            eval(['handles.data.PIV' num2str(e) '.gridres = get(handles.gridres,''String'');'])
            eval(['handles.data.PIV' num2str(e) '.gridbuf = get(handles.gridbuffer,''String'');'])
        end
    elseif get(hObject,'Value')>=5
        set(handles.velocityinterptype,'backgroundcolor',0.5*[1 1 1]);
        set(handles.imageinterptype,'backgroundcolor',0.5*[1 1 1]);
        set(handles.smoothingsize,'backgroundcolor',0.5*[1 1 1]);
        set(handles.framestep,'backgroundcolor',[1 1 1]);
        set(handles.PIVerror,'backgroundcolor',[1 1 1]);
    else
        set(handles.velocityinterptype,'backgroundcolor',[1 1 1]);
        set(handles.framestep,'backgroundcolor',0.5*[1 1 1]);
        set(handles.PIVerror,'backgroundcolor',0.5*[1 1 1]);
        if get(hObject,'Value')==3
            set(handles.imageinterptype,'backgroundcolor',[1 1 1]);
        else
            set(handles.imageinterptype,'backgroundcolor',0.5*[1 1 1]);
        end
        if get(handles.smoothingcheckbox,'Value')==1
            set(handles.smoothingsize,'backgroundcolor',[1 1 1]);
        else
            set(handles.smoothingsize,'backgroundcolor',0.5*[1 1 1]);
        end 
    end

    if get(hObject,'Value')==4 && strcmp(handles.data.masktype,'dynamic')
        errordlg('Dynamic Masking is not compatible with the Ensemble correlation.','Warning')
    end
    handles=set_PIVcontrols(handles);
    load_data(handles);
    guidata(hObject,handles)
end
function passtype_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- ? Button Next to Algorithm Drop-Down Menu ---
function algorithmhelp_Callback(hObject, eventdata, handles)
PIVhelp(5)

% --- Window Resolution Text Box ---
function windowres_Callback(hObject, eventdata, handles)
if str2double(handles.Njob)>0
    eval(['handles.data.PIV' handles.data.cpass '.winres = get(hObject,''String'');'])
    a=get(hObject,'String');
    wx=str2double(a(1:(strfind(a,',')-1)));
    wy=str2double(a((strfind(a,',')+1):end));
    if get(handles.autowinsizecheckbox,'Value')==1
        xbin = 2^(ceil(log(2*wx)/log(2)));
        ybin = 2^(ceil(log(2*wy)/log(2)));
        eval(['handles.data.PIV' handles.data.cpass '.winsize = [num2str(xbin) '','' num2str(ybin)];'])
        eval(['set(handles.windowsize,''String'',handles.data.PIV' handles.data.cpass '.winsize)'])
    end
    if wx*wy<256
        set(hObject,'backgroundcolor',[1 0.5 0]);
    else
        set(hObject,'backgroundcolor',[1 1 1]);
    end
    if get(handles.setgridresbutton,'Value')==1
        A=get(handles.gridres,'String');
        gx=str2double(A(1:(strfind(A,',')-1)));
        gy=str2double(A((strfind(A,',')+1):end));
        overX=(wx-gx)/wx*100;
        overY=(wy-gy)/wy*100;
        eval(['handles.data.PIV' handles.data.cpass '.winoverlap = [num2str(overX),'','',num2str(overY)];'])
    else
        A=get(handles.winoverlap,'String');
        overX=str2double(A(1:(strfind(A,',')-1)));
        overY=str2double(A((strfind(A,',')+1):end));
        gx=round(wx*(1-overX/100));
        gy=round(wy*(1-overY/100));
        eval(['handles.data.PIV' handles.data.cpass '.gridres = [num2str(gx),'','',num2str(gy)];'])
    end
    handles=set_PIVcontrols(handles);
    guidata(hObject,handles)
end
function windowres_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Window Size Text Box ---
function windowsize_Callback(hObject, eventdata, handles)
if str2double(handles.Njob)>0
    if get(handles.autowinsizecheckbox,'Value')==0
        A=get(hObject,'String');
        Wx=str2double(A(1:(strfind(A,',')-1)));
        Wy=str2double(A((strfind(A,',')+1):end));
        B=get(handles.windowres,'String');
        Rx=str2double(B(1:(strfind(B,',')-1)));
        Ry=str2double(B((strfind(B,',')+1):end));
        if Wx<Rx
            Wx=Rx;
            set(hObject,'String',[num2str(Wx) ',' num2str(Wy)]);
        end
        if Wy<Ry
            Wy=Ry;
            set(hObject,'String',[num2str(Wx) ',' num2str(Wy)]);
        end
        eval(['handles.data.PIV' handles.data.cpass '.winsize = get(hObject,''String'');'])
        handles=set_PIVcontrols(handles);
        guidata(hObject,handles)
    else
        eval(['set(hObject,''String'',handles.data.PIV' handles.data.cpass '.winsize);'])
    end
end
function windowsize_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Auto-Size Window Check Box ---
function autowinsizecheckbox_Callback(hObject, eventdata, handles)
if str2double(handles.Njob)>0
    eval(['handles.data.PIV' handles.data.cpass '.winauto = num2str(get(hObject,''Value''));'])
    handles=set_PIVcontrols(handles);
    guidata(hObject,handles)
end

% --- Set Grid Resolution Radio Button ---
function setgridresbutton_Callback(hObject, eventdata, handles)
set(hObject,'Value',1)
set(handles.setwinoverlapbutton,'Value',0)
if str2double(handles.Njob)>0
    if str2double(handles.data.method)==1
        for e=1:str2double(handles.data.passes)
            eval(['handles.data.PIV' num2str(e) '.gridtype = ''1'';']);
        end
    else
        eval(['handles.data.PIV' handles.data.cpass '.gridtype = ''1'';']);
    end
    handles=set_PIVcontrols(handles);
    guidata(hObject,handles)
end

% --- Grid Resolution Text Box ---
function gridres_Callback(hObject, eventdata, handles)
if str2double(handles.Njob)>0
    if get(handles.setgridresbutton,'Value')==1
        A=get(hObject,'String');
        gx=str2double(A(1:(strfind(A,',')-1)));
        gy=str2double(A((strfind(A,',')+1):end));
        if str2double(handles.data.method)==1
            for e=1:str2double(handles.data.passes)
                eval(['A=handles.data.PIV', num2str(e) '.winres;'])
                wx=str2double(A(1:(strfind(A,',')-1)));
                wy=str2double(A((strfind(A,',')+1):end));
                overX=(wx-gx)/wx*100;
                overY=(wy-gy)/wy*100;
                eval(['handles.data.PIV' num2str(e) '.gridres = get(hObject,''String'');'])
                eval(['handles.data.PIV' num2str(e) '.winoverlap = [num2str(overX),'','',num2str(overY)];'])
            end
        else
            A=get(handles.windowres,'String');
            wx=str2double(A(1:(strfind(A,',')-1)));
            wy=str2double(A((strfind(A,',')+1):end));
            overX=(wx-gx)/wx*100;
            overY=(wy-gy)/wy*100;
            eval(['handles.data.PIV' handles.data.cpass '.gridres = get(hObject,''String'');'])
            eval(['handles.data.PIV' handles.data.cpass '.winoverlap = [num2str(overX),'','',num2str(overY)];'])
        end
        handles=set_PIVcontrols(handles);
        guidata(hObject,handles)
    else
        eval(['set(hObject,''String'',handles.data.PIV' handles.data.cpass '.gridres);'])
    end
end
function gridres_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Set Window Overlap Button ---
function setwinoverlapbutton_Callback(hObject, eventdata, handles)
if str2double(handles.Njob)>0 && str2double(handles.data.method)~=1
    set(hObject,'Value',1)
    set(handles.setgridresbutton,'Value',0)
    eval(['handles.data.PIV' handles.data.cpass '.gridtype = ''2'';']);
    handles=set_PIVcontrols(handles);
    guidata(hObject,handles)
else
    set(hObject,'Value',0)
end

% --- Window Overlap Text Box ---
function winoverlap_Callback(hObject, eventdata, handles)
if str2double(handles.Njob)>0
    if get(handles.setwinoverlapbutton,'Value')==1
        A=get(handles.windowres,'String');
        wx=str2double(A(1:(strfind(A,',')-1)));
        wy=str2double(A((strfind(A,',')+1):end));
        A=get(hObject,'String');
        overX=str2double(A(1:(strfind(A,',')-1)));
        overY=str2double(A((strfind(A,',')+1):end));
        gx=round(wx*(1-overX/100));
        gy=round(wy*(1-overY/100));
        eval(['handles.data.PIV' handles.data.cpass '.winoverlap = get(hObject,''String'');'])
        eval(['handles.data.PIV' handles.data.cpass '.gridres = [num2str(gx),'','',num2str(gy)];'])
        handles=set_PIVcontrols(handles);
        guidata(hObject,handles)
    else
        eval(['set(hObject,''String'',handles.data.PIV' handles.data.cpass '.winoverlap);'])
    end
end
function winoverlap_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Grid Buffer Text Box ---
function gridbuffer_Callback(hObject, eventdata, handles)
if str2double(handles.Njob)>0
    if str2double(handles.data.method)==1
        for e=1:str2double(handles.data.passes)
            eval(['handles.data.PIV' num2str(e) '.gridbuf = get(hObject,''String'');'])
        end
    else
        eval(['handles.data.PIV' handles.data.cpass '.gridbuf = get(hObject,''String'');'])
    end
    guidata(hObject,handles)
end
function gridbuffer_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Bulk Window Offset Text Box ---
function bulkwinoffset_Callback(hObject, eventdata, handles)
if str2double(handles.Njob)>0
    if str2double(handles.data.cpass)==1
        eval(['handles.data.PIV' handles.data.cpass '.BWO = get(hObject,''String'');'])
        guidata(hObject,handles)
    else
        eval(['set(hObject,''String'',handles.data.PIV' handles.data.cpass '.BWO)'])
    end
end
function bulkwinoffset_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Correlation Type Drop-Down Menu ---
function correlationtype_Callback(hObject, eventdata, handles)
if str2double(handles.Njob)>0
    eval(['handles.data.PIV' handles.data.cpass '.corr = num2str(get(hObject,''Value''));'])
    guidata(hObject,handles)
    N=handles.data.cpass;
    A=eval(['handles.data.PIV' num2str(N)]);
    if str2double(A.corr)==1
        set(handles.rpcdiameter,'backgroundcolor',0.5*[1 1 1]);
    else
        set(handles.rpcdiameter,'backgroundcolor',[1 1 1]);
    end
end
function correlationtype_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- RPC Diameter Text Box ---
function rpcdiameter_Callback(hObject, eventdata, handles)
if str2double(handles.Njob)>0
    eval(['handles.data.PIV' handles.data.cpass '.RPCd = get(hObject,''String'');'])
    guidata(hObject,handles)
    if get(handles.correlationtype,'Value')==2
        if str2double(get(hObject,'String'))<2
            if str2double(get(hObject,'String'))==0
                set(hObject,'backgroundcolor','r');
            else
                set(hObject,'backgroundcolor',[1 0.5 0]);
            end

        else
            set(hObject,'backgroundcolor',[1 1 1]);
        end
    else
        set(hObject,'backgroundcolor',0.5*[1 1 1]);
    end
end
function rpcdiameter_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Smoothing Check Box ---
function smoothingcheckbox_Callback(hObject, eventdata, handles)
if str2double(handles.Njob)>0
    eval(['handles.data.PIV' handles.data.cpass '.velsmooth = num2str(get(hObject,''Value''));'])
    eval(['size=handles.data.PIV' handles.data.cpass '.velsmoothfilt;'])
    if get(hObject,'Value')==1 && get(handles.passtype,'Value')>1
        set(handles.smoothingsize,'string',size,'backgroundcolor',[1 1 1]);
    else
        set(handles.smoothingsize,'string',size,'backgroundcolor',0.5*[1 1 1]);
    end
    guidata(hObject,handles)
end

% --- Smoothing Size Text Box ---
function smoothingsize_Callback(hObject, eventdata, handles)
if str2double(handles.Njob)>0 && get(handles.smoothingcheckbox,'Value')==1
    eval(['handles.data.PIV' handles.data.cpass '.velsmoothfilt = num2str(get(hObject,''Value''));'])
    guidata(hObject,handles)
end
function smoothingsize_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Validation Check Box ---
function validatecheckbox_Callback(hObject, eventdata, handles)
if str2double(handles.Njob)>0
    eval(['handles.data.PIV' handles.data.cpass '.val = num2str(get(hObject,''Value''));'])
    handles=set_PIVcontrols(handles);
    guidata(hObject,handles)
end

% --- ? Button Next to Validation Options ---
function validationhelp_Callback(hObject, eventdata, handles)
PIVhelp(8)

% --- Universal Outlier Detection Checkbox ---
function uodcheckbox_Callback(hObject, eventdata, handles)
if str2double(handles.Njob)>0
    eval(['handles.data.PIV' handles.data.cpass '.uod = num2str(get(hObject,''Value''));'])
    handles=set_PIVcontrols(handles);
    guidata(hObject,handles)
end

% --- UOD Type Drop-Down Menu ---
function uod_type_Callback(hObject, eventdata, handles)
if str2double(handles.Njob)>0 && get(handles.validatecheckbox,'Value')==1
    eval(['handles.data.PIV' handles.data.cpass '.uod_type = num2str(get(hObject,''Value''));'])
    guidata(hObject,handles)
end
function uod_type_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- UOD Window Size Text Box ---
function uod_window_Callback(hObject, eventdata, handles)
if str2double(handles.Njob)>0 && get(handles.validatecheckbox,'Value')==1
    A=get(hObject,'String');
    if strcmp(A(end),';')
        A=A(1:end-1);
    end
    set(hObject,'String',A)
    eval(['handles.data.PIV' handles.data.cpass '.uod_window = get(hObject,''String'');'])
    guidata(hObject,handles)
end
function uod_window_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- UOD Threshold Text Box ---
function uod_thresh_Callback(hObject, eventdata, handles)
if str2double(handles.Njob)>0 && get(handles.validatecheckbox,'Value')==1
    eval(['handles.data.PIV' handles.data.cpass '.uod_thresh = get(hObject,''String'');'])
    guidata(hObject,handles)
end
function uod_thresh_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Bootstrapping Checkbox ---
function bootstrapcheckbox_Callback(hObject, eventdata, handles)
if str2double(handles.Njob)>0
    eval(['handles.data.PIV' handles.data.cpass '.bootstrap = num2str(get(hObject,''Value''));'])
    handles=set_PIVcontrols(handles);
    guidata(hObject,handles)
end

% --- Bootstrapping Percent Removed Text Box ---
function bootstrap_percentsampled_Callback(hObject, eventdata, handles)
if str2double(handles.Njob)>0  && get(handles.validatecheckbox,'Value')==1
    eval(['handles.data.PIV' handles.data.cpass '.bootstrap_percentsampled = get(hObject,''String'');'])
    guidata(hObject,handles)
end
function bootstrap_percentsampled_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Bootstrapping Interpolations per Frame Text Box ---
function bootstrap_iterations_Callback(hObject, eventdata, handles)
if str2double(handles.Njob)>0  && get(handles.validatecheckbox,'Value')==1
    eval(['handles.data.PIV' handles.data.cpass '.bootstrap_iterations = get(hObject,''String'');'])
    guidata(hObject,handles)
end
function bootstrap_iterations_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Bootstrapping Number of Passes Text Box ---
function bootstrap_passes_Callback(hObject, eventdata, handles)
if str2double(handles.Njob)>0  && get(handles.validatecheckbox,'Value')==1
    eval(['handles.data.PIV' handles.data.cpass '.bootstrap_passes = get(hObject,''String'');'])
    guidata(hObject,handles)
end
function bootstrap_passes_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Thresholding Check Box ---
function thresholdingcheckbox_Callback(hObject, eventdata, handles)
if str2double(handles.Njob)>0 
    eval(['handles.data.PIV' handles.data.cpass '.thresh = num2str(get(hObject,''Value''));'])
    handles=set_PIVcontrols(handles);
    guidata(hObject,handles)
end

% --- U Threshold Text Box ---
function thresh_U_Callback(hObject, eventdata, handles)
if str2double(handles.Njob)>0  && get(handles.validatecheckbox,'Value')==1
    eval(['handles.data.PIV' handles.data.cpass '.valuthresh = get(hObject,''String'');'])
    a=get(hObject,'String');
    tx=str2double(a(1:(strfind(a,',')-1)));
    ty=str2double(a((strfind(a,',')+1):end));
    if get(handles.validatecheckbox,'Value')+get(handles.thresholdingcheckbox,'Value')==2
        if tx>=ty
            set(hObject,'backgroundcolor','r');
        else
            set(hObject,'backgroundcolor',[1 1 1]);
        end
    else
        set(hObject,'backgroundcolor',0.5*[1 1 1]);
    end
    guidata(hObject,handles)
end
function thresh_U_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- V Threshold Text Box ---
function thresh_V_Callback(hObject, eventdata, handles)
if str2double(handles.Njob)>0  && get(handles.validatecheckbox,'Value')==1
    eval(['handles.data.PIV' handles.data.cpass '.valvthresh = get(hObject,''String'');'])
    a=get(hObject,'String');
    tx=str2double(a(1:(strfind(a,',')-1)));
    ty=str2double(a((strfind(a,',')+1):end));
    if get(handles.validatecheckbox,'Value')+get(handles.thresholdingcheckbox,'Value')==2
        if tx>=ty
            set(hObject,'backgroundcolor','r');
        else
            set(hObject,'backgroundcolor',[1 1 1]);
        end
    else
        set(hObject,'backgroundcolor',0.5*[1 1 1]);
    end
    guidata(hObject,handles)
end
function thresh_V_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Validate Extra Peaks if Initial Validation Fails Checkbox ---
function valextrapeaks_Callback(hObject, eventdata, handles)
if str2double(handles.Njob)>0
    eval(['handles.data.PIV' handles.data.cpass '.valextrapeaks = num2str(get(hObject,''Value''));'])
    guidata(hObject,handles)
end

% --- Write Output Check Box ---
function writeoutputcheckbox_Callback(hObject, eventdata, handles)
if str2double(handles.Njob)>0
    if get(hObject,'Value')==0
        set(handles.corrpeaknum,'BackgroundColor',0.5*[1 1 1])
    else
        set(handles.corrpeaknum,'BackgroundColor',[1 1 1])
    end
    eval(['handles.data.PIV' handles.data.cpass '.write = num2str(get(hObject,''Value''));'])
    handles=set_PIVcontrols(handles);
    guidata(hObject,handles)
end

% --- Save Additional Peak Information Checkbox ---
function savepeakinfo_Callback(hObject, eventdata, handles)
if str2double(handles.Njob)>0
    eval(['handles.data.PIV' handles.data.cpass '.savepeakinfo = num2str(get(hObject,''Value''));'])
    handles=set_PIVcontrols(handles);
    guidata(hObject,handles)
end

% --- Number of Correlation Peaks to Save Drop-Down Menu ---
function corrpeaknum_Callback(hObject, eventdata, handles)
if str2double(handles.Njob)>0
    eval(['handles.data.PIV' handles.data.cpass '.corrpeaknum = num2str(get(hObject,''Value''));'])
    guidata(hObject,handles)
end
function corrpeaknum_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Save Peak Magnitude Checkbox ---
function savepeakmag_Callback(hObject, eventdata, handles)
if str2double(handles.Njob)>0
    eval(['handles.data.PIV' handles.data.cpass '.savepeakmag = num2str(get(hObject,''Value''));'])
    guidata(hObject,handles)
end

% --- Save Resulting Velocity Checkbox ---
function savepeakvel_Callback(hObject, eventdata, handles)
if str2double(handles.Njob)>0
    eval(['handles.data.PIV' handles.data.cpass '.savepeakvel = num2str(get(hObject,''Value''));'])
    guidata(hObject,handles)
end

% --- Output Basename Text Box ---
function outputbasename_Callback(hObject, eventdata, handles)
if str2double(handles.Njob)>0
    eval(['handles.data.PIV' handles.data.cpass '.outbase = get(hObject,''String'');'])
    guidata(hObject,handles)
end
function outputbasename_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Laser Pulse Separation Text Box ---
function pulseseparation_Callback(hObject, eventdata, handles)
if str2double(handles.Njob)>0
    handles.data.wrsep = get(hObject,'String');
    guidata(hObject,handles)
    if str2double(get(hObject,'String'))<=0
        set(hObject,'backgroundcolor','r')
    else
        set(hObject,'backgroundcolor',[1 1 1])
    end
end
function pulseseparation_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Camera Frame Rate Textbox ---
function samplingrate_Callback(hObject, eventdata, handles)
if str2double(handles.Njob)>0
    handles.data.wrsamp = get(hObject,'String');
    guidata(hObject,handles)
    if str2double(get(hObject,'String'))<=0
        set(hObject,'backgroundcolor','r')
    else
        set(hObject,'backgroundcolor',[1 1 1])
    end
end
function samplingrate_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Image Resolution Textbox ---
function magnification_Callback(hObject, eventdata, handles)
if str2double(handles.Njob)>0
    handles.data.wrmag = get(hObject,'String');
    [handles]=load_data(handles);
    guidata(hObject,handles)
    if str2double(get(hObject,'String'))<=0
        set(hObject,'backgroundcolor','r')
    else
        set(hObject,'backgroundcolor',[1 1 1])
    end
end
function magnification_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Experiment Date Textbox ---
function exp_date_Callback(hObject, eventdata, handles)
if str2double(handles.Njob)>0
    handles.data.exp_date = get(hObject,'String');
    guidata(hObject,handles)
end
function exp_date_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Experiment Characteristic Length Textbox ---
function exp_L_Callback(hObject, eventdata, handles)
if str2double(handles.Njob)>0
    handles.data.exp_L = get(hObject,'String');
    [handles]=load_data(handles);
    guidata(hObject,handles)
end
function exp_L_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Experiment Characteristic Velocity Textbox ---
function exp_v0_Callback(hObject, eventdata, handles)
if str2double(handles.Njob)>0
    handles.data.exp_v0 = get(hObject,'String');
    [handles]=load_data(handles);
    guidata(hObject,handles)
end
function exp_v0_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Experiment Notes Textbox ---
function exp_notesbox_Callback(hObject, eventdata, handles)
if str2double(handles.Njob)>0
    handles.data.exp_notes = get(hObject,'String');
    guidata(hObject,handles)
end
function exp_notesbox_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Experiment Density Textbox ---
function exp_density_Callback(hObject, eventdata, handles)
if str2double(handles.Njob)>0
    handles.data.exp_density = get(hObject,'String');
    [handles]=load_data(handles);
    guidata(hObject,handles)
end
function exp_density_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Experiment Viscosity Textbox ---
function exp_viscosity_Callback(hObject, eventdata, handles)
if str2double(handles.Njob)>0
    handles.data.exp_viscosity = get(hObject,'String');
    [handles]=load_data(handles);
    guidata(hObject,handles)
end
function exp_viscosity_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Experiment Surface Tension Textbox ---
function exp_surfacetension_Callback(hObject, eventdata, handles)
if str2double(handles.Njob)>0
    handles.data.exp_surfacetension = get(hObject,'String');
    [handles]=load_data(handles);
    guidata(hObject,handles)
end
function exp_surfacetension_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Experiment Particle Diameter Textbox ---
function exp_partD_Callback(hObject, eventdata, handles)
if str2double(handles.Njob)>0
    handles.data.exp_partD = get(hObject,'String');
    guidata(hObject,handles)
end
function exp_partD_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Experiment Particle Density Textbox ---
function exp_partdensity_Callback(hObject, eventdata, handles)
if str2double(handles.Njob)>0
    handles.data.exp_partdensity = get(hObject,'String');
    guidata(hObject,handles)
end
function exp_partdensity_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Experiment Laser Wavelength Textbox ---
function exp_wavelength_Callback(hObject, eventdata, handles)
if str2double(handles.Njob)>0
    handles.data.exp_wavelength = get(hObject,'String');
    guidata(hObject,handles)
end
function exp_wavelength_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Physical Pixel Size Textbox ---
function exp_pixelsize_Callback(hObject, eventdata, handles)
if str2double(handles.Njob)>0
    handles.data.exp_pixelsize = get(hObject,'String');
    guidata(hObject,handles)
end
function exp_pixelsize_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Experiment Lens Focal Length Textbox ---
function exp_lensfocal_Callback(hObject, eventdata, handles)
if str2double(handles.Njob)>0
    handles.data.exp_lensfocal = get(hObject,'String');
    [handles]=load_data(handles);
    guidata(hObject,handles)
end
function exp_lensfocal_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Experiment Lens f# Textbox ---
function exp_lensfnum_Callback(hObject, eventdata, handles)
if str2double(handles.Njob)>0
    if ~str2double(handles.data.exp_micro)
        handles.data.exp_lensfnum = get(hObject,'String');
        [handles]=load_data(handles);
        guidata(hObject,handles)
    else
        set(hObject,'String',handles.data.exp_lensfnum)
    end
end
function exp_lensfnum_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Experiment Micro-PIV Checkbox ---
function exp_microcheckbox_Callback(hObject, eventdata, handles)
if str2double(handles.Njob)>0
    handles.data.exp_micro=num2str(get(hObject,'Value'));
    [handles]=load_data(handles);
    guidata(hObject,handles)
end

% --- Experiment Numerical Aperture Textbox ---
function exp_NA_Callback(hObject, eventdata, handles)
if str2double(handles.Njob)>0
    if str2double(handles.data.exp_micro)
        handles.data.exp_NA=get(hObject,'String');
        [handles]=load_data(handles);
        guidata(hObject,handles)
    else
        set(hObject,'String',handles.data.exp_NA)
    end
end
function exp_NA_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Experiment Index of Refraction Textbox ---
function exp_n_Callback(hObject, eventdata, handles)
if str2double(handles.Njob)>0
    if str2double(handles.data.exp_micro)
        handles.data.exp_n=get(hObject,'String');
        [handles]=load_data(handles);
        guidata(hObject,handles)
    else
        set(hObject,'String',handles.data.exp_n)
    end
end
function exp_n_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Experiment Reynolds Number Textbox ---
function exp_Re_Callback(hObject, eventdata, handles)
if str2double(handles.Njob)>0
    set(hObject,'String',handles.data.exp_Re);
    guidata(hObject,handles)
end
function exp_Re_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Experiment Particle Stokes Number Textbox ---
function exp_St_Callback(hObject, eventdata, handles)
if str2double(handles.Njob)>0
    set(hObject,'String',handles.data.exp_St);
    guidata(hObject,handles)
end
function exp_St_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Experiment Magnification Textbox ---
function exp_M_Callback(hObject, eventdata, handles)
if str2double(handles.Njob)>0
    set(hObject,'String',handles.data.exp_M);
    guidata(hObject,handles)
end
function exp_M_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Experiment Region of Interest (m) Textbox ---
function exp_ROI_Callback(hObject, eventdata, handles)
if str2double(handles.Njob)>0
    set(hObject,'String',handles.data.exp_ROI);
    guidata(hObject,handles)
end
function exp_ROI_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Experiment Diffraction Diameter Textbox ---
function exp_diffractiondiameter_Callback(hObject, eventdata, handles)
if str2double(handles.Njob)>0
    set(hObject,'String',handles.data.exp_diffractiondiameter);
    guidata(hObject,handles)
end
function exp_diffractiondiameter_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Experiment Depth of Focus Textbox ---
function exp_depthoffocus_Callback(hObject, eventdata, handles)
if str2double(handles.Njob)>0
    set(hObject,'String',handles.data.exp_depthoffocus);
    guidata(hObject,handles)
end
function exp_depthoffocus_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Help Topics List ---
function helplistbox_Callback(hObject, eventdata, handles)
update_helptext(handles)
function helplistbox_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Help Text Display ---
function helptextbox_Callback(hObject, eventdata, handles)
update_helptext(handles)
function helptextbox_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Update All Data ---
function handles=update_data(handles)
set(handles.exp_Re,'backgroundcolor',0.8*[1 1 1]);
set(handles.exp_St,'backgroundcolor',0.8*[1 1 1]);
set(handles.exp_M,'backgroundcolor',0.8*[1 1 1]);
set(handles.exp_ROI,'backgroundcolor',0.8*[1 1 1]);
set(handles.exp_diffractiondiameter,'backgroundcolor',0.8*[1 1 1]);
set(handles.exp_depthoffocus,'backgroundcolor',0.8*[1 1 1]);
if str2double(handles.Njob)==0
    set(handles.passlist,'String','','backgroundcolor',0.5*[1 1 1]);
    set(handles.joblist,'String','','backgroundcolor',0.5*[1 1 1]);
    set(handles.imagelist,'String','','backgroundcolor','r');
    set(handles.uod_window,'String','','backgroundcolor',0.5*[1 1 1]);
    set(handles.bootstrap_percentsampled,'String','','backgroundcolor',0.5*[1 1 1]);
    set(handles.bootstrap_iterations,'String','','backgroundcolor',0.5*[1 1 1]);
    set(handles.bootstrap_passes,'String','','backgroundcolor',0.5*[1 1 1]);
    set(handles.thresh_V,'String','','backgroundcolor',0.5*[1 1 1]);
    set(handles.thresh_U,'String','','backgroundcolor',0.5*[1 1 1]);
    set(handles.uod_thresh,'String','','backgroundcolor',0.5*[1 1 1]);
    set(handles.corrpeaknum,'backgroundcolor',0.5*[1 1 1]);
    set(handles.outputbasename,'String','','backgroundcolor',0.5*[1 1 1]);
    set(handles.outputdirectory,'String','','backgroundcolor',0.5*[1 1 1]);
    set(handles.imagedirectory,'String','','backgroundcolor',0.5*[1 1 1]);
    set(handles.imagebasename,'String','','backgroundcolor',0.5*[1 1 1]);
    set(handles.imagezeros,'String','','backgroundcolor',0.5*[1 1 1]);
    set(handles.imageextension,'String','','backgroundcolor',0.5*[1 1 1]);
    set(handles.imagecorrelationstep,'String','','backgroundcolor',0.5*[1 1 1]);
    set(handles.imageframestep,'String','','backgroundcolor',0.5*[1 1 1]);
    set(handles.imageframestart,'String','','backgroundcolor',0.5*[1 1 1]);
    set(handles.imageframeend,'String','','backgroundcolor',0.5*[1 1 1]);
    set(handles.magnification,'String','','backgroundcolor',0.5*[1 1 1]);
    set(handles.currentjobname,'String','','backgroundcolor',0.5*[1 1 1]);
    set(handles.pulseseparation,'String','','backgroundcolor',0.5*[1 1 1]);
    set(handles.samplingrate,'String','','backgroundcolor',0.5*[1 1 1]);
    set(handles.smoothingsize,'String','','backgroundcolor',0.5*[1 1 1]);
    set(handles.windowres,'String','','backgroundcolor',0.5*[1 1 1]);
    set(handles.gridres,'String','','backgroundcolor',0.5*[1 1 1]);
    set(handles.winoverlap,'String','','backgroundcolor',0.5*[1 1 1]);
    set(handles.gridbuffer,'String','','backgroundcolor',0.5*[1 1 1]);
    set(handles.windowsize,'String','','backgroundcolor',0.5*[1 1 1]);
    set(handles.rpcdiameter,'String','','backgroundcolor',0.5*[1 1 1]);
    set(handles.bulkwinoffset,'String','','backgroundcolor',0.5*[1 1 1]);
    set(handles.writeoutputcheckbox,'Value',0,'backgroundcolor',handles.syscolor);
    set(handles.validatecheckbox,'Value',0,'backgroundcolor',handles.syscolor);
    set(handles.writeoutputcheckbox,'Value',0,'backgroundcolor',handles.syscolor);
    set(handles.thresholdingcheckbox,'Value',0,'backgroundcolor',handles.syscolor);
    set(handles.smoothingcheckbox,'Value',0,'backgroundcolor',handles.syscolor);
    set(handles.imageinterptype,'Value',1,'backgroundcolor',0.5*[1 1 1]);
    set(handles.framestep,'String','','backgroundcolor',0.5*[1 1 1]);
    set(handles.PIVerror,'String','','backgroundcolor',0.5*[1 1 1]);
    set(handles.passtype,'Value',1,'backgroundcolor',0.5*[1 1 1]);
    set(handles.uod_type,'Value',1,'backgroundcolor',0.5*[1 1 1]);
    set(handles.correlationtype,'Value',1,'backgroundcolor',0.5*[1 1 1]);
    set(handles.velocityinterptype,'Value',1,'backgroundcolor',0.5*[1 1 1]);
    set(handles.staticmaskfile,'String','','backgroundcolor',0.5*[1 1 1]);
    set(handles.maskdirectory,'String','','backgroundcolor',0.5*[1 1 1]);
    set(handles.maskbasename,'String','','backgroundcolor',0.5*[1 1 1]);
    set(handles.maskzeros,'String','','backgroundcolor',0.5*[1 1 1]);
    set(handles.maskextension,'String','','backgroundcolor',0.5*[1 1 1]);
    set(handles.maskframestep,'String','','backgroundcolor',0.5*[1 1 1]);
    set(handles.maskframestart,'String','','backgroundcolor',0.5*[1 1 1]);
    set(handles.masklist,'String','','backgroundcolor',0.5*[1 1 1]);
    set(handles.exp_date,'String','','backgroundcolor',0.5*[1 1 1]);
    set(handles.exp_wavelength,'String','','backgroundcolor',0.5*[1 1 1]);
    set(handles.exp_pixelsize,'String','','backgroundcolor',0.5*[1 1 1]);
    set(handles.exp_lensfocal,'String','','backgroundcolor',0.5*[1 1 1]);
    set(handles.exp_lensfnum,'String','','backgroundcolor',0.5*[1 1 1]);
    set(handles.exp_partD,'String','','backgroundcolor',0.5*[1 1 1]);
    set(handles.exp_partdensity,'String','','backgroundcolor',0.5*[1 1 1]);
    set(handles.exp_viscosity,'String','','backgroundcolor',0.5*[1 1 1]);
    set(handles.exp_density,'String','','backgroundcolor',0.5*[1 1 1]);
    set(handles.exp_surfacetension,'String','','backgroundcolor',0.5*[1 1 1]);
    set(handles.exp_v0,'String','','backgroundcolor',0.5*[1 1 1]);
    set(handles.exp_L,'String','','backgroundcolor',0.5*[1 1 1]);
    set(handles.exp_notesbox,'String','','backgroundcolor',0.5*[1 1 1]);
    set(handles.exp_NA,'String','','backgroundcolor',0.5*[1 1 1]);
    set(handles.exp_n,'String','','backgroundcolor',0.5*[1 1 1]);
else
    a=get(handles.joblist,'String');
    eval(['handles.data=handles.' char(a(str2double(handles.Cjob),:)) ';']);

    set(handles.passlist,'String','','backgroundcolor',[1 1 1]);
    set(handles.joblist,'backgroundcolor',[1 1 1]);
    set(handles.imagelist,'backgroundcolor',[1 1 1]);
    set(handles.windowres,'String','','backgroundcolor',[1 1 1]);
    set(handles.windowsize,'String','','backgroundcolor',[1 1 1]);
    set(handles.gridres,'String','','backgroundcolor',[1 1 1]);
    set(handles.gridbuffer,'String','','backgroundcolor',[1 1 1]);
    set(handles.winoverlap,'String','','backgroundcolor',[1 1 1]);
    set(handles.rpcdiameter,'String','','backgroundcolor',[1 1 1]);
    set(handles.bulkwinoffset,'String','','backgroundcolor',[1 1 1]);
    set(handles.uod_window,'String','','backgroundcolor',[1 1 1]);
    set(handles.bootstrap_percentsampled,'String','','backgroundcolor',[1 1 1]);
    set(handles.bootstrap_iterations,'String','','backgroundcolor',[1 1 1]);
    set(handles.bootstrap_passes,'String','','backgroundcolor',[1 1 1]);
    set(handles.thresh_V,'String','','backgroundcolor',[1 1 1]);
    set(handles.thresh_U,'String','','backgroundcolor',[1 1 1]);
    set(handles.uod_thresh,'String','','backgroundcolor',[1 1 1]);
    set(handles.corrpeaknum,'backgroundcolor',[1 1 1]);
    set(handles.outputbasename,'String','','backgroundcolor',[1 1 1]);
    set(handles.imagedirectory,'String','','backgroundcolor',[1 1 1]);
    set(handles.outputdirectory,'String','','backgroundcolor',[1 1 1]);
    set(handles.imagebasename,'String','','backgroundcolor',[1 1 1]);
    set(handles.imagezeros,'String','','backgroundcolor',[1 1 1]);
    set(handles.imageextension,'String','','backgroundcolor',[1 1 1]);
    set(handles.imagecorrelationstep,'String','','backgroundcolor',[1 1 1]);
    set(handles.imageframestep,'String','','backgroundcolor',[1 1 1]);
    set(handles.imageframestart,'String','','backgroundcolor',[1 1 1]);
    set(handles.imageframeend,'String','','backgroundcolor',[1 1 1]);
    set(handles.magnification,'String','','backgroundcolor',[1 1 1]);
    set(handles.staticmaskfile,'String','','backgroundcolor',[1 1 1]);
    set(handles.currentjobname,'String','','backgroundcolor',[1 1 1]);
    set(handles.pulseseparation,'String','','backgroundcolor',[1 1 1]);
    set(handles.samplingrate,'String','','backgroundcolor',[1 1 1]);
    set(handles.smoothingsize,'String','','backgroundcolor',[1 1 1]);
    set(handles.imageinterptype,'Value',1,'backgroundcolor',[1 1 1]);
    set(handles.passtype,'Value',1,'backgroundcolor',[1 1 1]);
    set(handles.uod_type,'Value',1,'backgroundcolor',[1 1 1]);
    set(handles.correlationtype,'Value',1,'backgroundcolor',[1 1 1]);
    set(handles.velocityinterptype,'Value',1,'backgroundcolor',[1 1 1]);
    set(handles.framestep,'String','','backgroundcolor',[1 1 1]);
    set(handles.PIVerror,'String','','backgroundcolor',[1 1 1]);
    set(handles.exp_date,'String','','backgroundcolor',[1 1 1]);
    set(handles.exp_wavelength,'String','','backgroundcolor',[1 1 1]);
    set(handles.exp_lensfocal,'String','','backgroundcolor',[1 1 1]);
    set(handles.exp_pixelsize,'String','','backgroundcolor',[1 1 1]);
    set(handles.exp_partD,'String','','backgroundcolor',[1 1 1]);
    set(handles.exp_partdensity,'String','','backgroundcolor',[1 1 1]);
    set(handles.exp_viscosity,'String','','backgroundcolor',[1 1 1]);
    set(handles.exp_density,'String','','backgroundcolor',[1 1 1]);
    set(handles.exp_surfacetension,'String','','backgroundcolor',[1 1 1]);
    set(handles.exp_v0,'String','','backgroundcolor',[1 1 1]);
    set(handles.exp_L,'String','','backgroundcolor',[1 1 1]);
    set(handles.exp_notesbox,'String','','backgroundcolor',[1 1 1]);
    
    if strcmp(handles.data.masktype,'static')
        set(handles.staticmaskfile,'String','','backgroundcolor',[1 1 1]);
    else
        set(handles.staticmaskfile,'String','','backgroundcolor',0.5*[1 1 1]);
    end
    if strcmp(handles.data.masktype,'dynamic')
        set(handles.maskdirectory,'String','','backgroundcolor',[1 1 1]);
        set(handles.maskbasename,'String','','backgroundcolor',[1 1 1]);
        set(handles.maskzeros,'String','','backgroundcolor',[1 1 1]);
        set(handles.maskextension,'String','','backgroundcolor',[1 1 1]);
        set(handles.maskframestep,'String','','backgroundcolor',[1 1 1]);
        set(handles.maskframestart,'String','','backgroundcolor',[1 1 1]);
        set(handles.masklist,'String','','backgroundcolor',[1 1 1]);
        load_masklist(handles);
    else
        set(handles.maskdirectory,'String','','backgroundcolor',0.5*[1 1 1]);
        set(handles.maskbasename,'String','','backgroundcolor',0.5*[1 1 1]);
        set(handles.maskzeros,'String','','backgroundcolor',0.5*[1 1 1]);
        set(handles.maskextension,'String','','backgroundcolor',0.5*[1 1 1]);
        set(handles.maskframestep,'String','','backgroundcolor',0.5*[1 1 1]);
        set(handles.maskframestart,'String','','backgroundcolor',0.5*[1 1 1]);
        set(handles.masklist,'backgroundcolor',0.5*[1 1 1]);
    end
    if str2double(handles.data.datout)==1
        set(handles.datcheckbox,'Value',1)
    else
        set(handles.datcheckbox,'Value',0)
    end
    if str2double(handles.data.multiplematout)==1
        set(handles.multiplematcheckbox,'Value',1)
    else
        set(handles.multiplematcheckbox,'Value',0)
    end
        
    load_data(handles);
    load_PIVlist(handles);
    handles=set_PIVcontrols(handles);
    load_imlist(handles);
end

% --- Load Image List ---
function load_imlist(handles)
dir_struct = dir(handles.data.imdirec);
if isempty(dir_struct)
    set(handles.imagedirectory,'backgroundcolor','r');
else
    set(handles.imagedirectory,'backgroundcolor',[1 1 1]);
end
N=length(str2double(handles.data.imfstart):str2double(handles.data.imfstep):str2double(handles.data.imfend));
files = cell(N,2);
e=0;
for f=str2double(handles.data.imfstart):str2double(handles.data.imfstep):str2double(handles.data.imfend)
    e=e+1;
    files(e,1)={[handles.data.imbase sprintf(['%0.' handles.data.imzeros 'i.' handles.data.imext],f)]};
    files(e,2)={[handles.data.imbase sprintf(['%0.' handles.data.imzeros 'i.' handles.data.imext],f+str2double(handles.data.imcstep))]};
end

[sorted_names,sorted_index] = sortrows({dir_struct.name}');
[files1,id,id1] = intersect(sorted_names,files(:,1));
[files2,id,id2] = intersect(sorted_names,files(:,2));
[idf]=intersect(id1,id2);
if isempty(idf)
    set(handles.imagelist,'backgroundcolor','r');
else
    set(handles.imagelist,'backgroundcolor',[1 1 1]);
    if length(idf)~=size(files,1)
        set(handles.imagelist,'backgroundcolor','r');
    end
end
files=files(idf,:);
filesf=cell(length(files),1);
for e=1:size(files,1)
    filesf(e)={[char(files(e,1)) ' and ' char(files(e,2))]};
end
set(handles.imagelist,'String',filesf,'Value',1);
handles=load_data(handles);

% --- Load Mask List ---
function load_masklist(handles)
dir_struct = dir(handles.data.maskdirec);
if isempty(dir_struct)
    set(handles.maskdirectory,'backgroundcolor','r');
else
    set(handles.maskdirectory,'backgroundcolor',[1 1 1]);
end

maskfend=str2double(handles.data.maskfstart)+str2double(handles.data.maskfstep)*length(str2double(handles.data.imfstart):str2double(handles.data.imfstep):str2double(handles.data.imfend))-1;
N=length(str2double(handles.data.maskfstart):str2double(handles.data.maskfstep):maskfend);
files = cell(N,1);
e=0;
for f=str2double(handles.data.maskfstart):str2double(handles.data.maskfstep):maskfend
    e=e+1;
    files(e,1)={[handles.data.maskbase sprintf(['%0.' handles.data.maskzeros 'i.' handles.data.maskext],f)]};
end

[sorted_names,sorted_index] = sortrows({dir_struct.name}');
[files1,id,id1] = intersect(sorted_names,files(:,1));

if isempty(id1)
    set(handles.masklist,'backgroundcolor','r');
    set(handles.masklist,'UserData',[]);
else
    set(handles.masklist,'backgroundcolor',[1 1 1]);
    if length(id1)~=size(files,1)
        set(handles.masklist,'backgroundcolor','r');
        set(handles.masklist,'UserData',[]);
    end
end
files=files(id1,:);
filesf=cell(length(files),1);
for e=1:size(files,1)
    filesf(e)={[char(files(e,1))]};
end
set(handles.masklist,'String',filesf,'Value',1);

% --- Load PIV Pass List ---
function load_PIVlist(handles)
f=cell(str2double(handles.data.passes),1);
for e=1:str2double(handles.data.passes)
    f(e)={['Pass ' num2str(e)]};
end
set(handles.passlist,'String',f);

% --- Load PIV Data ---
function handles=set_PIVcontrols(handles)
N=get(handles.passlist,'Value');
A=eval(['handles.data.PIV' num2str(N)]);
set(handles.windowres,'string',A.winres);
set(handles.windowsize,'string',A.winsize);
set(handles.autowinsizecheckbox,'Value',str2double(A.winauto));
set(handles.gridres,'string',A.gridres);
set(handles.winoverlap,'string',A.winoverlap);
set(handles.gridbuffer,'string',A.gridbuf);
set(handles.bulkwinoffset,'string',A.BWO);
set(handles.correlationtype,'Value',str2double(A.corr));
set(handles.rpcdiameter,'string',str2double(A.RPCd));
set(handles.smoothingsize,'String',A.velsmoothfilt);
set(handles.smoothingcheckbox,'Value',str2double(A.velsmooth));
set(handles.validatecheckbox,'Value',str2double(A.val));
set(handles.uod_type,'Value',str2double(A.uod_type));
set(handles.thresholdingcheckbox,'Value',str2double(A.thresh));
set(handles.uod_window,'string',A.uod_window);
set(handles.uod_thresh,'string',A.uod_thresh);
set(handles.bootstrap_percentsampled,'String',A.bootstrap_percentsampled);
set(handles.bootstrap_iterations,'String',A.bootstrap_iterations);
set(handles.bootstrap_passes,'String',A.bootstrap_passes);
set(handles.thresh_U,'string',A.valuthresh);
set(handles.thresh_V,'string',A.valvthresh);
set(handles.valextrapeaks,'value',str2double(A.valextrapeaks));
set(handles.corrpeaknum,'value',str2double(A.corrpeaknum));
set(handles.savepeakinfo,'value',str2double(A.savepeakinfo));
set(handles.savepeakmag,'value',str2double(A.savepeakmag));
set(handles.savepeakvel,'value',str2double(A.savepeakvel));
set(handles.writeoutputcheckbox,'Value',str2double(A.write));
set(handles.outputbasename,'string',A.outbase);
handles.data.cpass=num2str(N);

if str2double(A.val)==1
    if strcmp(A.uod,'1')
        set(handles.uodcheckbox,'Value',1);
        set(handles.uod_window,'backgroundcolor',[1 1 1]);
        set(handles.uod_thresh,'backgroundcolor',[1 1 1]);
        set(handles.uod_type,'backgroundcolor',[1 1 1]);
    else
        set(handles.uodcheckbox,'Value',0);
        set(handles.uod_window,'backgroundcolor',0.5*[1 1 1]);
        set(handles.uod_thresh,'backgroundcolor',0.5*[1 1 1]);
        set(handles.uod_type,'backgroundcolor',0.5*[1 1 1]);
    end
    if strcmp(A.bootstrap,'1')
        set(handles.bootstrapcheckbox,'Value',1)
        set(handles.bootstrap_percentsampled,'backgroundcolor',[1 1 1]);
        set(handles.bootstrap_iterations,'backgroundcolor',[1 1 1]);
        set(handles.bootstrap_passes,'backgroundcolor',[1 1 1]);
    else
        set(handles.bootstrapcheckbox,'Value',0)
        set(handles.bootstrap_percentsampled,'backgroundcolor',0.5*[1 1 1]);
        set(handles.bootstrap_iterations,'backgroundcolor',0.5*[1 1 1]);
        set(handles.bootstrap_passes,'backgroundcolor',0.5*[1 1 1]);
    end
    if str2double(A.thresh)==1
        a=get(handles.thresh_U,'String');
        tx=str2double(a(1:(strfind(a,',')-1)));
        ty=str2double(a((strfind(a,',')+1):end));
        if tx>=ty
            set(handles.thresh_U,'backgroundcolor','r');
        else
            set(handles.thresh_U,'backgroundcolor',[1 1 1]);
        end
        a=get(handles.thresh_V,'String');
        tx=str2double(a(1:(strfind(a,',')-1)));
        ty=str2double(a((strfind(a,',')+1):end));
        if tx>=ty
            set(handles.thresh_V,'backgroundcolor','r');
        else
            set(handles.thresh_V,'backgroundcolor',[1 1 1]);
        end
    else
        set(handles.thresh_U,'backgroundcolor',0.5*[1 1 1]);
        set(handles.thresh_V,'backgroundcolor',0.5*[1 1 1]);
    end
else
    set(handles.bootstrap_percentsampled,'backgroundcolor',0.5*[1 1 1]);
    set(handles.bootstrap_iterations,'backgroundcolor',0.5*[1 1 1]);
    set(handles.bootstrap_passes,'backgroundcolor',0.5*[1 1 1]);
    set(handles.uod_window,'backgroundcolor',0.5*[1 1 1]);
    set(handles.uod_thresh,'backgroundcolor',0.5*[1 1 1]);
    set(handles.uod_type,'backgroundcolor',0.5*[1 1 1]);
    set(handles.thresh_U,'backgroundcolor',0.5*[1 1 1]);
    set(handles.thresh_V,'backgroundcolor',0.5*[1 1 1]);
end
if str2double(A.write)==1
    set(handles.outputbasename,'backgroundcolor',[1 1 1]);
    if str2double(A.savepeakinfo)==1
        set(handles.corrpeaknum,'backgroundcolor',[1 1 1]);
    else
        set(handles.corrpeaknum,'backgroundcolor',0.5*[1 1 1]);
    end
else
    set(handles.corrpeaknum,'backgroundcolor',0.5*[1 1 1]);
    set(handles.outputbasename,'backgroundcolor',0.5*[1 1 1]);
end
if str2double(A.winauto)==1
    B=get(handles.windowres,'String');
    Rx=str2double(B(1:(strfind(B,',')-1)));
    Ry=str2double(B((strfind(B,',')+1):end));
    Rx = 2^(ceil(log(2*Rx)/log(2)));
    Ry = 2^(ceil(log(2*Ry)/log(2)));
    set(handles.windowsize,'backgroundcolor',0.5*[1 1 1]);
    set(handles.windowsize,'String',[num2str(Rx) ',' num2str(Ry)]);
    eval(['handles.data.PIV' handles.data.cpass '.winsize = get(handles.windowsize,''String'');'])
else
    set(handles.windowsize,'backgroundcolor',[1 1 1]);
end

if strcmp(A.gridtype,'2') && str2double(handles.data.method)~=1
    set(handles.setgridresbutton,'Value',0)
    set(handles.setwinoverlapbutton,'Value',1)
    set(handles.gridres,'backgroundcolor',0.5*[1 1 1]);
    set(handles.winoverlap,'backgroundcolor',[1 1 1])
else
    eval(['handles.data.PIV' handles.data.cpass '.gridtype = ''1'';']);
    set(handles.setgridresbutton,'Value',1)
    set(handles.setwinoverlapbutton,'Value',0)
    set(handles.gridres,'backgroundcolor',[1 1 1]);
    set(handles.winoverlap,'backgroundcolor',0.5*[1 1 1]);
end
if get(handles.correlationtype,'Value')==2
    if str2double(get(handles.rpcdiameter,'String'))<2
        if str2double(get(handles.rpcdiameter,'String'))==0
            set(handles.rpcdiameter,'backgroundcolor','r');
        else
            set(handles.rpcdiameter,'backgroundcolor',[1 0.5 0]);
        end
    else
        set(handles.rpcdiameter,'backgroundcolor',[1 1 1]);
    end
else
    set(handles.rpcdiameter,'backgroundcolor',0.5*[1 1 1]);
end
a=get(handles.windowres,'String');
wx=str2double(a(1:(strfind(a,',')-1)));
wy=str2double(a((strfind(a,',')+1):end));
if wx*wy<256
    set(handles.windowres,'backgroundcolor',[1 0.5 0]);
else
    set(handles.windowres,'backgroundcolor',[1 1 1]);
end
if get(handles.writeoutputcheckbox,'Value')==0 && get(handles.passlist,'Value')==str2double(handles.data.passes)
    set(handles.writeoutputcheckbox,'backgroundcolor',[1 0.5 0]);
else
    set(handles.writeoutputcheckbox,'backgroundcolor',handles.syscolor);
end

% --- Load Extra Data ---
function [handles]=load_data(handles)
set(handles.imagedirectory,'String',handles.data.imdirec);
set(handles.imagebasename,'String',handles.data.imbase);
set(handles.imagezeros,'String',handles.data.imzeros);
set(handles.imageextension,'String',handles.data.imext);
set(handles.imagecorrelationstep,'String',handles.data.imcstep);
set(handles.imageframestep,'String',handles.data.imfstep);
set(handles.imageframestart,'String',handles.data.imfstart);
set(handles.imageframeend,'String',handles.data.imfend);

set(handles.staticmaskfile,'String',handles.data.staticmaskname)
set(handles.maskdirectory,'String',handles.data.maskdirec);
set(handles.maskbasename,'String',handles.data.maskbase);
set(handles.maskzeros,'String',handles.data.maskzeros);
set(handles.maskextension,'String',handles.data.maskext);
set(handles.maskframestep,'String',handles.data.maskfstep);
set(handles.maskframestart,'String',handles.data.maskfstart);
if strcmp(handles.data.masktype,'none')
    set(handles.nomaskbutton,'Value',1)
    set(handles.staticmaskbutton,'Value',0)
    set(handles.dynamicmaskbutton,'Value',0)
elseif strcmp(handles.data.masktype,'static')
    set(handles.nomaskbutton,'Value',0)
    set(handles.staticmaskbutton,'Value',1)
    set(handles.dynamicmaskbutton,'Value',0)
elseif strcmp(handles.data.masktype,'dynamic')
    set(handles.nomaskbutton,'Value',0)
    set(handles.staticmaskbutton,'Value',0)
    set(handles.dynamicmaskbutton,'Value',1)
end

set(handles.magnification,'String',handles.data.wrmag);
set(handles.pulseseparation,'String',handles.data.wrsep);
set(handles.samplingrate,'String',handles.data.wrsamp);
set(handles.currentjobname,'String',handles.data.batchname);
set(handles.outputdirectory,'String',handles.data.outdirec);

set(handles.exp_date,'String',handles.data.exp_date);
set(handles.exp_wavelength,'String',handles.data.exp_wavelength);
set(handles.exp_pixelsize,'String',handles.data.exp_pixelsize);
set(handles.exp_lensfocal,'String',handles.data.exp_lensfocal);
set(handles.exp_lensfnum,'String',handles.data.exp_lensfnum);
set(handles.exp_partD,'String',handles.data.exp_partD);
set(handles.exp_partdensity,'String',handles.data.exp_partdensity);
set(handles.exp_viscosity,'String',handles.data.exp_viscosity);
set(handles.exp_density,'String',handles.data.exp_density);
set(handles.exp_surfacetension,'String',handles.data.exp_surfacetension);
set(handles.exp_L,'String',handles.data.exp_L);
set(handles.exp_v0,'String',handles.data.exp_v0);
set(handles.exp_notesbox,'String',handles.data.exp_notes);
set(handles.exp_Re,'String',handles.data.exp_Re);
set(handles.exp_St,'String',handles.data.exp_St);
set(handles.exp_M,'String',handles.data.exp_M);
set(handles.exp_ROI,'String',handles.data.exp_ROI);
set(handles.exp_diffractiondiameter,'String',handles.data.exp_diffractiondiameter);
set(handles.exp_depthoffocus,'String',handles.data.exp_depthoffocus);
set(handles.exp_microcheckbox,'Value',str2double(handles.data.exp_micro));
set(handles.exp_NA,'String',handles.data.exp_NA);
set(handles.exp_n,'String',handles.data.exp_n);

if str2double(handles.data.exp_micro)
    set(handles.exp_NA,'backgroundcolor',[1 1 1]);
    set(handles.exp_n,'backgroundcolor',[1 1 1]);
    set(handles.exp_lensfnum,'backgroundcolor',0.5*[1 1 1]);
else
    set(handles.exp_NA,'backgroundcolor',0.5*[1 1 1]);
    set(handles.exp_n,'backgroundcolor',0.5*[1 1 1]);
    set(handles.exp_lensfnum,'backgroundcolor',[1 1 1]);
end

set(handles.passtype,'Value',str2double(handles.data.method));
set(handles.velocityinterptype,'Value',str2double(handles.data.velinterp));
set(handles.imageinterptype,'Value',str2double(handles.data.iminterp));
set(handles.framestep,'String',handles.data.framestep);
set(handles.PIVerror,'String',handles.data.PIVerror);

if get(handles.passtype,'Value')>=5
    set(handles.velocityinterptype,'backgroundcolor',0.5*[1 1 1]);
    set(handles.imageinterptype,'backgroundcolor',0.5*[1 1 1]);
    set(handles.smoothingsize,'backgroundcolor',0.5*[1 1 1]);
    set(handles.framestep,'backgroundcolor',[1 1 1]);
    set(handles.PIVerror,'backgroundcolor',[1 1 1]);
elseif get(handles.passtype,'Value')>1
    set(handles.framestep,'backgroundcolor',0.5*[1 1 1]);
    set(handles.PIVerror,'backgroundcolor',0.5*[1 1 1]);
    set(handles.velocityinterptype,'backgroundcolor',[1 1 1]);
    if get(handles.passtype,'Value')==3
        set(handles.imageinterptype,'backgroundcolor',[1 1 1]);
    else
        set(handles.imageinterptype,'backgroundcolor',0.5*[1 1 1]);
    end
    if get(handles.smoothingcheckbox,'Value')==1
        set(handles.smoothingsize,'backgroundcolor',[1 1 1]);
    else
        set(handles.smoothingsize,'backgroundcolor',0.5*[1 1 1]);
    end
else
    set(handles.velocityinterptype,'backgroundcolor',0.5*[1 1 1]);
    set(handles.imageinterptype,'backgroundcolor',0.5*[1 1 1]);
    set(handles.smoothingsize,'backgroundcolor',0.5*[1 1 1]);
    set(handles.framestep,'backgroundcolor',0.5*[1 1 1]);
    set(handles.PIVerror,'backgroundcolor',0.5*[1 1 1]);
end

if strcmp(handles.data.imext(1),'.')
    set(handles.imageextension,'backgroundcolor','r');
else
    set(handles.imageextension,'backgroundcolor',[1 1 1]);
end

if isempty(dir(handles.data.imdirec))
    set(handles.imagedirectory,'backgroundcolor','r');
else
    set(handles.imagedirectory,'backgroundcolor',[1 1 1]);
end

if str2double(get(handles.magnification,'String'))<=0
    set(handles.magnification,'backgroundcolor','r')
else
    set(handles.magnification,'backgroundcolor',[1 1 1])
end

if str2double(get(handles.pulseseparation,'String'))<=0
    set(handles.pulseseparation,'backgroundcolor','r')
else
    set(handles.pulseseparation,'backgroundcolor',[1 1 1])
end

if str2double(get(handles.samplingrate,'String'))<=0
    set(handles.samplingrate,'backgroundcolor','r')
else
    set(handles.samplingrate,'backgroundcolor',[1 1 1])
end

if ~(isempty(handles.data.exp_viscosity) || isempty(handles.data.exp_density) || isempty(handles.data.exp_L) || isempty(handles.data.exp_v0))
    try
        Re=str2double(handles.data.exp_density)*str2double(handles.data.exp_v0)*str2double(handles.data.exp_L)/str2double(handles.data.exp_viscosity);
        Re=num2str(Re);
        handles.data.exp_Re=Re;
        set(handles.exp_Re,'String',Re);
    catch
        handles.data.exp_Re='';
        set(handles.exp_Re,'String','')
    end
else
    handles.data.exp_Re='';
    set(handles.exp_Re,'String','')
end

if ~(isempty(handles.data.exp_partD) || isempty(handles.data.exp_partdensity) || isempty(handles.data.exp_viscosity) || isempty(handles.data.exp_L) || isempty(handles.data.exp_v0))
    try
        St=str2double(handles.data.exp_partdensity)*(str2double(handles.data.exp_partD)*10^-6)^2*str2double(handles.data.exp_v0)/18/str2double(handles.data.exp_viscosity)/str2double(handles.data.exp_L);
        St=num2str(St);
        handles.data.exp_St=St;
        set(handles.exp_St,'String',St);
    catch
        handles.data.exp_St='';
        set(handles.exp_St,'String','')
    end
else
    handles.data.exp_St='';
    set(handles.exp_St,'String','')
end

if ~isempty(handles.data.wrmag)
    try
        im1=double(imread([handles.data.imdirec '\' handles.data.imbase sprintf(['%0.' handles.data.imzeros 'i.' handles.data.imext],str2double(handles.data.imfstart))]));
        L=size(im1)*str2double(handles.data.wrmag)*10^-6;
        handles.data.exp_ROI=[num2str(L(2)),',',num2str(L(1))];

        if ~isempty(handles.data.exp_pixelsize)
            try
                Sx=str2double(handles.data.exp_pixelsize)*size(im1,2)*10^-6;
                Sy=str2double(handles.data.exp_pixelsize)*size(im1,1)*10^-6;
                S=sqrt(Sx^2+Sy^2);
                M=S/sqrt(L(1)^2+L(2)^2);
                
                if isnan(M)
                    handles.data.exp_M='';
                    handles.data.exp_diffractiondiameter='';
                    handles.data.exp_depthoffocus='';
                else
                    handles.data.exp_M=num2str(M);
                    if ~(isempty(handles.data.exp_lensfnum) || isempty(handles.data.exp_wavelength)) && ~str2double(handles.data.exp_micro)
                        try
                            d=2.44*str2double(handles.data.exp_lensfnum)*(M+1)*str2double(handles.data.exp_wavelength);
                            dz=2*str2double(handles.data.exp_lensfnum)*d*(M+1)/M^2;

                            if isnan(d) || isnan(dz)
                                handles.data.exp_diffractiondaimeter='';
                                handles.data.exp_depthoffocus='';
                            else
                                handles.data.exp_diffractiondiameter=num2str(d);
                                handles.data.exp_depthoffocus=num2str(dz);
                            end
                        catch
                            handles.data.exp_diffractiondaimeter='';
                            handles.data.exp_depthoffocus='';
                        end
                    elseif ~(isempty(handles.data.exp_NA) || isempty(handles.data.exp_n) || isempty(handles.data.exp_wavelength)) && str2double(handles.data.exp_micro)
                        try
                            fnum=0.5*sqrt((str2double(handles.data.exp_n)/str2double(handles.data.exp_NA))^2-1);
                            d=2.44*fnum*(M+1)*str2double(handles.data.exp_wavelength);
                            dz=2*fnum*d*(M+1)/M^2;

                            if isnan(d) || isnan(dz)
                                handles.data.exp_diffractiondaimeter='';
                                handles.data.exp_depthoffocus='';
                            else
                                handles.data.exp_diffractiondiameter=num2str(d);
                                handles.data.exp_depthoffocus=num2str(dz);
                            end
                        catch
                            handles.data.exp_diffractiondaimeter='';
                            handles.data.exp_depthoffocus='';
                        end
                    else
                        handles.data.exp_diffractiondiameter='';
                        handles.data.exp_depthoffocus='';
                    end
                end
            catch
                handles.data.exp_M='';
                handles.data.exp_diffractiondiameter='';
                handles.data.exp_depthoffocus='';
            end
        end
    catch
        handles.data.exp_ROI='';
        handles.data.exp_M='';
        handles.data.exp_diffractiondiameter='';
        handles.data.exp_depthoffocus='';
    end
else
    handles.data.exp_ROI='';
    handles.data.exp_M='';
    handles.data.exp_diffractiondiameter='';
    handles.data.exp_depthoffocus='';
end
set(handles.exp_M,'String',handles.data.exp_M)
set(handles.exp_ROI,'String',handles.data.exp_ROI)
set(handles.exp_diffractiondiameter,'String',handles.data.exp_diffractiondiameter)
set(handles.exp_depthoffocus,'String',handles.data.exp_depthoffocus)

function write_expsummary(Data,handles)
if ispc
    fname=[Data.outdirec,'\ExpSummary_',Data.batchname,'.txt'];
else
    fname=[Data.outdirec,'/ExpSummary_',Data.batchname,'.txt'];
end
fid=fopen(fname,'w');
if fid==-1
    error(['error writing experiment summary ',fname])
end

fprintf(fid,['Summary for PIV Job: ',Data.batchname,'\n']);
fprintf(fid,'\n--------------------Experiment Parameters--------------------\n');
fprintf(fid,['Date of Experiment:            ',Data.exp_date,'\n']);
fprintf(fid,['Laser Wavelength (um):         ',Data.exp_wavelength,'\n']);
fprintf(fid,['Laser Pulse Separation (us):   ',Data.wrsep,'\n']);
fprintf(fid,['Camera Frame Rate (Hz):        ',Data.wrsamp,'\n']);
fprintf(fid,['Pixel Size (um):               ',Data.exp_pixelsize,'\n']);
fprintf(fid,['Image Resolution (um/pix):     ',Data.wrmag,'\n']);
fprintf(fid,['Lens Focal Length (mm):        ',Data.exp_lensfocal,'\n']);
if ~str2double(Data.exp_micro)
    fprintf(fid,['Lens f#:                       ',Data.exp_lensfnum,'\n']);
else
    fprintf(fid,['Numerical Aperture:            ',Data.exp_NA,'\n']);
    fprintf(fid,['Index of Refraction:           ',Data.exp_n,'\n']);
end
fprintf(fid,['Particle Diameter (um):        ',Data.exp_partD,'\n']);
fprintf(fid,['Particle Density (kg/m^3):     ',Data.exp_partdensity,'\n']);
fprintf(fid,['Fluid Density (kg/m^3):        ',Data.exp_density,'\n']);
fprintf(fid,['Dynamic Viscosity (Pa*s):      ',Data.exp_viscosity,'\n']);
fprintf(fid,['Surface Tension (N/m):         ',Data.exp_surfacetension,'\n']);
fprintf(fid,['Characteristic Length (m):     ',Data.exp_L,'\n']);
fprintf(fid,['Characteristic Velocity (m/s): ',Data.exp_v0,'\n']);
fprintf(fid,['Reynolds Number:               ',Data.exp_Re,'\n']);
fprintf(fid,['Particle Stokes Number:        ',Data.exp_St,'\n']);
fprintf(fid,['Magnification Factor:          ',Data.exp_M,'\n']);
fprintf(fid,['Region of Interest (m):        ',Data.exp_ROI,'\n']);
fprintf(fid,['Diffraction Diameter (um):     ',Data.exp_diffractiondiameter,'\n']);
fprintf(fid,['Depth of Focus (um):           ',Data.exp_depthoffocus,'\n']);
for i=1:size(Data.exp_notes,1)
    fprintf(fid,[Data.exp_notes{i},'\n']);
end

methods={'Multipass','Multigrid','Deform','Ensemble','Multiframe - Persoons Method','Multiframe - Exp. Method','Direct-Phase'};
velinterp={'Nearest Neighbor','Bilinear','Cubic'};
iminterp={'Cardinal Function','Cardinal Function w/ Blackman Filter'};
fprintf(fid,'\n-----------------------PIV Processing------------------------\n');
    fprintf(fid,['Algorithm:                     ',methods{str2double(Data.method)},'\n']);
if str2double(Data.method)~=1 && str2double(Data.method)~=5
    fprintf(fid,['Velocity Interp Function:      ',velinterp{str2double(Data.velinterp)},'\n']);
end
if str2double(Data.method)==4
    fprintf(fid,['Image Interpolation Function:  ',iminterp{str2double(Data.iminterp)},'\n']);
end
if str2double(Data.method)>=5
    fprintf(fid,['PIV Error:                     ',Data.PIVerror,'\n']);
    fprintf(fid,['Maximum Framestep:             ',Data.framestep,'\n']);
end

for i=1:str2double(Data.passes)
    corr={'SCC','RPC'};
    y_n={'No','Yes'};
    A=eval(['Data.PIV' num2str(i)]);
    fprintf(fid,['\n------------------------Pass ',num2str(i),' Setup-------------------------\n']);
    fprintf(fid,['Window Resolution (pix):       ',A.winres,'\n']);
    fprintf(fid,['Window Size (pix):             ',A.winsize,'\n']);
    fprintf(fid,['Grid Resolution (pix):         ',A.gridres,'\n']);
    fprintf(fid,['Window Overlap Percentage:     ',A.winoverlap,'\n']);
    fprintf(fid,['Grid Buffer (pix):             ',A.gridbuf,'\n']);
    fprintf(fid,['Bulk Window Offset (pix):      ',A.BWO,'\n']);
    fprintf(fid,['Correlation:                   ',corr{str2double(A.corr)},'\n']);
    if str2double(Data.method)~=1 && str2double(Data.method)~=5
        fprintf(fid,['Smoothing?:                    ',y_n{str2double(A.velsmooth)+1},'\n']);
        if str2double(A.velsmooth)
            fprintf(fid,['Smoothing Size:                ',A.velsmoothfilt,'\n']);
        end
    end
    fprintf(fid,'Validation Type(s):            ');
    if str2double(A.val)
        if str2double(A.thresh)
            fprintf(fid,'Thresholding ');
        end
        if str2double(A.uod)
            fprintf(fid,'UOD ');
        end
        if str2double(A.bootstrap)
            fprintf(fid,'Bootstrapping');
        end
        fprintf(fid,'\n');
        if str2double(A.thresh)
            fprintf(fid,['Umin, Umax (pix):              ',A.valuthresh,'\n']);
            fprintf(fid,['Vmin, Vmax (pix):              ',A.valvthresh,'\n']);
        end
        if str2double(A.uod)
            uod_type={'Mean','Median'};
            fprintf(fid,['UOD Type:                      ',uod_type{str2double(A.uod_type)},'\n']);
            fprintf(fid,['UOD Window Sizes:              ',A.uod_window,'\n']);
            fprintf(fid,['UOD Thresholds:                ',A.uod_thresh,'\n']);
        end
        if str2double(A.bootstrap)
            fprintf(fid,['Bootstrap Percent Sampled:     ',A.bootstrap_percentsampled,'\n']);
            fprintf(fid,['Bootstrap Iterations:          ',A.bootstrap_iterations,'\n']);
            fprintf(fid,['Bootstrap Passes:              ',A.bootstrap_passes,'\n']);
        end
    else
        fprintf(fid,'None\n');
    end
    fprintf(fid,['Write Output?:                 ',y_n{str2double(A.write)+1},'\n']);
    if str2double(A.write)
        fprintf(fid,['Output Basename:               ',A.outbase,'\n']);
        fprintf(fid,['Save Add. Peak Info?:          ',y_n{str2double(A.savepeakinfo)+1},'\n']);
        if str2double(A.savepeakinfo)
            peaks={'1','1,2','1,2,3'};
            fprintf(fid,['Save Data for Peaks:           ',peaks{str2double(A.corrpeaknum)},'\n']);
            fprintf(fid,['Save Peak Magnitude?:          ',y_n{str2double(A.savepeakmag)+1},'\n']);
            fprintf(fid,['Save Resulting Vel.?:          ',y_n{str2double(A.savepeakvel)+1},'\n']);
        end
    end
end

fprintf(fid,'\n----------------------Images and Masking---------------------\n');
% fprintf(fid,['Image Directory: ',Data.imdirec,'\n']);
fprintf(fid,'Masking Type: ');
if strcmp(Data.masktype,'static')
    fprintf(fid,['Static\nMask File: ',Data.staticmaskname,'\n']);
elseif strcmp(Data.masktype,'dynamic')
    fprintf(fid,'Dynamic\n');
%     fprintf(fid,['Mask Directory: ',Data.maskdirec,'\n']);
else
    fprintf(fid,'None\n');
end
fprintf(fid,'Correlation Order:\n');
imagelist=get(handles.imagelist,'String');
if strcmp(Data.masktype,'dynamic')
    masklist=get(handles.masklist,'String');
    for i=1:size(imagelist,1)
        if ~isempty(imagelist{i})
            fprintf(fid,[imagelist{i},' --- ',masklist{i},'\n']);
        end
    end
else
    for i=1:size(imagelist,1)
        fprintf(fid,[imagelist{i},'\n']);
    end
end

fclose(fid);
