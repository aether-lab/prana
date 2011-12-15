function varargout = prana(varargin)
% HELP TEXT for running prana from the command line goes here.

%     This file is part of prana, an open-source GUI-driven program for
%     calculating velocity fields using PIV or PTV.
%     Copyright (C) 2010  Virginia Polytechnic Institute and State
%     University
% 
%     prana is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     (at your option) any later version.
% 
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
% 
%     You should have received a copy of the GNU General Public License
%     along with this program.  If not, see <http://www.gnu.org/licenses/>.

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @prana_OpeningFcn, ...
                   'gui_OutputFcn',  @prana_OutputFcn, ...
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
function prana_OpeningFcn(hObject, eventdata, handles, varargin)
warning off
handles.syscolor=get(hObject,'color');

verstr=version('-release');
if str2double(verstr(1:4))<2009
    errordlg('Prana requires Matlab R2009 or later.', 'prana')
end

try
    load('defaultsettings.mat')
catch
    defaultdata.clientversion='2.0';
    defaultdata.version='2.0';
    defaultdata.imbase='Img_';
    defaultdata.imzeros='6';
    defaultdata.imext='tif';
    defaultdata.imcstep='1';
    defaultdata.imfstep='2';
    defaultdata.imfstart='1';
    defaultdata.imfend='1';

    defaultdata.wrmag='1';
    defaultdata.wrsamp='1';
    defaultdata.wrsep='1';
    defaultdata.batchname='Proc1';
    defaultdata.datout='1';
    defaultdata.multiplematout='0';

    defaultdata.exp_date='';
    defaultdata.exp_L='';
    defaultdata.exp_v0='';
    defaultdata.exp_notes={'Camera Description:' '' 'Lens Description:' '' 'Notes:' ''};
    defaultdata.exp_density='1000';
    defaultdata.exp_viscosity='1.308e-3';
    defaultdata.exp_surfacetension='0.07197';
    defaultdata.exp_partD='';
    defaultdata.exp_partdensity='';
    defaultdata.exp_wavelength='.532';
    defaultdata.exp_pixelsize='';
    defaultdata.exp_lensfocal='';
    defaultdata.exp_lensfnum='';
    defaultdata.exp_micro='0';
    defaultdata.exp_NA='';
    defaultdata.exp_n='';
    defaultdata.exp_Re='';
    defaultdata.exp_St='';
    defaultdata.exp_M='';
    defaultdata.exp_ROI='';
    defaultdata.exp_diffractiondiameter='';
    defaultdata.exp_depthoffocus='';

    defaultdata.masktype='none';
    defaultdata.staticmaskname='';
    defaultdata.maskbase='maskfor_Img_';
    defaultdata.maskzeros='6';
    defaultdata.maskext='tif';
    defaultdata.maskfstep='1';
    defaultdata.maskfstart='1';

    defaultdata.PIV0.winres='32,32; 32,32';
    defaultdata.PIV0.winres1='32,32';
    defaultdata.PIV0.winres2='32,32';
    defaultdata.PIV0.winsize='64,64';
    defaultdata.PIV0.winauto='1';
    defaultdata.PIV0.gridres='8,8';
    defaultdata.PIV0.winoverlap='75,75';
    defaultdata.PIV0.gridtype='1';
    defaultdata.PIV0.gridbuf='8,8';
    defaultdata.PIV0.BWO='0,0';
    defaultdata.PIV0.corr='2';
    defaultdata.PIV0.RPCd='2.8';
    defaultdata.PIV0.zeromean='0';
    defaultdata.PIV0.peaklocator='1';
    defaultdata.PIV0.velsmooth='0';
    defaultdata.PIV0.velsmoothfilt='2';
    defaultdata.PIV0.val='0';
    defaultdata.PIV0.uod='1';
    defaultdata.PIV0.bootstrap='0';
    defaultdata.PIV0.thresh='0';
    defaultdata.PIV0.uod_type='2';
    defaultdata.PIV0.uod_window='3,3;3,3';
    defaultdata.PIV0.uod_thresh='3,2';
    defaultdata.PIV0.bootstrap_percentsampled='15';
    defaultdata.PIV0.bootstrap_iterations='700';
    defaultdata.PIV0.bootstrap_passes='12';
    defaultdata.PIV0.valuthresh='-16,16';
    defaultdata.PIV0.valvthresh='-16,16';
    defaultdata.PIV0.valextrapeaks='0';
    defaultdata.PIV0.savepeakinfo='0';
    defaultdata.PIV0.corrpeaknum='1';
    defaultdata.PIV0.savepeakmag='0';
    defaultdata.PIV0.savepeakvel='0';
    defaultdata.PIV0.outbase='PIV_';
    defaultdata.PIV0.write='1';

    defaultdata.PIV1=defaultdata.PIV0;
    defaultdata.PIV2=defaultdata.PIV0;

    defaultdata.passes='2';
    defaultdata.method='1';
    defaultdata.velinterp='3';
    defaultdata.iminterp='1';
    defaultdata.framestep='3';
    defaultdata.PIVerror='0.1';
    defaultdata.channel = '1';
    
    
    if ispc
%         defaultdata.loaddirec=[pwd '\'];
        defaultdata.ID.save_dir        = [pwd,'\ID\'];
        defaultdata.Size.save_dir      = [pwd,'\Size\'];
        defaultdata.Track.save_dir     = [pwd,'\Track\'];
        defaultdata.Track.PIVprops.load_dir= [pwd,'\'];
    else
%         defaultdata.loaddirec=[pwd '/'];
        defaultdata.ID.save_dir        = [pwd,'/ID/'];
        defaultdata.Size.save_dir      = [pwd,'/Size/'];
        defaultdata.Track.save_dir     = [pwd,'/Track/'];
        defaultdata.Track.PIVprops.load_dir= [pwd,'/'];
    end
    
    defaultdata.splash='1';
end

if ~isfield(defaultdata,'outputpassbasename')
    defaultdata.outputpassbase = 'PIV_';
    defaultdata.PIV0.outbase = [defaultdata.outputpassbase 'pass0_'];
    defaultdata.PIV1.outbase = [defaultdata.outputpassbase 'pass1_'];
    defaultdata.PIV2.outbase = [defaultdata.outputpassbase 'pass2_'];
end

if str2double(defaultdata.splash)==1 || str2double(defaultdata.clientversion)<2.0
    splash=splashdlg({...
         'What''s new in Prana v2.0?' ...
         '' ...
         'Particle Identification and Sizing' ...
         'Particle Tracking Velocimetry' ...
         },...
         'Prana v2.0','Ok','Don''t show this anymore','Ok');
    if strcmp(splash,'Don''t show this anymore')
        defaultdata.splash='0';
    end
    defaultdata.clientversion='0.99';
    defaultdata.version='2.0';
    save('defaultsettings.mat','defaultdata')
end
defaultdata.version=pranaPIVcode('version');
handles.data=defaultdata;
pranadir=which('prana');
addpath([pranadir(1:end-7),'documentation']);

if ispc
    handles.loaddirec=[pwd '\'];
    handles.data.ID.save_dir        = [pwd,'\ID\'];
    handles.data.Size.save_dir      = [pwd,'\Size\'];
    handles.data.Track.save_dir     = [pwd,'\Track\'];
    handles.data.Track.PIVprops.load_dir= [pwd,'\'];
else
    handles.data.ID.save_dir        = [pwd,'/ID/'];
    handles.data.Size.save_dir      = [pwd,'/Size/'];
    handles.data.Track.save_dir     = [pwd,'/Track/'];
    handles.data.Track.PIVprops.load_dir= [pwd,'/'];
    handles.loaddirec=[pwd '/'];
end

handles.data.par='0';
try
    compinfo=findResource('scheduler','configuration','local');
    handles.data.parprocessors=num2str(compinfo.clustersize);
catch
    handles.data.parprocessors='1';
end

handles.data.imdirec=pwd;
handles.data.maskdirec=pwd;
handles.data.outdirec=pwd;

try
    windowdiagram=imread(fullfile(pranadir(1:end-8),'documentation','windowdiagram.tif'),'tif');
catch
    windowdiagram=zeros(564,531);
end
set(gca,'children',imshow(windowdiagram));
axis off;

handles.data.cpass=num2str(get(handles.passlist,'Value'));
handles.data0=handles.data;
handles.Njob=num2str(size(get(handles.joblist,'String'),1));
handles.Cjob=num2str(get(handles.joblist,'String'));

handles=rmfield(handles,'data');
handles=update_data(handles);
handles.output = hObject;
guidata(hObject, handles);


% --- Outputs from this function are returned to the command line ---
function varargout = prana_OutputFcn(hObject, eventdata, handles) 
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
msgbox({'Prana','',['Client Version ',num2str(handles.data0.clientversion)],'',...
    ['Data Version ',num2str(handles.data0.clientversion)],'',...
    'Authors: Brady Drew, Adric Eckstein, John Charonko, Sam Raben, and the rest of the AEThER Lab','',...
    'Copyright 2010 - Virginia Polytechnic Institute and State University','',...
    'This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version. This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details. You should have received a copy of the GNU General Public License along with this program. If not, see www.gnu.org/licenses'},'About')

% --- Help Menu -> Getting Started ---
function gettingstarted_Callback(hObject, eventdata, handles)
web('gettingstarted.htm')

% --- Help Menu -> Help Topics ---
function helpmenu_helptopics_Callback(hObject, eventdata, handles)
PIVhelp

% --- Experiment Parameters Tab ---
function exptoggle_Callback(hObject, eventdata, handles)
set(hObject,'Value',1)
set(handles.imagetoggle,'Value',0)
set(handles.processingtoggle,'Value',0)
set(handles.particletoggle,'Value',0)
set(handles.exppanel,'Visible','on')
set(handles.imagepanel,'Visible','off')
set(handles.processingpanel,'Visible','off')
set(handles.particlepanel,'Visible','off')

% --- Image and Data I / O Tab ---
function imagetoggle_Callback(hObject, eventdata, handles)
set(hObject,'Value',1)
set(handles.exptoggle,'Value',0)
set(handles.processingtoggle,'Value',0)
set(handles.particletoggle,'Value',0)
set(handles.exppanel,'Visible','off')
set(handles.imagepanel,'Visible','on')
set(handles.processingpanel,'Visible','off')
set(handles.particlepanel,'Visible','off')

% --- PIV Processing Tab ---
function processingtoggle_Callback(hObject, eventdata, handles)
set(hObject,'Value',1)
set(handles.exptoggle,'Value',0)
set(handles.imagetoggle,'Value',0)
set(handles.particletoggle,'Value',0)
set(handles.exppanel,'Visible','off')
set(handles.imagepanel,'Visible','off')
set(handles.processingpanel,'Visible','on')
set(handles.particlepanel,'Visible','off')

% --- ID, Sizing & Tracking Tab ---
function particletoggle_Callback(hObject, eventdata, handles)
set(hObject,'Value',1)
set(handles.exptoggle,'Value',0)
set(handles.imagetoggle,'Value',0)
set(handles.processingtoggle,'Value',0)
set(handles.exppanel,'Visible','off')
set(handles.imagepanel,'Visible','off')
set(handles.processingpanel,'Visible','off')
set(handles.particlepanel,'Visible','on')

% --- Grid and Correlation Setup Tab ---
function gridsetuptoggle_Callback(hObject, eventdata, handles)
set(hObject,'Value',1)
set(handles.validationtoggle,'Value',0)
set(handles.validationpanel,'Visible','off')
set(handles.gridsetuppanel,'Visible','on')

% --- Validation and Output Tab ---
function validationtoggle_Callback(hObject, eventdata, handles)
set(hObject,'Value',1)
set(handles.gridsetuptoggle,'Value',0)
set(handles.gridsetuppanel,'Visible','off')
set(handles.validationpanel,'Visible','on')

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

% --- Correlate Image Pairs in Parallel Checkbox ---
function parcheckbox_Callback(hObject, eventdata, handles)
if str2double(handles.Njob)>0
    handles.data.par=num2str(get(hObject,'Value'));
    handles=load_data(handles);
    guidata(hObject,handles)
end
    

% --- Number of Processors to Use for Parallel Processing ---
function parprocessors_Callback(hObject, eventdata, handles)
if str2double(handles.Njob)>0
    handles.data.parprocessors=get(hObject,'String');
    handles=load_data(handles);
    guidata(hObject,handles)
end

function parprocessors_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Run Current Job Button ---
function runcurrent_Callback(hObject, eventdata, handles)
if str2double(handles.Njob)>0
    Data=handles.data;
    write_expsummary(Data,handles);
    pranaPIVcode(Data);
end

% --- Run All Jobs Button ---
function runall_Callback(hObject, eventdata, handles)
if str2double(handles.Njob)>0
    Jlist=char(get(handles.joblist,'String'));
    eval(['handles.' Jlist(str2double(handles.Cjob),:) '=handles.data;']);
    
    for e=1:size(Jlist,1)
        Data=eval(['handles.' Jlist(e,:)]);
        write_expsummary(Data,handles);
        pranaPIVcode(Data);
    end
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
    Jlist = char(get(handles.joblist,'String'));
    eval(['handles.' Jlist(str2double(handles.Cjob),:) '=handles.data;']);
end
handlesP = handles;

% Select job to load
%  f is the name of the saved job file. This variable should be renamed
%  something more descriptive.
[f, handles.loaddirec]=uigetfile('*.mat','LOAD JOB',handles.loaddirec,'multiselect','on');

if ischar(f) == 1
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
                    %Attempt to make backwards-compatible with older
                    %versions of prana
                    if ~isfield(Data,'version')
                        if ~isfield(Data,'par')
                            Data.par='0';
                            Data.parprocessors='1';
                            Data.version='1.5';
                            
                        else
                            handles.data.version='1.9';
                        end
                    end
                    if ~isfield(Data.PIV0,'zeromean')
                        for pass=0:str2double(Data.passes)
                            eval(['Data.PIV',num2str(pass),'.zeromean=''0'';']);
                            eval(['Data.PIV',num2str(pass),'.peaklocator=''1'';']);
                        end
                    end
                    
                    % This performs a check to see if the job files
                    % contains the field 'outputpassbase' if not then it
                    % used the output name from the final pass.
                    if ~isfield(Data,'outputpassbase')
                        ll=1; bb = 0;                        
                        while bb == 0 
                            if isfield(Data,['PIV' num2str(ll)])
                                ll = ll+1;
                            else
                                bb = 1;
                            end
                        end
                        eval(['Data.outputpassbase = Data.PIV' num2str(ll-1) '.outbase;']);
                    end

                    if ~isfield(Data,'ID')
                        Data.runPIV = '1';
                        
                        load defaultsettings.mat defaultdata
                        Data.ID=defaultdata.ID;
                        Data.Size=defaultdata.Size;
                        Data.Track=defaultdata.Track;
                        
                        if ispc
                            Data.ID.save_dir        = [Data.outdirec,'\ID\'];
                            Data.Size.save_dir      = [Data.outdirec,'\Size\'];
                            Data.Track.save_dir     = [Data.outdirec,'\Track\'];
                            Data.Track.PIVprops.load_dir       = [Data.outdirec,'\'];
                        else
                            Data.ID.save_dir        = [Data.outdirec,'/ID/'];
                            Data.Size.save_dir      = [Data.outdirec,'/Size/'];
                            Data.Track.save_dir     = [Data.outdirec,'/Track/'];
                            Data.Track.PIVprops.load_dir       = [Data.outdirec,'/'];
                        end
                    end
                                        
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

                if vn ~= -1
                    if str2double(handles.Njob) > 0
                        Jlist = char(get(handles.joblist,'String'));
                        Jlist = {Jlist;Data.batchname};
                    else
                        Jlist = {Data.batchname};
                    end
                    handles.Njob = num2str(str2double(handles.Njob)+1);
                    handles.Cjob = handles.Njob;
                    set(handles.joblist,'String',Jlist,'Value',str2double(handles.Cjob));
                    
                    
                    handles = update_data(handles);
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
    Data = handles.data;

       uisave('Data',[handles.data.imdirec handles.data.batchname '.mat']);
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
%     if strcmp(handles.data.masktype,'dynamic')
%         load_masklist(handles)
%     end
    load_imlist(handles);
    guidata(hObject,handles)
end

% --- Image Directory Text Box ---
function imagedirectory_Callback(hObject, eventdata, handles)
if str2double(handles.Njob)>0
    handles.data.imdirec = get(hObject,'String');
%     if strcmp(handles.data.masktype,'dynamic')
%         load_masklist(handles)
%     end
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
%     if strcmp(handles.data.masktype,'dynamic')
%         load_masklist(handles)
%     end
    load_imlist(handles);
    guidata(hObject,handles)
end

% --- Compute IMmin Button ---
function immin_Callback(hObject, eventdata, handles)
if str2double(handles.Njob)>0
    try
        h = waitbar(0,'1','Name','Computing Minimum of Test Images','CreateCancelBtn','setappdata(gcbf,''canceling'',1)');setappdata(h,'canceling',0);
        e=0;
        for im=str2double(handles.data.imfstart):str2double(handles.data.imfend)
            if mod(im,2)
                waitbar((im-str2double(handles.data.imfstart))/length(str2double(handles.data.imfstart):str2double(handles.data.imfend)),h,'Analyzing images...')
            end
            if getappdata(h,'canceling')
                e=-1;
                break
            end
            im_read=double(imread(fullfile(handles.data.imdirec, [handles.data.imbase, sprintf(['%0.' handles.data.imzeros 'i.' handles.data.imext],im)])));
            if im==str2double(handles.data.imfstart)
                min_im=ones(size(im_read)) .* 255;
            end
            min_im(im_read<min_im)=im_read(im_read<min_im);
        end
        delete(h)
    catch
        msgbox('Error Loading Images');
        e=-1;
    end
    if e==0
        try 
            imwrite(uint8(min_im),fullfile(handles.data.imdirec,'IMmin.tif'),'tif');
        catch
            e=-1;
            msgbox('Error Writing IMmin.tif')
        end
    end
    
    if e==0
        figure
        if gcf~=gcbf
            imagesc(min_im,[0 255]),colormap gray,axis off
        end
    end
end

% --- Image Basename Text Box ---
function imagebasename_Callback(hObject, eventdata, handles)
if str2double(handles.Njob)>0
    handles.data.imbase = get(hObject,'String');
%     if strcmp(handles.data.masktype,'dynamic')
%         load_masklist(handles)
%     end
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
%     if strcmp(handles.data.masktype,'dynamic')
%         load_masklist(handles)
%     end
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
%     if strcmp(handles.data.masktype,'dynamic')
%         load_masklist(handles)
%     end
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
%     if strcmp(handles.data.masktype,'dynamic')
%         load_masklist(handles)
%     end
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
%     if strcmp(handles.data.masktype,'dynamic')
%         load_masklist(handles)
%     end
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
%     if strcmp(handles.data.masktype,'dynamic')
%         load_masklist(handles)
%     end
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
%     if strcmp(handles.data.masktype,'dynamic')
%         load_masklist(handles)
%     end
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
%     load_masklist(handles);
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

    handles.data.maskdirec = uigetdir(handles.data.imdirec);
    if handles.data.maskdirec==0
        handles.data.maskdirec = D;
    end
    set(handles.maskdirectory,'string',handles.data.maskdirec);
%     load_masklist(handles);
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
%     load_masklist(handles);
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
%     load_masklist(handles);
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
%     load_masklist(handles);
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
%     load_masklist(handles);
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
%     load_masklist(handles);
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

    [a,b] = uigetfile([handles.data.imdirec '/*.*'], 'Select static mask file...');
    handles.data.staticmaskname = [b a];
    if handles.data.staticmaskname==0
        handles.data.staticmaskname = D;
    end
    guidata(hObject,handles)
    set(handles.staticmaskfile,'string',handles.data.staticmaskname);
end

% --- Static Mask Tool Button ---
function masktool_Callback(hObject, eventdata, handles)
if str2double(handles.Njob)>0
    try
        im = double(imread(fullfile(handles.data.imdirec, [handles.data.imbase, sprintf(['%0.' handles.data.imzeros 'i.' handles.data.imext],str2double(handles.data.imfstart))])));
        
        channel = str2double(handles.data.channel);
        
         if size(im, 3) > 1
%              Extract only red channel
             if channel == 1;
                im = im(:,:,1);
%                 Extract only green channel
             elseif channel == 2;
                im = im(:,:,2);
%                 Extract only blue channel
             elseif channel == 3;
                im = im(:,:,3);
%                 Weighted average of channels (see rgb2gray for
%                 explanation of weighting factors)
             elseif channel == 4;
                im = 0.2989 * im(:, :, 1) + 0.5870 * im(:, :, 2) + 0.1140 * im(:, :, 3);
%                 Evenly weighted mean of channels
             elseif channel == 5;
                im = (im(:,:,1) + im(:,:,2) + im(:,:,3))/3;
             elseif channel == 6
                 im = im;%#ok
             end
         else
%             Take only red channel
            im =im(:,:,1);
         end
         
        mask=ones(size(im,1),size(im,2));
        roiwindow = CROIEditor(im./max(im(:)));
        while isempty(roiwindow.labels)
        addlistener(roiwindow,'MaskDefined',@your_roi_defined_callback);
        drawnow
        end
        
        mask(roiwindow.labels==0) = 0;
        handles.data.staticmaskname=fullfile(handles.data.imdirec,'staticmask.tif');
        imwrite(mask,handles.data.staticmaskname,'tif')
        set(handles.staticmaskfile,'String',handles.data.staticmaskname);
        handles.data.masktype='static';
        set(handles.staticmaskbutton,'Value',1)
        set(handles.dynamicmaskbutton,'Value',0)
        set(handles.nomaskbutton,'Value',0)
        Jlist=char(get(handles.joblist,'String'));
        eval(['handles.' Jlist(str2double(handles.Cjob),:) '=handles.data;']);
        handles=update_data(handles);
        guidata(hObject,handles)

        close('Analyzer - ROI Editor')
        
%         stillmasking='Yes';mask=ones(size(im,1),size(im,2));h=figure;
%         while strcmp(stillmasking,'Yes')
%             figure(h),imshow(im.*repmat(mask,[1 1 size(im,3)]),[0 255])
%             masktemp=roipoly;
%             mask(masktemp==1)=0;
%             figure(h),imshow(im.*repmat(mask,[1 1 size(im,3)]),[0 255])
%             stillmasking=questdlg('Mask another region?','Region Masked','Yes','No','No');
%         end
%         close(h)
%         if ~isempty(stillmasking)
%             handles.data.staticmaskname=fullfile(handles.data.imdirec,'staticmask.tif');
%             inc_check = get(handles.mask_inclusion,'Value');
%             imwrite(abs(inc_check-mask),handles.data.staticmaskname,'tif')
%             set(handles.staticmaskfile,'String',handles.data.staticmaskname);
%             handles.data.masktype='static';
%             set(handles.staticmaskbutton,'Value',1)
%             set(handles.dynamicmaskbutton,'Value',0)
%             set(handles.nomaskbutton,'Value',0)
%             Jlist=char(get(handles.joblist,'String'));
%             eval(['handles.' Jlist(str2double(handles.Cjob),:) '=handles.data;']);
%             handles=update_data(handles);
%             guidata(hObject,handles)
%         else
%             msgbox('Masking Cancelled')
%         end
    catch ME
        msgbox('Image Frame Not Found');
        e=-1;
    end

end
function your_roi_defined_callback(h,e)
[mask, labels, n] = roiwindow.getROIData;
delete(roiwindow);

% --- Preview Image + Mask Button ---
function impreview_Callback(hObject, eventdata, handles)
if str2double(handles.Njob)>0
    
    % Read color channel
    channel = str2double(handles.data.channel);
    
    try
        im = double(imread(fullfile(handles.data.imdirec, [handles.data.imbase, sprintf(['%0.' handles.data.imzeros 'i.' handles.data.imext],str2double(handles.data.imfstart))])));
        if size(im, 3) > 1
            %   Extract only red channel
            if channel == 1;
                im = im(:,:,1);
                %   Extract only green channel
            elseif channel == 2;
                im = im(:,:,2);
                %	Extract only blue channel
            elseif channel == 3;
                im = im(:,:,3);
                % 	Weighted average of channels (see rgb2gray for
                %	explanation of weighting factors)
            elseif channel == 4;
                im = 0.2989 * im(:, :, 1) + 0.5870 * im(:, :, 2) + 0.1140 * im(:, :, 3);
                %	Evenly weighted mean of channels
            elseif channel == 5;
                im = (im(:,:,1) + im(:,:,2) + im(:,:,3))/3;
            elseif channel == 6;
                im = im;%#ok
            end
        else
            %   Take only red channel
            im =im(:,:,1);
        end

    %   Flip and normalize image
    im1 = im(end:-1:1,:,:)./255;

        try
            if strcmp(handles.data.masktype,'static')
                mask = double(imread(handles.data.staticmaskname));
                mask = flipud(mask);
            elseif strcmp(handles.data.masktype,'dynamic')
                mask = double(imread(fullfile(handles.data.maskdirec, [handles.data.maskbase, sprintf(['%0.' handles.data.maskzeros 'i.' handles.data.maskext],str2double(handles.data.maskfstart))])));
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
        h=figure;hold on

        L=size(im1);
        Ps=get(0,'screensize');
        if (L(1)/Ps(4))>(L(2)/Ps(3))
            dy=0.8*Ps(4);
            dx=dy*(L(2)/L(1));
        else
            dx=0.8*Ps(3);
            dy=dx*(L(1)/L(2));
        end
        figure(h)
        set(gcf,'position',[(Ps(3)-dx)/2 (Ps(4)-dy)/2 dx dy])
        set(gcf,'color',0.5*[1 1 1])
        imagesc(im1,[0 1]),axis image,colormap gray,axis off,set(gca,'position',[0 0 1 1])

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

    handles.data.outdirec = uigetdir(handles.data.imdirec);
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
    handles.data.passes=num2str(N+1);
    if str2double(handles.data.method)==1
        eval(['handles.data=setfield(handles.data,''PIV' num2str(N+1) ''',handles.data.PIV1);']);
    else
        eval(['handles.data=setfield(handles.data,''PIV' num2str(N+1) ''',handles.data.PIV0);']);
    end
    passbase = get(handles.outputpassbasename,'String');
    eval(['handles.data.PIV' num2str(N+1) '.outbase=[''' passbase 'pass' num2str(N+1) '_''];']);
    load_PIVlist(handles);
    handles=set_PIVcontrols(handles);
    guidata(hObject,handles)
end

% --- Copy Pass Button ---
function copypassbutton_Callback(hObject, eventdata, handles)
if str2double(handles.Njob)>0
    N=str2double(handles.data.passes);
    handles.data.passes=num2str(N+1);
    cpass = get(handles.passlist,'Value'); % This grabs the currently selected pass number
    passbase = get(handles.outputpassbasename,'String');
    for i = N+1:-1:cpass+1
        eval(['handles.data=setfield(handles.data,''PIV' num2str(i) ''',handles.data.PIV' num2str(i-1) ');']);
        eval(['handles.data.PIV' num2str(i) '.outbase=[''' passbase 'pass'  num2str(i) '_''];']);
    end
    eval(['handles.data=setfield(handles.data,''PIV' num2str(cpass+1) ''',handles.data.PIV' num2str(cpass) ');']);
    eval(['handles.data.PIV' num2str(cpass+1) '.outbase=[''' passbase 'pass'  num2str(cpass+1) '_''];']);
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
        passbase = get(handles.outputpassbasename,'String');
        for e=(p+1):str2double(handles.data.passes)
            eval(['handles.data=setfield(handles.data,''PIV' num2str(e-1) ''',handles.data.PIV' num2str(e) ');']);
            eval(['handles.data.PIV' num2str(e-1) '.outbase=[''' passbase 'pass' num2str(e-1) '_''];']);
            if e==str2double(handles.data.passes)
                eval(['handles.data=rmfield(handles.data,''PIV' num2str(e) ''');']);
            end
        end
        handles.data.passes=num2str(str2double(handles.data.passes)-1);
        set(handles.passlist,'Value',1);
        load_PIVlist(handles)
        guidata(hObject,handles)
    end
    if str2double(handles.data.passes) < 2
        set(handles.writeoutputcheckbox, 'Value', 1);
        handles.data.PIV1.write = '1';
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
        set(handles.velocityinterptype,'backgroundcolor',[1 1 1]);
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
PIVhelp(8)

% --- Window Resolution Text Box ---
function windowres_Callback(hObject, eventdata, handles)
%  If a job has been created or loaded...
if str2double(handles.Njob) > 0
        
%    Read the string in the window resolution textbox
    A = get(hObject,'String');

% Parse string in window resolution textbox to determine the resolutions of each window
    [wx1 wy1 wx2 wy2] = parseNum(A);
    
% Set "window resolution" field in the "PIV pass" data structure to the string in the winres textbox
    eval(['handles.data.PIV' handles.data.cpass '.winres = [num2str(wx1) '','' num2str(wy1) '';'' num2str(wx2) '','' num2str(wy2)];'])

% Determine the largest window dimensions for the pass
    wxMax = max([wx1 wx2]);
    wyMax = max([wy1 wy2]);
    
% Check whether actual window size should be determined automatically
    if get(handles.autowinsizecheckbox, 'Value') == 1
        
% Calculate x- and y- window sizes (pixels)
        [xbin ybin] = roiSize(A);
        
% Set "window size" fields in "PIV pass" data structure to the window sizes calculated above
        eval(['handles.data.PIV' handles.data.cpass '.winsize = [num2str(xbin) '','' num2str(ybin)];'])
        eval(['set(handles.windowsize,''String'', handles.data.PIV' handles.data.cpass '.winsize)']) 
    end
    
% Alert user if window size is smaller than 256 pixels. Not sure why this is important.
    if wxMax * wyMax < 256
        set(hObject,'backgroundcolor',[1 0.5 0]);
    else
        set(hObject,'backgroundcolor',[1 1 1]);
    end
    
% If specifying grid resolution...
    if get(handles.setgridresbutton,'Value') == 1
        
% Read string in "grid resolution" text box 
        A = get(handles.gridres,'String');
        
% Parse string in "grid resolution" text box to determine x- and y- grid resolutions
         [gx gy] = parseNum(A);
        
% Calculate window overlaps (%)
        overX=(wxMax - gx) / wxMax * 100;
        overY=(wyMax - gy) / wyMax * 100;
        
%  Set overlaps to zero if overlaps are calculated to be less than zero
        overX = overX * (overX > 0);
        overY = overY * (overY > 0);
        
% Set "overlap" fields in "PIV pass" data structure to the overlaps calculated above
        eval(['handles.data.PIV' handles.data.cpass '.winoverlap = [num2str(overX),'','',num2str(overY)];'])
        
% Otherwise, if specifying window overlap....
    else
        
% Read string in "window overlap" textbox
        A = get(handles.winoverlap,'String');
        
% Parse string in "window overlap" textbox to determine the x- and y- window overlaps (%)
        [overX overY] = parseNum(A);
        
% Calculate x- and y- grid resolutions from overlaps 
        gx = round(wxMax * (1-overX/100));
        gy = round(wyMax * (1-overY/100));
        
% Update "grid resolution" field in "PIV pass" data structure to the grid resolutions calculated above
        eval(['handles.data.PIV' handles.data.cpass '.gridres = [num2str(gx),'','',num2str(gy)];'])
    end
    
%  Update "window resolution" fields in "PIV pass" data structure to the
%  window resolutions calculated above
    eval(['handles.data.PIV' handles.data.cpass '.winres1 = [num2str(wx1) '','' num2str(wy1)];']); 
    eval(['handles.data.PIV' handles.data.cpass '.winres2 = [num2str(wx2) '','' num2str(wy2)];']); 
    
% Update GUI
    handles=set_PIVcontrols(handles);
    guidata(hObject,handles)
    
end


function windowres_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Window Size Text Box ---
function windowsize_Callback(hObject, eventdata, handles)
% If a job has been loaded or created....
if str2double(handles.Njob ) > 0

% If the "auto window size" checkbox is NOT checked (i.e. if the window size is specified explicitly)
    if get(handles.autowinsizecheckbox,'Value') == 0
        
% Read the string in the "Actual Window Size" text box
        A = get(hObject,'String');
        
% Determine the dimensions of the interrogation region by parsing the string in the "Actual Window Size" text box
        [Rx Ry] = parseNum(A);

% Read the string in the "Window Resolution" text box
        B = get(handles.windowres,'String');
        
% Determine the window resolutions by parsing the string in the "Window Resolution" text box
        [wx1 wy1 wx2 wy2] = parseNum(B);        

% Determine largest window resolutions
        wxMax = max(wx1, wx2);
        wyMax = max(wy1, wy2);

%  If the x-dimension of the interrogation region is smaller than that of the effective window resolution,
%  re-size the x-dimension of the interrogation regions to equal that of the effective window resolution
        if Rx < wxMax
            Rx = wxMax;
            set(hObject,'String',[num2str(Rx) ',' num2str(Ry)]);
        end
        
%  If the x-dimension of the interrogation region is smaller than that of the effective window resolution,
%  re-size the x-dimension of the interrogation regions to equal that of the effective window resolution
        if Ry < wyMax
            Ry = wyMax;
            set(hObject,'String',[num2str(Rx) ',' num2str(Ry)]);
        end
        
% Update the "window size" field of the "PIV Pass" data structure to the value of the 
% string in the "Actual Window Size" text box
eval(['handles.data.PIV' handles.data.cpass '.winsize = [num2str(Rx) '','' num2str(Ry)];']);
 
% If the "auto window size" checkbox IS checked (i.e. if the window size is calculated automatically)
    else
        
% Set the text in the "Actual Window Size" text box to the value stored in
% the "PIV Pass" data structure
        eval(['set(hObject,''String'',handles.data.PIV' handles.data.cpass '.winsize);'])
       
    end
    % Update the GUI
        handles = set_PIVcontrols(handles);
        guidata(hObject,handles)
end

function windowsize_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Auto-Size Window Check Box ---
function autowinsizecheckbox_Callback(hObject, eventdata, handles)
% If a job has been created or specified...
if str2double(handles.Njob) > 0
  
    A = get(handles.windowres, 'String');
    [xROI yROI] = roiSize(A);
    
   eval(['handles.data.PIV' handles.data.cpass '.winsize = [num2str(xROI) '','' num2str(yROI)];']); 
    
    eval(['handles.data.PIV' handles.data.cpass '.winauto = num2str(get(hObject,''Value''));']);
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

% If a job has been specified or loaded....
if str2double(handles.Njob) > 0
    
%  If performing a multi-pass job...
    if get(handles.setgridresbutton,'Value') == 1
        
%  Read string in grid resolution texbox
        A = get(hObject,'String');
        
% Parse string in "grid resolution" text box to determine x- and y- grid resolutions
        [gx gy] = parseNum(A);

%  If performing a multipass run....
        if str2double(handles.data.method) == 1
            
% Then calculate the window overlaps for each pass...
            for e=1:str2double(handles.data.passes)
                
% Read window resolutions from data structure
                eval(['A=handles.data.PIV', num2str(e) '.winres1;'])
                eval(['B=handles.data.PIV', num2str(e) '.winres2;'])
                
%  Parse the strings containing the window resolution information for each image
                [wx1 wy1] = parseNum(A);
                [wx2 wy2] = parseNum(B);

%  Determine the size of the largest window resolution in each dimension
                wxMax = max(wx1, wx2);
                wyMax = max(wy1, wy2);
                
% Calculate the x- and y- window overlaps
                overX = (wxMax - gx) / wxMax * 100;
                overY = (wyMax - gy) / wyMax * 100;
                
%  Set overlaps to zero if "overlap" is calculated to be less than zero
            overX = overX * (overX > 0);
            overY = overY * (overY > 0);
                
%  Update the "grid resolution" and "window overlap" fields in the "pass" data structure
                eval(['handles.data.PIV' num2str(e) '.gridres = [num2str(gx) '','' num2str(gy)];'])
                eval(['handles.data.PIV' num2str(e) '.winoverlap = [num2str(overX),'','',num2str(overY)];'])
                
%  End of "multipass" case...
        end
            
%  Otherwise, If NOT performing a multipass run...  
        else
       
% Read string in the "Grid Resolution" text box;
        A = get(hObject, 'String');
        [gx gy] = parseNum(A);
            
% Read string in window resolution textbox
            B = get(handles.windowres,'String');
            
% Parse string in window resolution textbox to determine the resolutions of each window
            [wx1 wy1 wx2 wy2] = parseNum(B);
            
% Determine the size of the largest window resolution in each dimension
            wxMax = max(wx1, wx2);
            wyMax = max(wy1, wy2);
            
% Calculate the x- and y- window overlaps
            overX = (wxMax - gx) / wxMax * 100;
            overY = (wyMax - gy) / wyMax * 100;
            
%  Set overlaps to zero if "overlap" is calculated to be less than zero
            overX = overX * (overX > 0);
            overY = overY * (overY > 0);
            
%  Update the "grid resolution" and "window overlap" fields in the "pass" data structure   
            eval(['handles.data.PIV' handles.data.cpass '.gridres = [num2str(gx) '','' num2str(gy)];'])
            eval(['handles.data.PIV' handles.data.cpass '.winoverlap = [num2str(overX),'','',num2str(overY)];'])
        end
        
% Update GUI
        handles=set_PIVcontrols(handles);
        guidata(hObject,handles)
    else
        
%  If specifying window overlap rather than grid resolution, just read the
%  grid resolution. 
        eval(['set(hObject,''String'',handles.data.PIV' handles.data.cpass '.gridres);'])
    
    end
    
end


function gridres_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Set Window Overlap Button ---
function setwinoverlapbutton_Callback(hObject, eventdata, handles)
if str2double(handles.Njob) > 0 && str2double(handles.data.method)~=1
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

% If a job has been loaded or created...
if str2double(handles.Njob) > 0

% If specifying window overlap...
    if get(handles.setwinoverlapbutton,'Value') == 1

% Read string in "window resolution" text box
        A=get(handles.windowres,'String');

% Parse string in "window resolution" text box to determine the desired window resolutions
        [wx1 wy1 wx2 wy2] = parseNum(A);
        
% Determine sizes of largest window dimensions
        wxMax = max(wx1, wx2);
        wyMax = max(wy1, wy2);

% Read string in "Window Overlap" text box
        B = get(hObject,'String');
        
% Parse string in "Window Overlap" text box to determine the desired window overlaps
        [overX overY] = parseNum(B);
        
% Calculate grid resolution from specified overlaps
        gx = round(wxMax * (1 - overX / 100));
        gy = round(wyMax * (1 - overY / 100));
        
% Update "Window Overlap" field in "Piv Pass" data structure
        eval(['handles.data.PIV' handles.data.cpass '.winoverlap = get(hObject,''String'');'])
        
% Update "Grid Resolution" field in "Piv Pass" data structure
        eval(['handles.data.PIV' handles.data.cpass '.gridres = [num2str(gx),'','',num2str(gy)];'])
        
% Update GUI
        handles=set_PIVcontrols(handles);
        guidata(hObject,handles)
        
% Otherwise, if specifying grid resolution... 
    else
        
% Update the "Window Overlap" field in the "PIV Pass" data structure to the value calculated automatically
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

% --- Subpixel Correlation Peak Location Drop-down Menu ---
function subpixelinterp_Callback(hObject, eventdata, handles)
if str2double(handles.Njob)>0
    eval(['handles.data.PIV' handles.data.cpass '.peaklocator = num2str(get(hObject,''Value''));'])
    guidata(hObject,handles)
end

function subpixelinterp_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Zero-Mean Windows Checkbox ---
function zeromeancheckbox_Callback(hObject, eventdata, handles)
if str2double(handles.Njob)>0
    eval(['handles.data.PIV' handles.data.cpass '.zeromean = num2str(get(hObject,''Value''));'])
    guidata(hObject,handles)
end

% --- Smoothing Check Box ---
function smoothingcheckbox_Callback(hObject, eventdata, handles)
if str2double(handles.Njob)>0
    eval(['handles.data.PIV' handles.data.cpass '.velsmooth = num2str(get(hObject,''Value''));'])
    eval(['smoothingsize=handles.data.PIV' handles.data.cpass '.velsmoothfilt;'])
    if get(hObject,'Value')==1 && get(handles.passtype,'Value')>1
        set(handles.smoothingsize,'string',smoothingsize,'backgroundcolor',[1 1 1]);
    else
        set(handles.smoothingsize,'string',smoothingsize,'backgroundcolor',0.5*[1 1 1]);
    end
    guidata(hObject,handles)
end

% --- Smoothing Size Text Box ---
function smoothingsize_Callback(hObject, eventdata, handles)
if str2double(handles.Njob)>0 && get(handles.smoothingcheckbox,'Value')==1
    eval(['handles.data.PIV' handles.data.cpass '.velsmoothfilt = ''' get(hObject,'String') ''';'])
    set(hObject,'String',eval(['handles.data.PIV' num2str(handles.data.cpass) '.velsmoothfilt']));
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
PIVhelp(11)

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
        if strcmp(handles.data.passes, handles.data.cpass);
             set(hObject, 'Value', 1);
        end
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
    [handles]=load_data(handles);
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
    [handles]=load_data(handles);
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
    [handles]=load_data(handles);
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
    [handles]=load_data(handles);
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
function handles = update_data(handles)
set(handles.exp_Re,'backgroundcolor',0.8*[1 1 1]);
set(handles.exp_St,'backgroundcolor',0.8*[1 1 1]);
set(handles.exp_M,'backgroundcolor',0.8*[1 1 1]);
set(handles.exp_ROI,'backgroundcolor',0.8*[1 1 1]);
set(handles.exp_diffractiondiameter,'backgroundcolor',0.8*[1 1 1]);
set(handles.exp_depthoffocus,'backgroundcolor',0.8*[1 1 1]);
vernum = pranaPIVcode('version');
set(handles.version_box,'String',vernum,'backgroundcolor',0.5*[1 1 1]);
if str2double(handles.Njob) == 0
    set(handles.passlist,'String','','backgroundcolor',0.5*[1 1 1]);
    set(handles.joblist,'String','','backgroundcolor',0.5*[1 1 1]);
    set(handles.parprocessors,'backgroundcolor',0.5*[1 1 1]);
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
    set(handles.outputpassbasename,'String','','backgroundcolor',0.5*[1 1 1]);
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
%     set(handles.currentjobname,'String','','backgroundcolor',0.5*[1 1 1]);
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
    set(handles.subpixelinterp,'Value',1,'backgroundcolor',0.5*[1 1 1]);
    set(handles.velocityinterptype,'Value',1,'backgroundcolor',0.5*[1 1 1]);    
    set(handles.staticmaskfile,'String','','backgroundcolor',0.5*[1 1 1]);
    set(handles.maskdirectory,'String','','backgroundcolor',0.5*[1 1 1]);
    set(handles.maskbasename,'String','','backgroundcolor',0.5*[1 1 1]);
    set(handles.maskzeros,'String','','backgroundcolor',0.5*[1 1 1]);
    set(handles.maskextension,'String','','backgroundcolor',0.5*[1 1 1]);
%     set(handles.maskframestep,'String','','backgroundcolor',0.5*[1 1 1]);
%     set(handles.maskframestart,'String','','backgroundcolor',0.5*[1 1 1]);
%     set(handles.masklist,'String','','backgroundcolor',0.5*[1 1 1]);
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
    set(handles.colorchannel_popupMenu,'Value',1,'backgroundcolor',0.5*[1 1 1]);
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
	set(handles.outputpassbasename,'String','','backgroundcolor',[1 1 1]);
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
%     set(handles.currentjobname,'String','','backgroundcolor',[1 1 1]);
    set(handles.pulseseparation,'String','','backgroundcolor',[1 1 1]);
    set(handles.samplingrate,'String','','backgroundcolor',[1 1 1]);
    set(handles.smoothingsize,'String','','backgroundcolor',[1 1 1]);
    set(handles.imageinterptype,'Value',1,'backgroundcolor',[1 1 1]);
    set(handles.passtype,'Value',1,'backgroundcolor',[1 1 1]);
    set(handles.uod_type,'Value',1,'backgroundcolor',[1 1 1]);
    set(handles.correlationtype,'Value',1,'backgroundcolor',[1 1 1]);
    set(handles.subpixelinterp,'Value',1,'backgroundcolor',[1 1 1]);  
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
    set(handles.colorchannel_popupMenu,'Value',1,'backgroundcolor',[1 1 1]);

    if str2double(handles.data.par)==1
        set(handles.parprocessors,'string','','backgroundcolor',[1 1 1])
    else
        set(handles.parprocessors,'string','','backgroundcolor',0.5*[1 1 1])
    end
    
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
%         set(handles.maskframestep,'String','','backgroundcolor',[1 1 1]);
%         set(handles.maskframestart,'String','','backgroundcolor',[1 1 1]);
%         set(handles.masklist,'String','','backgroundcolor',[1 1 1]);
%         load_masklist(handles);
    else
        set(handles.maskdirectory,'String','','backgroundcolor',0.5*[1 1 1]);
        set(handles.maskbasename,'String','','backgroundcolor',0.5*[1 1 1]);
        set(handles.maskzeros,'String','','backgroundcolor',0.5*[1 1 1]);
        set(handles.maskextension,'String','','backgroundcolor',0.5*[1 1 1]);
%         set(handles.maskframestep,'String','','backgroundcolor',0.5*[1 1 1]);
%         set(handles.maskframestart,'String','','backgroundcolor',0.5*[1 1 1]);
%         set(handles.masklist,'backgroundcolor',0.5*[1 1 1]);
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

% % --- Load Mask List ---
% function load_masklist(handles)
% dir_struct = dir(handles.data.maskdirec);
% if isempty(dir_struct)
%     set(handles.maskdirectory,'backgroundcolor','r');
% else
%     set(handles.maskdirectory,'backgroundcolor',[1 1 1]);
% end
% 
% maskfend=str2double(handles.data.maskfstart)+str2double(handles.data.maskfstep)*length(str2double(handles.data.imfstart):str2double(handles.data.imfstep):str2double(handles.data.imfend))-1;
% N=length(str2double(handles.data.maskfstart):str2double(handles.data.maskfstep):maskfend);
% files = cell(N,1);
% e=0;
% for f=str2double(handles.data.maskfstart):str2double(handles.data.maskfstep):maskfend
%     e=e+1;
%     files(e,1)={[handles.data.maskbase sprintf(['%0.' handles.data.maskzeros 'i.' handles.data.maskext],f)]};
% end
% 
% [sorted_names,sorted_index] = sortrows({dir_struct.name}');
% [files1,id,id1] = intersect(sorted_names,files(:,1));
% 
% if isempty(id1)
%      set(handles.masklist,'backgroundcolor','r');
%     set(handles.masklist,'UserData',[]);
% else
%     set(handles.masklist,'backgroundcolor',[1 1 1]);
%     if length(id1)~=size(files,1)
%         set(handles.masklist,'backgroundcolor','r');
%         set(handles.masklist,'UserData',[]);
%     end
% end
% files=files(id1,:);
% filesf=cell(length(files),1);
% for e=1:size(files,1)
%     filesf(e)={[char(files(e,1))]};
% end
% set(handles.masklist,'String',filesf,'Value',1);

% --- Load PIV Pass List ---
function load_PIVlist(handles)
f=cell(str2double(handles.data.passes),1);
for e=1:str2double(handles.data.passes)
    f(e)={['Pass ' num2str(e)]};
end
set(handles.passlist,'String',f);
if get(handles.passlist,'Value')>str2double(handles.data.passes)
    set(handles.passlist,'Value',str2double(handles.data.passes))
end

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
set(handles.subpixelinterp,'Value',str2double(A.peaklocator));
set(handles.zeromeancheckbox,'Value',str2double(A.zeromean));
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

if N>1
    set(handles.bulkwinoffset,'backgroundcolor',0.5*[1 1 1]);
else
    set(handles.bulkwinoffset,'backgroundcolor',[1 1 1]);
end

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

if str2double(A.winauto) == 1
     B=get(handles.windowsize,'String');
     [Rx Ry] = parseNum(B);
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
if get(handles.correlationtype,'Value')>=2
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
[wx1 wy1 wx2 wy2] = parseNum(a);
wxMax = max(wx1, wx2);
wyMax = max(wy1, wy2);
if wxMax * wyMax<256
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
set(handles.parprocessors,'String',handles.data.parprocessors);
set(handles.parcheckbox,'Value',str2double(handles.data.par));
if str2double(handles.data.par)==1
    set(handles.parprocessors,'backgroundcolor',[1 1 1])
else
    set(handles.parprocessors,'backgroundcolor',0.5*[1 1 1])
end

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
% set(handles.maskframestep,'String',handles.data.maskfstep);
% set(handles.maskframestart,'String',handles.data.maskfstart);
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
% set(handles.currentjobname,'String',handles.data.batchname);
set(handles.outputdirectory,'String',handles.data.outdirec);
set(handles.outputpassbasename,'String',handles.data.outputpassbase);

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
set(handles.colorchannel_popupMenu,'Value',str2double(handles.data.channel));
set(handles.framestep,'String',handles.data.framestep);
set(handles.PIVerror,'String',handles.data.PIVerror);

if get(handles.passtype,'Value')>=5
    set(handles.velocityinterptype,'backgroundcolor',[1 1 1]);
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
        im1=double(imread(fullfile(handles.data.imdirec, [handles.data.imbase sprintf(['%0.' handles.data.imzeros 'i.' handles.data.imext],str2double(handles.data.imfstart))])));
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
                                handles.data.exp_diffractiondiameter='';
                                handles.data.exp_depthoffocus='';
                            else
                                handles.data.exp_diffractiondiameter=num2str(d);
                                handles.data.exp_depthoffocus=num2str(dz);
                            end
                        catch
                            handles.data.exp_diffractiondiameter='';
                            handles.data.exp_depthoffocus='';
                        end
                    elseif ~(isempty(handles.data.exp_NA) || isempty(handles.data.exp_n) || isempty(handles.data.exp_wavelength)) && str2double(handles.data.exp_micro)
                        try
                            fnum=0.5*sqrt((str2double(handles.data.exp_n)/str2double(handles.data.exp_NA))^2-1);
                            d=2.44*fnum*(M+1)*str2double(handles.data.exp_wavelength);
                            dz=2*fnum*d*(M+1)/M^2;

                            if isnan(d) || isnan(dz)
                                handles.data.exp_diffractiondiameter='';
                                handles.data.exp_depthoffocus='';
                            else
                                handles.data.exp_diffractiondiameter=num2str(d);
                                handles.data.exp_depthoffocus=num2str(dz);
                            end
                        catch
                            handles.data.exp_diffractiondiameter='';
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
        else
            handles.data.exp_M='';
            handles.data.exp_diffractiondiameter='';
            handles.data.exp_depthoffocus='';
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
    try
        mkdir(Data.outdirec)
        fid=fopen(fname,'w');
        if fid==-1
            error(['error writing experiment summary ',fname])
        end
    catch
        error(['error writing experiment summary ',fname])
    end
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

methods={'Multipass - DWO','Multigrid - DWO','Multigrid - Deform (DWO)','Multigrid - Ensemble (DWO)','Multigrid - Multiframe (DWO)'};
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
    corr={'SCC','RPC','SPC'};
    peak={'Three-Point Gaussian','Four-Point Gaussian','Gaussian Least Squares'};
    y_n={'No','Yes'};
    A=eval(['Data.PIV' num2str(i)]);
    fprintf(fid,['\n------------------------Pass ',num2str(i),' Setup-------------------------\n']);
    fprintf(fid,['Window Resolution (first image) (pix):       ',A.winres1,'\n']);
    fprintf(fid,['Window Resolution (second image) (pix):       ',A.winres2,'\n']);
    fprintf(fid,['Window Size (pix):             ',A.winsize,'\n']);
    fprintf(fid,['Grid Resolution (pix):         ',A.gridres,'\n']);
    fprintf(fid,['Window Overlap Percentage:     ',A.winoverlap,'\n']);
    fprintf(fid,['Grid Buffer (pix):             ',A.gridbuf,'\n']);
    fprintf(fid,['Bulk Window Offset (pix):      ',A.BWO,'\n']);
    fprintf(fid,['Correlation:                   ',corr{str2double(A.corr)},'\n']);
    fprintf(fid,['Zero-Mean Image Windows:       ',y_n{str2double(A.zeromean)+1},'\n']);
    fprintf(fid,['Subpixel Peak Location Method: ',peak{str2double(A.peaklocator)},'\n']);
    if str2double(Data.method)~=1 && str2double(Data.method)~=5
        fprintf(fid,['Smoothing:                     ',y_n{str2double(A.velsmooth)+1},'\n']);
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
    fprintf(fid,['Write Output:                  ',y_n{str2double(A.write)+1},'\n']);
    if str2double(A.write)
        fprintf(fid,['Output Basename:               ',A.outbase,'\n']);
        fprintf(fid,['Save Add. Peak Info:           ',y_n{str2double(A.savepeakinfo)+1},'\n']);
        if str2double(A.savepeakinfo)
            peaks={'1','1,2','1,2,3'};
            fprintf(fid,['Save Data for Peaks:           ',peaks{str2double(A.corrpeaknum)},'\n']);
            fprintf(fid,['Save Peak Magnitude:           ',y_n{str2double(A.savepeakmag)+1},'\n']);
            fprintf(fid,['Save Resulting Vel.:           ',y_n{str2double(A.savepeakvel)+1},'\n']);
        end
    end
end

fprintf(fid,'\n----------------------Images and Masking---------------------\n');
% fprintf(fid,['Image Directory:               ',Data.imdirec,'\n']);
fprintf(fid,['Image Basename:                ',Data.imbase,'\n']);
fprintf(fid,['Image Zeros:                   ',Data.imbase,'\n']);
fprintf(fid,['Image Extension:               ',Data.imbase,'\n']);
fprintf(fid,['Image Frame Start:             ',Data.imbase,'\n']);
fprintf(fid,['Image Frame Step:              ',Data.imbase,'\n']);
fprintf(fid,['Image Frame End:               ',Data.imbase,'\n']);
fprintf(fid,['Image Correlation Step:        ',Data.imbase,'\n']);
fprintf(fid, ['Color Channel:                   ', Data.channel,'\n']);
fprintf(fid,'Masking Type: ');
if strcmp(Data.masktype,'static')
    if ispc
        slshind=strfind(Data.staticmaskname,'\');
    else
        slshind=strfind(Data.staticmaskname,'/');
    end
    fprintf(fid,['Static\nMask File: ',Data.staticmaskname(slshind(end)+1:end),'\n']);
elseif strcmp(Data.masktype,'dynamic')
    fprintf(fid,'Dynamic\n');
%     fprintf(fid,['Mask Directory:                ',Data.maskdirec,'\n']);
    fprintf(fid,['Mask Basename:                 ',Data.maskbase,'\n']);
    fprintf(fid,['Mask Zeros:                    ',Data.maskzeros,'\n']);
    fprintf(fid,['Mask Extension:                ',Data.maskext,'\n']);
    fprintf(fid,['Mask Frame Start:              ',Data.maskfstart,'\n']);
    fprintf(fid,['Mask Frame Step:               ',Data.maskfstep,'\n']);
else
    fprintf(fid,'None\n');
end

fclose(fid);

function ButtonName=splashdlg(Question,Title,Btn1,Btn2,Btn3,Default)
%QUESTDLG Question dialog box.
%  ButtonName = QUESTDLG(Question) creates a modal dialog box that
%  automatically wraps the cell array or string (vector or matrix)
%  Question to fit an appropriately sized window.  The name of the
%  button that is pressed is returned in ButtonName.  The Title of
%  the figure may be specified by adding a second string argument:
%
%    ButtonName = questdlg(Question, Title)
%
%  Question will be interpreted as a normal string.
%
%  QUESTDLG uses UIWAIT to suspend execution until the user responds.
%
%  The default set of buttons names for QUESTDLG are 'Yes','No' and
%  'Cancel'.  The default answer for the above calling syntax is 'Yes'.
%  This can be changed by adding a third argument which specifies the
%  default Button:
%
%    ButtonName = questdlg(Question, Title, 'No')
%
%  Up to 3 custom button names may be specified by entering
%  the button string name(s) as additional arguments to the function
%  call.  If custom button names are entered, the default button
%  must be specified by adding an extra argument, DEFAULT, and
%  setting DEFAULT to the same string name as the button you want
%  to use as the default button:
%
%    ButtonName = questdlg(Question, Title, Btn1, Btn2, DEFAULT);
%
%  where DEFAULT is set to Btn1.  This makes Btn1 the default answer.
%  If the DEFAULT string does not match any of the button string names,
%  a warning message is displayed.
%
%  To use TeX interpretation for the Question string, a data
%  structure must be used for the last argument, i.e.
%
%    ButtonName = questdlg(Question, Title, Btn1, Btn2, OPTIONS);
%
%  The OPTIONS structure must include the fields Default and Interpreter.
%  Interpreter may be 'none' or 'tex' and Default is the default button
%  name to be used.
%
%  If the dialog is closed without a valid selection, the return value
%  is empty.
%
%  Example:
%
%  ButtonName = questdlg('What is your favorite color?', ...
%                        'Color Question', ...
%                        'Red', 'Green', 'Blue', 'Green');
%  switch ButtonName,
%    case 'Red',
%     disp('Your favorite color is Red');
%    case 'Blue',
%     disp('Your favorite color is Blue.')
%     case 'Green',
%      disp('Your favorite color is Green.');
%  end % switch
%
%  See also DIALOG, ERRORDLG, HELPDLG, INPUTDLG, LISTDLG,
%    MSGBOX, WARNDLG, FIGURE, TEXTWRAP, UIWAIT, UIRESUME.


%  Copyright 1984-2007 The MathWorks, Inc.
%  $Revision: 5.55.4.14 $


if nargin<1
  error('MATLAB:questdlg:TooFewArguments', 'Too few arguments for QUESTDLG');
end

Interpreter='none';
if ~iscell(Question),Question=cellstr(Question);end

%%%%%%%%%%%%%%%%%%%%%
%%% General Info. %%%
%%%%%%%%%%%%%%%%%%%%%
Black      =[0       0        0      ]/255;
% LightGray  =[192     192      192    ]/255;
% LightGray2 =[160     160      164    ]/255;
% MediumGray =[128     128      128    ]/255;
% White      =[255     255      255    ]/255;

%%%%%%%%%%%%%%%%%%%%
%%% Nargin Check %%%
%%%%%%%%%%%%%%%%%%%%
if nargout>1
  error('MATLAB:questdlg:WrongNumberOutputs', 'Wrong number of output arguments for QUESTDLG');
end
if nargin==1,Title=' ';end
if nargin<=2, Default='Yes';end
if nargin==3, Default=Btn1 ;end
if nargin<=3, Btn1='Yes'; Btn2='No'; Btn3='Cancel';NumButtons=3;end
if nargin==4, Default=Btn2;Btn2=[];Btn3=[];NumButtons=1;end
if nargin==5, Default=Btn3;Btn3=[];NumButtons=2;end
if nargin==6, NumButtons=3;end
if nargin>6
  error('MATLAB:questdlg:TooManyInputs', 'Too many input arguments');NumButtons=3; %#ok
end

if isstruct(Default),
  Interpreter=Default.Interpreter;
  Default=Default.Default;
end


%%%%%%%%%%%%%%%%%%%%%%%
%%% Create QuestFig %%%
%%%%%%%%%%%%%%%%%%%%%%%
FigPos    = get(0,'DefaultFigurePosition');
FigPos(3) = 300;
FigPos(4) =  90;
FigPos    = getnicedialoglocation(FigPos, get(0,'DefaultFigureUnits'));

QuestFig=dialog(                                    ...
  'Visible'         ,'off'                      , ...
  'Name'            ,Title                      , ...
  'Pointer'         ,'arrow'                    , ...
  'Position'        ,FigPos                     , ...
  'KeyPressFcn'     ,@doFigureKeyPress          , ...
  'IntegerHandle'   ,'off'                      , ...
  'WindowStyle'     ,'normal'                   , ...
  'HandleVisibility','callback'                 , ...
  'CloseRequestFcn' ,@delete                  , ...
  'Tag'             ,Title                      , ...
  'Color'           ,[1 1 1]                      ...
  );

%%%%%%%%%%%%%%%%%%%%%
%%% Set Positions %%%
%%%%%%%%%%%%%%%%%%%%%
DefOffset  =20;

IconWidth  =300;
IconHeight =150;
IconXOffset=DefOffset;
IconYOffset=FigPos(4)-DefOffset-IconHeight;  %#ok
IconCMap=[Black;get(QuestFig,'Color')];  %#ok

DefBtnWidth =56;
BtnHeight   =22;

BtnYOffset=DefOffset;

BtnWidth=DefBtnWidth;

ExtControl=uicontrol(QuestFig   , ...
  'Style'    ,'pushbutton', ...
  'String'   ,' '          ...
  );

btnMargin=1.4;
set(ExtControl,'String',Btn1);
BtnExtent=get(ExtControl,'Extent');
BtnWidth=max(BtnWidth,BtnExtent(3)+20);
if NumButtons > 1
  set(ExtControl,'String',Btn2);
  BtnExtent=get(ExtControl,'Extent');
  BtnWidth=max(BtnWidth,BtnExtent(3)+20);
  if NumButtons > 2
    set(ExtControl,'String',Btn3);
    BtnExtent=get(ExtControl,'Extent');
    BtnWidth=max(BtnWidth,BtnExtent(3)*btnMargin);
  end
end
BtnHeight = max(BtnHeight,BtnExtent(4)*btnMargin);

delete(ExtControl);

MsgTxtXOffset=IconXOffset+IconWidth;

FigPos(3)=max(FigPos(3),MsgTxtXOffset+NumButtons*(BtnWidth+2*DefOffset));
set(QuestFig,'Position',FigPos);

BtnXOffset=zeros(NumButtons,1);

if NumButtons==1,
  BtnXOffset=(FigPos(3)-BtnWidth)/2;
elseif NumButtons==2,
  BtnXOffset=[MsgTxtXOffset
    FigPos(3)-DefOffset-BtnWidth];
elseif NumButtons==3,
  BtnXOffset=[MsgTxtXOffset
    0
    FigPos(3)-DefOffset-BtnWidth];
  BtnXOffset(2)=(BtnXOffset(1)+BtnXOffset(3))/2;
end

MsgTxtYOffset=DefOffset+BtnYOffset+BtnHeight;
MsgTxtWidth=FigPos(3)-DefOffset-MsgTxtXOffset-IconWidth;
MsgTxtHeight=FigPos(4)-DefOffset-MsgTxtYOffset;
MsgTxtForeClr=Black;
MsgTxtBackClr=[1 1 1];%get(QuestFig,'Color');

CBString='uiresume(gcbf)';
DefaultValid = false;
DefaultWasPressed = false;
BtnHandle = [];
DefaultButton = 0;

for i = 1:NumButtons
  switch i
    case 1
      ButtonString=Btn1;
      ButtonTag='Btn1';
      ButtonType='pushbutton';
      if strcmp(ButtonString, Default)
        DefaultValid = true;
        DefaultButton = 1;
      end

    case 2
      ButtonString=Btn2;
      ButtonTag='Btn2';
      ButtonType='checkbox';
      if strcmp(ButtonString, Default)
        DefaultValid = true;
        DefaultButton = 2;
      end
  end

  BtnHandle(end+1)=uicontrol(QuestFig            , ...
    'Style'              ,ButtonType, ...
    'Position'           ,[ BtnXOffset(1) BtnYOffset BtnWidth BtnHeight ]           , ...
    'KeyPressFcn'        ,@doControlKeyPress , ...
    'BackgroundColor'    ,[1 1 1]     , ...
    'CallBack'           ,CBString    , ...
    'String'             ,ButtonString, ...
    'HorizontalAlignment','center'    , ...
     'Tag'                ,ButtonTag     ...
    );
end

if ~DefaultValid
  warnstate = warning('backtrace','off');
  warning('MATLAB:QUESTDLG:stringMismatch','Default string does not match any button string name.');
  warning(warnstate);
end

MsgHandle=uicontrol(QuestFig            , ...
  'Style'              ,'text'         , ...
  'Position'           ,[MsgTxtXOffset MsgTxtYOffset .95*MsgTxtWidth MsgTxtHeight ]              , ...
  'String'             ,{' '}          , ...
  'Tag'                ,'Question'     , ...
  'HorizontalAlignment','left'         , ...
  'FontWeight'         ,'bold'         , ...
  'BackgroundColor'    ,MsgTxtBackClr  , ...
  'ForegroundColor'    ,MsgTxtForeClr    ...
  );

[WrapString,NewMsgTxtPos]=textwrap(MsgHandle,Question,75);

% NumLines=size(WrapString,1);

AxesHandle=axes('Parent',QuestFig,'Position',[0 0 1 1],'Visible','off');

texthandle=text( ...  
    'Parent'              ,AxesHandle                      , ...
    'Units'               ,'pixels'                        , ...
    'Color'               ,get(BtnHandle(1),'ForegroundColor')   , ...
    'HorizontalAlignment' ,'left'                          , ...
    'FontName'            ,get(BtnHandle(1),'FontName')    , ...
    'FontSize'            ,14.0                            , ...
    'FontWeight'          ,'normal'                          , ...
    'VerticalAlignment'   ,'bottom'                        , ...
    'String'              ,WrapString                      , ...
    'Interpreter'         ,Interpreter                     , ...
    'Tag'                 ,'Question'                        ...
    );  %#ok

textExtent = get(texthandle, 'extent');

% (g357851)textExtent and extent from uicontrol are not the same. For window, extent from uicontrol is larger
%than textExtent. But on Mac, it is reverse. Pick the max value.
MsgTxtWidth=max([MsgTxtWidth NewMsgTxtPos(3)+2 textExtent(3)]);
MsgTxtHeight=max([MsgTxtHeight NewMsgTxtPos(4)+2 textExtent(4)]);

MsgTxtXOffset=IconXOffset+IconWidth+DefOffset;
FigPos(3)=max(NumButtons*(BtnWidth+DefOffset)+DefOffset, ...
  MsgTxtXOffset+MsgTxtWidth+DefOffset);


% Center Vertically around icon
if IconHeight>MsgTxtHeight,
  IconYOffset=BtnYOffset+BtnHeight+DefOffset;
  MsgTxtYOffset=IconYOffset+(IconHeight-MsgTxtHeight)/2;
  FigPos(4)=IconYOffset+IconHeight+DefOffset;
  % center around text
else
  MsgTxtYOffset=BtnYOffset+BtnHeight+DefOffset;
  IconYOffset=MsgTxtYOffset+(MsgTxtHeight-IconHeight)/2;
  FigPos(4)=MsgTxtYOffset+MsgTxtHeight+DefOffset;
end

if NumButtons==1,
  BtnXOffset=(FigPos(3)-BtnWidth)/2;
elseif NumButtons==2,
  BtnXOffset=[(FigPos(3)-DefOffset)/2-BtnWidth
    (FigPos(3)+DefOffset)/2
    ];

elseif NumButtons==3,
  BtnXOffset(2)=(FigPos(3)-BtnWidth)/2;
  BtnXOffset=[BtnXOffset(2)-DefOffset-BtnWidth
    BtnXOffset(2)
    BtnXOffset(2)+BtnWidth+DefOffset
    ];
end

set(QuestFig ,'Position',getnicedialoglocation(FigPos, get(QuestFig,'Units')));

BtnPos=get(BtnHandle,{'Position'});
BtnPos=cat(1,BtnPos{:});
BtnPos(:,1)=BtnXOffset;
BtnPos=num2cell(BtnPos,2);
set(BtnHandle,{'Position'},BtnPos);

if DefaultValid
  setdefaultbutton(QuestFig, BtnHandle(DefaultButton));
end

delete(MsgHandle);


set(texthandle, 'Position',[MsgTxtXOffset MsgTxtYOffset 0]);


IconAxes=axes(                                      ...
  'Parent'      ,QuestFig              , ...
  'Units'       ,'Pixels'              , ...
  'Position'    ,[IconXOffset IconYOffset IconWidth IconHeight], ...
  'NextPlot'    ,'replace'             , ...
  'Tag'         ,'IconAxes'              ...
  );

set(QuestFig ,'NextPlot','add');

pranadir=which('prana');
try
    logo=imread(fullfile(pranadir(1:end-8),'documentation','logo.tif'),'tif');
catch
    logo=zeros(500,1000);
end
Img=image('Cdata',logo,'Parent',IconAxes);

set(IconAxes, ...
  'Visible','off'           , ...
  'YDir'   ,'reverse'       , ...
  'XLim'   ,get(Img,'XData'), ...
  'YLim'   ,get(Img,'YData')  ...
  );

% make sure we are on screen
movegui(QuestFig)


set(QuestFig ,'WindowStyle','modal','Visible','on');
drawnow;

if DefaultButton ~= 0
  uicontrol(BtnHandle(DefaultButton));
end

if ishghandle(QuestFig)
  % Go into uiwait if the figure handle is still valid.
  % This is mostly the case during regular use.
  uiwait(QuestFig);
end

% Check handle validity again since we may be out of uiwait because the
% figure was deleted.
if ishghandle(QuestFig)
  if DefaultWasPressed
    ButtonName=Default;
  else
    ButtonName=get(get(QuestFig,'CurrentObject'),'String');
  end
    delete(QuestFig);
else
  ButtonName='';
end

function doFigureKeyPress(obj, evd)  %#ok
switch(evd.Key)
case {'return','space'}
% if DefaultValid
  DefaultWasPressed = true;
  uiresume(gcbf);
% end
case 'escape'
    delete(QuestFig);
end

function doControlKeyPress(obj, evd)  %#ok
switch(evd.Key)
case {'return'}
% if DefaultValid
  DefaultWasPressed = true;
  uiresume(gcbf);
% end
case 'escape'
    delete(QuestFig);
end

function doDelete(varargin)  %#ok
delete(QuestFig);

function figure_size = getnicedialoglocation(figure_size, figure_units)
% adjust the specified figure position to fig nicely over GCBF
% or into the upper 3rd of the screen

%  Copyright 1999-2006 The MathWorks, Inc.
%  $Revision: 1.1.6.3 $

parentHandle = gcbf;
propName = 'Position';
if isempty(parentHandle)
    parentHandle = 0;
    propName = 'ScreenSize';
end

old_u = get(parentHandle,'Units');
set(parentHandle,'Units',figure_units);
container_size=get(parentHandle,propName);
set(parentHandle,'Units',old_u);

figure_size(1) = container_size(1)  + 1/2*(container_size(3) - figure_size(3));
figure_size(2) = container_size(2)  + 2/3*(container_size(4) - figure_size(4));

function setdefaultbutton(figHandle, btnHandle)
% WARNING: This feature is not supported in MATLAB and the API and
% functionality may change in a future release.

%SETDEFAULTBUTTON Set default button for a figure.
%  SETDEFAULTBUTTON(BTNHANDLE) sets the button passed in to be the default button
%  (the button and callback used when the user hits "enter" or "return"
%  when in a dialog box.
%
%  This function is used by inputdlg.m, msgbox.m, questdlg.m and
%  uigetpref.m.
%
%  Example:
%
%  f = figure;
%  b1 = uicontrol('style', 'pushbutton', 'string', 'first', ...
%       'position', [100 100 50 20]);
%  b2 = uicontrol('style', 'pushbutton', 'string', 'second', ...
%       'position', [200 100 50 20]);
%  b3 = uicontrol('style', 'pushbutton', 'string', 'third', ...
%       'position', [300 100 50 20]);
%  setdefaultbutton(b2);
%

%  Copyright 2005-2007 The MathWorks, Inc.

% Nargin Check
if nargin<1, error('MATLAB:setdefaultbutton:InvalidNumberOfArguments','Too few arguments for setdefaultbutton'); end
if nargin>2, error('MATLAB:setdefaultbutton:InvalidNumberOfArguments','Too many arguments for setdefaultbutton'); end

if (usejava('awt') == 1)
    % We are running with Java Figures
    useJavaDefaultButton(figHandle, btnHandle)
else
    % We are running with Native Figures
    useHGDefaultButton(figHandle, btnHandle);
end

function useJavaDefaultButton(figH, btnH)
% Get a UDD handle for the figure.
fh = handle(figH);
% Call the setDefaultButton method on the figure handle
fh.setDefaultButton(btnH);

function useHGDefaultButton(figHandle, btnHandle)
% First get the position of the button.
btnPos = getpixelposition(btnHandle);

% Next calculate offsets.
leftOffset   = btnPos(1) - 1;
bottomOffset = btnPos(2) - 2;
widthOffset  = btnPos(3) + 3;
heightOffset = btnPos(4) + 3;

% Create the default button look with a uipanel.
% Use black border color even on Mac or Windows-XP (XP scheme) since
% this is in natve figures which uses the Win2K style buttons on Windows
% and Motif buttons on the Mac.
h1 = uipanel(get(btnHandle, 'Parent'), 'HighlightColor', 'black', ...
    'BorderType', 'etchedout', 'units', 'pixels', ...
    'Position', [leftOffset bottomOffset widthOffset heightOffset]);

% Make sure it is stacked on the bottom.
uistack(h1, 'bottom');


% --------------------------------------------------------------------





% --- Executes on button press in runidcheckbox.
function runidcheckbox_Callback(hObject, eventdata, handles)


% --- Executes on selection change in idmethod.
function idmethod_Callback(hObject, eventdata, handles)
function idmethod_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function idimthresh_Callback(hObject, eventdata, handles)
function idimthresh_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function idsaveloc_Callback(hObject, eventdata, handles)
function idsaveloc_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in loadidsaveloc.
function loadidsaveloc_Callback(hObject, eventdata, handles)



function idsavebase_Callback(hObject, eventdata, handles)
function idsavebase_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in runsizingcheckbox.
function runsizingcheckbox_Callback(hObject, eventdata, handles)

% --- Executes on selection change in sizingmethod.
function sizingmethod_Callback(hObject, eventdata, handles)
function sizingmethod_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function sizingstd_Callback(hObject, eventdata, handles)
function sizingstd_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function sizingsaveloc_Callback(hObject, eventdata, handles)
function sizingsaveloc_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in loadsizingsaveloc.
function loadsizingsaveloc_Callback(hObject, eventdata, handles)



function sizingsavebase_Callback(hObject, eventdata, handles)
function sizingsavebase_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in runtrackingcheckbox.
function runtrackingcheckbox_Callback(hObject, eventdata, handles)


% --- Executes on selection change in trackingmethod.
function trackingmethod_Callback(hObject, eventdata, handles)
function trackingmethod_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function trackingPIVweight_Callback(hObject, eventdata, handles)
function trackingPIVweight_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function trackingsaveloc_Callback(hObject, eventdata, handles)
function trackingsaveloc_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function [num1 num2 num3 num4] = parseNum(str)
% parseNUM parses a numerical string containing delimiters (, ; : and blank space).
% Note that the allowed delimiters are specified here as comma, semicolon,
% colon, and blank space. Also note that each delimiter here is functionally equivalent;
% i.e., inputs of str = '32,32:16,16' and str = '32:32,16;16' will yield identical results.
% 
% INPUTS
%   str = String containing numbers and delimiters (string)
%
% OUTPUTS
%   [num1 num2 num3 num4] = Numbers parsed from input string (double precision or integer)
% 
% EXAMPLE
%   str = '32, 32 ; 64, 64';
%   [num1 num2 num3 num4] = parseNum(str);
% 
% ans =
% 
%     32    32    64    64
% 
% SEE ALSO
%   strread
% 

% If an empty string is input...
if strcmp(str, '')
%  Prompt user for input
    fprintf(1, 'Enter dimensions...\n\n');
%  Set numbers to default values
    num1 = 32;
    num2 = 32;
    num3 = 32;
    num4 = 32;
else

    
 %  Read numerical values from string and ignore delimiters (, ; :)
    nums = strread(str, '%d', -1, 'delimiter', ',;:');
  
% Replicate array 4 times (max number of inputs in any texbox in this GUI)
% to allow simplified input for homogeneous window sizes.
nums = repmat(nums, 4, 1);

% Specify output numbers
num1 = nums(1);
num2 = nums(2);
num3 = nums(3);
num4 = nums(4);

end

function [xROI yROI] = roiSize(winRes)
% roiSize determines the dimensions (in pixels) of a region of interest via
% the largest window dimension in each direction. winRes is a
% numerical string containing delimiters (, ; : and blank space).
% Note that the allowed delimiters are specified here as comma, semicolon,
% colon, and blank space. Also note that each delimiter here is functionally equivalent;
% i.e., inputs of "winRes = 32,32:16,16" and "winRes = 32:32,16;16"
% will yield identical results.
% 
% INPUTS
%   winRes =  String containing window resolutions and delimiters (string)
% 
% OUTPUTS
%   xROI = x-dimension of region of interest (pixels)
%   yROI = y-dimension of region of interest (pixels)
% 
% EXAMPLE
%   winRes = '32, 64 ; 64, 32';
%   [xROI yROI] = roiSize(winRes);
% 
% ans =
% 
%     64    64
% 
% SEE ALSO
%   parseNums, strread
% 

% Parse string containing window size information
[wx1 wy1 wx2 wy2] = parseNum(winRes);

% Calculate max dimensions of effective window resolutions
wxMax = max(wx1, wx2);
wyMax = max(wy1, wy2);

% Calculate x- and y- sizes of the interrogation region (pixels)
xROI = 2^(ceil(log(2*wxMax)/log(2)));
yROI = 2^(ceil(log(2*wyMax)/log(2)));


% --------------------------------------------------------------------
function tools_menuItem_Callback(hObject, eventdata, handles)
% hObject    handle to tools_menuItem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --------------------------------------------------------------------
function saveworkspace_menuItem_Callback(hObject, eventdata, handles)
% hObject    handle to saveworkspace_menuItem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Use image directory as default directory for UIGET function
imageDir = get(handles.imagedirectory, 'String');

% Specify directory in which to save workspace
workspaceDir = uigetdir(imageDir);

% Exit save if the user selects 'cancel' in UIGET
if workspaceDir ~= 0;

%     Read the names of jobs in the present workspace
    jobs = get(handles.joblist, 'String');

%     Save each job to the directory specified by workspaceDir
    for n = 1:length(jobs);
                eval( ['Data = handles.' char(jobs(n)) ';'] );
                save([workspaceDir '/' char(jobs(n)) '.mat'], 'Data');
    end

else
    
    return
    
end


% --------------------------------------------------------------------

function trackingestweight_Callback(hObject, eventdata, handles)
function trackingestweight_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function trackingiterations_Callback(hObject, eventdata, handles)

function trackingiterations_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Tracking Validation Checkbox ---
function trackingvalcheckbox_Callback(hObject, eventdata, handles)
if str2double(handles.Njob)>0
    handles.data.Track.valprops.run=num2str(get(hObject,'Value'));
    update_PTV(handles)
    guidata(hObject,handles)
end

% --- Tracking Validation Coefficient ---
function trackingvalcoefficient_Callback(hObject, eventdata, handles)
if str2double(handles.Njob)>0
    handles.data.Track.valprops.C_cutoff=get(hObject,'String');
end
function trackingvalcoefficient_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Tracking Validation Radius ---
function trackingvalradius_Callback(hObject, eventdata, handles)
if str2double(handles.Njob)>0
    handles.data.Track.valprops.s_radius=get(hObject,'String');
end
function trackingvalradius_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Tracking Validation U Threshold ---
function trackingvalUthresh_Callback(hObject, eventdata, handles)
if str2double(handles.Njob)>0
    handles.data.Track.valprops.MAD_U=get(hObject,'String');
end
function trackingvalUthresh_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Tracking Validation V Threshold ---
function trackingvalVthresh_Callback(hObject, eventdata, handles)
if str2double(handles.Njob)>0
    handles.data.Track.valprops.MAD_V=get(hObject,'String');
end
function trackingvalVthresh_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Run PIV Checkbox ---
function runPIVcheckbox_Callback(hObject, eventdata, handles)
if str2double(handles.Njob)>0
    handles.data.runPIV=num2str(get(hObject,'Value'));
end


% --- Executes on button press in loadtrackingsaveloc.
function loadtrackingsaveloc_Callback(hObject, eventdata, handles)



function trackingsavebase_Callback(hObject, eventdata, handles)
function trackingsavebase_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function trackingradius_Callback(hObject, eventdata, handles)
function trackingradius_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function trackingdistweight_Callback(hObject, eventdata, handles)
function trackingdistweight_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function trackingsizeweight_Callback(hObject, eventdata, handles)
function trackingsizeweight_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function trackingintensityweight_Callback(hObject, eventdata, handles)
function trackingintensityweight_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function trackingestradius_Callback(hObject, eventdata, handles)
function trackingestradius_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function trackingvectors_Callback(hObject, eventdata, handles)
function trackingvectors_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in colorchannel_popupMenu.
function colorchannel_popupMenu_Callback(hObject, eventdata, handles)
% hObject    handle to colorchannel_popupMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if str2double(handles.Njob)>0
    handles.data.channel = num2str(get(hObject,'Value'));
    guidata(hObject,handles)
end


% --- Executes during object creation, after setting all properties.
function colorchannel_popupMenu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to colorchannel_popupMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called



% --- Executes on selection change in trackingprediction.
function trackingprediction_Callback(hObject, eventdata, handles)
function trackingprediction_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function verion_box_Callback(hObject, eventdata, handles)
% hObject    handle to verion_box (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of verion_box as text
%        str2double(get(hObject,'String')) returns contents of verion_box as a double


% --- Executes during object creation, after setting all properties.
function verion_box_CreateFcn(hObject, eventdata, handles)
% hObject    handle to verion_box (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function version_box_Callback(hObject, eventdata, handles)
% hObject    handle to version_box (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of version_box as text
%        str2double(get(hObject,'String')) returns contents of version_box as a double


% --- Executes during object creation, after setting all properties.
function version_box_CreateFcn(hObject, eventdata, handles)
% hObject    handle to version_box (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function outputpassbasename_Callback(hObject, eventdata, handles)
% hObject    handle to outputpassbasename (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of outputpassbasename as text
%        str2double(get(hObject,'String')) returns contents of outputpassbasename as a double
if str2double(handles.Njob)>0
    N=str2double(handles.data.passes); % Number of passes
    passbase = get(handles.outputpassbasename,'String'); % Read text in outbasename textbox
    for i = 1:N
        eval(['handles.data.PIV' num2str(i) '.outbase=[''' passbase 'pass' num2str(i) '_''];']);
    end
    cpass = get(handles.passlist,'Value'); % This grabs the currently selected pass number
    set(handles.outputbasename,'String',eval(['handles.data.PIV' num2str(cpass) '.outbase']));
    handles.data.outputpassbase = passbase;
    guidata(hObject,handles)
end


% --- Executes during object creation, after setting all properties.
function outputpassbasename_CreateFcn(hObject, eventdata, handles)
% hObject    handle to outputpassbasename (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
