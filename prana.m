function varargout = prana(varargin)
% HELP TEXT for running prana from the command line goes here.
%
%     This file is part of prana, an open-source GUI-driven program for
%     calculating velocity fields using PIV or PTV.
%
%     Copyright (C) 2010-2014  Virginia Polytechnic Institute and State
%     University
% 
%     Copyright 2014-2015.  Los Alamos National Security, LLC. This material was
%     produced under U.S. Government contract DE-AC52-06NA25396 for Los 
%     Alamos National Laboratory (LANL), which is operated by Los Alamos 
%     National Security, LLC for the U.S. Department of Energy. The U.S. 
%     Government has rights to use, reproduce, and distribute this software.
%     NEITHER THE GOVERNMENT NOR LOS ALAMOS NATIONAL SECURITY, LLC MAKES ANY
%     WARRANTY, EXPRESS OR IMPLIED, OR ASSUMES ANY LIABILITY FOR THE USE OF
%     THIS SOFTWARE.  If software is modified to produce derivative works,
%     such modified software should be clearly marked, so as not to confuse
%     it with the version available from LANL.
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
                   'gui_LayoutFcn',  @prana_LayoutFcn, ...
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
    %JJC 2011-12-22:
    % this can cause problems if defaultsettings.mat exists, but is from a previous version
    % we should implement a version check before trying to proceed and use the values stored in it
    %load('defaultsettings.mat')
    
    %for now, just always create the default from scratch
    error('prana:pranaGUI:defaultsettingsVersionError','defaultsettings.mat doesn''t match the current version')
    
    %Are we still saving defaultsettings.mat anywhere?  I couldn't find it.
catch ERR  %#ok<NASGU>
    
    % Now we call a function that will build a default prana Job.  This is
    % done in part to allow other (exterieur) codes to build default jobs.
    [defaultdata] = buildDefaultPranaJob();
end



if str2double(defaultdata.splash)==0 || str2double(defaultdata.clientversion)<2.0
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
    
    %JJC: Why are these here?  We define them above - redefining them here overwrites them,
    % but ONLY if the splash screen is displayed.  Also, the version numbers are wrong,
    % so the splash screen will always show.  
    % Finally, The save should be outside this splash screen check, or else the default will
    % be overwritten each time it runs if the user doesn't turn off the splash screen, but 
    % never saved if it's already disabled.  
    % It needs to be AFTER the pranaPIVcode('version') check, at any rate.
    
    %defaultdata.clientversion='0.99';
    %defaultdata.version='2.0';
    
    %save('defaultsettings.mat','defaultdata')
end

%JJC: for now, disable this save - running prana drops defaultsettings.mat files in whatever your current working directory is - ANNOYING!
%save('defaultsettings.mat','defaultdata')

handles.data=defaultdata;
pranadir=which('prana');
addpath([pranadir(1:end-7),'documentation']);

if ispc
    handles.loaddirec=[pwd '\'];
else
    handles.loaddirec=[pwd '/'];
end

handles.data.par='0';
try
    %findResource was removed in recent version (>2014b ?)
    if str2double(verstr(1:4))<2014
        compinfo=findResource('scheduler','configuration','local');
        handles.data.parprocessors=num2str(compinfo.clustersize);
    else %use parcluster instead
        %get information about the parallel options and settings
        compinfo = parcluster('local');
        %use them to set the default clustersize for parallel computing
        handles.data.parprocessors=num2str(compinfo.NumWorkers);
    end
catch
    handles.data.parprocessors='1';
end

handles.data.imdirec=pwd;
handles.data.imdirec2=pwd;
handles.data.maskdirec=pwd;
handles.data.outdirec=pwd;

try
    windowdiagram=imread(fullfile(pranadir(1:end-8),'documentation','windowdiagram.tif'),'tif');
catch
    windowdiagram=zeros(564,531);
end
set(gca,'children',imshow(windowdiagram));
axis off;

% the following commented lines were the original commands used by the
% GUIDE version of the GUI, and needed to be set to a default value due to
% a change in execution order when the GUI changed to text-only.
% handles.data.cpass=num2str(get(handles.passlist,'Value'));
handles.data.cpass='0';
handles.data0=handles.data;
% handles.Njob=num2str(size(get(handles.joblist,'String'),1));
% handles.Cjob=num2str(get(handles.joblist,'String'));
handles.Njob='0';
handles.Cjob='';

handles=rmfield(handles,'data');
% handles=update_data(handles);  %same problem, can't update yet because items don't exist yet
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

% --- Execute Menu -> Compile codes ---
function jobmenu_compile_codes_Callback(hObject, eventdata, handles)
compile_codes_Callback(hObject, eventdata, handles)

% --- Help Menu ---
function helpmenu_Callback(hObject, eventdata, handles)

% --- Help Menu -> About ---
function helpmenu_about_Callback(hObject, eventdata, handles)
msgbox( ...
    {...
    'prana','',['Client Version ',pranaPIVcode('version')],...
    ['Data Version ',num2str(handles.data0.clientversion)],'',...
    ['Authors: Sayantan Bhattacharya, Matt Giarra, Nick Cardwell, Brady Drew, ',...
    'Adric Eckstein, John Charonko, Sam Raben, and the rest of the AEThER Lab'],'',...
    'Copyright 2010-2015 - Virginia Polytechnic Institute and State University', ...
    'Copyright 2014-2015.  Los Alamos National Security, LLC.','',...
    ['This program is free software: you can redistribute it and/or modify ',...
    'it under the terms of the GNU General Public License as published by ',...
    'the Free Software Foundation, either version 3 of the License, or (at ',...
    'your option) any later version. This program is distributed in the hope ',...
    'that it will be useful, but WITHOUT ANY WARRANTY; without even the ',...
    'implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. ',...
    'See the GNU General Public License for more details. You should have ',...
    'received a copy of the GNU General Public License along with this program. ',...
    'If not, see www.gnu.org/licenses'],'',...
    ['Portions of this program are copyright John D''Errico and Jonas Reber ',...
    'under the following license:'],'',...
    ['Redistribution and use in source and binary forms, with or without ',...
     'modification, are permitted provided that the following conditions are ',...
     'met:'],...
     '    * Redistributions of source code must retain the above copyright',...
     '      notice, this list of conditions and the following disclaimer.',...
     '    * Redistributions in binary form must reproduce the above copyright',...
     '      notice, this list of conditions and the following disclaimer in',...
     '      the documentation and/or other materials provided with the',...
     '      distribution',...
     ['THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" ',...
     'AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE ',...
     'IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ',...
     'ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE ',...
     'LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR ',...
     'CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF ',...
     'SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS ',...
     'INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN ',...
     'CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ',...
     'ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE ',...
     'POSSIBILITY OF SUCH DAMAGE.']
    },...
    'About')

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
	
	% Check for compiled codes, and compile them
	% if they don't already exist.
	check_compiled_codes(hObject, eventdata, handles);
	
    Data=handles.data;
    if str2double(Data.runPIV)
        pranaPIVcode(Data);
    end
    if str2double(Data.ID.runid) || str2double(Data.Size.runsize) || str2double(Data.Track.runtrack)
        pranaPTVcode(Data)
    end
end

% --- Run All Jobs Button ---
function runall_Callback(hObject, eventdata, handles)
if str2double(handles.Njob)>0
	
	% Check for compiled codes, and compile them
	% if they don't already exist.
	check_compiled_codes(hObject, eventdata, handles);
    
    Jlist=char(get(handles.joblist,'String'));
    eval(['handles.' Jlist(str2double(handles.Cjob),:) '=handles.data;']);
    
    for e=1:size(Jlist,1)
        Data=eval(['handles.' Jlist(e,:)]);
        if get(handles.runPIVcheckbox,'value')
            pranaPIVcode(Data);
        end
        if get(handles.runidcheckbox,'value') || get(handles.runsizingcheckbox,'value') || get(handles.runtrackingcheckbox,'value')
            pranaPTVcode(Data)
        end
    end
end

function check_compiled_codes(hObject, eventdata, handles)
	% File name for compiled codes
	% mexext is a built-in matlab command.
	compiled_code_file_name = ['whittaker_blackman.' mexext];
	
	% Compile the compileable codes if their compiled versions
	% don't already exist.
	if ~exist(compiled_code_file_name, 'file')
		compile_codes_Callback(hObject, eventdata, handles)
	end

% --- Run Current Job Button ---
function compile_codes_Callback(hObject, eventdata, handles)

% Inform the user.
fprintf('Compiling whittaker_blackman.c to mex file....\n');

% Compile the codes using the mac command
% Case for linux machines
try
	if isunix && ~ismac
	    mex -O CFLAGS="\$CFLAGS -std=c99" whittaker_blackman.c
	else
		% This command should work with both mac and windows.
		mex whittaker_blackman.c;
	end
	
	% Inform the user
	fprintf(['Compiled whittaker_blackman.c to whittaker_blackman.' mexext '\n']);

catch
	fprintf('Error compiling codes.\n');
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
            load(fullfile(handles.loaddirec,char(f(pp))));
            if exist('Data','var')~=0
                vn=0;
                while vn==0

                    [Data] = jobfile_validator(Data);

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
    uisave('Data',fullfile(handles.data.outdirec,[handles.data.batchname '.mat']));
end

% --- Copy Job Button ---
function copyjobbutton_Callback(hObject, eventdata, handles)
if str2double(handles.Njob)>0
    Jlist = char(get(handles.joblist,'String'));
    eval(['handles.' Jlist(str2double(handles.Cjob),:) '=handles.data;']);
end
Data = handles.data;
vn = 0;
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
    C = handles.data.imcstep;
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
        handles.data.imext=imname((i(end)+1):end);
        % This finds everything that is not a number
        j = regexp(imname(1:i(end)-1),'\D');
        % Now looking at the first thing that is not a number excluding the
        % extension to determine the number of zeros.
        zeros=i(end)-1-j(end);
        handles.data.imbase=imname(1:(i(end)-1-zeros));
        handles.data.imzeros=num2str(zeros);
        fstart=str2double(imname((i(end)-zeros):(i(end)-1)));
        handles.data.imfstart=num2str(fstart);
        dirinfo = dir([handles.data.imdirec handles.data.imbase '*.' handles.data.imext]);
        handles.data.imfend=num2str(str2double(dirinfo(end).name(i(end)-zeros:i(end)-1))-str2double(C));
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
    %load_imlist(handles);
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
%     load_imlist(handles);
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
    if any(get(handles.passtype,'Value')==[4 5]) && get(hObject,'Value')==1
        errordlg('Dynamic Masking is not compatible with the Ensemble correlation.','Warning')
    end
    handles=update_data(handles);
    guidata(hObject,handles)
end

% --- Static Mask File Text Box ---
function staticmaskfile_Callback(hObject, eventdata, handles)
if str2double(handles.Njob)>0 && get(handles.staticmaskbutton,'Value')==1
    handles.data.staticmaskname = get(hObject,'String');
    if ~isempty(dir(handles.data.staticmaskname))
        set(handles.staticmaskfile,'Backgroundcolor',[1 1 1])
    end
    load_imlist(handles);
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
        createStaticMask(uint8(im), handles.data.staticmaskname);
         
%          mask=ones(size(im,1),size(im,2));
%         roiwindow = CROIEditor(im./max(im(:)));
%         while isempty(roiwindow.labels)
%         addlistener(roiwindow,'MaskDefined',@your_roi_defined_callback);
%         drawnow
%         end
%         
%         mask(roiwindow.labels==0) = 0;
%         
%         
%         
%         
%         
%         imwrite(mask,handles.data.staticmaskname,'tif')

        set(handles.staticmaskfile,'String',handles.data.staticmaskname);
        handles.data.masktype='static';
        set(handles.staticmaskbutton,'Value',1)
        set(handles.dynamicmaskbutton,'Value',0)
        set(handles.nomaskbutton,'Value',0)
        Jlist=char(get(handles.joblist,'String'));
        eval(['handles.' Jlist(str2double(handles.Cjob),:) '=handles.data;']);
        handles=update_data(handles);
        guidata(hObject,handles)

%         close('Analyzer - ROI Editor')
        
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
        disp(ME.message);
        msgbox('Image Frame Not Found');
        e=-1;
    end

end

% function your_roi_defined_callback(h,e)
% [mask, labels, n] = roiwindow.getROIData;
% delete(roiwindow);

% --- Preview Image + Mask Button ---
function impreview_Callback(hObject, eventdata, handles)
if str2double(handles.Njob)>0
    
    % Read color channel
    channel = str2double(handles.data.channel);
    
    try
        fileloc = fullfile(handles.data.imdirec, [handles.data.imbase, sprintf(['%0.' handles.data.imzeros 'i.' handles.data.imext],str2double(handles.data.imfstart))]);
        % Get information about the file with the most important being the
        % MaxSampleValue which will be need to normalizing the image.
        imageInfo = imfinfo(fileloc);
        im = double(imread(fileloc));
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
    im1 = im(end:-1:1,:,:) ./ max(im(:));
%     im1 = im(end:-1:1,:,:)./max(imageInfo.MaxSampleValue);

        try
            if strcmp(handles.data.masktype,'static')
                mask = double(imread(handles.data.staticmaskname));
                mask = flipud(mask);
            elseif strcmp(handles.data.masktype,'dynamic')
                mask = double(imread(fullfile(handles.data.maskdirec, [handles.data.maskbase, sprintf(['%0.' handles.data.maskzeros 'i.' handles.data.maskext],str2double(handles.data.imfstart))])));
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
    catch ER
        disp(ER.message);
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
        set(gcf,'name',sprintf('Image + Mask, %0.0f total vectors locations shown for pass %s',length(X),handles.data.cpass))
    end
end

% --- Output Directory Text Box ---
function outputdirectory_Callback(hObject, eventdata, handles)
if str2double(handles.Njob)>0
    handles.data.outdirec = get(hObject,'String');
    if ~isempty(dir(handles.data.outdirec))
        set(handles.outputdirectory,'Backgroundcolor',[1 1 1])
    end
    load_imlist(handles);
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

    handles.data.outdirec = uigetdir(D,'Location for PIV Output');
    if handles.data.outdirec==0
        handles.data.outdirec = D;
    end
    set(handles.outputdirectory,'string',handles.data.outdirec);
    if ~isempty(dir(handles.data.outdirec))
        set(handles.outputdirectory,'Backgroundcolor',[1 1 1])
    end
    load_imlist(handles);
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

% --- Min Number of Iterations for Deformation ---
function deform_min_iter_Callback(hObject, eventdata, handles)
if str2double(handles.Njob)>0
    eval(['handles.data.PIV' handles.data.cpass '.deform_min=get(hObject,''String'');']);
    if str2double(get(handles.deform_max_iter,'String'))>1 && eval(['str2double(handles.data.PIV' handles.data.cpass '.deform_min) == 1'])
        eval(['handles.data.PIV' handles.data.cpass '.deform_min = ''2'';'])
        set(handles.deform_min_iter,'String','2')
    end
    guidata(hObject,handles)
end
function deform_min_iter_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Max Number of Iterations for Deformation ---
function deform_max_iter_Callback(hObject, eventdata, handles)
if str2double(handles.Njob)>0
    eval(['handles.data.PIV' handles.data.cpass '.deform_max=get(hObject,''String'');']);
    if str2double(get(hObject,'String'))>1 && eval(['str2double(handles.data.PIV' handles.data.cpass '.deform_min) == 1'])
        eval(['handles.data.PIV' handles.data.cpass '.deform_min = ''2'';'])
        set(handles.deform_min_iter,'String','2')
    end
    guidata(hObject,handles)
end
function deform_max_iter_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Convergence value of Iterations for Deformation ---
function deform_conv_Callback(hObject, eventdata, handles)
if str2double(handles.Njob)>0
    eval(['handles.data.PIV' handles.data.cpass '.deform_conv=get(hObject,''String'');']);
    guidata(hObject,handles)
end
function deform_conv_CreateFcn(hObject, eventdata, handles)
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
        set(handles.deform_min_iter,'backgroundcolor',0.5*[1 1 1]);
        set(handles.deform_max_iter,'backgroundcolor',0.5*[1 1 1]);
        set(handles.deform_conv,'backgroundcolor',0.5*[1 1 1]);
        set(handles.smoothingsize,'backgroundcolor',0.5*[1 1 1]);
        set(handles.framestep,'backgroundcolor',0.5*[1 1 1]);
        set(handles.PIVerror,'backgroundcolor',0.5*[1 1 1]);
        for e=1:str2double(handles.data.passes)
            eval(['handles.data.PIV' num2str(e) '.gridres = get(handles.gridres,''String'');'])
            eval(['handles.data.PIV' num2str(e) '.gridbuf = get(handles.gridbuffer,''String'');'])
        end
    elseif get(hObject,'Value')==6
        set(handles.velocityinterptype,'backgroundcolor',[1 1 1]);
        set(handles.imageinterptype,'backgroundcolor',0.5*[1 1 1]);
        set(handles.deform_min_iter,'backgroundcolor',0.5*[1 1 1]);
        set(handles.deform_max_iter,'backgroundcolor',0.5*[1 1 1]);
        set(handles.deform_conv,'backgroundcolor',0.5*[1 1 1]);
        set(handles.smoothingsize,'backgroundcolor',0.5*[1 1 1]);
        set(handles.framestep,'backgroundcolor',[1 1 1]);
        set(handles.PIVerror,'backgroundcolor',[1 1 1]);
    else
        set(handles.velocityinterptype,'backgroundcolor',[1 1 1]);
        set(handles.framestep,'backgroundcolor',0.5*[1 1 1]);
        set(handles.PIVerror,'backgroundcolor',0.5*[1 1 1]);
        if any(get(hObject,'Value')==[3 5 7])
            set(handles.imageinterptype,'backgroundcolor',[1 1 1]);
            set(handles.deform_min_iter,'backgroundcolor',[1 1 1]);
            set(handles.deform_max_iter,'backgroundcolor',[1 1 1]);
            set(handles.deform_conv,'backgroundcolor',[1 1 1]);
        else
            set(handles.imageinterptype,'backgroundcolor',0.5*[1 1 1]);
            set(handles.deform_min_iter,'backgroundcolor',0.5*[1 1 1]);
            set(handles.deform_max_iter,'backgroundcolor',0.5*[1 1 1]);
            set(handles.deform_conv,'backgroundcolor',0.5*[1 1 1]);
        end
        if get(handles.smoothingcheckbox,'Value')==1
            set(handles.smoothingsize,'backgroundcolor',[1 1 1]);
        else
            set(handles.smoothingsize,'backgroundcolor',0.5*[1 1 1]);
        end 
    end

    if any(get(hObject,'Value')==[4 5]) && strcmp(handles.data.masktype,'dynamic')
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
% If a job has been created or loaded...
if str2double(handles.Njob) > 0
        
    % Read the string in the window resolution textbox
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
        [xbin ybin] = roiSize(A); %#ok<NASGU,ASGLU>
        
    % Set "window size" fields in "PIV pass" data structure to the window sizes calculated above
        eval(['handles.data.PIV' handles.data.cpass '.winsize = [num2str(xbin) '','' num2str(ybin)];'])
        eval(['set(handles.windowsize,''String'', handles.data.PIV' handles.data.cpass '.winsize)']) 
    end
    
    if get(handles.correlationtype,'Value') == 7 %if DCC
        %Update "Actual Window Size" text box based on "Window Resolution"
        Rx = wx1 + wx2 - 1;
        Ry = wy1 + wy2 - 1;
        %update the text on the GUI
        set(handles.windowsize,'String',[num2str(Rx) ',' num2str(Ry)]);
        % Update the "window size" field of the "PIV Pass" data structure to the value of the 
        % string in the "Actual Window Size" text box
        eval(['handles.data.PIV' handles.data.cpass '.winsize = [num2str(Rx) '','' num2str(Ry)];']);    
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
        
        % Set overlaps to zero if overlaps are calculated to be less than zero
        overX = overX * (overX > 0); %#ok<NASGU>
        overY = overY * (overY > 0); %#ok<NASGU>
        
        % Set "overlap" fields in "PIV pass" data structure to the overlaps calculated above
        eval(['handles.data.PIV' handles.data.cpass '.winoverlap = [num2str(overX),'','',num2str(overY)];'])
        
    % Otherwise, if specifying window overlap....
    else
        
        % Read string in "window overlap" textbox
        A = get(handles.winoverlap,'String');
        
        % Parse string in "window overlap" textbox to determine the x- and y- window overlaps (%)
        [overX overY] = parseNum(A);
        
        % Calculate x- and y- grid resolutions from overlaps 
        gx = round(wxMax * (1-overX/100)); %#ok<NASGU>
        gy = round(wyMax * (1-overY/100)); %#ok<NASGU>
        
        % Update "grid resolution" field in "PIV pass" data structure to the grid resolutions calculated above
        eval(['handles.data.PIV' handles.data.cpass '.gridres = [num2str(gx),'','',num2str(gy)];'])
    end
    
    % Update "window resolution" fields in "PIV pass" data structure to the
    % window resolutions calculated above
    % eval(['handles.data.PIV' handles.data.cpass '.winres1 = [num2str(wx1) '','' num2str(wy1)];']); 
    % eval(['handles.data.PIV' handles.data.cpass '.winres2 = [num2str(wx2) '','' num2str(wy2)];']); 
    
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

        % If the x-dimension of the interrogation region is smaller than that of the effective window resolution,
        % re-size the x-dimension of the interrogation regions to equal that of the effective window resolution
        if Rx < wxMax
            Rx = wxMax;
            set(hObject,'String',[num2str(Rx) ',' num2str(Ry)]);
        end
        
        % If the y-dimension of the interrogation region is smaller than that of the effective window resolution,
        % re-size the y-dimension of the interrogation regions to equal that of the effective window resolution
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
    [xROI yROI] = roiSize(A); %#ok<NASGU,ASGLU>
    
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
    
    % If performing a multi-pass job...
    if get(handles.setgridresbutton,'Value') == 1
        
        % Read string in grid resolution texbox
        A = get(hObject,'String');
        
        % Parse string in "grid resolution" text box to determine x- and y- grid resolutions
        [gx gy] = parseNum(A);

        % If performing a multipass run....
        if str2double(handles.data.method) == 1
            
            % Then calculate the window overlaps for each pass...
            for e=1:str2double(handles.data.passes)
                
                % Read window resolutions from data structure
                % eval(['A=handles.data.PIV', num2str(e) '.winres1;'])
                % eval(['B=handles.data.PIV', num2str(e) '.winres2;'])
                eval(['A=handles.data.PIV', num2str(e) '.winres;'])
                
                % Parse the strings containing the window resolution information for each image
                % [wx1 wy1] = parseNum(A);
                % [wx2 wy2] = parseNum(B);
                [wx1 wy1 wx2 wy2] = parseNum(A);

                % Determine the size of the largest window resolution in each dimension
                wxMax = max(wx1, wx2);
                wyMax = max(wy1, wy2);
                
                % Calculate the x- and y- window overlaps
                overX = (wxMax - gx) / wxMax * 100;
                overY = (wyMax - gy) / wyMax * 100;
                
                % Set overlaps to zero if "overlap" is calculated to be less than zero
                overX = overX * (overX > 0); %#ok<NASGU>
                overY = overY * (overY > 0); %#ok<NASGU>
                
                % Update the "grid resolution" and "window overlap" fields in the "pass" data structure
                eval(['handles.data.PIV' num2str(e) '.gridres = [num2str(gx) '','' num2str(gy)];'])
                eval(['handles.data.PIV' num2str(e) '.winoverlap = [num2str(overX),'','',num2str(overY)];'])
                
            % End of "multipass" case...
            end
            
        % Otherwise, If NOT performing a multipass run...  
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
            overX = overX * (overX > 0); %#ok<NASGU>
            overY = overY * (overY > 0); %#ok<NASGU>
            
            % Update the "grid resolution" and "window overlap" fields in the "pass" data structure   
            eval(['handles.data.PIV' handles.data.cpass '.gridres = [num2str(gx) '','' num2str(gy)];'])
            eval(['handles.data.PIV' handles.data.cpass '.winoverlap = [num2str(overX),'','',num2str(overY)];'])
        end
        
        % Update GUI
        handles=set_PIVcontrols(handles);
        guidata(hObject,handles)
    else
        
        % If specifying window overlap rather than grid resolution, just read the
        % grid resolution. 
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
        gx = round(wxMax * (1 - overX / 100)); %#ok<NASGU>
        gy = round(wyMax * (1 - overY / 100)); %#ok<NASGU>
        
        % Update "Window Overlap" field in "Piv Pass" data structure
        eval(['handles.data.PIV' handles.data.cpass '.winoverlap = [num2str(overX) '','' num2str(overY)];'])
        
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
    [Bufx Bufy] = parseNum(get(hObject,'string')); 
    if str2double(handles.data.method)==1
        for e=1:str2double(handles.data.passes)
            eval(['handles.data.PIV' num2str(e) '.gridbuf = [num2str(Bufx) '','' num2str(Bufy)];'])
        end
    else
        eval(['handles.data.PIV' handles.data.cpass '.gridbuf = [num2str(Bufx) '','' num2str(Bufy)];'])
    end
    handles=set_PIVcontrols(handles);
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
        [BWOx BWOy] = parseNum(get(hObject,'string')); %#ok<ASGLU,NASGU> 
        eval(['handles.data.PIV' handles.data.cpass '.BWO = [num2str(BWOx) '','' num2str(BWOy)];'])
        handles=set_PIVcontrols(handles);
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
    ctype = {'SCC','RPC','DRPC','GCC','FWC','SPC','DCC'};
    eval(['handles.data.PIV' handles.data.cpass '.corr = ctype{get(hObject,''Value'')};'])
    N=handles.data.cpass;
    A=eval(['handles.data.PIV' num2str(N)]);
    if any(strcmpi(A.corr,{'SCC','GCC','DRPC'}))%any(str2double(A.corr)== [1 3]) %2 and 5 are RPC and SPC, both need rpcdiameter
        set(handles.rpcdiameter,'backgroundcolor',0.5*[1 1 1]);
    else
        set(handles.rpcdiameter,'backgroundcolor',[1 1 1]);
    end
    if strcmpi(A.corr,'FWC')
        set(handles.frac_filter_weight,'backgroundcolor',[1 1 1]);
        set(handles.rpcdiameter,'backgroundcolor',0.5*[1 1 1]);
    else
        set(handles.frac_filter_weight,'backgroundcolor',0.5*[1 1 1]);
    end
    %For DCC, we don't allow windowing, and the actual window size is fixed based on the resolution
    if strcmpi(A.corr,{'DCC'})
        %disable windowing
        set(handles.autowinsizecheckbox,'enable','off');
        set(handles.autowinsizecheckbox,'value',0.0);
        eval(['handles.data.PIV' handles.data.cpass '.winauto = num2str(get(handles.autowinsizecheckbox,''Value''));']);
        %Update "Actual Window Size" text box based on "Window Resolution"
        % Read the string in the window resolution textbox
        A = get(handles.windowres,'String');
        % Parse string in window resolution textbox to determine the resolutions of each window
        [wx1 wy1 wx2 wy2] = parseNum(A);
        Rx = wx1 + wx2 - 1;
        Ry = wy1 + wy2 - 1;
        %update the text on the GUI
        set(handles.windowsize,'String',[num2str(Rx) ',' num2str(Ry)]);
        % Update the "window size" field of the "PIV Pass" data structure to the value of the 
        % string in the "Actual Window Size" text box
        eval(['handles.data.PIV' handles.data.cpass '.winsize = [num2str(Rx) '','' num2str(Ry)];']);    
    else
        set(handles.autowinsizecheckbox,'enable','on');
    end
    set(handles.correlationtype,'backgroundcolor',[1 1 1]);
    handles=set_PIVcontrols(handles);
    guidata(hObject,handles)
end
function correlationtype_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- RPC Diameter Text Box ---
function rpcdiameter_Callback(hObject, eventdata, handles)
if str2double(handles.Njob)>0
        
    % Parse the string to find the differnt RPC diameter for the X and Y
    % directions.
    [RPCdx RPCdy] = parseNum(get(hObject,'string'));
    
    % Put the new values into the job file.
    eval(['handles.data.PIV' handles.data.cpass '.RPCd = [num2str(RPCdx) '','' num2str(RPCdy)];'])
    
    % Update the boxes in the GUI
    handles=set_PIVcontrols(handles);
    
    guidata(hObject,handles)
    if any(get(handles.correlationtype,'Value')==[2 5]) %need to check for RPC, SPC
        if RPCdx < 2 || RPCdy < 2%str2double(get(hObject,'String'))<2
            if RPCdx <= 0 || RPCdy <= 0%str2double(get(hObject,'String'))==0
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

% --- Fractional Filter Text Box ---
function frac_filter_weight_Callback(hObject, eventdata, handles)
if str2double(handles.Njob)>0
    eval(['handles.data.PIV' handles.data.cpass '.frac_filt = get(hObject,''String'');'])
    guidata(hObject,handles)
    if get(handles.correlationtype,'Value')==4
        if str2double(get(hObject,'String'))==0
            set(hObject,'backgroundcolor',[1 0.5 0]);
        else
            set(hObject,'backgroundcolor',[1 1 1]);
        end
    else
        set(hObject,'backgroundcolor',0.5*[1 1 1]);
    end
end
function frac_filter_weight_CreateFcn(hObject, eventdata, handles)
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

% --- Correlation Peak Threshold Check Box ---
function corrpeakthreshold_checkbox_Callback(hObject, eventdata, handles)
if str2double(handles.Njob)>0
    eval(['handles.data.PIV' handles.data.cpass '.corrpeaktest = num2str(get(hObject,''Value''));'])
    handles=set_PIVcontrols(handles);
    guidata(hObject,handles)
end

% --- Correlation Peak Absolute Threshold Text Box ---
function corrpeak_absthresh_Callback(hObject, eventdata, handles)
if str2double(handles.Njob)>0 && get(handles.validatecheckbox,'Value')==1
    eval(['handles.data.PIV' handles.data.cpass '.corrpeak_absthresh = get(hObject,''String'');'])
    guidata(hObject,handles)
end

% --- Correlation Peak Ratio Threshold Text Box ---
function corrpeak_ratiothresh_Callback(hObject, eventdata, handles)
if str2double(handles.Njob)>0 && get(handles.validatecheckbox,'Value')==1
    eval(['handles.data.PIV' handles.data.cpass '.corrpeak_ratiothresh = get(hObject,''String'');'])
    guidata(hObject,handles)
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

% --- Executes on button press in saveplane.
function saveplane_Callback(hObject, eventdata, handles)
if str2double(handles.Njob)>0
    eval(['handles.data.PIV' handles.data.cpass '.saveplane = num2str(get(hObject,''Value''));'])
    if get(hObject,'Value') == 1
        splash=splashdlg_planes({...
            'You have selected to save all of the correlation planes.' ...
            '' ...
            'CAUTION: This option is extremely memory and storage intensive.'...
            'Selecting this option can freeze your computer and/or fill your hard drive'...
            ''...
            'If you use this option, it is advised that you select only a '...
            'small number of points and frames to limit the number of correlation planes that will be exported'
            },...
            'Saving Correlation Planes','Ok','Ok');
    end
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
    if str2double(get(hObject,'String'))<=0 || isnan(str2double(get(hObject,'String')))
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
    if str2double(get(hObject,'String'))<=0 || isnan(str2double(get(hObject,'String')))
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
    if str2double(get(hObject,'String'))<=0 || isnan(str2double(get(hObject,'String')))
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
    set(handles.corrpeak_absthresh,'String','','backgroundcolor',0.5*[1 1 1]);
    set(handles.corrpeak_ratiothresh,'String','','backgroundcolor',0.5*[1 1 1]);
    set(handles.corrpeaknum,'backgroundcolor',0.5*[1 1 1]);
    set(handles.outputbasename,'String','','backgroundcolor',0.5*[1 1 1]);
    set(handles.outputpassbasename,'String','','backgroundcolor',0.5*[1 1 1]);
    set(handles.outputdirectory,'String','','backgroundcolor',0.5*[1 1 1]);
    set(handles.imagedirectory,'String','','backgroundcolor',0.5*[1 1 1]);
    set(handles.imagebasename,'String','','backgroundcolor',0.5*[1 1 1]);
    set(handles.imagedirectory2,'String','','backgroundcolor',0.5*[1 1 1]);
    set(handles.imagebasename2,'String','','backgroundcolor',0.5*[1 1 1]);
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
    set(handles.frac_filter_weight,'String','','backgroundcolor',0.5*[1 1 1]);
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
    set(handles.input_velocity,'String','','backgroundcolor',0.5*[1 1 1]);
    set(handles.input_veldirec,'String','','backgroundcolor',0.5*[1 1 1]);
    set(handles.input_velbase,'String','','backgroundcolor',0.5*[1 1 1]);
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
    set(handles.cameranumber_popupMenu,'Value',1,'backgroundcolor',0.5*[1 1 1]);
    % --- Tracking ---
    set(handles.idmethod,'backgroundcolor',0.5*[1 1 1]);
    set(handles.idimthresh,'String','','backgroundcolor',0.5*[1 1 1]);
    set(handles.idsavebase,'String','','backgroundcolor',0.5*[1 1 1]);
    set(handles.idsaveloc,'String','','backgroundcolor',0.5*[1 1 1]);
    set(handles.sizingmethod,'backgroundcolor',0.5*[1 1 1]);
    set(handles.sizingstd,'String','','backgroundcolor',0.5*[1 1 1]);
    set(handles.sizing_min_area,'String','','backgroundcolor',0.5*[1 1 1]);
    set(handles.sizingsavebase,'String','','backgroundcolor',0.5*[1 1 1]);
    set(handles.sizingsaveloc,'String','','backgroundcolor',0.5*[1 1 1]);
    set(handles.trackingmethod,'backgroundcolor',0.5*[1 1 1]);
    set(handles.trackingprediction,'backgroundcolor',0.5*[1 1 1]);
    set(handles.trackingPIVweight,'String','','backgroundcolor',0.5*[1 1 1]);
    set(handles.trackingradius,'String','','backgroundcolor',0.5*[1 1 1]);
    set(handles.trackingdistweight,'String','','backgroundcolor',0.5*[1 1 1]);
    set(handles.trackingsizeweight,'String','','backgroundcolor',0.5*[1 1 1]);
    set(handles.trackingintensityweight,'String','','backgroundcolor',0.5*[1 1 1]);
    set(handles.trackingestradius,'String','','backgroundcolor',0.5*[1 1 1]);
    set(handles.trackingestweight,'String','','backgroundcolor',0.5*[1 1 1]);
    set(handles.trackingvectors,'String','','backgroundcolor',0.5*[1 1 1]);
    set(handles.trackingiterations,'String','','backgroundcolor',0.5*[1 1 1]);
    set(handles.trackingsavebase,'String','','backgroundcolor',0.5*[1 1 1]);
    set(handles.trackingsaveloc,'String','','backgroundcolor',0.5*[1 1 1]);
    set(handles.trackingvalcoefficient,'String','','backgroundcolor',0.5*[1 1 1]);
    set(handles.trackingvalradius,'String','','backgroundcolor',0.5*[1 1 1]);
    set(handles.trackingvalUthresh,'String','','backgroundcolor',0.5*[1 1 1]);
    set(handles.trackingvalVthresh,'String','','backgroundcolor',0.5*[1 1 1]);
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
    set(handles.corrpeak_absthresh,'String','','backgroundcolor',[1 1 1]);
    set(handles.corrpeak_ratiothresh,'String','','backgroundcolor',[1 1 1]);
    set(handles.uod_thresh,'String','','backgroundcolor',[1 1 1]);
    set(handles.corrpeaknum,'backgroundcolor',[1 1 1]);
    set(handles.outputbasename,'String','','backgroundcolor',[1 1 1]);
	set(handles.outputpassbasename,'String','','backgroundcolor',[1 1 1]);
    set(handles.imagedirectory,'String','','backgroundcolor',[1 1 1]);
    set(handles.imagedirectory2,'String','','backgroundcolor',[1 1 1]);
    set(handles.outputdirectory,'String','','backgroundcolor',[1 1 1]);
    set(handles.imagebasename,'String','','backgroundcolor',[1 1 1]);
    set(handles.imagebasename2,'String','','backgroundcolor',[1 1 1]);
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
    set(handles.deform_min_iter,'String','','backgroundcolor',[1 1 1]);
    set(handles.deform_max_iter,'String','','backgroundcolor',[1 1 1]);
    set(handles.deform_conv,'String','','backgroundcolor',[1 1 1]);
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
    set(handles.cameranumber_popupMenu,'Value',1,'backgroundcolor',[1 1 1]);
    % --- Tracking ---
    set(handles.idmethod,'Value',1,'backgroundcolor',[1 1 1]);
    set(handles.idimthresh,'String','','backgroundcolor',[1 1 1]);
    set(handles.idsavebase,'String','','backgroundcolor',[1 1 1]);
    set(handles.idsaveloc,'String','','backgroundcolor',[1 1 1]);
    set(handles.sizingmethod,'Value',1,'backgroundcolor',[1 1 1]);
    set(handles.sizingstd,'String','','backgroundcolor',[1 1 1]);
    set(handles.sizing_min_area,'String','','backgroundcolor',[1 1 1]);
    set(handles.sizingsavebase,'String','','backgroundcolor',[1 1 1]);
    set(handles.sizingsaveloc,'String','','backgroundcolor',[1 1 1]);
    set(handles.trackingmethod,'Value',1,'backgroundcolor',[1 1 1]);
    set(handles.trackingprediction,'Value',1,'backgroundcolor',[1 1 1]);
    set(handles.trackingPIVweight,'String','','backgroundcolor',[1 1 1]);
    set(handles.trackingradius,'String','','backgroundcolor',[1 1 1]);
    set(handles.trackingdistweight,'String','','backgroundcolor',[1 1 1]);
    set(handles.trackingsizeweight,'String','','backgroundcolor',[1 1 1]);
    set(handles.trackingintensityweight,'String','','backgroundcolor',[1 1 1]);
    set(handles.trackingestradius,'String','','backgroundcolor',[1 1 1]);
    set(handles.trackingestweight,'String','','backgroundcolor',[1 1 1]);
    set(handles.trackingvectors,'String','','backgroundcolor',[1 1 1]);
    set(handles.trackingiterations,'String','','backgroundcolor',[1 1 1]);
    set(handles.trackingsavebase,'String','','backgroundcolor',[1 1 1]);
    set(handles.trackingsaveloc,'String','','backgroundcolor',[1 1 1]);
    set(handles.trackingvalcoefficient,'String','','backgroundcolor',[1 1 1]);
    set(handles.trackingvalradius,'String','','backgroundcolor',[1 1 1]);
    set(handles.trackingvalUthresh,'String','','backgroundcolor',[1 1 1]);
    set(handles.trackingvalVthresh,'String','','backgroundcolor',[1 1 1]);

    if str2double(handles.data.par)==1
        set(handles.parprocessors,'string','','backgroundcolor',[1 1 1])
    else
        set(handles.parprocessors,'string','','backgroundcolor',0.5*[1 1 1])
    end
    
    if strcmpi(handles.data.input_vel_type,'static')
        set(handles.input_velocity,'String','','backgroundcolor',[1 1 1]);
        set(handles.input_veldirec,'String','','backgroundcolor',0.5*[1 1 1]);
        set(handles.input_velbase,'String','','backgroundcolor',0.5*[1 1 1]);
    elseif strcmpi(handles.data.input_vel_type,'dynamic')
        set(handles.input_velocity,'String','','backgroundcolor',0.5*[1 1 1]);
        set(handles.input_veldirec,'String','','backgroundcolor',[1 1 1]);
        set(handles.input_velbase,'String','','backgroundcolor',[1 1 1]);
    else %'none'
        set(handles.input_velocity,'String','','backgroundcolor',0.5*[1 1 1]);
        set(handles.input_veldirec,'String','','backgroundcolor',0.5*[1 1 1]);
        set(handles.input_velbase,'String','','backgroundcolor',0.5*[1 1 1]);
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
dir_struct = dir(fullfile(handles.data.imdirec,['*.' handles.data.imext]));
dir_struct2 = dir(fullfile(handles.data.imdirec2,['*.' handles.data.imext]));
if isempty(dir_struct)
    set(handles.imagedirectory,'backgroundcolor','r');
    %set(handles.imagedirectory2,'backgroundcolor','r');
else
    set(handles.imagedirectory,'backgroundcolor',[1 1 1]);
    %set(handles.imagedirectory2,'backgroundcolor',[1 1 1]);
end
if isempty(dir_struct2)
    %set(handles.imagedirectory,'backgroundcolor','r');
    set(handles.imagedirectory2,'backgroundcolor','r');
else
    %set(handles.imagedirectory,'backgroundcolor',[1 1 1]);
    set(handles.imagedirectory2,'backgroundcolor',[1 1 1]);
end
N=length(str2double(handles.data.imfstart):str2double(handles.data.imfstep):str2double(handles.data.imfend));
files = cell(N,2);
e=0;
for f=str2double(handles.data.imfstart):str2double(handles.data.imfstep):str2double(handles.data.imfend)
    e=e+1;
    files(e,1)={[handles.data.imbase sprintf(['%0.' handles.data.imzeros 'i.' handles.data.imext],f)]};
    files(e,2)={[handles.data.imbase2 sprintf(['%0.' handles.data.imzeros 'i.' handles.data.imext],f+str2double(handles.data.imcstep))]};
end

% %JJC: as far as I can tell, pranaprocessing correlates imdir/imageNNN1 to
% imdir/imageNNN2, NOT all of imdir/imageNNN1 to imdir/imageNNN2 and then
% imdir2/imageNNN2 to imdir2/imageNNN2
% for f=str2double(handles.data.imfstart):str2double(handles.data.imfstep):str2double(handles.data.imfend)
%     e=e+1;
%     files(e,1)={[handles.data.imbase2 sprintf(['%0.' handles.data.imzeros 'i.' handles.data.imext],f)]};
%     files(e,2)={[handles.data.imbase2 sprintf(['%0.' handles.data.imzeros 'i.' handles.data.imext],f+str2double(handles.data.imcstep))]};
% end

[sorted_names,sorted_index] = sortrows({dir_struct.name}');
[sorted_names2,~] = sortrows({dir_struct2.name}');
% sorted_names=[sorted_names;sorted_names2];
[files1,id,id1] = intersect(sorted_names ,files(:,1));
[files2,id,id2] = intersect(sorted_names2,files(:,2));
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
if strcmpi(A.corr,'SCC');
    corr = 1;
elseif strcmpi(A.corr,'RPC');
    corr = 2;
elseif strcmpi(A.corr,'DRPC');
    corr = 3;
elseif strcmpi(A.corr,'GCC');
    corr = 4;
elseif strcmpi(A.corr,'FWC');
    corr = 5;
elseif strcmpi(A.corr,'SPC');
    corr = 6;
elseif strcmpi(A.corr,'DCC');
    corr = 7;
end
set(handles.windowres,'string',A.winres);
set(handles.windowsize,'string',A.winsize);
set(handles.autowinsizecheckbox,'Value',str2double(A.winauto));
set(handles.gridres,'string',A.gridres);
set(handles.winoverlap,'string',A.winoverlap);
set(handles.gridbuffer,'string',A.gridbuf);
set(handles.bulkwinoffset,'string',A.BWO);
set(handles.correlationtype,'Value',corr);
set(handles.subpixelinterp,'Value',str2double(A.peaklocator));
set(handles.zeromeancheckbox,'Value',str2double(A.zeromean));
set(handles.rpcdiameter,'string',A.RPCd);
set(handles.frac_filter_weight,'string',str2double(A.frac_filt));
set(handles.deform_min_iter,'string',str2double(A.deform_min));
set(handles.deform_max_iter,'string',str2double(A.deform_max));
set(handles.deform_conv,'string',str2double(A.deform_conv));
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
set(handles.corrpeak_absthresh,'string',A.corrpeak_absthresh);
set(handles.corrpeak_ratiothresh,'string',A.corrpeak_ratiothresh);
set(handles.valextrapeaks,'value',str2double(A.valextrapeaks));
set(handles.corrpeaknum,'value',str2double(A.corrpeaknum));
set(handles.savepeakinfo,'value',str2double(A.savepeakinfo));
set(handles.savepeakmag,'value',str2double(A.savepeakmag));
set(handles.savepeakvel,'value',str2double(A.savepeakvel));
set(handles.saveplane,'value',str2double(A.saveplane));
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
    if str2double(A.corrpeaktest)==1
        set(handles.corrpeakthreshold_checkbox,'Value',1);
        set(handles.corrpeak_absthresh,'backgroundcolor',[1 1 1]);
        set(handles.corrpeak_ratiothresh,'backgroundcolor',[1 1 1]);
    else
        set(handles.corrpeakthreshold_checkbox,'Value',0);
        set(handles.corrpeak_absthresh,'backgroundcolor',0.5*[1 1 1]);
        set(handles.corrpeak_ratiothresh,'backgroundcolor',0.5*[1 1 1]);
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
    set(handles.corrpeak_absthresh,'backgroundcolor',0.5*[1 1 1]);
    set(handles.corrpeak_ratiothresh,'backgroundcolor',0.5*[1 1 1]);
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
if any(get(handles.correlationtype,'Value')==[2 6]) %check diameter if RPC, SPC
    if str2double(get(handles.rpcdiameter,'String'))<2
        if str2double(get(handles.rpcdiameter,'String'))==0
            set(handles.rpcdiameter,'backgroundcolor','r');
        else
            set(handles.rpcdiameter,'backgroundcolor',[1 0.5 0]);
        end
    else
        set(handles.rpcdiameter,'backgroundcolor',[1 1 1]);
    end
    set(handles.frac_filter_weight,'backgroundcolor',0.5.*[1 1 1]);
else
    if get(handles.correlationtype,'Value')==5
        set(handles.frac_filter_weight,'backgroundcolor',[1 1 1]);
    else
        set(handles.frac_filter_weight,'backgroundcolor',0.5.*[1 1 1]);
    end
    set(handles.rpcdiameter,'backgroundcolor',0.5*[1 1 1]);
end
%May need to check for DCC and disable checkbox, force windowres?
if get(handles.correlationtype,'Value')==7 % check if DCC
    set(handles.autowinsizecheckbox,'enable','off');
    set(handles.windowsize,'backgroundcolor',0.5*[1 1 1]); %override color, because it's really off
    %string value is set above
else
    %if it's not DCC, make sure the autowinsize checkbox works
    set(handles.autowinsizecheckbox,'enable','on');
    %state of winsize control is updated above based on value of autowin
end

% Grays out the smoothing filt size box when smoothing is not selected when
% switching between passes.
if str2double(A.velsmooth)  == 0 || any(get(handles.passtype,'Value')==[1 6])
    set(handles.smoothingsize,'backgroundcolor',0.5*[1 1 1]);
else
    set(handles.smoothingsize,'backgroundcolor',[1 1 1]);
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


function [handles]=update_PTV(handles)
if str2double(handles.data.ID.runid)
    set(handles.runidcheckbox,'Value',str2double(handles.data.ID.runid));
    set(handles.idmethod,'Value',str2double(handles.data.ID.method),'backgroundcolor',[1 1 1]);
    set(handles.idimthresh,'String',handles.data.ID.imthresh,'backgroundcolor',[1 1 1]);
    set(handles.idsavebase,'String',handles.data.ID.savebase,'backgroundcolor',[1 1 1]);
%     set(handles.idsaveloc,'String',handles.data.ID.save_dir,'backgroundcolor',[1 1 1]);
    if exist(handles.data.ID.save_dir,'dir')
        set(handles.idsaveloc,'String',handles.data.ID.save_dir,'backgroundcolor',[1 1 1]);
    else
        set(handles.idsaveloc,'String',handles.data.ID.save_dir,'backgroundcolor','r');
    end
else
    set(handles.idmethod,'backgroundcolor',0.5*[1 1 1]);
    set(handles.idimthresh,'backgroundcolor',0.5*[1 1 1]);
    set(handles.idsavebase,'backgroundcolor',0.5*[1 1 1]);
    set(handles.idsaveloc,'backgroundcolor',0.5*[1 1 1]);
end
if str2double(handles.data.Size.runsize)
    set(handles.runsizingcheckbox,'Value',str2double(handles.data.Size.runsize));
    set(handles.sizingmethod,'Value',sizeMethodLookup(handles.data.Size.method),'backgroundcolor',[1 1 1]);
    set(handles.sizing_min_area,'String',handles.data.Size.min_area,'backgroundcolor',[1 1 1]);
    set(handles.sizingstd,'String',handles.data.Size.std,'backgroundcolor',[1 1 1]);
    set(handles.sizingsavebase,'String',handles.data.Size.savebase,'backgroundcolor',[1 1 1]);
    set(handles.sizingsaveloc,'String',handles.data.Size.save_dir,'backgroundcolor',[1 1 1]);
    if exist(handles.data.Size.save_dir,'dir')
        set(handles.sizingsaveloc,'String',handles.data.Size.save_dir,'backgroundcolor',[1 1 1]);
    else
        set(handles.sizingsaveloc,'String',handles.data.Size.save_dir,'backgroundcolor','r');
    end
else
    set(handles.sizingmethod,'backgroundcolor',0.5*[1 1 1]);
    set(handles.sizing_min_area,'backgroundcolor',0.5*[1 1 1]);
    set(handles.sizingstd,'backgroundcolor',0.5*[1 1 1]);
    set(handles.sizingsavebase,'backgroundcolor',0.5*[1 1 1]);
    set(handles.sizingsaveloc,'backgroundcolor',0.5*[1 1 1]);
end
if str2double(handles.data.Track.runtrack)
    set(handles.runtrackingcheckbox,'Value',str2double(handles.data.Track.runtrack));
    set(handles.trackingmethod,'Value',str2double(handles.data.Track.method),'backgroundcolor',[1 1 1]);
    set(handles.trackingprediction,'Value',str2double(handles.data.Track.prediction),'backgroundcolor',[1 1 1]);
    set(handles.trackingPIVweight,'String',handles.data.Track.PIVweight,'backgroundcolor',[1 1 1]);
    set(handles.trackingradius,'String',handles.data.Track.radius,'backgroundcolor',[1 1 1]);
    set(handles.trackingdistweight,'String',handles.data.Track.disweight,'backgroundcolor',[1 1 1]);
    set(handles.trackingsizeweight,'String',handles.data.Track.sizeweight,'backgroundcolor',[1 1 1]);
    set(handles.trackingintensityweight,'String',handles.data.Track.intensityweight,'backgroundcolor',[1 1 1]);
    set(handles.trackingestradius,'String',handles.data.Track.estradius,'backgroundcolor',[1 1 1]);
    set(handles.trackingestweight,'String',handles.data.Track.estweight,'backgroundcolor',[1 1 1]);
    set(handles.trackingvectors,'String',handles.data.Track.vectors,'backgroundcolor',[1 1 1]);
    set(handles.trackingiterations,'String',handles.data.Track.iterations,'backgroundcolor',[1 1 1]);
    set(handles.trackingsavebase,'String',handles.data.Track.savebase,'backgroundcolor',[1 1 1]);
    set(handles.trackingsaveloc,'String',handles.data.Track.save_dir,'backgroundcolor',[1 1 1]);
    if exist(handles.data.Track.save_dir,'dir')
        set(handles.trackingsaveloc,'String',handles.data.Track.save_dir,'backgroundcolor',[1 1 1]);
    else
        set(handles.trackingsaveloc,'String',handles.data.Track.save_dir,'backgroundcolor','r');
    end
    if str2double(handles.data.Track.valprops.run)
        set(handles.trackingvalcheckbox,'Value',str2double(handles.data.Track.valprops.run));
        set(handles.trackingvalcoefficient,'String',handles.data.Track.valprops.valcoef,'backgroundcolor',[1 1 1]);
        set(handles.trackingvalradius,'String',handles.data.Track.valprops.valrad,'backgroundcolor',[1 1 1]);
        set(handles.trackingvalUthresh,'String',handles.data.Track.valprops.MAD_U,'backgroundcolor',[1 1 1]);
        set(handles.trackingvalVthresh,'String',handles.data.Track.valprops.MAD_V,'backgroundcolor',[1 1 1]);
    else
        set(handles.trackingvalcoefficient,'backgroundcolor',0.5*[1 1 1],'String','');
        set(handles.trackingvalradius,'backgroundcolor',0.5*[1 1 1],'String','');
        set(handles.trackingvalUthresh,'backgroundcolor',0.5*[1 1 1],'String','');
        set(handles.trackingvalVthresh,'backgroundcolor',0.5*[1 1 1],'String','');
    end
    set(handles.trackingoutputdat,'Value',str2double(handles.data.Track.trackdat));
    set(handles.trackingoutputmat,'Value',str2double(handles.data.Track.trackmat));

else
    set(handles.trackingmethod,'backgroundcolor',0.5*[1 1 1]);
    set(handles.trackingprediction,'backgroundcolor',0.5*[1 1 1]);
    set(handles.trackingPIVweight,'backgroundcolor',0.5*[1 1 1]);
    set(handles.trackingradius,'backgroundcolor',0.5*[1 1 1]);
    set(handles.trackingdistweight,'backgroundcolor',0.5*[1 1 1]);
    set(handles.trackingsizeweight,'backgroundcolor',0.5*[1 1 1]);
    set(handles.trackingintensityweight,'backgroundcolor',0.5*[1 1 1]);
    set(handles.trackingestradius,'backgroundcolor',0.5*[1 1 1]);
    set(handles.trackingestweight,'backgroundcolor',0.5*[1 1 1]);
    set(handles.trackingvectors,'backgroundcolor',0.5*[1 1 1]);
    set(handles.trackingiterations,'backgroundcolor',0.5*[1 1 1]);
    set(handles.trackingsavebase,'backgroundcolor',0.5*[1 1 1]);
    set(handles.trackingsaveloc,'backgroundcolor',0.5*[1 1 1]);
    set(handles.trackingvalcoefficient,'backgroundcolor',0.5*[1 1 1]);
    set(handles.trackingvalradius,'backgroundcolor',0.5*[1 1 1]);
    set(handles.trackingvalUthresh,'backgroundcolor',0.5*[1 1 1]);
    set(handles.trackingvalVthresh,'backgroundcolor',0.5*[1 1 1]);
    set(handles.trackingoutputdat,'Value',0);
    set(handles.trackingoutputmat,'Value',0);
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
set(handles.imagedirectory2,'String',handles.data.imdirec2);
set(handles.imagebasename2,'String',handles.data.imbase2);
set(handles.imagezeros,'String',handles.data.imzeros);
set(handles.imageextension,'String',handles.data.imext);
set(handles.imagecorrelationstep,'String',handles.data.imcstep);
set(handles.imageframestep,'String',handles.data.imfstep);
set(handles.imageframestart,'String',handles.data.imfstart);
set(handles.imageframeend,'String',handles.data.imfend);

set(handles.input_velocity,'String',handles.data.input_velocity)
set(handles.input_veldirec,'String',handles.data.input_veldirec)
set(handles.input_velbase,'String',handles.data.input_velbase)

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

if strcmpi(handles.data.input_vel_type,'none')
    set(handles.no_input_vel_button,'Value',1)
    set(handles.static_input_vel_button,'Value',0)
    set(handles.dynamic_input_vel_button,'Value',0)
    set(handles.input_velocity,'backgroundcolor',0.5*[1 1 1])
    set(handles.input_veldirec,'backgroundcolor',0.5*[1 1 1])
    set(handles.input_velbase,'backgroundcolor',0.5*[1 1 1])
elseif strcmpi(handles.data.input_vel_type,'static')
    set(handles.no_input_vel_button,'Value',0)
    set(handles.static_input_vel_button,'Value',1)
    set(handles.dynamic_input_vel_button,'Value',0)
    set(handles.input_velocity,'string',handles.data.input_velocity,'backgroundcolor',[1 1 1]);
    set(handles.input_veldirec,'backgroundcolor',0.5*[1 1 1])
    set(handles.input_velbase,'backgroundcolor',0.5*[1 1 1])
elseif strcmpi(handles.data.input_vel_type,'dynamic')
    set(handles.no_input_vel_button,'Value',0)
    set(handles.static_input_vel_button,'Value',0)
    set(handles.dynamic_input_vel_button,'Value',1)
    set(handles.input_velocity,'backgroundcolor',0.5*[1 1 1])
    set(handles.input_veldirec,'string',handles.data.input_veldirec,'backgroundcolor',[1 1 1]);
    set(handles.input_velbase,'string',handles.data.input_velbase,'backgroundcolor',[1 1 1]);
else
    error('Unknown input velocity type')
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

set(handles.runPIVcheckbox,'Value',str2double(handles.data.runPIV));

set(handles.passtype,'Value',str2double(handles.data.method));
set(handles.velocityinterptype,'Value',str2double(handles.data.velinterp));
set(handles.imageinterptype,'Value',str2double(handles.data.iminterp));
set(handles.colorchannel_popupMenu,'Value',str2double(handles.data.channel));
set(handles.cameranumber_popupMenu,'Value',str2double(handles.data.numcams));
set(handles.framestep,'String',handles.data.framestep);
set(handles.PIVerror,'String',handles.data.PIVerror);

if get(handles.passtype,'Value')==6
    set(handles.velocityinterptype,'backgroundcolor',[1 1 1]);
    set(handles.imageinterptype,'backgroundcolor',0.5*[1 1 1]);
    set(handles.deform_min_iter,'backgroundcolor',0.5*[1 1 1]);
    set(handles.deform_max_iter,'backgroundcolor',0.5*[1 1 1]);
    set(handles.deform_conv,'backgroundcolor',0.5*[1 1 1]);    
    set(handles.smoothingsize,'backgroundcolor',0.5*[1 1 1]);
    set(handles.framestep,'backgroundcolor',[1 1 1]);
    set(handles.PIVerror,'backgroundcolor',[1 1 1]);
elseif get(handles.passtype,'Value')>1
    set(handles.framestep,'backgroundcolor',0.5*[1 1 1]);
    set(handles.PIVerror,'backgroundcolor',0.5*[1 1 1]);
    set(handles.velocityinterptype,'backgroundcolor',[1 1 1]);
    if any(get(handles.passtype,'Value')==[3 5 7])
        set(handles.imageinterptype,'backgroundcolor',[1 1 1]);
        set(handles.deform_min_iter,'backgroundcolor',[1 1 1]);
        set(handles.deform_max_iter,'backgroundcolor',[1 1 1]);
        set(handles.deform_conv,'backgroundcolor',[1 1 1]);
    else
        set(handles.imageinterptype,'backgroundcolor',0.5*[1 1 1]);
        set(handles.deform_min_iter,'backgroundcolor',0.5*[1 1 1]);
        set(handles.deform_max_iter,'backgroundcolor',0.5*[1 1 1]);
        set(handles.deform_conv,'backgroundcolor',0.5*[1 1 1]);
    end
    if get(handles.smoothingcheckbox,'Value')==1
        set(handles.smoothingsize,'backgroundcolor',[1 1 1]);
    else
        set(handles.smoothingsize,'backgroundcolor',0.5*[1 1 1]);
    end
else
    set(handles.velocityinterptype,'backgroundcolor',0.5*[1 1 1]);
    set(handles.imageinterptype,'backgroundcolor',0.5*[1 1 1]);
    set(handles.deform_min_iter,'backgroundcolor',0.5*[1 1 1]);
    set(handles.deform_max_iter,'backgroundcolor',0.5*[1 1 1]);
    set(handles.deform_conv,'backgroundcolor',0.5*[1 1 1]);
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
if isempty(dir(handles.data.imdirec2))
    set(handles.imagedirectory2,'backgroundcolor','r');
else
    set(handles.imagedirectory2,'backgroundcolor',[1 1 1]);
end

if ~exist(handles.data.outdirec,'dir')
    set(handles.outputdirectory,'backgroundcolor','r');
else
    set(handles.outputdirectory,'backgroundcolor',[1 1 1]);
end

if strcmpi(handles.data.masktype,'static')
    if isempty(dir(handles.data.staticmaskname))
        set(handles.staticmaskfile,'backgroundcolor','r');
    else
        set(handles.staticmaskfile,'backgroundcolor',[1 1 1]);
    end
elseif strcmpi(handles.data.masktype,'dynamic')
    if isempty(dir(handles.data.maskdirec))
        set(handles.maskdirectory,'backgroundcolor','r');
    else
        set(handles.maskdirectory,'backgroundcolor',[1 1 1]);
    end
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

% --- Tracking ---
update_PTV(handles);

    



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

function ButtonName=splashdlg_planes(Question,Title,Btn1,Btn2,Btn3,Default)
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
Orange     =[243     132       0     ]/255;
Black      =[  0       0        0    ]/255;
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
if nargin==1, Title=' ';end
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
  'Color'           ,Orange                      ...
  );

%%%%%%%%%%%%%%%%%%%%%
%%% Set Positions %%%
%%%%%%%%%%%%%%%%%%%%%
DefOffset  =10;

IconWidth  =0;
IconHeight =150;
IconXOffset=DefOffset;
IconYOffset=FigPos(4)-DefOffset-IconHeight;  %#ok
IconCMap=[Orange;get(QuestFig,'Color')];  %#ok

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


% IconAxes=axes(                                      ...
%   'Parent'      ,QuestFig              , ...
%   'Units'       ,'Pixels'              , ...
%   'Position'    ,[IconXOffset IconYOffset IconWidth IconHeight], ...
%   'NextPlot'    ,'replace'             , ...
%   'Tag'         ,'IconAxes'              ...
%   );
% 
% set(QuestFig ,'NextPlot','add');
% 
% pranadir=which('prana');
% try
%     logo=imread(fullfile(pranadir(1:end-8),'documentation','logo.tif'),'tif');
% catch
%     logo=zeros(500,1000);
% end
% Img=image('Cdata',logo,'Parent',IconAxes);
% 
% set(IconAxes, ...
%   'Visible','off'           , ...
%   'YDir'   ,'reverse'       , ...
%   'XLim'   ,get(Img,'XData'), ...
%   'YLim'   ,get(Img,'YData')  ...
%   );

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
% --- ID ---
% --------------------------------------------------------------------
% --- Executes on button press in runidcheckbox.
function runidcheckbox_Callback(hObject, eventdata, handles)
if str2double(handles.Njob)>0
    handles.data.ID.runid = num2str(get(hObject,'Value'));
    update_PTV(handles);
    guidata(hObject,handles)
end

% --- Executes on selection change in idmethod.
function idmethod_Callback(hObject, eventdata, handles)
if str2double(handles.Njob)>0
    handles.data.ID.method = num2str(get(hObject,'Value'));
    update_PTV(handles);
    guidata(hObject,handles)
end

function idmethod_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function idimthresh_Callback(hObject, eventdata, handles)
if str2double(handles.Njob)>0
    handles.data.ID.imthresh = get(hObject,'String');
    update_PTV(handles);
    guidata(hObject,handles)
end

function idimthresh_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function idsavebase_Callback(hObject, eventdata, handles)
if str2double(handles.Njob)>0
    handles.data.ID.savebase = get(hObject,'String');
    update_PTV(handles);
    guidata(hObject,handles)
end
function idsavebase_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function idsaveloc_Callback(hObject, eventdata, handles)
if str2double(handles.Njob)>0
    handles.data.ID.save_dir = get(hObject,'String');
    update_PTV(handles);
    guidata(hObject,handles)
end
function idsaveloc_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in loadidsaveloc.
function loadidsaveloc_Callback(hObject, eventdata, handles)
if str2double(handles.Njob)>0
    D = handles.data.ID.save_dir;
    handles.data.ID.save_dir = uigetdir(handles.data.ID.save_dir);
    if handles.data.ID.save_dir==0
        handles.data.ID.save_dir = D;
    end
    set(handles.idsaveloc,'string',handles.data.ID.save_dir);
    update_PTV(handles);
    guidata(hObject,handles)
end

% --------------------------------------------------------------------
% --- Sizing ---
% --------------------------------------------------------------------
% --- Executes on button press in runsizingcheckbox.
function runsizingcheckbox_Callback(hObject, eventdata, handles)
if str2double(handles.Njob)>0
    handles.data.Size.runsize = num2str(get(hObject,'Value'));
    update_PTV(handles);
    guidata(hObject,handles);
end

% --- Executes on selection change in sizingmethod.
function sizingmethod_Callback(hObject, eventdata, handles)
if str2double(handles.Njob)>0
    
    handles.data.Size.method = sizeMethodLookup(get(hObject,'Value'));
%     handles.data.Size.method = num2str(get(hObject,'Value'));
    update_PTV(handles);
    guidata(hObject,handles)
end
function sizingmethod_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function sizingstd_Callback(hObject, eventdata, handles)
if str2double(handles.Njob)>0
    handles.data.Size.std = get(hObject,'String');
    update_PTV(handles);
    guidata(hObject,handles)
end
function sizingstd_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- This is the area filter for the sizing code ---
function sizing_min_area_Callback(hObject, eventdata, handles)
if str2double(handles.Njob)>0
    handles.data.Size.min_area = get(hObject,'String');
    update_PTV(handles);
    guidata(hObject,handles)
end

% --- Executes during object creation, after setting all properties.
function sizing_min_area_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function sizingsavebase_Callback(hObject, eventdata, handles)
if str2double(handles.Njob)>0
    handles.data.Size.savebase = get(hObject,'String');
    update_PTV(handles);
    guidata(hObject,handles)
end
function sizingsavebase_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function sizingsaveloc_Callback(hObject, eventdata, handles)
if str2double(handles.Njob)>0
    handles.data.Size.save_dir = get(hObject,'String');
    update_PTV(handles);
    guidata(hObject,handles)
end
function sizingsaveloc_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in loadsizingsaveloc.
function loadsizingsaveloc_Callback(hObject, eventdata, handles)
if str2double(handles.Njob)>0
    D = handles.data.Size.save_dir;
    handles.data.Size.save_dir = uigetdir(handles.data.Size.save_dir);
    if handles.data.Size.save_dir==0
        handles.data.Size.save_dir = D;
    end
    set(handles.sizingsaveloc,'string',handles.data.Size.save_dir);
    update_PTV(handles);
    guidata(hObject,handles)
end

% --------------------------------------------------------------------
% --- Tracking ---
% --------------------------------------------------------------------
% --- Executes on button press in runtrackingcheckbox.
function runtrackingcheckbox_Callback(hObject, eventdata, handles)
if str2double(handles.Njob)>0
    handles.data.Track.runtrack = num2str(get(hObject,'Value'));
    update_PTV(handles);
    guidata(hObject,handles);
end

% --- Executes on selection change in trackingmethod.
function trackingmethod_Callback(hObject, eventdata, handles)
if str2double(handles.Njob)>0
    handles.data.Track.method = num2str(get(hObject,'Value'));
    update_PTV(handles);
    guidata(hObject,handles)
end
function trackingmethod_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
% --- Executes on selection change in trackingprediction.
function trackingprediction_Callback(hObject, eventdata, handles)
if str2double(handles.Njob)>0
    handles.data.Track.prediction = num2str(get(hObject,'Value'));
    update_PTV(handles);
    guidata(hObject,handles)
end
function trackingprediction_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function trackingPIVweight_Callback(hObject, eventdata, handles)
if str2double(handles.Njob)>0
    handles.data.Track.PIVweight = get(hObject,'String');
    update_PTV(handles);
    guidata(hObject,handles)
end
function trackingPIVweight_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function trackingradius_Callback(hObject, eventdata, handles)
if str2double(handles.Njob)>0
    handles.data.Track.radius = get(hObject,'String');
    update_PTV(handles);
    guidata(hObject,handles)
end
function trackingradius_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function trackingdistweight_Callback(hObject, eventdata, handles)
if str2double(handles.Njob)>0
    handles.data.Track.disweight = get(hObject,'String');
    if str2double(handles.data.Track.disweight)>1
        handles.data.Track.disweight = '1';        
    end
    update_PTV(handles);
    guidata(hObject,handles)
end
function trackingdistweight_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function trackingsizeweight_Callback(hObject, eventdata, handles)
if str2double(handles.Njob)>0
    handles.data.Track.sizeweight = get(hObject,'String');
    if str2double(handles.data.Track.sizeweight)>1
        handles.data.Track.sizeweight = '1';        
    end
    update_PTV(handles);
    guidata(hObject,handles)
end
function trackingsizeweight_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function trackingintensityweight_Callback(hObject, eventdata, handles)
if str2double(handles.Njob)>0
    handles.data.Track.intensityweight = get(hObject,'String');
    if str2double(handles.data.Track.intensityweight)>1
        handles.data.Track.intensityweight = '1';        
    end
    update_PTV(handles);
    guidata(hObject,handles)
end
function trackingintensityweight_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function trackingestradius_Callback(hObject, eventdata, handles)
if str2double(handles.Njob)>0
    handles.data.Track.estradius = get(hObject,'String');
    update_PTV(handles);
    guidata(hObject,handles)
end
function trackingestradius_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function trackingestweight_Callback(hObject, eventdata, handles)
if str2double(handles.Njob)>0
    handles.data.Track.estweight = get(hObject,'String');
    if str2double(handles.data.Track.estweight)>1
        handles.data.Track.estweight = '1';
    elseif str2double(handles.data.Track.estweight)<0
        handles.data.Track.estweight = '0';
    end
    update_PTV(handles);
    guidata(hObject,handles)
end
function trackingestweight_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function trackingvectors_Callback(hObject, eventdata, handles)
if str2double(handles.Njob)>0
    handles.data.Track.vectors = get(hObject,'String');
    update_PTV(handles);
    guidata(hObject,handles)
end
function trackingvectors_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function trackingiterations_Callback(hObject, eventdata, handles)
if str2double(handles.Njob)>0
    handles.data.Track.iterations = get(hObject,'String');
    update_PTV(handles);
    guidata(hObject,handles)
end
function trackingiterations_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function trackingsavebase_Callback(hObject, eventdata, handles)
if str2double(handles.Njob)>0
    handles.data.Track.savebase = get(hObject,'String');
    update_PTV(handles);
    guidata(hObject,handles)
end
function trackingsavebase_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function trackingsaveloc_Callback(hObject, eventdata, handles)
if str2double(handles.Njob)>0
    handles.data.Track.save_dir = get(hObject,'String');
    update_PTV(handles);
    guidata(hObject,handles)
end
function trackingsaveloc_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in trackingoutputdat.
function trackingoutputdat_Callback(hObject, eventdata, handles)
if str2double(handles.Njob)>0
    handles.data.Track.trackdat = num2str(get(hObject,'Value'));
    update_PTV(handles);
    guidata(hObject,handles);
end

% --- Executes on button press in trackingoutputmat.
function trackingoutputmat_Callback(hObject, eventdata, handles)
if str2double(handles.Njob)>0
    handles.data.Track.trackmat = num2str(get(hObject,'Value'));
    update_PTV(handles);
    guidata(hObject,handles);
end


% --- Executes on button press in loadtrackingsaveloc.
function loadtrackingsaveloc_Callback(hObject, eventdata, handles)
if str2double(handles.Njob)>0
    D = handles.data.Track.save_dir;
    handles.data.Track.save_dir = uigetdir(handles.data.Track.save_dir);
    if handles.data.Track.save_dir==0
        handles.data.Track.save_dir = D;
    end
    set(handles.trackingsaveloc,'string',handles.data.Track.save_dir);
    update_PTV(handles);
    guidata(hObject,handles)
end
%%

% --- Tracking Validation ---
% --- Tracking Validation Checkbox ---
function trackingvalcheckbox_Callback(hObject, eventdata, handles)
if str2double(handles.Njob)>0
    handles.data.Track.valprops.run=num2str(get(hObject,'Value'));
    update_PTV(handles);
    guidata(hObject,handles)
end

% --- Tracking Validation Coefficient ---
function trackingvalcoefficient_Callback(hObject, eventdata, handles)
if str2double(handles.Njob)>0
    handles.data.Track.valprops.valcoef=get(hObject,'String');
    update_PTV(handles);
    guidata(hObject,handles)
end
function trackingvalcoefficient_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Tracking Validation Radius ---
function trackingvalradius_Callback(hObject, eventdata, handles)
if str2double(handles.Njob)>0
    handles.data.Track.valprops.valrad=get(hObject,'String');
    update_PTV(handles);
    guidata(hObject,handles)
end
function trackingvalradius_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Tracking Validation U Threshold ---
function trackingvalUthresh_Callback(hObject, eventdata, handles)
if str2double(handles.Njob)>0
    handles.data.Track.valprops.MAD_U=get(hObject,'String');
    update_PTV(handles);
    guidata(hObject,handles)
end
function trackingvalUthresh_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Tracking Validation V Threshold ---
function trackingvalVthresh_Callback(hObject, eventdata, handles)
if str2double(handles.Njob)>0
    handles.data.Track.valprops.MAD_V=get(hObject,'String');
    update_PTV(handles);
    guidata(hObject,handles)
end
function trackingvalVthresh_CreateFcn(hObject, eventdata, handles)
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
    nums = strread(str, '%f', -1, 'delimiter', ',;:');
  
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
[workspacepath] = uigetdir(imageDir,'Location to Save Workspace');

% Exit save if the user selects 'cancel' in UIGET
if workspacepath ~= 0;
    
    %     Read the names of jobs in the present workspace
    jobs = get(handles.joblist, 'String');

    %     Save each job to the directory specified by workspaceDir
    for n = 1:length(jobs(:,1));
        endtest = 1;jobtemp=char(jobs(n,:));
        while endtest
            if strcmpi(jobtemp(end),' ')
                jobtemp=jobtemp(1:end-1);
            else
                endtest = 0;
            end
        end
        eval( ['Data = handles.' char(jobtemp) ';'] );
        save([workspacepath '/' char(jobtemp) '.mat'], 'Data');
        clear jobtemp
    end
else
    
    return
    
end

% --------------------------------------------------------------------
% --- Run PIV Checkbox ---
function runPIVcheckbox_Callback(hObject, eventdata, handles)
if str2double(handles.Njob)>0
    handles.data.runPIV=num2str(get(hObject,'Value'));
    set_PIVcontrols(handles);
    guidata(hObject,handles)
end

% --- Executes on selection change in colorchannel_popupMenu.
function colorchannel_popupMenu_Callback(hObject, eventdata, handles)
if str2double(handles.Njob)>0
    handles.data.channel = num2str(get(hObject,'Value'));
    guidata(hObject,handles)
end

% --- Executes during object creation, after setting all properties.
function colorchannel_popupMenu_CreateFcn(hObject, eventdata, handles)

function version_box_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function version_box_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function outputpassbasename_Callback(hObject, eventdata, handles)
if str2double(handles.Njob)>0
    N=str2double(handles.data.passes); % Number of passes
    passbase = get(handles.outputpassbasename,'String'); % Read text in outbasename textbox
    for i = 1:N
        eval(['handles.data.PIV' num2str(i) '.outbase=[''' passbase 'pass' num2str(i) '_''];']);
    end
    % Rename the ID Size and Track base names as well.
    handles.data.ID.savebase    = [passbase 'ID_'];
    handles.data.Size.savebase  = [passbase 'Size_'];
    handles.data.Track.savebase = [passbase 'Track_'];
    
    cpass = get(handles.passlist,'Value'); % This grabs the currently selected pass number
    set(handles.outputbasename,'String',eval(['handles.data.PIV' num2str(cpass) '.outbase']));
    handles.data.outputpassbase = passbase;
    guidata(hObject,handles)
end

% --- Executes during object creation, after setting all properties.
function outputpassbasename_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in renamejob_pushButton.
function renamejob_pushButton_Callback(hObject, eventdata, handles)
if str2double(handles.Njob) > 0 % If more than one job exists
    vn = 0; % Initialize some flag?
    selectedJob = str2double(handles.Cjob); % Determine which job is selected
    Jlist = char(get(handles.joblist,'String')); % read the names of the jobs

    % Ask user for the new name
    newJobName=char(inputdlg('Rename job as                  ','Rename Job',1,{strtrim(Jlist(str2double(handles.Cjob),:))})); % Input dialog
    if isempty(newJobName) % If no name was input, exit renamer.
        vn=-1; % Set valid name flag to -1 
    end

% Check for conditions that could cause errors
    while vn == 0 % If a valid name was input...
        if isfield(handles, newJobName) % If the job already exists...
            newJobName = char(inputdlg('Job already exists, rename?','NEW JOB',1,{strtrim(Jlist(str2double(handles.Cjob),:))})); % Inform user that job exists and ask for a different name
            if isempty(newJobName) % If an empty string was input for the new name, exit the renamer.
                vn=-1; % Set valid name flag to -1
            end
        else % If a valid string was input for the job name ...
            vn=1; % Set valid name flag to 1
        end
    end

if vn ~= -1 % If a valid name was input...
        JlistStruct = cellstr(Jlist); % Convert joblist to structure
        oldJobName = JlistStruct{selectedJob}; % Determine the old name of the job
        JlistStruct{selectedJob} = newJobName; % Change the job name in the GUI list box (but not in the data structure)
        Jlist = char(JlistStruct); % Convert job list structure to characters
        handles.(oldJobName) = handles.data; %pulls in all the changes that have been made to this job in the GUI
        [handles.(newJobName)] = handles.(oldJobName); % Copy the old job data to the new job (this line uses dynamic fields)
        handles = rmfield(handles, oldJobName); % Remove the old job from the data structure
        eval(['handles.' newJobName '.batchname = ''' newJobName ''';'])
end

    set(handles.joblist, 'String', Jlist, 'Value', str2double(handles.Cjob)); % Update joblist text box
    handles=update_data(handles); % Update handles
    guidata(hObject,handles)
    
end

function S = sizeMethodLookup(M)
% This function switch between numerials and strings for the sizing method.
% If a string is entered it will pass back the corresponding number and if
% a number is entered it will pass back the correct string.

size_str = {'GEO','IWC','TPG','FTG','CFPG','LSG','CLSG'};

if ischar(M)
    T = strfind(size_str,M);
    for i = 1:length(T)
        if ~isempty(T{i})
            S=i;
        end
    end
elseif isnumeric(M)
    S = size_str{M};
elseif isnan(M)
    S = 'GEO';
end


% --- Executes on button press in disp_exp_summary.
function disp_exp_summary_Callback(hObject, eventdata, handles)
if str2double(handles.Njob) > 0
    expsummary = write_expsummary(handles.data);
    fprintf(expsummary);
    fprintf('\n\n');
end



function input_velocity_Callback(hObject, eventdata, handles)
% hObject    handle to input_velocity (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if str2double(handles.Njob)>0  && get(handles.static_input_vel_button,'Value')==1
    handles.data.input_velocity = get(handles.input_velocity,'String');
    guidata(hObject,handles)
else
    set(handles.input_velocity,'String',handles.data.input_velocity);
end

% --- Executes during object creation, after setting all properties.
function input_velocity_CreateFcn(hObject, eventdata, handles)
% hObject    handle to input_velocity (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in browse_input_vel.
function browse_input_vel_Callback(hObject, eventdata, handles)
% hObject    handle to browse_input_vel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


if str2double(handles.Njob)>0 && get(handles.static_input_vel_button,'Value')==1
    D = handles.data.input_velocity;
    
    [a,b] = uigetfile([handles.data.imdirec '/*.*'], 'Select static mask file...');
    handles.data.input_velocity = [b a];
    if handles.data.input_velocity==0
        handles.data.input_velocity = D;
    end
    guidata(hObject,handles)
    set(handles.input_velocity,'string',handles.data.input_velocity);
end

% --- Executes on button press in no_input_vel_button.
function no_input_vel_button_Callback(hObject, eventdata, handles)
% hObject    handle to no_input_vel_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of no_input_vel_button

handles.data.input_vel_type='none';
set(hObject,'Value',1)
set(handles.static_input_vel_button,'Value',0)
set(handles.dynamic_input_vel_button,'Value',0)
if str2double(handles.Njob)>0
    Jlist=char(get(handles.joblist,'String'));
    eval(['handles.' Jlist(str2double(handles.Cjob),:) '=handles.data;']);
    handles=update_data(handles);
    guidata(hObject,handles)
end

% --- Executes on button press in static_input_vel_button.
function static_input_vel_button_Callback(hObject, eventdata, handles)
% hObject    handle to static_input_vel_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.data.input_vel_type='static';
set(hObject,'Value',1)
set(handles.no_input_vel_button,'Value',0)
set(handles.dynamic_input_vel_button,'Value',0)
if str2double(handles.Njob)>0
    Jlist=char(get(handles.joblist,'String'));
    eval(['handles.' Jlist(str2double(handles.Cjob),:) '=handles.data;']);
    handles=update_data(handles);
    guidata(hObject,handles)
end

% --- Executes on button press in dynamic_input_vel_button.
function dynamic_input_vel_button_Callback(hObject, eventdata, handles)
% hObject    handle to dynamic_input_vel_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.data.input_vel_type='dynamic';
set(hObject,'Value',1)
set(handles.no_input_vel_button,'Value',0)
set(handles.static_input_vel_button,'Value',0)
if str2double(handles.Njob)>0
    Jlist=char(get(handles.joblist,'String'));
    eval(['handles.' Jlist(str2double(handles.Cjob),:) '=handles.data;']);
    handles=update_data(handles);
    guidata(hObject,handles)
end



function input_veldirec_Callback(hObject, eventdata, handles)
% hObject    handle to input_veldirec (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if str2double(handles.Njob)>0 && get(handles.dynamic_input_vel_button,'Value')==1
    handles.data.input_veldirec = get(hObject,'String');
    guidata(hObject,handles)
else
    set(handles.input_veldirec,'String',handles.data.input_veldirec);
end


% --- Executes during object creation, after setting all properties.
function input_veldirec_CreateFcn(hObject, eventdata, handles)
% hObject    handle to input_veldirec (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in browse_input_veldirec.
function browse_input_veldirec_Callback(hObject, eventdata, handles)
% hObject    handle to browse_input_veldirec (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if str2double(handles.Njob)>0  && get(handles.dynamic_input_vel_button,'Value')==1
    D = handles.data.input_veldirec;

    handles.data.input_veldirec = uigetdir(D,'Directory for Input Velocities');
    if handles.data.input_veldirec == 0
        handles.data.input_veldirec = D;
    end
    set(handles.input_veldirec,'string',handles.data.input_veldirec);
    if ~isempty(dir(handles.data.input_veldirec))
        set(handles.input_veldirec,'Backgroundcolor',[1 1 1])
    end
    %load_imlist(handles);
    guidata(hObject,handles)
end



function input_velbase_Callback(hObject, eventdata, handles)
% hObject    handle to input_velbase (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if str2double(handles.Njob)>0 && get(handles.dynamic_input_vel_button,'Value')==1
    handles.data.input_velbase = get(hObject,'String');
    guidata(hObject,handles)
else
    set(handles.input_velbase,'String',handles.data.input_velbase);
end


% --- Executes during object creation, after setting all properties.
function input_velbase_CreateFcn(hObject, eventdata, handles)
% hObject    handle to input_velbase (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in cameranumber_popupMenu.
function cameranumber_popupMenu_Callback(hObject, eventdata, handles)
% hObject    handle to cameranumber_popupMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if str2double(handles.Njob)>0
    handles.data.numcams = num2str(get(hObject,'Value'));
    guidata(hObject,handles)
end


% --- Executes during object creation, after setting all properties.
function cameranumber_popupMenu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to cameranumber_popupMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in loadcam2firstimage.
function loadcam2firstimage_Callback(hObject, eventdata, handles)
if str2double(handles.Njob)>0
    D = handles.data.imdirec2;
    B = handles.data.imbase2;
%     Z = handles.data.imzeros;
%     E = handles.data.imext;
%     S = handles.data.imfstart;
%     F = handles.data.imfend;
%     C = handles.data.imcstep;
    [imname,handles.data.imdirec2] = uigetfile('*.*','Select an image file',handles.data.imdirec);
    if handles.data.imdirec2==0
        handles.data.imdirec2 = D;
        handles.data.imbase2 = B;
%         handles.data.imzeros = Z;
%         handles.data.imext = E;
%         handles.data.imfstart = S;
%         handles.data.imfend = F;
    else
        i=strfind(imname,'.');
        handles.data.imext=imname((i(end)+1):end);
        % This finds everything that is not a number
        j = regexp(imname(1:i(end)-1),'\D');
        % Now looking at the first thing that is not a number excluding the
        % extension to determine the number of zeros.
        zeros=i(end)-1-j(end);
        handles.data.imbase2=imname(1:(i(end)-1-zeros));
        handles.data.imzeros=num2str(zeros);
%         fstart=str2double(imname((i(end)-zeros):(i(end)-1)));
%         handles.data.imfstart=num2str(fstart);
%         dirinfo = dir([handles.data.imdirec handles.data.imbase '*.' handles.data.imext]);
%         handles.data.imfend=num2str(str2double(dirinfo(end).name(i(end)-zeros:i(end)-1))-str2double(C));
    end
    set(handles.imagedirectory2,'string',handles.data.imdirec2);
    set(handles.imagebasename2,'string',handles.data.imbase2);
%     set(handles.imagezeros,'string',handles.data.imzeros);
%     set(handles.imageextension,'string',handles.data.imext);
%     set(handles.imageframestart,'string',handles.data.imfstart);
%     set(handles.imageframeend,'string',handles.data.imfend);
%     if strcmp(handles.data.masktype,'dynamic')
%         load_masklist(handles)
%     end
    load_imlist(handles);
    guidata(hObject,handles)
end

function imagedirectory2_Callback(hObject, eventdata, handles)
if str2double(handles.Njob)>0
    D = handles.data.imdirec2;
    handles.data.imdirec2 = uigetdir(handles.data.imdirec2);
    if handles.data.imdirec2==0
        handles.data.imdirec2 = D;
    end
    set(handles.imagedirectory2,'string',handles.data.imdirec2);
%     if strcmp(handles.data.masktype,'dynamic')
%         load_masklist(handles)
%     end
    load_imlist(handles);
    guidata(hObject,handles)
end

% --- Executes during object creation, after setting all properties.
function imagedirectory2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to imagedirectory2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in loadimagedirectory2button.
function loadimagedirectory2button_Callback(hObject, eventdata, handles)
% hObject    handle to loadimagedirectory2button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if str2double(handles.Njob)>0
    D = handles.data.imdirec2;
    handles.data.imdirec2 = uigetdir(handles.data.imdirec2);
    if handles.data.imdirec2==0
        handles.data.imdirec2 = D;
    end
    set(handles.imagedirectory2,'string',handles.data.imdirec2);
%     if strcmp(handles.data.masktype,'dynamic')
%         load_masklist(handles)
%     end
%     load_imlist(handles);
    guidata(hObject,handles)
end


function imagebasename2_Callback(hObject, eventdata, handles)
if str2double(handles.Njob)>0
    handles.data.imbase2 = get(hObject,'String');
%     if strcmp(handles.data.masktype,'dynamic')
%         load_masklist(handles)
%     end
    load_imlist(handles);
    guidata(hObject,handles)
end

% --- Executes during object creation, after setting all properties.
function imagebasename2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to imagebasename2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end





% --- Handles default GUIDE GUI creation and callback dispatch
function varargout = gui_mainfcn(gui_State, varargin)

gui_StateFields =  {'gui_Name'
    'gui_Singleton'
    'gui_OpeningFcn'
    'gui_OutputFcn'
    'gui_LayoutFcn'
    'gui_Callback'};
gui_Mfile = '';
for i=1:length(gui_StateFields)
    if ~isfield(gui_State, gui_StateFields{i})
        error(message('MATLAB:guide:StateFieldNotFound', gui_StateFields{ i }, gui_Mfile));
    elseif isequal(gui_StateFields{i}, 'gui_Name')
        gui_Mfile = [gui_State.(gui_StateFields{i}), '.m'];
    end
end

numargin = length(varargin);

if numargin == 0
    % PRANA
    % create the GUI only if we are not in the process of loading it
    % already
    gui_Create = true;
elseif local_isInvokeActiveXCallback(gui_State, varargin{:})
    % PRANA(ACTIVEX,...)
    vin{1} = gui_State.gui_Name;
    vin{2} = [get(varargin{1}.Peer, 'Tag'), '_', varargin{end}];
    vin{3} = varargin{1};
    vin{4} = varargin{end-1};
    vin{5} = guidata(varargin{1}.Peer);
    feval(vin{:});
    return;
elseif local_isInvokeHGCallback(gui_State, varargin{:})
    % PRANA('CALLBACK',hObject,eventData,handles,...)
    gui_Create = false;
else
    % PRANA(...)
    % create the GUI and hand varargin to the openingfcn
    gui_Create = true;
end

if ~gui_Create
    % In design time, we need to mark all components possibly created in
    % the coming callback evaluation as non-serializable. This way, they
    % will not be brought into GUIDE and not be saved in the figure file
    % when running/saving the GUI from GUIDE.
    designEval = false;
    if (numargin>1 && ishghandle(varargin{2}))
        fig = varargin{2};
        while ~isempty(fig) && ~ishghandle(fig,'figure')
            fig = get(fig,'parent');
        end
        
        designEval = isappdata(0,'CreatingGUIDEFigure') || (isscalar(fig)&&isprop(fig,'GUIDEFigure'));
    end
        
    if designEval
        beforeChildren = findall(fig);
    end
    
    % evaluate the callback now
    varargin{1} = gui_State.gui_Callback;
    if nargout
        [varargout{1:nargout}] = feval(varargin{:});
    else       
        feval(varargin{:});
    end
    
    % Set serializable of objects created in the above callback to off in
    % design time. Need to check whether figure handle is still valid in
    % case the figure is deleted during the callback dispatching.
    if designEval && ishghandle(fig)
        set(setdiff(findall(fig),beforeChildren), 'Serializable','off');
    end
else
    if gui_State.gui_Singleton
        gui_SingletonOpt = 'reuse';
    else
        gui_SingletonOpt = 'new';
    end

    % Check user passing 'visible' P/V pair first so that its value can be
    % used by oepnfig to prevent flickering
    gui_Visible = 'auto';
    gui_VisibleInput = '';
    for index=1:2:length(varargin)
        if length(varargin) == index || ~ischar(varargin{index})
            break;
        end

        % Recognize 'visible' P/V pair
        len1 = min(length('visible'),length(varargin{index}));
        len2 = min(length('off'),length(varargin{index+1}));
        if ischar(varargin{index+1}) && strncmpi(varargin{index},'visible',len1) && len2 > 1
            if strncmpi(varargin{index+1},'off',len2)
                gui_Visible = 'invisible';
                gui_VisibleInput = 'off';
            elseif strncmpi(varargin{index+1},'on',len2)
                gui_Visible = 'visible';
                gui_VisibleInput = 'on';
            end
        end
    end
    
    % Open fig file with stored settings.  Note: This executes all component
    % specific CreateFunctions with an empty HANDLES structure.

    
    % Do feval on layout code in m-file if it exists
    gui_Exported = ~isempty(gui_State.gui_LayoutFcn);
    % this application data is used to indicate the running mode of a GUIDE
    % GUI to distinguish it from the design mode of the GUI in GUIDE. it is
    % only used by actxproxy at this time.   
    setappdata(0,genvarname(['OpenGuiWhenRunning_', gui_State.gui_Name]),1);
    if gui_Exported
        gui_hFigure = feval(gui_State.gui_LayoutFcn, gui_SingletonOpt);

        % make figure invisible here so that the visibility of figure is
        % consistent in OpeningFcn in the exported GUI case
        if isempty(gui_VisibleInput)
            gui_VisibleInput = get(gui_hFigure,'Visible');
        end
        set(gui_hFigure,'Visible','off')

        % openfig (called by local_openfig below) does this for guis without
        % the LayoutFcn. Be sure to do it here so guis show up on screen.
        movegui(gui_hFigure,'onscreen');
    else
        gui_hFigure = local_openfig(gui_State.gui_Name, gui_SingletonOpt, gui_Visible);
        % If the figure has InGUIInitialization it was not completely created
        % on the last pass.  Delete this handle and try again.
        if isappdata(gui_hFigure, 'InGUIInitialization')
            delete(gui_hFigure);
            gui_hFigure = local_openfig(gui_State.gui_Name, gui_SingletonOpt, gui_Visible);
        end
    end
    if isappdata(0, genvarname(['OpenGuiWhenRunning_', gui_State.gui_Name]))
        rmappdata(0,genvarname(['OpenGuiWhenRunning_', gui_State.gui_Name]));
    end

    % Set flag to indicate starting GUI initialization
    setappdata(gui_hFigure,'InGUIInitialization',1);

    % Fetch GUIDE Application options
    gui_Options = getappdata(gui_hFigure,'GUIDEOptions');
    % Singleton setting in the GUI M-file takes priority if different
    gui_Options.singleton = gui_State.gui_Singleton;

    if ~isappdata(gui_hFigure,'GUIOnScreen')
        % Adjust background color
        if gui_Options.syscolorfig
            set(gui_hFigure,'Color', get(0,'DefaultUicontrolBackgroundColor'));
        end

        % Generate HANDLES structure and store with GUIDATA. If there is
        % user set GUI data already, keep that also.
        data = guidata(gui_hFigure);
        handles = guihandles(gui_hFigure);
        if ~isempty(handles)
            if isempty(data)
                data = handles;
            else
                names = fieldnames(handles);
                for k=1:length(names)
                    data.(char(names(k)))=handles.(char(names(k)));
                end
            end
        end
        guidata(gui_hFigure, data);
    end

    % Apply input P/V pairs other than 'visible'
    for index=1:2:length(varargin)
        if length(varargin) == index || ~ischar(varargin{index})
            break;
        end

        len1 = min(length('visible'),length(varargin{index}));
        if ~strncmpi(varargin{index},'visible',len1)
            try set(gui_hFigure, varargin{index}, varargin{index+1}), catch break, end
        end
    end

    % If handle visibility is set to 'callback', turn it on until finished
    % with OpeningFcn
    gui_HandleVisibility = get(gui_hFigure,'HandleVisibility');
    if strcmp(gui_HandleVisibility, 'callback')
        set(gui_hFigure,'HandleVisibility', 'on');
    end

    feval(gui_State.gui_OpeningFcn, gui_hFigure, [], guidata(gui_hFigure), varargin{:});

    if isscalar(gui_hFigure) && ishghandle(gui_hFigure)
        % Handle the default callbacks of predefined toolbar tools in this
        % GUI, if any
        guidemfile('restoreToolbarToolPredefinedCallback',gui_hFigure); 
        
        % Update handle visibility
        set(gui_hFigure,'HandleVisibility', gui_HandleVisibility);

        % Call openfig again to pick up the saved visibility or apply the
        % one passed in from the P/V pairs
        if ~gui_Exported
            gui_hFigure = local_openfig(gui_State.gui_Name, 'reuse',gui_Visible);
        elseif ~isempty(gui_VisibleInput)
            set(gui_hFigure,'Visible',gui_VisibleInput);
        end
        if strcmpi(get(gui_hFigure, 'Visible'), 'on')
            figure(gui_hFigure);
            
            if gui_Options.singleton
                setappdata(gui_hFigure,'GUIOnScreen', 1);
            end
        end

        % Done with GUI initialization
        if isappdata(gui_hFigure,'InGUIInitialization')
            rmappdata(gui_hFigure,'InGUIInitialization');
        end

        % If handle visibility is set to 'callback', turn it on until
        % finished with OutputFcn
        gui_HandleVisibility = get(gui_hFigure,'HandleVisibility');
        if strcmp(gui_HandleVisibility, 'callback')
            set(gui_hFigure,'HandleVisibility', 'on');
        end
        gui_Handles = guidata(gui_hFigure);
    else
        gui_Handles = [];
    end

    if nargout
        [varargout{1:nargout}] = feval(gui_State.gui_OutputFcn, gui_hFigure, [], gui_Handles);
    else
        feval(gui_State.gui_OutputFcn, gui_hFigure, [], gui_Handles);
    end

    if isscalar(gui_hFigure) && ishghandle(gui_hFigure)
        set(gui_hFigure,'HandleVisibility', gui_HandleVisibility);
    end
end

function gui_hFigure = local_openfig(name, singleton, visible)

% openfig with three arguments was new from R13. Try to call that first, if
% failed, try the old openfig.
if nargin('openfig') == 2
    % OPENFIG did not accept 3rd input argument until R13,
    % toggle default figure visible to prevent the figure
    % from showing up too soon.
    gui_OldDefaultVisible = get(0,'defaultFigureVisible');
    set(0,'defaultFigureVisible','off');
    gui_hFigure = matlab.hg.internal.openfigLegacy(name, singleton);
    set(0,'defaultFigureVisible',gui_OldDefaultVisible);
else
    % Call version of openfig that accepts 'auto' option"
    gui_hFigure = matlab.hg.internal.openfigLegacy(name, singleton, visible);  
%     %workaround for CreateFcn not called to create ActiveX
%     if feature('HGUsingMATLABClasses')
%         peers=findobj(findall(allchild(gui_hFigure)),'type','uicontrol','style','text');    
%         for i=1:length(peers)
%             if isappdata(peers(i),'Control')
%                 actxproxy(peers(i));
%             end            
%         end
%     end
end

function result = local_isInvokeActiveXCallback(gui_State, varargin)

try
    result = ispc && iscom(varargin{1}) ...
             && isequal(varargin{1},gcbo);
catch
    result = false;
end

function result = local_isInvokeHGCallback(gui_State, varargin)

try
    fhandle = functions(gui_State.gui_Callback);
    result = ~isempty(findstr(gui_State.gui_Name,fhandle.file)) || ...
             (ischar(varargin{1}) ...
             && isequal(ishghandle(varargin{2}), 1) ...
             && (~isempty(strfind(varargin{1},[get(varargin{2}, 'Tag'), '_'])) || ...
                ~isempty(strfind(varargin{1}, '_CreateFcn'))) );
catch
    result = false;
end


