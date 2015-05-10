function varargout = prana_SPIV(varargin)
% PRANA_SPIV M-file for prana_SPIV.fig
%      PRANA_SPIV, by itself, creates a new PRANA_SPIV or raises the existing
%      singleton*.
%
%      H = PRANA_SPIV returns the handle to a new PRANA_SPIV or the handle to
%      the existing singleton*.
%
%      PRANA_SPIV('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in PRANA_SPIV.M with the given input arguments.
%
%      PRANA_SPIV('Property','Value',...) creates a new PRANA_SPIV or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before stereoreconstruct_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to prana_SPIV_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help prana_SPIV

% Last Modified by GUIDE v2.5 09-Apr-2015 19:14:58

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @prana_SPIV_OpeningFcn, ...
    'gui_OutputFcn',  @prana_SPIV_OutputFcn, ...
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

end

function prana_SPIV_OpeningFcn(hObject, eventdata, handles, varargin)
% Choose default command line output for prana_SPIV
handles.output = hObject;

% NOTE: This section loads in/pre-allocates all settings and input handles
% for the algorithm.

guidata(hObject);
guiprops = guidata(hObject);

guiprops.markerdiam    = str2double(get(handles.markerdiam,'String'));
guiprops.xplatelevelshift =str2double(get(handles.xplatelevelshift,'String'));
guiprops.yplatelevelshift =str2double(get(handles.yplatelevelshift,'String'));
guiprops.targetsidesbox = get(handles.targetsidesbox,'Value'); % cam1 orientation
guiprops.targetsidesbox2 = get(handles.targetsidesbox2,'Value'); % cam2 orientation
guiprops.vertmarkerspacing  = str2double(get(handles.vertmarkerspacing,'String'));
guiprops.hormarkerspacing   = str2double(get(handles.hormarkerspacing,'String'));
guiprops.targetcolor = get(handles.targetcolor,'Value');
guiprops.zplanenumber    = str2double(get(handles.zplanenumber,'String'));
guiprops.zstart = str2double(get(handles.zstart,'String'));
guiprops.planespacing = str2num(get(handles.planespacing,'String')); % takes an array of plane_spacing if that is the input
guiprops.channeldepth  = str2num(get(handles.channeldepth,'String'));% takes an array of channel depth like [0 1]
guiprops.levelnumber = str2double(get(handles.levelnumber,'String'));
guiprops.platethickness= str2double(get(handles.platethickness,'String'));
guiprops.modeltypemenu    = get(handles.modeltypemenu,'Value');   %1 is for cubic xy and linear z, 2 is for cubic xy and quadratic z, 3 is for DLT
% 

contents= cellstr(get(handles.multilevelpopup,'String'));
guiprops.multilevelpopup=contents{get(handles.multilevelpopup,'Value')};


if guiprops.targetcolor==1
    guiprops.caljob.grid_points_bright='true';
elseif guiprops.targcolor==2
    guiprops.caljob.grid_points_bright='false';
end
guiprops.caljob.targetsidecam1=guiprops.targetsidesbox;
guiprops.caljob.targetsidecam2=guiprops.targetsidesbox2;
guiprops.caljob.grid_point_diameter=guiprops.markerdiam;

guiprops.caljob.x_grid_spacing=guiprops.hormarkerspacing;
guiprops.caljob.y_grid_spacing=guiprops.vertmarkerspacing;

guiprops.caljob.zplanenumber=guiprops.zplanenumber;
guiprops.caljob.z_grid_start=guiprops.zstart;
guiprops.caljob.plane_spacing=guiprops.planespacing;
guiprops.caljob.plane_numbers_list=1:1:guiprops.caljob.zplanenumber;

%guiprops.caljob.multiple_level_plate ='true';
if strcmp(guiprops.multilevelpopup,'True')
    guiprops.caljob.multiple_level_plate ='true';
    set(handles.multilevelpanel,'Visible','on');
    % set(handles.levelnumber,'String','3');
     guiprops.levelnumber=2;
    guiprops.caljob.plate_level_number=guiprops.levelnumber;
elseif strcmp(guiprops.multilevelpopup,'False')
    guiprops.caljob.multiple_level_plate ='false';
    set(handles.multilevelpanel,'Visible','off');
    guiprops.levelnumber = 1;
    guiprops.caljob.plate_level_number=guiprops.levelnumber;
end

guiprops.caljob.x_plate_level_shift=[0,guiprops.xplatelevelshift];
guiprops.caljob.y_plate_level_shift=[0,guiprops.yplatelevelshift];
guiprops.caljob.plate_level_spacing = [0,-guiprops.channeldepth];
guiprops.caljob.plate_level_number=guiprops.levelnumber;
guiprops.caljob.front_to_back_plate_thickness =guiprops.platethickness;

guiprops.caljob.y_pixel_number=256;
guiprops.caljob.x_pixel_number=256;
guiprops.caljob.length_unit='mm';
guiprops.caljob.image_extension='.tif';
guiprops.caljob.camnumber=[1,2];
guiprops.caljob.calimagelist={};

guiprops.caljob.JOBFILE={};
guiprops.caljob.calibration_data={};
guiprops.caljob.calibration_plane_data={};

guiprops.caljob.allx1data=[];
guiprops.caljob.allx2data=[];
guiprops.caljob.allX1data=[];
guiprops.caljob.allX2data=[];

guiprops.caljob.modeltype    = get(handles.modeltypemenu,'Value'); 
guiprops.caljob.optionsls=optimset('MaxIter',30000,'MaxFunEvals',30000,'TolX',1e-11,'TolFun',1e-7,...
    'LargeScale','off','Display','off');
guiprops.caljob.aXcam1             = [];   % polynomial coefficients for the least squares fit, mapping matrix for image X (camera 1) onto object x,y,z
guiprops.caljob.aYcam1             = [];
guiprops.caljob.aXcam2             = [];
guiprops.caljob.aYcam2             = [];
guiprops.caljob.a_cam1             = [];   % coeff. of fund. matrix of camera pinhole model for cam. 1
guiprops.caljob.a_cam2             = [];
guiprops.caljob.convergemessage = {''};

guiprops.selfcaljob=[];
guiprops.planarjob=[];
guiprops.rectype='';
guiprops.scaling=[];
guiprops.propertiessavename='StereoJob.mat';
guiprops.propertiessavepath='';
guiprops.stereorecdirlist=[];

guidata(hObject,guiprops);
%keyboard;
end

function varargout = prana_SPIV_OutputFcn(hObject, eventdata, handles)

end

function loadpoints_Callback(hObject, eventdata, handles)

% NOTE: This section allows for the loading of all critical data saved in
% .MAT files
 %guiprops = guidata(hObject);
[testptsname,testptspath] = uigetfile('*.mat','Select .mat file');

if ~isequal(testptsname,0)                                                  % if the filename is valid
    load([testptspath testptsname]);
    %keyboard;
%     handles.output = hObject;
    guiprops = guidata(hObject);
     %set(handles.markerdiam,'String',num2str(guiprops.markerdiam));
    
    guiprops.propertiessavename = testptsname;
    guiprops.propertiessavepath = testptspath;
    
    namefields = fieldnames(datasave);                                      % load saved variables into memory
    
    for j=1:length(namefields)   
        %set(handles.outputdirectorybox,'String',guiprops.outputdirectory);
        guiprops.(namefields{j}) = datasave.(namefields{j});
    end
    %keyboard;
    %%%SOMEHOW HANDLES ARE NOT BEING RESET AFTER LOADING A JOB NOT SURE WHY
% %      set(handles.markerdiam,'String',num2str(guiprops.markerdiam));
% %      set(handles.xplatelevelshift,'String',num2str(guiprops.xplatelevelshift));
% %      set(handles.yplatelevelshift,'String',num2str(guiprops.yplatelevelshift));
% %      set(handles.targetsidesbox,'Value',guiprops.targetsidesbox);
% %      set(handles.hormarkerspacing,'String',num2str(guiprops.hormarkerspacing));
% %      set(handles.vertmarkerspacing,'String',num2str(guiprops.vertmarkerspacing));
% %      set(handles.targetcolor,'Value',guiprops.targetcolor);
% %      set(handles.multilevelpopup,'String',guiprops.multilevelpopup);
% %      
% %      set(handles.zplanenumber,'String',num2str(guiprops.zplanenumber));
% %      set(handles.zstart,'String',num2str(guiprops.zstart));
% %      set(handles.planespacing,'String',num2str(guiprops.planespacing));
% %      set(handles.channeldepth,'String',num2str(guiprops.channeldepth));
% %      set(handles.levelnumber,'String',num2str(guiprops.levelnumber));
% %      set(handles.platethickness,'String',num2str(guiprops.platethickness));
% %    
% %     
% %     if strcmp(guiprops.multilevelpopup,'True')
% %         %guiprops.caljob.multiple_level_plate ='true';
% %         set(handles.multilevelpanel,'Visible','on');
% %         % set(handles.levelnumber,'String','3');
% % %         guiprops.levelnumber=2;
% % %         guiprops.caljob.plate_level_number=guiprops.levelnumber;
% %     elseif strcmp(guiprops.multilevelpopup,'False')
% %         %guiprops.caljob.multiple_level_plate ='false';
% %         set(handles.multilevelpanel,'Visible','off');
% % %         guiprops.levelnumber = 1;
% % %         guiprops.caljob.plate_level_number=guiprops.levelnumber;
% %     end
% %     set(handles.modeltypemenu,      'Value' ,guiprops.modeltypemenu);
% %     set(handles.convergebox,        'String',guiprops.convergemessage);
    
    
%     set(handles.markerdiam,         'String',num2str(guiprops.markerd));
%     set(handles.vertmarkerspacing,  'String',num2str(guiprops.verspace));
%     set(handles.hormarkerspacing,   'String',num2str(guiprops.horspace));
%     set(handles.targetcolor,        'Value' ,guiprops.targetcolor);
%     set(handles.channeldepth,       'String',num2str(guiprops.chandepth));
   
%    
    
%keyboard;
    guidata(hObject,guiprops);
end
end

function savepts_Callback(hObject, eventdata, handles)

guiprops = guidata(hObject);

%keyboard;

datasave.caljob             =guiprops.caljob;
datasave.selfcaljob         =guiprops.selfcaljob;
datasave.planarjob          = guiprops.planarjob;
datasave.modeltypemenu      = guiprops.modeltypemenu;
datasave.convergemessage    = guiprops.convergemessage;
datasave.rectype            = guiprops.rectype;
datasave.scaling            = guiprops.scaling;

datasave.markerdiam =guiprops.markerdiam ;
datasave.xplatelevelshift=guiprops.xplatelevelshift;
datasave.yplatelevelshift=guiprops.yplatelevelshift;
datasave.targetsidesboxcam1=guiprops.targetsidesbox;
datasave.targetsidesboxcam2=guiprops.targetsidesbox2;
datasave.vertmarkerspacing=guiprops.vertmarkerspacing;
datasave.hormarkerspacing=guiprops.hormarkerspacing;
datasave.targetcolor=guiprops.targetcolor;
datasave.zplanenumber=guiprops.zplanenumber;
datasave.zstart=guiprops.zstart;
datasave.planespacing=guiprops.planespacing;
datasave.channeldepth=guiprops.channeldepth;
datasave.levelnumber=guiprops.levelnumber;
datasave.platethickness=guiprops.platethickness;
datasave.modeltypemenu=guiprops.modeltypemenu;   %1 is for cubic xy and linear z, 2 is for cubic xy and quadratic z, 3 is for DLT
% 
datasave.multilevelpopup=guiprops.multilevelpopup;
datasave.stereorecdirlist=guiprops.stereorecdirlist;

datasave.propertiessavename = guiprops.propertiessavename;
datasave.propertiessavepath = guiprops.propertiessavepath;

[file1,path1] = uiputfile(guiprops.propertiessavename,'Save file');
if ~isequal(file1,0)
    save([path1 file1],'datasave');
    guiprops.propertiessavename=file1;
    guiprops.propertiessavepath=path1;
end
guidata(hObject,guiprops);
end

function quitbutton_Callback(hObject, eventdata, handles)
close all;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% CALIBRATION SECTION OF CODE                                      %%%%%
%%%% All handles contained within this section pertain to Calibration %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function markerdiam_Callback(hObject, eventdata, handles)

guiprops            = guidata(hObject);
guiprops.markerdiam    = str2double(get(hObject,'String'));
guiprops.caljob.grid_point_diameter=guiprops.markerdiam;
guidata(hObject,guiprops);

end

function markerdiam_CreateFcn(hObject, eventdata, handles)
end

function targetsidesbox_Callback(hObject, eventdata, handles)

guiprops            = guidata(hObject);
guiprops.targetsidesbox = get(hObject,'Value');
guiprops.caljob.targetsidecam1=guiprops.targetsidesbox;
guidata(hObject,guiprops);

end

function targetsidesbox_CreateFcn(hObject, eventdata, handles)
end

function targetsidesbox2_Callback(hObject, eventdata, handles)

guiprops            = guidata(hObject);
guiprops.targetsidesbox2 = get(hObject,'Value');
guiprops.caljob.targetsidecam2=guiprops.targetsidesbox2;
guidata(hObject,guiprops);

end

% --- Executes during object creation, after setting all properties.
function targetsidesbox2_CreateFcn(hObject, eventdata, handles)
end

function vertmarkerspacing_Callback(hObject, eventdata, handles)

guiprops            = guidata(hObject);
guiprops.vertmarkerspacing  = str2double(get(hObject,'String'));
guiprops.caljob.y_grid_spacing=guiprops.vertmarkerspacing;
guidata(hObject,guiprops);

end

function vertmarkerspacing_CreateFcn(hObject, eventdata, handles)
end

function hormarkerspacing_Callback(hObject, eventdata, handles)

guiprops            = guidata(hObject);
guiprops.hormarkerspacing   = str2double(get(hObject,'String'));
guiprops.caljob.x_grid_spacing=guiprops.hormarkerspacing;
guidata(hObject,guiprops);

end

function hormarkerspacing_CreateFcn(hObject, eventdata, handles)
end



function targetcolor_Callback(hObject, eventdata, handles)

guiprops           = guidata(hObject);
guiprops.targetcolor = get(hObject,'Value');
if guiprops.targetcolor==1
    guiprops.caljob.grid_points_bright='true';
elseif guiprops.targetcolor==2
    guiprops.caljob.grid_points_bright='false';
end

%keyboard;
guidata(hObject,guiprops);

end

function targetcolor_CreateFcn(hObject, eventdata, handles)

guiprops            = guidata(hObject);
guiprops.targetcolor  = get(hObject,'Value');
guidata(hObject,guiprops);

end


function xplatelevelshift_Callback(hObject, eventdata, handles)
guiprops            = guidata(hObject);
guiprops.xplatelevelshift = str2double(get(hObject,'String'));
guiprops.caljob.x_plate_level_shift=[0,guiprops.xplatelevelshift];
guidata(hObject,guiprops);
end

% --- Executes during object creation, after setting all properties.
function xplatelevelshift_CreateFcn(hObject, eventdata, handles)
% hObject    handle to xplatelevelshift (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end


function yplatelevelshift_Callback(hObject, eventdata, handles)
guiprops            = guidata(hObject);
guiprops.yplatelevelshift = str2double(get(hObject,'String'));
guiprops.caljob.y_plate_level_shift=[0,guiprops.yplatelevelshift];
guidata(hObject,guiprops);
end

% --- Executes during object creation, after setting all properties.
function yplatelevelshift_CreateFcn(hObject, eventdata, handles)
% hObject    handle to yplatelevelshift (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end

% --- Executes on selection change in multilevelpopup.
function multilevelpopup_Callback(hObject, eventdata, handles)
guiprops            = guidata(hObject);
contents= cellstr(get(hObject,'String'));
guiprops.multilevelpopup=contents{get(hObject,'Value')};
if strcmp(guiprops.multilevelpopup,'True')
    guiprops.caljob.multiple_level_plate ='true';
    set(handles.multilevelpanel,'Visible','on');
    % set(handles.levelnumber,'String','3');
     guiprops.levelnumber=2;
    guiprops.caljob.plate_level_number=guiprops.levelnumber;
elseif strcmp(guiprops.multilevelpopup,'False')
    guiprops.caljob.multiple_level_plate ='false';
    set(handles.multilevelpanel,'Visible','off');
    guiprops.levelnumber = 1;
    guiprops.caljob.plate_level_number=guiprops.levelnumber;
end
guidata(hObject,guiprops);

% Hints: contents = cellstr(get(hObject,'String')) returns multilevelpopup contents as cell array
%        contents{get(hObject,'Value')} returns selected item from multilevelpopup
end

% --- Executes during object creation, after setting all properties.
function multilevelpopup_CreateFcn(hObject, eventdata, handles)
% hObject    handle to multilevelpopup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end


function zplanenumber_Callback(hObject, eventdata, handles)
guiprops            = guidata(hObject);
guiprops.zplanenumber    = str2double(get(hObject,'String'));
guiprops.caljob.plane_numbers_list=1:1:guiprops.zplanenumber;
guiprops.caljob.zplanenumber=guiprops.zplanenumber;
guidata(hObject,guiprops);
end

% --- Executes during object creation, after setting all properties.
function zplanenumber_CreateFcn(hObject, eventdata, handles)
% hObject    handle to zplanenumber (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end


function zstart_Callback(hObject, eventdata, handles)
guiprops            = guidata(hObject);
guiprops.zstart = str2double(get(hObject,'String'));
guiprops.caljob.z_grid_start=guiprops.zstart;
guidata(hObject,guiprops);
end

% --- Executes during object creation, after setting all properties.
function zstart_CreateFcn(hObject, eventdata, handles)
% hObject    handle to zstart (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end


function planespacing_Callback(hObject, eventdata, handles)
guiprops            = guidata(hObject);
% STR2NUM is used because the plane spacing can vary in which case it will
% be a vector.Str2double only operates on scalars.
guiprops.planespacing = str2num(get(hObject,'String'));
guiprops.caljob.plane_spacing=guiprops.planespacing;
guidata(hObject,guiprops);
end

% --- Executes during object creation, after setting all properties.
function planespacing_CreateFcn(hObject, eventdata, handles)
% hObject    handle to planespacing (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end
function channeldepth_Callback(hObject, eventdata, handles)

guiprops            = guidata(hObject);
guiprops.channeldepth  = str2num(get(hObject,'String'));
guiprops.caljob.plate_level_spacing =guiprops.channeldepth;
guidata(hObject,guiprops);
end

function channeldepth_CreateFcn(hObject, eventdata, handles)
end

function levelnumber_Callback(hObject, eventdata, handles)
guiprops            = guidata(hObject);
guiprops.levelnumber = str2double(get(hObject,'String'));
guiprops.caljob.plate_level_number=guiprops.levelnumber;
guidata(hObject,guiprops);
end

% --- Executes during object creation, after setting all properties.
function levelnumber_CreateFcn(hObject, eventdata, handles)
% hObject    handle to levelnumber (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end


function platethickness_Callback(hObject, eventdata, handles)
guiprops            = guidata(hObject);
guiprops.platethickness= str2double(get(hObject,'String'));
guiprops.caljob.front_to_back_plate_thickness =guiprops.platethickness;
guidata(hObject,guiprops);
end

% --- Executes during object creation, after setting all properties.
function platethickness_CreateFcn(hObject, eventdata, handles)
% hObject    handle to platethickness (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end




function getcal1_Callback(hObject, eventdata, handles)
guiprops            = guidata(hObject);
[filename1,pathname1]   = uigetfile('*.*','Select All The Calibration Image Planes(which you want)','Multiselect','on');
guiprops.selectcam1file = pathname1;

    cn1= input('Please Enter First Camera Number:');
    if ~isempty(cn1)
        guiprops.caljob.camnumber(1)=cn1;
    else
        guiprops.caljob.camnumber(1)=1;
    end
% for multiple planes selected filename1 becomes a cell and so it checks
% for that to see if multiple calibration plane or a single plane is
% selected.
if iscell(filename1) 
    for t=1:length(filename1)
        guiprops.caljob.calimagelist{1,t}=[pathname1,filename1{t}];
    end
else
    guiprops.caljob.calimagelist{1,1}=[pathname1,filename1];
end

guidata(hObject,guiprops);
end

function getcal2_Callback(hObject, eventdata, handles)
guiprops            = guidata(hObject);
[filename2,pathname2]   = uigetfile('*.*','Select All The Calibration Image Planes(which you want)',guiprops.selectcam1file ,'Multiselect','on');
guiprops.selectcam2file = pathname2;

    cn2= input('Please Enter Second Camera Number:');
    if ~isempty(cn2)
        guiprops.caljob.camnumber(2)=cn2;
    else
        guiprops.caljob.camnumber(2)=2;
    end
% for multiple planes selected filename1 becomes a cell and so it checks
% for that to see if multiple calibration plane or a single plane is
% selected.

if iscell(filename2)
    for t=1:length(filename2)
        guiprops.caljob.calimagelist{2,t}=[pathname2,filename2{t}];
    end
else
    guiprops.caljob.calimagelist{2,1}=[pathname2,filename2];
end
guidata(hObject,guiprops);
end

function controlptselect_Callback(hObject, eventdata, handles)

guiprops    = guidata(hObject);
%guiprops.caljob.calimagelist{1,1}
%calibrationjob=guiprops.caljob.;
A=imread(guiprops.caljob.calimagelist{1,1});
[guiprops.caljob.y_pixel_number,guiprops.caljob.x_pixel_number]=size(A);
guiprops.caljob.length_unit='mm';
ext= regexp(guiprops.caljob.calimagelist{1,1},'\d','split');
guiprops.caljob.image_extension=ext{end};
%keyboard;
JOBFILE.Camera_Numbers_List					= guiprops.caljob.camnumber;% [1,2];%guiprops.caljob.camera_numbers_list;
JOBFILE.Plane_Numbers_List 					= guiprops.caljob.plane_numbers_list;
JOBFILE.Plane_Spacing 						= guiprops.caljob.plane_spacing;
JOBFILE.Z_Grid_Start   						= guiprops.caljob.z_grid_start;
JOBFILE.X_Grid_Spacing 						= guiprops.caljob.x_grid_spacing;
JOBFILE.Y_Grid_Spacing 						= guiprops.caljob.y_grid_spacing;
JOBFILE.Grid_Point_Diameter 				= guiprops.caljob.grid_point_diameter;
JOBFILE.Multiple_Level_Plate 				= guiprops.caljob.multiple_level_plate;
JOBFILE.Plate_Level_Number 					= guiprops.caljob.plate_level_number;
JOBFILE.Plate_Level_Spacing 				= guiprops.caljob.plate_level_spacing;
JOBFILE.X_Plate_Level_Shift 				= guiprops.caljob.x_plate_level_shift;
JOBFILE.Y_Plate_Level_Shift 				= guiprops.caljob.y_plate_level_shift;
JOBFILE.Front_To_Back_Plate_Thickness 		= guiprops.caljob.front_to_back_plate_thickness;
JOBFILE.Grid_Points_Bright 					= guiprops.caljob.grid_points_bright;
JOBFILE.Length_Unit 						= guiprops.caljob.length_unit;
JOBFILE.X_Pixel_Number 						= guiprops.caljob.x_pixel_number;
JOBFILE.Y_Pixel_Number 						= guiprops.caljob.y_pixel_number;
JOBFILE.ImageExtension						= guiprops.caljob.image_extension;
JOBFILE.CalImageList                        = guiprops.caljob.calimagelist;
JOBFILE.targetsidecam1                      = guiprops.caljob.targetsidecam1;
JOBFILE.targetsidecam2                      = guiprops.caljob.targetsidecam2;
guiprops.caljob.JOBFILE=JOBFILE;

[calibration_data,calibration_plane_data]=camera_calibration_new(JOBFILE);

guiprops.caljob.calibration_data=calibration_data;
guiprops.caljob.calibration_plane_data=calibration_plane_data;

guidata(hObject,guiprops);
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% FIT MODEL TYPE SECTION                                           %%%%%
%%%% All handles contained within this section pertain to Model Fit   %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function modeltypemenu_Callback(hObject, eventdata, handles)

guiprops            = guidata(hObject);
guiprops.modeltypemenu  = get(hObject,'Value');
guiprops.caljob.modeltype=guiprops.modeltypemenu;
guidata(hObject,guiprops);

end

function modeltypemenu_CreateFcn(hObject, eventdata, handles)
end
function convergebox_Callback(hObject, eventdata, handles)

end

function convergebox_CreateFcn(hObject, eventdata, handles)

end


function fitmodelbutton_Callback(hObject, eventdata, handles)

% first convert from image coords to regular coords, then perform fit.
guiprops = guidata(hObject);

set(handles.convergebox,'String','fitting...');
pause(0.001);


allx1data(:,1)   = guiprops.caljob.calibration_data.x_world_full{1};  % contains all x,y,z data for camera 1
allx1data(:,2)   = guiprops.caljob.calibration_data.y_world_full{1};
allx1data(:,3)   = guiprops.caljob.calibration_data.z_world_full{1};

allx2data(:,1)   = guiprops.caljob.calibration_data.x_world_full{2};       % contains all x,y,z data for camera 2
allx2data(:,2)   = guiprops.caljob.calibration_data.y_world_full{2};
allx2data(:,3)   = guiprops.caljob.calibration_data.z_world_full{2};

allX1data(:,1)   = guiprops.caljob.calibration_data.x_image_full{1};       % contains all X,Y data for camera 1
allX1data(:,2)   = guiprops.caljob.calibration_data.y_image_full{1};

allX2data(:,1)   = guiprops.caljob.calibration_data.x_image_full{2};
allX2data(:,2)   = guiprops.caljob.calibration_data.y_image_full{2};
%keyboard;
% [~,in1]=sort(allx1data(:,3),1,'ascend');
% [~,in2]=sort(allx2data(:,3),1,'ascend');
% allx1data=allx1data(in1,:);
% allX1data=allX1data(in1,:);
% allx2data=allx2data(in2,:);
% allX2data=allX2data(in2,:);

% A1          = guiprops.imagemat{ind1(1)};                          % read the first image for cam 1 to get the size of the image to convert units
% A2          = guiprops.imagemat{ind2(1)};
% [rA1, cA1]  = size(A1);
% [rA2 ,cA2]  = size(A2);

rA1=guiprops.caljob.y_pixel_number;

% %JJC: I think these are unnecessary now that I treat the image coordinates
% %everywhere else in the code as pixel centered.  But is it better to just
% %fix them here and be done with it?
% allX1data(:,1)  = allX1data(:,1)-0.5;      % convert from image coords to regular coords ((0,0) at bottom left corner)
% allX1data(:,2)  = rA1-allX1data(:,2)+0.5;
% allX2data(:,1)  = allX2data(:,1)-0.5;      % convert from image coords to regular coords ((0,0) at bottom left corner)
% allX2data(:,2)  = rA1-allX2data(:,2)+0.5;

%JJC: Still need to flip the coordinate system?
%I think so since these positions are referenced to the upper left, and 
%everywhere else the image coordinates will be reference to lower left.
allX1data(:,2)  = rA1-allX1data(:,2) + 1; %plus 1 is because the first index is at 1, not 0
allX2data(:,2)  = rA1-allX2data(:,2) + 1;

modeltype   = guiprops.caljob.modeltype;
optionsls   = guiprops.caljob.optionsls;

guiprops.caljob.allx1data=allx1data;
guiprops.caljob.allx2data=allx2data;
guiprops.caljob.allX1data=allX1data;
guiprops.caljob.allX2data=allX2data;

% fprintf('saving initial world and image coordinates....\n');
% save('calxworld.mat','allx1data');save('calyworld.mat','allx2data');save('calximage.mat','allX1data');save('calyimage.mat','allX2data');

[a_cam1, a_cam2, aXcam1, aYcam1, aXcam2, aYcam2, convergemessage] = fitmodels(allx1data,...
    allx2data,allX1data,allX2data,modeltype,optionsls);

%inselfcal(allx1data,allx2data,allX1data,alX2data);
guiprops.caljob.aXcam1 = aXcam1;     % save the mapping coefficients
guiprops.caljob.aYcam1 = aYcam1;
guiprops.caljob.aXcam2 = aXcam2;
guiprops.caljob.aYcam2 = aYcam2;
%[aXcam1 aYcam1 aXcam2 aYcam2]

%guiprops.scal   = 0;

guiprops.caljob.a_cam1 = a_cam1;
guiprops.caljob.a_cam2 = a_cam2;

% Do I need this section set?
%%%                                                                         %
% % zero_plane_cam1 = median(ind1);                                             %
% % zero_plane_cam2 = median(ind2);                                             %
% %                                                                             %
% % guiprops.ind1   = ind1;                                                     %
% % guiprops.ind2   = ind2;                                                     %
% % guiprops.zpc1   = zero_plane_cam1;                                          %
% % guiprops.zpc2   = zero_plane_cam2;                                          %
                                                                            %
guiprops.convergemessage = convergemessage;                                 %
guiprops.caljob.convergemessage = convergemessage; 
%
set(handles.convergebox,'String',guiprops.convergemessage);                 %
%keyboard;
%
guidata(hObject,guiprops);                                                  %
%%%                                                                         %
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%SELF CALIBRATION

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --- Executes on button press in setupselfcaljob.
function setupselfcaljob_Callback(hObject, eventdata, handles)
selfcalprana;
end

% --- Executes on button press in loadselfcaljob.
function loadselfcaljob_Callback(hObject, eventdata, handles)
guiprops        = guidata(hObject);
[f1,p1]=uigetfile('*.*','Load Self Calibration Job File');
load([p1,f1]);
guiprops.selfcaljob=Data;
clear Data;
guidata(hObject,guiprops);
end

% --- Executes on button press in selfcalbutton.
function selfcalbutton_Callback(hObject, eventdata, handles)

guiprops        = guidata(hObject);
caldata=guiprops.caljob;
selfcaljob=guiprops.selfcaljob;

reftrue=1;
 while(reftrue~=0)
     
 [caldatamod]=selfcalibration_main(caldata,selfcaljob);
 % Updating calibration data after each iteration of self calibration
 caldata.allx1data=caldatamod.allx1data;
 caldata.allx2data=caldatamod.allx2data;
 caldata.aXcam1=caldatamod.aXcam1;
 caldata.aYcam1=caldatamod.aYcam1;
 caldata.aXcam2=caldatamod.aXcam2;
 caldata.aYcam2=caldatamod.aYcam2;
 caldata.convergemessage = caldatamod.convergemessage; 
 clear caldatamod;
 %[caldatamod]=selfcalibration_v1(caldata,selfcaljob);
 
  refine= input('Do you want to Refine? (Y/N):','s');
  
  if strcmpi(refine,'Y')
      
      reftrue=1;
      
  else
      
      reftrue=0;
  end
 
 end
 guiprops.caljob=caldata;
guidata(hObject,guiprops);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%STEREO RECONSTRUCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% --- Executes on button press in pranajobbutton.
function pranajobbutton_Callback(hObject, eventdata, handles)
prana;
end

% --- Executes on button press in loadpranajob.
function loadpranajob_Callback(hObject, eventdata, handles)
guiprops        = guidata(hObject);
[f2,p2]=uigetfile('*.*','Load Prana Job File');
load([p2,f2]);
guiprops.planarjob=Data;
clear Data;
guidata(hObject,guiprops);
end

function willertreconstruction_Callback(hObject, eventdata, handles)
guiprops = guidata(hObject);
if (get(hObject,'Value') == get(hObject,'Max'))
	guiprops.rectype{1}='Willert';
else
	guiprops.rectype{1}='';
end
guidata(hObject,guiprops);

end
function soloffreconstruction_Callback(hObject, eventdata, handles)

guiprops = guidata(hObject);

if (get(hObject,'Value') == get(hObject,'Max'))
	guiprops.rectype{2}='Soloff';
else
	guiprops.rectype{2}='';
end

guidata(hObject,guiprops);

end

% --- Executes on button press in recbutton.
function recbutton_Callback(hObject, eventdata, handles)
guiprops = guidata(hObject);
rectype=guiprops.rectype;
planarjob=guiprops.planarjob;
caldata=guiprops.caljob;
[outputdirlist,scaling]=stereoreconstruction(caldata,planarjob,rectype);

guiprops.stereorecdirlist=outputdirlist;
guiprops.scaling=scaling;
%keyboard;

guidata(hObject,guiprops);
end
