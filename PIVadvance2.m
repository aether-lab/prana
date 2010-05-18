function varargout = PIVadvance2(varargin)


% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @PIVadvance2_OpeningFcn, ...
    'gui_OutputFcn',  @PIVadvance2_OutputFcn, ...
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


%% opening funciton for figure
function PIVadvance2_OpeningFcn(hObject, eventdata, handles, varargin)

handles.syscolor=get(hObject,'color');
handles.loaddirec=[pwd '\'];

handles.data.imdirec=pwd;
handles.data.imbase='ImgA';
handles.data.imzeros='6';
handles.data.imext='tif';
handles.data.imcstep='2';
handles.data.imfstep='1';
handles.data.imfstart='1';
handles.data.imfend='1';

handles.data.wrmag='1';
handles.data.wrsamp='1';
handles.data.wrsep='1';
handles.data.maskname='';
handles.data.batchname='Proc1';
handles.data.outdirec=pwd;

handles.data.PIV0.winres='32,32';
handles.data.PIV0.winsize='64,64';
handles.data.PIV0.winauto='1';
handles.data.PIV0.gridres='8,8';
handles.data.PIV0.gridbuf='8,8';
handles.data.PIV0.BWO='0,0';
handles.data.PIV0.corr='2';
handles.data.PIV0.RPCd='2.8';
handles.data.PIV0.val='0';
handles.data.PIV0.valtype='2';
handles.data.PIV0.thresh='0';
handles.data.PIV0.valsize='7,7;7,7';
handles.data.PIV0.valthresh='3,2';
handles.data.PIV0.valuthresh='-16,16';
handles.data.PIV0.valvthresh='-16,16';
handles.data.PIV0.outbase='PIV_';
handles.data.PIV0.write='1';

handles.data.PIV1=handles.data.PIV0;
handles.data.PIV2=handles.data.PIV0;

handles.data.passes='2';
handles.data.cpass=num2str(get(handles.listbox3,'Value'));
handles.data.method='1';
handles.data.velsmooth='0';
handles.data.velsmoothfilt='2';
handles.data.velinterp='3';
handles.data.iminterp='1';

handles.data0=handles.data;
handles.Njob=num2str(size(get(handles.listbox5,'String'),1));
handles.Cjob=num2str(get(handles.listbox5,'String'));
handles=rmfield(handles,'data');

handles=update_data(handles);

handles.output = hObject;
guidata(hObject, handles);


%output
function varargout = PIVadvance2_OutputFcn(hObject, eventdata, handles)
varargout{1} = handles.output;


%PIV pass window
function listbox3_Callback(hObject, eventdata, handles)
if str2double(handles.Njob)>0
    handles=set_PIVcontrols(handles);
    guidata(hObject,handles)
end

function listbox3_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


%% PIV method
function popupmenu3_Callback(hObject, eventdata, handles)
if str2double(handles.Njob)>0
    handles.data.method=num2str(get(hObject,'Value'));
    if get(hObject,'Value')>1
        set(handles.popupmenu8,'backgroundcolor',[1 1 1]);
        if get(hObject,'Value')==3
            set(handles.popupmenu6,'backgroundcolor',[1 1 1]);
        else
            set(handles.popupmenu6,'backgroundcolor',0.5*[1 1 1]);
        end
        if get(handles.checkbox10,'Value')==1
            set(handles.edit45,'backgroundcolor',[1 1 1]);
        else
            set(handles.edit45,'backgroundcolor',0.5*[1 1 1]);
        end
    else
        set(handles.popupmenu8,'backgroundcolor',0.5*[1 1 1]);
        set(handles.popupmenu6,'backgroundcolor',0.5*[1 1 1]);
        set(handles.edit45,'backgroundcolor',0.5*[1 1 1]);
    end
    for e=1:str2double(handles.data.passes)
        eval(['handles.data.PIV' num2str(e) '.gridres = get(handles.edit24,''String'');'])
        eval(['handles.data.PIV' num2str(e) '.gridbuf = get(handles.edit46,''String'');'])
    end
    load_data(handles);
    guidata(hObject,handles)
end

function popupmenu3_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


%% add pass
function pushbutton7_Callback(hObject, eventdata, handles)
if str2double(handles.Njob)>0
    N=str2double(handles.data.passes);
    eval(['handles.data=setfield(handles.data,''PIV' num2str(N+1) ''',handles.data.PIV0);']);
    eval(['handles.data.PIV' num2str(N+1) '.outbase=[''Pass'' num2str(N+1) ''_''];']);
    handles.data.passes=num2str(N+1);
    load_PIVlist(handles);
    handles=set_PIVcontrols(handles);
    guidata(hObject,handles)
end


%% delete pass
function pushbutton8_Callback(hObject, eventdata, handles)
if str2double(handles.Njob)>0
    if str2double(handles.data.passes)>1
        p=get(handles.listbox3,'Value');
        eval(['handles.data=rmfield(handles.data,''PIV' num2str(p) ''');']);
        for e=(p+1):str2double(handles.data.passes)
            eval(['handles.data=setfield(handles.data,''PIV' num2str(e-1) ''',handles.data.PIV' num2str(e) ');']);
            eval(['handles.data.PIV' num2str(e-1) '.outbase=[''Pass'' num2str(e-1) ''_''];']);
            if e==str2double(handles.data.passes)
                eval(['handles.data=rmfield(handles.data,''PIV' num2str(e) ''');']);
            end
        end
        handles.data.passes=num2str(str2double(handles.data.passes)-1);
        set(handles.listbox3,'Value',1);
        load_PIVlist(handles)
        guidata(hObject,handles)
    end
    handles=set_PIVcontrols(handles);
    guidata(hObject,handles)
end

%% window resolution
function edit23_Callback(hObject, eventdata, handles)
if str2double(handles.Njob)>0
    eval(['handles.data.PIV' handles.data.cpass '.winres = get(hObject,''String'');'])
    a=get(hObject,'String');
    wx=str2double(a(1:(strfind(a,',')-1)));
    wy=str2double(a((strfind(a,',')+1):end));
    if get(handles.checkbox5,'Value')==1
        xbin = 2^(ceil(log(2*wx)/log(2)));
        ybin = 2^(ceil(log(2*wy)/log(2)));
        eval(['handles.data.PIV' handles.data.cpass '.winsize = [num2str(xbin) '','' num2str(ybin)];'])
        eval(['set(handles.edit25,''String'',handles.data.PIV' handles.data.cpass '.winsize)'])
    end
    if wx*wy<256
        set(hObject,'backgroundcolor',[1 0.5 0]);
    else
        set(hObject,'backgroundcolor',[1 1 1]);
    end
    guidata(hObject,handles)
end

function edit23_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


%% grid resolution
function edit24_Callback(hObject, eventdata, handles)
if str2double(handles.Njob)>0
    if str2double(handles.data.method)==1
        for e=1:str2double(handles.data.passes)
            eval(['handles.data.PIV' num2str(e) '.gridres = get(hObject,''String'');'])
        end
    else
        eval(['handles.data.PIV' handles.data.cpass '.gridres = get(hObject,''String'');'])
    end
    guidata(hObject,handles)
end

function edit24_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


%% window size
function edit25_Callback(hObject, eventdata, handles)
if str2double(handles.Njob)>0
    if get(handles.checkbox5,'Value')==0
        A=get(hObject,'String');
        Wx=str2double(A(1:(strfind(A,',')-1)));
        Wy=str2double(A((strfind(A,',')+1):end));
        B=get(handles.edit23,'String');
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
        guidata(hObject,handles)
    else
        eval(['set(hObject,''String'',handles.data.PIV' handles.data.cpass '.winsize);'])
    end
end

function edit25_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


%%  autowindow switch
function checkbox5_Callback(hObject, eventdata, handles)
if str2double(handles.Njob)>0
    eval(['handles.data.PIV' handles.data.cpass '.winauto = num2str(get(hObject,''Value''));'])
    handles=set_PIVcontrols(handles);
    guidata(hObject,handles)
end


%% correlation technique
function popupmenu4_Callback(hObject, eventdata, handles)
if str2double(handles.Njob)>0
    eval(['handles.data.PIV' handles.data.cpass '.corr = num2str(get(hObject,''Value''));'])
    guidata(hObject,handles)
    N=handles.data.cpass;
    A=eval(['handles.data.PIV' num2str(N)]);
    if str2double(A.corr)==1
        set(handles.edit26,'backgroundcolor',0.5*[1 1 1]);
    else
        set(handles.edit26,'backgroundcolor',[1 1 1]);
    end
end

function popupmenu4_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


%% RPC diameter
function edit26_Callback(hObject, eventdata, handles)
if str2double(handles.Njob)>0
    eval(['handles.data.PIV' handles.data.cpass '.RPCd = get(hObject,''String'');'])
    guidata(hObject,handles)
    if get(handles.popupmenu4,'Value')==2
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

function edit26_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


%% validation window size
function edit27_Callback(hObject, eventdata, handles)
if str2double(handles.Njob)>0
    A=get(hObject,'String');
    if strcmp(A(end),';')
        A=A(1:end-1);
    end
    set(hObject,'String',A)
    eval(['handles.data.PIV' handles.data.cpass '.valsize = get(hObject,''String'');'])
    guidata(hObject,handles)
end

function edit27_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


%% validation v threshold
function edit28_Callback(hObject, eventdata, handles)
if str2double(handles.Njob)>0
    eval(['handles.data.PIV' handles.data.cpass '.valvthresh = get(hObject,''String'');'])
    a=get(hObject,'String');
    tx=str2double(a(1:(strfind(a,',')-1)));
    ty=str2double(a((strfind(a,',')+1):end));
    if get(handles.checkbox6,'Value')+get(handles.checkbox8,'Value')==2
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

function edit28_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


%% validation switch
function checkbox6_Callback(hObject, eventdata, handles)
if str2double(handles.Njob)>0
    eval(['handles.data.PIV' handles.data.cpass '.val = num2str(get(hObject,''Value''));'])
    handles=set_PIVcontrols(handles);
    guidata(hObject,handles)
end


%% write plt switch
function checkbox7_Callback(hObject, eventdata, handles)
if str2double(handles.Njob)>0
    eval(['handles.data.PIV' handles.data.cpass '.write = num2str(get(hObject,''Value''));'])
    handles=set_PIVcontrols(handles);
    guidata(hObject,handles)
end


%% validation u threhsold
function edit29_Callback(hObject, eventdata, handles)
if str2double(handles.Njob)>0
    eval(['handles.data.PIV' handles.data.cpass '.valuthresh = get(hObject,''String'');'])
    a=get(hObject,'String');
    tx=str2double(a(1:(strfind(a,',')-1)));
    ty=str2double(a((strfind(a,',')+1):end));
    if get(handles.checkbox6,'Value')+get(handles.checkbox8,'Value')==2
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

function edit29_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


%% validation threshold
function edit30_Callback(hObject, eventdata, handles)
if str2double(handles.Njob)>0
    eval(['handles.data.PIV' handles.data.cpass '.valthresh = get(hObject,''String'');'])
    guidata(hObject,handles)
end

function edit30_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


%% validation threshold switch
function checkbox8_Callback(hObject, eventdata, handles)
if str2double(handles.Njob)>0
    eval(['handles.data.PIV' handles.data.cpass '.thresh = num2str(get(hObject,''Value''));'])
    handles=set_PIVcontrols(handles);
    guidata(hObject,handles)
end

%% write plt basename
function edit31_Callback(hObject, eventdata, handles)
if str2double(handles.Njob)>0
    eval(['handles.data.PIV' handles.data.cpass '.outbase = get(hObject,''String'');'])
    guidata(hObject,handles)
end

function edit31_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


%% image directory
function edit32_Callback(hObject, eventdata, handles)
if str2double(handles.Njob)>0
    handles.data.imdirec = get(hObject,'String');
    load_imlist(handles);
    guidata(hObject,handles)
end

function edit32_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%% image directory browse
function pushbutton9_Callback(hObject, eventdata, handles)
if str2double(handles.Njob)>0
    D = handles.data.imdirec;
    handles.data.imdirec = uigetdir(handles.data.imdirec);
    if handles.data.imdirec==0
        handles.data.imdirec = D;
    end
    set(handles.edit32,'string',handles.data.imdirec);
    load_imlist(handles);
    guidata(hObject,handles)
end


%% output directory
function edit33_Callback(hObject, eventdata, handles)
if str2double(handles.Njob)>0
    handles.data.outdirec = get(hObject,'String');
    guidata(hObject,handles)
end

function edit33_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


%% browse output directory
function pushbutton10_Callback(hObject, eventdata, handles)
if str2double(handles.Njob)>0
    D = handles.data.outdirec;
    handles.data.outdirec = uigetdir(handles.data.outdirec);
    if handles.data.outdirec==0
        handles.data.outdirec = D;
    end
    set(handles.edit33,'string',handles.data.outdirec);
    guidata(hObject,handles)
end


%% image list
function listbox4_Callback(hObject, eventdata, handles)

function listbox4_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


%% image base name
function edit34_Callback(hObject, eventdata, handles)
if str2double(handles.Njob)>0
    handles.data.imbase = get(hObject,'String');
    load_imlist(handles);
    guidata(hObject,handles)
end

function edit34_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


%% image extension
function edit35_Callback(hObject, eventdata, handles)
if str2double(handles.Njob)>0
    handles.data.imext = get(hObject,'String');
    if strcmp(handles.data.imext(1),'.')
        set(hObject,'backgroundcolor','r');
    else
        set(hObject,'backgroundcolor',[1 1 1]);
    end
    load_imlist(handles);
    guidata(hObject,handles)
end

function edit35_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


%% image pair step
function edit36_Callback(hObject, eventdata, handles)
if str2double(handles.Njob)>0
    handles.data.imcstep = get(hObject,'String');
    load_imlist(handles);
    guidata(hObject,handles)
end

function edit36_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


%% image frame step
function edit37_Callback(hObject, eventdata, handles)
if str2double(handles.Njob)>0
    handles.data.imfstep = get(hObject,'String');
    load_imlist(handles);
    guidata(hObject,handles)
end

function edit37_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


%% image frame start
function edit38_Callback(hObject, eventdata, handles)
if str2double(handles.Njob)>0
    handles.data.imfstart = get(hObject,'String');
    load_imlist(handles);
    guidata(hObject,handles)
end

function edit38_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


%% image frame end
function edit39_Callback(hObject, eventdata, handles)
if str2double(handles.Njob)>0
    handles.data.imfend = get(hObject,'String');
    load_imlist(handles);
    guidata(hObject,handles)
end

function edit39_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


%% magnification factor
function edit40_Callback(hObject, eventdata, handles)
if str2double(handles.Njob)>0
    handles.data.wrmag = get(hObject,'String');
    guidata(hObject,handles)
    if str2double(get(hObject,'String'))<=0
        set(hObject,'backgroundcolor','r')
    else
        set(hObject,'backgroundcolor',[1 1 1])
    end
end

function edit40_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


%% image mask name
function edit41_Callback(hObject, eventdata, handles)
if str2double(handles.Njob)>0
    handles.data.maskname = get(hObject,'String');
    guidata(hObject,handles)
end

function edit41_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


%% browse image mask
function pushbutton11_Callback(hObject, eventdata, handles)
if str2double(handles.Njob)>0
    D = handles.data.maskname;
    [a,b] = uigetfile('*.*');
    handles.data.maskname = [b a];
    if handles.data.maskname==0
        handles.data.maskname = D;
    end
    set(handles.edit41,'string',handles.data.maskname);
    guidata(hObject,handles)
end


%% batch processing filename
function edit42_Callback(hObject, eventdata, handles)
if str2double(handles.Njob)>0
    Jlist0=char(get(handles.listbox5,'String'));
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
    set(handles.listbox5,'String',char(Jlist));
    guidata(hObject,handles)
end

function edit42_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%% pulse separation
function edit43_Callback(hObject, eventdata, handles)
if str2double(handles.Njob)>0
    handles.data.wrsep = get(hObject,'String');
    guidata(hObject,handles)
    if str2double(get(hObject,'String'))<=0
        set(hObject,'backgroundcolor','r')
    else
        set(hObject,'backgroundcolor',[1 1 1])
    end
end

function edit43_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


%% image file zeros
function edit44_Callback(hObject, eventdata, handles)
if str2double(handles.Njob)>0
    handles.data.imzeros = get(hObject,'String');
    load_imlist(handles);
    guidata(hObject,handles)
end

function edit44_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%% validation method
function popupmenu5_Callback(hObject, eventdata, handles)
if str2double(handles.Njob)>0
    eval(['handles.data.PIV' handles.data.cpass '.valtype = num2str(get(hObject,''Value''));'])
    guidata(hObject,handles)
end

function popupmenu5_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


%% velocity smoothing Gaussian width
function edit45_Callback(hObject, eventdata, handles)
if str2double(handles.Njob)>0
    handles.data.velsmoothfilt=get(hObject,'String');
    guidata(hObject,handles)
end

function edit45_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


%% image interpolation method
function popupmenu6_Callback(hObject, eventdata, handles)
if str2double(handles.Njob)>0
    handles.data.iminterp=num2str(get(hObject,'Value'));
    guidata(hObject,handles)
end

function popupmenu6_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


%% velocity smoothing switch
function checkbox10_Callback(hObject, eventdata, handles)
if str2double(handles.Njob)>0
    handles.data.velsmooth = num2str(get(hObject,'Value'));
    if get(hObject,'Value')==1 && get(handles.popupmenu3,'Value')>1
        set(handles.edit45,'backgroundcolor',[1 1 1]);
    else
        set(handles.edit45,'backgroundcolor',0.5*[1 1 1]);
    end
    guidata(hObject,handles)
end


%% grid buffer
function edit46_Callback(hObject, eventdata, handles)
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

function edit46_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


%% bulk window offset
function edit47_Callback(hObject, eventdata, handles)
if str2double(handles.Njob)>0
    if str2double(handles.data.cpass)==1
        eval(['handles.data.PIV' handles.data.cpass '.BWO = get(hObject,''String'');'])
        guidata(hObject,handles)
    else
        eval(['set(hObject,''String'',handles.data.PIV' handles.data.cpass '.BWO)'])
    end
end

function edit47_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


%% velocity interpolation method
function popupmenu8_Callback(hObject, eventdata, handles)
if str2double(handles.Njob)>0
    handles.data.velinterp=num2str(get(hObject,'Value'));
    guidata(hObject,handles)
end

function popupmenu8_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


%% load IM list
function load_imlist(handles)
dir_struct = dir(handles.data.imdirec);
if isempty(dir_struct)
    set(handles.edit32,'backgroundcolor','r');
else
    set(handles.edit32,'backgroundcolor',[1 1 1]);
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
    set(handles.listbox4,'backgroundcolor','r');
else
    set(handles.listbox4,'backgroundcolor',[1 1 1]);
    if length(idf)~=size(files,1)
        set(handles.listbox4,'backgroundcolor','r');
    end
end
files=files(idf,:);
filesf=cell(length(files),1);
for e=1:size(files,1)
    filesf(e)={[char(files(e,1)) ' and ' char(files(e,2))]};
end
set(handles.listbox4,'String',filesf,'Value',1)


%% load PIV list
function load_PIVlist(handles)
f=cell(str2double(handles.data.passes),1);
for e=1:str2double(handles.data.passes)
    f(e)={['Pass ' num2str(e)]};
end
set(handles.listbox3,'String',f);


%% load PIV data
function handles=set_PIVcontrols(handles)
N=get(handles.listbox3,'Value');
A=eval(['handles.data.PIV' num2str(N)]);
set(handles.edit23,'string',A.winres);
set(handles.edit25,'string',A.winsize);
set(handles.checkbox5,'Value',str2double(A.winauto));
set(handles.edit24,'string',A.gridres);
set(handles.edit46,'string',A.gridbuf);
set(handles.edit47,'string',A.BWO);
set(handles.popupmenu4,'Value',str2double(A.corr));
set(handles.edit26,'string',str2double(A.RPCd));
set(handles.checkbox6,'Value',str2double(A.val));
set(handles.popupmenu5,'Value',str2double(A.valtype));
set(handles.checkbox8,'Value',str2double(A.thresh));
set(handles.edit27,'string',A.valsize);
set(handles.edit30,'string',A.valthresh);
set(handles.edit29,'string',A.valuthresh);
set(handles.edit28,'string',A.valvthresh);
set(handles.checkbox7,'Value',str2double(A.write));
set(handles.edit31,'string',A.outbase);
handles.data.cpass=num2str(N);

if str2double(A.val)==1
    set(handles.edit27,'backgroundcolor',[1 1 1]);
    set(handles.edit30,'backgroundcolor',[1 1 1]);
    set(handles.popupmenu5,'backgroundcolor',[1 1 1]);
    if str2double(A.thresh)==1
        a=get(handles.edit29,'String');
        tx=str2double(a(1:(strfind(a,',')-1)));
        ty=str2double(a((strfind(a,',')+1):end));
        if tx>=ty
            set(handles.edit29,'backgroundcolor','r');
        else
            set(handles.edit29,'backgroundcolor',[1 1 1]);
        end
        a=get(handles.edit28,'String');
        tx=str2double(a(1:(strfind(a,',')-1)));
        ty=str2double(a((strfind(a,',')+1):end));
        if tx>=ty
            set(handles.edit28,'backgroundcolor','r');
        else
            set(handles.edit28,'backgroundcolor',[1 1 1]);
        end
    else
        set(handles.edit29,'backgroundcolor',0.5*[1 1 1]);
        set(handles.edit28,'backgroundcolor',0.5*[1 1 1]);
    end
else
    set(handles.edit27,'backgroundcolor',0.5*[1 1 1]);
    set(handles.edit30,'backgroundcolor',0.5*[1 1 1]);
    set(handles.popupmenu5,'backgroundcolor',0.5*[1 1 1]);
    set(handles.edit29,'backgroundcolor',0.5*[1 1 1]);
    set(handles.edit28,'backgroundcolor',0.5*[1 1 1]);
end

if str2double(A.write)==1
    set(handles.edit31,'backgroundcolor',[1 1 1]);
else
    set(handles.edit31,'backgroundcolor',0.5*[1 1 1]);
end

if str2double(A.winauto)==1
    B=get(handles.edit23,'String');
    Rx=str2double(B(1:(strfind(B,',')-1)));
    Ry=str2double(B((strfind(B,',')+1):end));
    Rx = 2^(ceil(log(2*Rx)/log(2)));
    Ry = 2^(ceil(log(2*Ry)/log(2)));
    set(handles.edit25,'backgroundcolor',0.5*[1 1 1]);
    set(handles.edit25,'String',[num2str(Rx) ',' num2str(Ry)]);
    eval(['handles.data.PIV' handles.data.cpass '.winsize = get(handles.edit25,''String'');'])
else
    set(handles.edit25,'backgroundcolor',[1 1 1]);
end

if get(handles.popupmenu4,'Value')==2
    if str2double(get(handles.edit26,'String'))<2
        if str2double(get(handles.edit26,'String'))==0
            set(handles.edit26,'backgroundcolor','r');
        else
            set(handles.edit26,'backgroundcolor',[1 0.5 0]);
        end
    else
        set(handles.edit26,'backgroundcolor',[1 1 1]);
    end
else
    set(handles.edit26,'backgroundcolor',0.5*[1 1 1]);
end

a=get(handles.edit23,'String');
wx=str2double(a(1:(strfind(a,',')-1)));
wy=str2double(a((strfind(a,',')+1):end));
if wx*wy<256
    set(handles.edit23,'backgroundcolor',[1 0.5 0]);
else
    set(handles.edit23,'backgroundcolor',[1 1 1]);
end

if get(handles.checkbox7,'Value')==0 && get(handles.listbox3,'Value')==str2double(handles.data.passes)
    set(handles.checkbox7,'backgroundcolor',[1 0.5 0]);
else
    set(handles.checkbox7,'backgroundcolor',handles.syscolor);
end



%% load extra data
function [handles]=load_data(handles)
set(handles.edit32,'String',handles.data.imdirec);
set(handles.edit34,'String',handles.data.imbase);
set(handles.edit44,'String',handles.data.imzeros);
set(handles.edit35,'String',handles.data.imext);
set(handles.edit36,'String',handles.data.imcstep);
set(handles.edit37,'String',handles.data.imfstep);
set(handles.edit38,'String',handles.data.imfstart);
set(handles.edit39,'String',handles.data.imfend);

set(handles.edit40,'String',handles.data.wrmag);
set(handles.edit43,'String',handles.data.wrsep);
set(handles.edit48,'String',handles.data.wrsamp);
set(handles.edit41,'String',handles.data.maskname);
set(handles.edit42,'String',handles.data.batchname);
set(handles.edit33,'String',handles.data.outdirec);

set(handles.popupmenu3,'Value',str2double(handles.data.method));
set(handles.popupmenu8,'Value',str2double(handles.data.velinterp));
set(handles.popupmenu6,'Value',str2double(handles.data.iminterp));
set(handles.edit45,'String',handles.data.velsmoothfilt);
set(handles.checkbox10,'Value',str2double(handles.data.velsmooth));

if get(handles.popupmenu3,'Value')>1
    set(handles.popupmenu8,'backgroundcolor',[1 1 1]);
    if get(handles.popupmenu3,'Value')==3
        set(handles.popupmenu6,'backgroundcolor',[1 1 1]);
    else
        set(handles.popupmenu6,'backgroundcolor',0.5*[1 1 1]);
    end
    if get(handles.checkbox10,'Value')==1
        set(handles.edit45,'backgroundcolor',[1 1 1]);
    else
        set(handles.edit45,'backgroundcolor',0.5*[1 1 1]);
    end
else
    set(handles.popupmenu8,'backgroundcolor',0.5*[1 1 1]);
    set(handles.popupmenu6,'backgroundcolor',0.5*[1 1 1]);
    set(handles.edit45,'backgroundcolor',0.5*[1 1 1]);
end

if strcmp(handles.data.imext(1),'.')
    set(handles.edit35,'backgroundcolor','r');
else
    set(handles.edit35,'backgroundcolor',[1 1 1]);
end

if isempty(dir(handles.data.imdirec))
    set(handles.edit32,'backgroundcolor','r');
else
    set(handles.edit32,'backgroundcolor',[1 1 1]);
end

if str2double(get(handles.edit40,'String'))<=0
    set(handles.edit40,'backgroundcolor','r')
else
    set(handles.edit40,'backgroundcolor',[1 1 1])
end

if str2double(get(handles.edit43,'String'))<=0
    set(handles.edit43,'backgroundcolor','r')
else
    set(handles.edit43,'backgroundcolor',[1 1 1])
end

if str2double(get(handles.edit48,'String'))<=0
    set(handles.edit48,'backgroundcolor','r')
else
    set(handles.edit48,'backgroundcolor',[1 1 1])
end


%% Update all data
function handles=update_data(handles)
if str2double(handles.Njob)==0
    set(handles.listbox3,'String','','backgroundcolor',0.5*[1 1 1]);
    set(handles.listbox4,'String','','backgroundcolor',0.5*[1 1 1]);
    set(handles.listbox5,'String','','backgroundcolor','r');
    set(handles.edit23,'String','','backgroundcolor',0.5*[1 1 1]);
    set(handles.edit24,'String','','backgroundcolor',0.5*[1 1 1]);
    set(handles.edit25,'String','','backgroundcolor',0.5*[1 1 1]);
    set(handles.edit26,'String','','backgroundcolor',0.5*[1 1 1]);
    set(handles.edit27,'String','','backgroundcolor',0.5*[1 1 1]);
    set(handles.edit28,'String','','backgroundcolor',0.5*[1 1 1]);
    set(handles.edit29,'String','','backgroundcolor',0.5*[1 1 1]);
    set(handles.edit30,'String','','backgroundcolor',0.5*[1 1 1]);
    set(handles.edit31,'String','','backgroundcolor',0.5*[1 1 1]);
    set(handles.edit32,'String','','backgroundcolor',0.5*[1 1 1]);
    set(handles.edit33,'String','','backgroundcolor',0.5*[1 1 1]);
    set(handles.edit34,'String','','backgroundcolor',0.5*[1 1 1]);
    set(handles.edit35,'String','','backgroundcolor',0.5*[1 1 1]);
    set(handles.edit36,'String','','backgroundcolor',0.5*[1 1 1]);
    set(handles.edit37,'String','','backgroundcolor',0.5*[1 1 1]);
    set(handles.edit38,'String','','backgroundcolor',0.5*[1 1 1]);
    set(handles.edit39,'String','','backgroundcolor',0.5*[1 1 1]);
    set(handles.edit40,'String','','backgroundcolor',0.5*[1 1 1]);
    set(handles.edit41,'String','','backgroundcolor',0.5*[1 1 1]);
    set(handles.edit42,'String','','backgroundcolor',0.5*[1 1 1]);
    set(handles.edit43,'String','','backgroundcolor',0.5*[1 1 1]);
    set(handles.edit44,'String','','backgroundcolor',0.5*[1 1 1]);
    set(handles.edit45,'String','','backgroundcolor',0.5*[1 1 1]);
    set(handles.edit46,'String','','backgroundcolor',0.5*[1 1 1]);
    set(handles.edit47,'String','','backgroundcolor',0.5*[1 1 1]);
    set(handles.edit48,'String','','backgroundcolor',0.5*[1 1 1]);
    set(handles.checkbox5,'Value',0,'backgroundcolor',handles.syscolor);
    set(handles.checkbox6,'Value',0,'backgroundcolor',handles.syscolor);
    set(handles.checkbox7,'Value',0,'backgroundcolor',handles.syscolor);
    set(handles.checkbox8,'Value',0,'backgroundcolor',handles.syscolor);
    set(handles.checkbox10,'Value',0,'backgroundcolor',handles.syscolor);
    set(handles.popupmenu3,'Value',1,'backgroundcolor',0.5*[1 1 1]);
    set(handles.popupmenu4,'Value',1,'backgroundcolor',0.5*[1 1 1]);
    set(handles.popupmenu5,'Value',1,'backgroundcolor',0.5*[1 1 1]);
    set(handles.popupmenu6,'Value',1,'backgroundcolor',0.5*[1 1 1]);
    set(handles.popupmenu8,'Value',1,'backgroundcolor',0.5*[1 1 1]);
else
    a=get(handles.listbox5,'String');
    eval(['handles.data=handles.' char(a(str2double(handles.Cjob),:)) ';'])

    set(handles.listbox3,'String','','backgroundcolor',[1 1 1]);
    set(handles.listbox4,'String','','backgroundcolor',[1 1 1]);
    set(handles.listbox5,'backgroundcolor',[1 1 1]);
    set(handles.edit23,'String','','backgroundcolor',[1 1 1]);
    set(handles.edit24,'String','','backgroundcolor',[1 1 1]);
    set(handles.edit25,'String','','backgroundcolor',[1 1 1]);
    set(handles.edit26,'String','','backgroundcolor',[1 1 1]);
    set(handles.edit27,'String','','backgroundcolor',[1 1 1]);
    set(handles.edit28,'String','','backgroundcolor',[1 1 1]);
    set(handles.edit29,'String','','backgroundcolor',[1 1 1]);
    set(handles.edit30,'String','','backgroundcolor',[1 1 1]);
    set(handles.edit31,'String','','backgroundcolor',[1 1 1]);
    set(handles.edit32,'String','','backgroundcolor',[1 1 1]);
    set(handles.edit33,'String','','backgroundcolor',[1 1 1]);
    set(handles.edit34,'String','','backgroundcolor',[1 1 1]);
    set(handles.edit35,'String','','backgroundcolor',[1 1 1]);
    set(handles.edit36,'String','','backgroundcolor',[1 1 1]);
    set(handles.edit37,'String','','backgroundcolor',[1 1 1]);
    set(handles.edit38,'String','','backgroundcolor',[1 1 1]);
    set(handles.edit39,'String','','backgroundcolor',[1 1 1]);
    set(handles.edit40,'String','','backgroundcolor',[1 1 1]);
    set(handles.edit41,'String','','backgroundcolor',[1 1 1]);
    set(handles.edit42,'String','','backgroundcolor',[1 1 1]);
    set(handles.edit43,'String','','backgroundcolor',[1 1 1]);
    set(handles.edit44,'String','','backgroundcolor',[1 1 1]);
    set(handles.edit45,'String','','backgroundcolor',[1 1 1]);
    set(handles.edit46,'String','','backgroundcolor',[1 1 1]);
    set(handles.edit47,'String','','backgroundcolor',[1 1 1]);
    set(handles.edit48,'String','','backgroundcolor',[1 1 1]);
    set(handles.popupmenu3,'Value',1,'backgroundcolor',[1 1 1]);
    set(handles.popupmenu4,'Value',1,'backgroundcolor',[1 1 1]);
    set(handles.popupmenu5,'Value',1,'backgroundcolor',[1 1 1]);
    set(handles.popupmenu6,'Value',1,'backgroundcolor',[1 1 1]);
    set(handles.popupmenu8,'Value',1,'backgroundcolor',[1 1 1]);

    load_imlist(handles);
    load_PIVlist(handles);
    load_data(handles);
    handles=set_PIVcontrols(handles);
end


%% select job file
function listbox5_Callback(hObject, eventdata, handles)
if str2double(handles.Njob)>0
    Jlist=char(get(handles.listbox5,'String'));
    eval(['handles.' Jlist(str2double(handles.Cjob),:) '=handles.data;']);
    handles.Cjob=num2str(get(hObject,'Value'));
    handles=update_data(handles);
    guidata(hObject,handles)
end

function listbox5_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


%% sampling frequency
function edit48_Callback(hObject, eventdata, handles)
if str2double(handles.Njob)>0
    handles.data.wrsamp = get(hObject,'String');
    guidata(hObject,handles)
    if str2double(get(hObject,'String'))<=0
        set(hObject,'backgroundcolor','r')
    else
        set(hObject,'backgroundcolor',[1 1 1])
    end
end

function edit48_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end





%% mask preview button
function pushbutton21_Callback(hObject, eventdata, handles)
if str2double(handles.Njob)>0
    try
        if length(handles.data.maskname)>0
            mask = double(imread(handles.data.maskname));
            mask=flipud(mask);
        else
            mask = 1+0*double(imread([handles.data.imdirec '\' handles.data.imbase sprintf(['%0.' handles.data.imzeros 'i.' handles.data.imext],str2double(handles.data.imfstart))]));
        end
        e=0;
    catch
        msgbox('Mask or Image Frame not found');
        e=-1;
    end
    if e==0
        mask(mask>0)=1;
        figure,hold on

        L=size(mask);
        Ps=get(0,'screensize');
        if (L(1)/Ps(4))>(L(2)/Ps(3))
            dy=0.8*Ps(4);
            dx=dy*(L(2)/L(1));
        else
            dx=0.8*Ps(3);
            dy=dx*(L(1)/L(2));
        end
        set(gcf,'position',[(Ps(3)-dx)/2 (Ps(4)-dy)/2 dx dy])
        set(gcf,'color',0.5*[1 1 1])

        %     im1=double(imread([handles.data.imdirec '\' handles.data.imbase sprintf(['%0.' handles.data.imzeros 'i.' handles.data.imext],str2double(handles.data.imfstart))]));
        %     im1=flipud(im1)/255;
        %     im1(mask==0)=0.5*im1(mask==0);
        imagesc(mask,[0 1]),axis image,colormap gray,axis off,set(gca,'position',[0 0 1 1])

        A=get(handles.edit46,'String');
        G=[str2double(A(1:(strfind(A,',')-1))) str2double(A((strfind(A,',')+1):end))];
        A=get(handles.edit24,'String');
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

        plot(X,Y,'r.')
        
    end
end




%% File Menu
function Untitled_1_Callback(hObject, eventdata, handles)
if str2double(handles.Njob)>0
    set(handles.Untitled_5,'Enable','on')
    set(handles.Untitled_6,'Enable','on')
    set(handles.Untitled_12,'Enable','on')
else
    set(handles.Untitled_5,'Enable','off')
    set(handles.Untitled_6,'Enable','off')
    set(handles.Untitled_12,'Enable','off')
end


%% New Job
function Untitled_2_Callback(hObject, eventdata, handles)
if str2double(handles.Njob)>0
    Jlist=char(get(handles.listbox5,'String'));
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
        Jlist=char(get(handles.listbox5,'String'));
        Jlist={Jlist;Data.batchname};
    else
        Jlist={Data.batchname};
    end
    handles.Njob=num2str(str2double(handles.Njob)+1);
    handles.Cjob=handles.Njob;
    set(handles.listbox5,'String',Jlist,'Value',str2double(handles.Cjob));
    handles=update_data(handles);
    guidata(hObject,handles)
end


%% Load Job
function Untitled_4_Callback(hObject, eventdata, handles)
if str2double(handles.Njob)>0
    Jlist=char(get(handles.listbox5,'String'));
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
                        Jlist=char(get(handles.listbox5,'String'));
                        Jlist={Jlist;Data.batchname};
                    else
                        Jlist={Data.batchname};
                    end
                    handles.Njob=num2str(str2double(handles.Njob)+1);
                    handles.Cjob=handles.Njob;
                    set(handles.listbox5,'String',Jlist,'Value',str2double(handles.Cjob));

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

%% Save Job
function Untitled_5_Callback(hObject, eventdata, handles)
if str2double(handles.Njob)>0
    Data=handles.data;
    uisave('Data',[handles.loaddirec handles.data.batchname '.mat']);
end


%% Delete Job
function Untitled_6_Callback(hObject, eventdata, handles)
if str2double(handles.Njob)>0
    del=questdlg('Are You Sure?','Delete Job File','yes','no','yes');
    if strcmp(del,'yes')==1
        Jlist=char(get(handles.listbox5,'String'));
        handles=rmfield(handles,Jlist(str2double(handles.Cjob),:));
        id=setdiff(1:str2double(handles.Njob),str2double(handles.Cjob));
        Jlist=Jlist(id,:);
        set(handles.listbox5,'String',Jlist,'Value',1);
        handles.Cjob='1';
        handles.Njob=num2str(str2double(handles.Njob)-1);
        handles=update_data(handles);
        guidata(hObject,handles)
    end
end


%% Run Job
function Untitled_10_Callback(hObject, eventdata, handles)
if str2double(handles.Njob)>0
    Data=handles.data;
    PIVadvance2code(Data);
end


%% Run All
function Untitled_11_Callback(hObject, eventdata, handles)
if str2double(handles.Njob)>0
    Jlist=char(get(handles.listbox5,'String'));
    for e=1:size(Jlist,1)
        Data=eval(['handles.' Jlist(e,:)]);
        PIVadvance2code(Data);
    end
end


%% Execute Menu
function Untitled_8_Callback(hObject, eventdata, handles)
if str2double(handles.Njob)>0
    set(handles.Untitled_10,'Enable','on')
    set(handles.Untitled_11,'Enable','on')
else
    set(handles.Untitled_10,'Enable','off')
    set(handles.Untitled_11,'Enable','off')
end

%% Help Menu
function Untitled_9_Callback(hObject, eventdata, handles)


%% Copy Job
function Untitled_12_Callback(hObject, eventdata, handles)
if str2double(handles.Njob)>0
    Jlist=char(get(handles.listbox5,'String'));
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
        Jlist=char(get(handles.listbox5,'String'));
        Jlist={Jlist;Data.batchname};
    else
        Jlist={Data.batchname};
    end
    handles.Njob=num2str(str2double(handles.Njob)+1);
    handles.Cjob=handles.Njob;
    set(handles.listbox5,'String',Jlist,'Value',str2double(handles.Cjob));
    handles=update_data(handles);
    guidata(hObject,handles)
end

%% Help PIV Controls > Processing
function Untitled_18_Callback(hObject, eventdata, handles)
PIVhelp(5);

%% Help PIV Controls > Validation
function Untitled_19_Callback(hObject, eventdata, handles)
PIVhelp(6);

%% Help Image Controls
function Untitled_13_Callback(hObject, eventdata, handles)
PIVhelp(2);

%% Help Write Controls
function Untitled_14_Callback(hObject, eventdata, handles)
PIVhelp(3);

%% Help PIV Controls
function Untitled_15_Callback(hObject, eventdata, handles)

%% Help Interpolation
function Untitled_16_Callback(hObject, eventdata, handles)
PIVhelp(8);

%% Help Job List
function Untitled_17_Callback(hObject, eventdata, handles)
PIVhelp(9);


%% Help About
function Untitled_20_Callback(hObject, eventdata, handles)


%% Help PIVadvance2
function Untitled_21_Callback(hObject, eventdata, handles)
PIVhelp(1);


%% Help PIV Controls Methods
function Untitled_22_Callback(hObject, eventdata, handles)
PIVhelp(4);

%% Help PIVcontrols Output
function Untitled_27_Callback(hObject, eventdata, handles)
PIVhelp(7)


%% Load Images
function Untitled_24_Callback(hObject, eventdata, handles)
if str2double(handles.Njob)>0
    [A,d]=uigetfile('*.tif;*.tiff;*.bmp;*.jpg;*.jpeg','LOAD PLT FILES',[handles.data.outdirec '\'],'Multiselect','on');
else
    [A,d]=uigetfile('*.tif;*.tiff;*.bmp;*.jpg;*.jpeg','LOAD PLT FILES',[pwd '\'],'Multiselect','on');
end
for e=1:length(A)
    im=imread([d A{e}]);
    assignin('base',A{e}(1:(findstr(A{e},'.')-1)),im);
end

%% Load Plt
function Untitled_25_Callback(hObject, eventdata, handles)
if str2double(handles.Njob)>0
    [A,d]=uigetfile('*.plt','LOAD PLT FILES',[handles.data.outdirec '\'],'Multiselect','on');
else
    [A,d]=uigetfile('*.plt','LOAD PLT FILES',[pwd '\'],'Multiselect','on');
end
A=sort(A);
for e=1:length(A)
    openplt([d A{e}]);
end

%% Load Menu
function Untitled_23_Callback(hObject, eventdata, handles)



%% Generic Open plt function
function openplt(name)

% if nargin==1
    varlist='';
% end
    
%opens file
fid = fopen(name);
if fid==-1
    error('cannot open file')
end
status = 1;

%default
dataformat='point';

while status==1
    
    %read a line from file
    mark=ftell(fid);
    temp=fgetl(fid);
    c=sscanf(temp,'%s');
    
    %header lines
    if isstrprop(c(1),'alpha')==1

        %reads variable names
        if length(strfind(lower(temp),'variables'))~=0
            try
                ind=strfind(temp,'"');
                nvar=length(ind)/2;
                q=0;
                var=cell(nvar,1);
                for n=1:2:2*nvar
                    q=q+1;
                    name=temp(ind(n)+1:ind(n+1)-1);
                    var(q)={name(isspace(name)==0)};
                end
            catch
                error('Cannot read variable names')
            end
            %replaces with new variable names
            if length(varlist)~=0
                if length(varlist)~=length(var)
                    error('Input variable list does not match file list')
                end
                var=varlist;
            end
        end

        %reads zone lengths
        if length(strfind(lower(temp),'zone '))~=0
            indi=(strfind(temp,'I=')+2):((strfind(temp,'I=')+2)+strfind(temp(strfind(temp,'I=')+2:end),' ')-2);
            indj=(strfind(temp,'J=')+2):((strfind(temp,'J=')+2)+strfind(temp(strfind(temp,'J=')+2:end),' ')-2);
            try
                zstr=[strtrim(temp(indi)) ',' strtrim(temp(indj))];
            catch
                error('Cannot read zone lengths')
            end
            
        end
        
        %reads solution time (needs to be generalized)
        if length(strfind(lower(temp),'solutiontime '))~=0
            inds=(strfind(temp,'SOLUTIONTIME = ')+15):length(temp);
            try
                time=str2double(strtrim(temp(inds)));
            catch
            end
        end

        %reads data format
        if length(strfind(lower(temp),'datapacking'))~=0
            ind=strfind(temp,'=');
            try
                dataformat = temp(ind+1:ind+5);
            catch
                error('Unrecognized datapacking format')
            end
        end

    else
        %numeric data lines        
        if length(temp)~=0
            
            %reads point format
            if strcmp(dataformat,'point')
                fseek(fid,mark,-1);
                try
                    variable=(fscanf(fid,'%g',[nvar inf]))';
                catch
                    error('Zone length does not match available data')
                end
                status=fclose(fid);
                for q=1:nvar
                    eval([char(var(q)) '=reshape(variable(:,q),' char(zstr) ');'])
                end

            %reads block format
            elseif strcmp(dataformat,'block')
                fseek(fid,mark,-1);
                try
                    for q=1:nvar
                        eval([char(var(q)) '=' '(fscanf(fid,''%g'',[' zstr ']))'';'])
                    end
                catch
                    error('Zone length does not match available data')
                end
                status=fclose(fid);

            else
                error('Unrecognized datapacking format')
            end
            
        end
    end
end


%check for existing variables in the workspace
evars=0;
for q=1:nvar
    tstr = ['exist(''' char(var(q)) ''',''var'')'];
    Wans = evalc('evalin(''base'',tstr)');
    evars = evars + str2num(sscanf(Wans(7:end),'%s'));
end

%for appending data to workspace variables
if evars~=0
    
    %find existing variable dimensions
    tstr=['ndims(' char(var(1)) ')'];
    Wans=evalc('evalin(''base'',tstr)');
    D   = str2num(sscanf(Wans(7:end),'%s'));
    eval(['Din=ndims(' char(var(1)) ');']);

    %when loading multiple frames
    if Din==D-1
        try
            tstr=['size(' char(var(1)) ',' num2str(D) ')'];  
            Wans=evalc('evalin(''base'',tstr)');
            ldmax = str2num(sscanf(Wans(7:end),'%s'));
            for qq=1:nvar
                dstr='(';
                for b=1:D-1
                    dstr=[dstr ':,'];
                end
                dstr=[dstr num2str(ldmax+1) ')'];
                assignin('base',['loaded_var_' num2str(qq)],eval(char(var(qq))));
                evalin('base',[char(var(qq)) dstr '=' 'loaded_var_' num2str(qq) ';'])
                evalin('base',['clear ' 'loaded_var_' num2str(qq) ';']);
            end
        catch
            keyboard
            error('problem appending data in multiple frames')
        end

        %when loading second frame
    elseif Din==D
        try
            for qq=1:nvar
                dstr='(';
                for b=1:D
                    dstr=[dstr ':,'];
                end
                dstr=[dstr '2)'];
                assignin('base',['loaded_var_' num2str(qq)],eval(char(var(qq))));
                evalin('base',[char(var(qq)) dstr '=' 'loaded_var_' num2str(qq) ';'])
            end
        catch
            error('problem appending data to initial frame')
        end
        
    %when dimensions do not agree
    else
        rem=input('workspace variable mismatch, clear conflicting variables (y/n)? >> ');
        if strcmp(rem,'y')==1
            for q=1:length(var)
                evalin('base',['clear ' char(var(q)) ';']);
            end
        else
            error('workspace variable mismatch, terminating loader')
        end
    end
end

%for replacing/creating workspace variables
if evars==0 
    for qq=1:nvar
        assignin('base',char(var(qq)),eval(char(var(qq))));
    end
end

%load solution time
if evalin('base','exist(''T'',''var'')')==1
    assignin('base','loaded_solution_time',time);
    evalin('base','T(length(T)+1)=loaded_solution_time;');
    evalin('base','clear loaded_solution_time');
else
    assignin('base','T',time);
end

evalin('base','clear ans');


