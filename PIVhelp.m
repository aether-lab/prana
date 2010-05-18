function varargout = PIVhelp(varargin)
% PIVHELP M-file for PIVhelp.fig
%      PIVHELP, by itself, creates a new PIVHELP or raises the existing
%      singleton*.
%
%      H = PIVHELP returns the handle to a new PIVHELP or the handle to
%      the existing singleton*.
%
%      PIVHELP('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in PIVHELP.M with the given input arguments.
%
%      PIVHELP('Property','Value',...) creates a new PIVHELP or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before PIVhelp_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to PIVhelp_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help PIVhelp

% Last Modified by GUIDE v2.5 10-Mar-2008 21:51:17

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @PIVhelp_OpeningFcn, ...
                   'gui_OutputFcn',  @PIVhelp_OutputFcn, ...
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


% --- Executes just before PIVhelp is made visible.
function PIVhelp_OpeningFcn(hObject, eventdata, handles, varargin)
if length(varargin)>0
    set(handles.listbox1,'Value',varargin{1});
end
update_text(handles);

handles.output = hObject;
guidata(hObject, handles);

% UIWAIT makes PIVhelp wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = PIVhelp_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on selection change in listbox1.
function listbox1_Callback(hObject, eventdata, handles)
update_text(handles);

function listbox1_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit1_Callback(hObject, eventdata, handles)

function edit1_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




function update_text(handles)

A=get(handles.listbox1,'Value');
%PIVadvance
if A==1
    set(handles.edit1,'String',...
       {'PIVadvance:' 'This program is designed to read in a set of single exposure images, which are then evaluated used one of the cross-correlation algorithms to estimate the displacement of local regions within the image. [1]' ''...
        'References' ... 
        '[1] J Westerweel, "Fundamentals of digital particle image velocimetry." Meas Sci Tech, 1997. 8(12): p. 1379-1392.' ''});
end
%read controls
if A==2
    set(handles.edit1,'String',...
       {'General:' 'This section controls the image loading and processing order.  Images must be numerically indexed at end of filename and have a valid file extension.  If images do not show up in listbox, no files will process.' ''...
        'Image Directory:' 'Select the directory in which the image files are located.  Use the browse button to search for files. (string)' '' ...
        'Image Basename:' 'Input the image file basename excluding any number indicies or file extension (string)' '' ...
        'Zeros:' 'Input number of zeros in the filename indexing (integer)' '' ...
        'Extension:' 'Image file extension, must be from available formats in imread.m (string)' '' ...
        'Image Correlation Step:' 'Input the separation between correlated image pair (integer)' '' ...
        'Frame Step:' 'Input the separation between velocity fields (integer)' '' ...
        'Frame Start:' 'Input start index of the PIV processing (integer)' '' ...
        'Frame End:' 'Input end index of the PIV processing (integer)' '' ...
        'Image Processing Order:' 'This listbox displays the processing order of the image set using the specified input parameters.  When this listbox is red, some or all of the images cannot be found and the job will not run.' ''});
end
%write controls
if A==3
    set(handles.edit1,'String',...
       {'General:' 'This section controls the plt file controls including experimental parameters, processing mask, and job names.  Defaults are in pixel units.' ''...
        'Pulse Separation:' 'Input the pulse separation in microseconds (positive real, us)' '' ...
        'Sampling:' 'Input the sampling frequency in Hz between sequential vector fields (positive real, Hz)' ''...
        'Magnification' 'Input the magnification in microns/pixel (positive real, um/pix' '' ...
        'Processing Mask' 'Input the complete filename for the image mask: 1 mask / data set. (string) Use the browse button to search for masks.  The filename should have a valid image extension.  Preview will display the current image mask and the PIV gridpoints using the current processing grid.  This mask does not filter images, It only indicates where vectors will be evaluated.' '' ...
        'Job Name' 'Input the job name, which used when saving and loading jobs only.  Name cannot begin with numbers or orther invalid characters since the name is used to create a structure. (string)' '' ...
        'Output Directory' 'Input the directory where plt vector text files will be placed. Use the browse button to search for files. (string)' ''});
end
%PIV controls
if A==4
    set(handles.edit1,'String',...
       {'PIV Method:' 'This section controls the PIV evaluation technique used to process images.  The listbox indicates the passes that will be used to process each image pair.' ''...
        'Multipass (DWO):' 'Evaluation grid is fixed for all passes.  Different Window sizes, Resolution, and Validation may be applied for each pass.  After each pass, the integer portion of the velocity field is used in a second order Discrete Window Offset (DWO) to increase the correlation strength of subsequent passes [1].  This method can be very computationally expensive, especially for initial passes, where the velocity field is highly oversampled.' '' ...
        'Multigrid (DWO):' 'Evaluation grid is variable between passes.  Different Window sizes, Resolution, and Validation may be applied for each pass.  After each pass, the velocity field is interpolated onto a new grid, which is then used in a second order Discrete Window Offset (DWO) to increase the correlation strength of subsequent passes (activates velocity interpolation function and velocity smoothing) [1].  This method greatly decreases computational costs and results in easier validation (since the velocity fields are not as oversampled).  However, the velocity estimate is limited to the accuracy of the interpolation scheme.  This error is often substantially below the pixel discretization used in the DWO.' '' ...
        'Multigrid (Deform):' 'Evaluation grid is variable between passes.  Different Window sizes, Resolution, and Validation may be applied for each pass.  After each pass, the velocity field is interpolated onto the image grid, which is then used in a second order continuous offset of each pixel.  The image is then interpolated back onto a rectilinear grid (activates velocity interpolation function, velocity smoothing, and image interpolation menu).  This method is a higher order iterative multigrid approach, where higher computational times are required to interpolate the images.  In addition, the accuracy of the interpolation scheme and velocity field plays a critical role in the sucess of the method.  However, this method often shows substantial increases in accuracy, especially for high shear flows [2].' '' ...
        'Multigrid (Ensemble DWO):' 'Correlates each image pair using a Multigrid DWO.  However, prior to subpixel interpolation, the correlations are averaged over the entire data set.  The subpixel estimator is then applied to the ensemble correlations to generate the ensemble averaged velocity field.  Will only output 1 plt.  This method can be very useful for low SNR images, where the correlation strength is low.  Also, very high spatial resolutions can be achieved due to the increased correlation in each frame pair.  However, this will only provide an average flowfield [3].' '' ...
        'References:' ...
        '[1] J Westerweel, D Dabiri, and M Gharib, "The effect of a discrete window offset on the accuracy of cross-correlation analysis of digital PIV recordings" Exp Fluids, 1997. 23(1): p. 20-28.' ...
        '[2] F Scarano, Iterative image deformation methods in PIV" Meas Sci Tech, 2002. 13: p. R1-R19.' ...
        '[3] CD Meinhart, ST Wereley, and JG Santiago, "A PIV algorithm for estimating time-averaged velocity fields" Journal of Fluids Engineering, 2000. 122: p. 285-289.' });
end
%Processing
if A==5
    set(handles.edit1,'String',...
       {'General:' 'This section sets the processing controls for the PIV pass displayed in the PIV controls listbox' ''...
        'Grid Resolution: "Gx,Gy"' 'Input the grid resolution in the x and y directions (pixel units)' '' ...
        'Grid Buffer: "Bx,By"' 'Input the buffer for the grid to be taken around the image (pixel units).  With grid buffer set to zero, buffer defaults to half the grid resolution.  General practice is to set the buffer to half the window size in order to prevent clipping of the window.' '' ...
        'Bulk Window Offset: "U0,V0"' 'Input the Bulk window offset (pixel units).  Can only be performed on the first pass.  Applies a Discrete Window Offset (DWO) [1] prior to the first pass based upon the uniform image shift given by the bulk window offset.  Use this feature to account for images with large displacement, but be cautious of any near wall or low flow regions.' '' ...
        'Window Resolution: "Rx,Ry"' 'Input the desired window resolution (pixel units).  At least 5 particles should be within this rectangular size on average.  In addition, the window size generally produces poor results when the area of this region is below 256 pixels^2.  The displacement should not exceed 1/4 of the first pass window size.  Several different resolutions should be examined for each data set.  This is the most important parameter in the processing!!' '' ...
        'Window Size: "Wx,Wy"' 'Input the window size (pixel units).  This setting is the actual rectangular window region taken from the image.  This setting, in conjunction with the window resolution controls the spatial window mask applied to each image region.  The rectangular window region is tapered with a Gaussian window, where the taper width is controlled by the spatial resolution.  The resolution is equal to the average width of the Gaussian.  Selecting "auto" will automatically select the nearest power of 2 window size that is at least a 50% Gaussian taper to minimize the effects of wraparound aliasing and spectral leakage [1].  When in doubt, leave on auto.' '' ...
        'Correlation:' 'SCC (Standard Fourier based Cross-Correlation) [2]' 'RPC (Robust Phase Correlation) [2]' '' ...
        'Diameter:' 'Input RPC diameter (pixel units).  Should be "roughly" equivalent to the particle-image diameter. [2]' '' ...
        'References:' ...
        '[1] A Eckstein, P Vlachos, 2008, �Improved DPIV Accuracy Using Advanced Windowing Techniques� Proc. ASME Fluids Engineering Conference, FEDSM2008-55152, Jacksonville, FL.' ...
        '[2] A Eckstein, P Vlachos, 2007, "A robust phase correlation DPIV processing algorithm for time resolved measurements" in PIV symposium. Rome, Italy.' ''});
end
%Validation
if A==6
    set(handles.edit1,'String',...
       {'General:' 'This section sets the validation controls for the PIV pass displayed in the PIV controls listbox.  The validation and thresholding checkboxes can be used to activate or deactivate.' '' ...
        'Validation Windows: "Wx1,Wy1;Wx2,Wy2..."' 'Input the size of the windows used to collect local statistics (vector units).  Number of passes is defined by this field and the threshold field (must be equivalent).' '' ...
        'Validation Threshold: "T1,T2..."' 'Input the threshold used in each validation pass (positive real).  The number of validation thresholds must be equivalent to number of validation windows.  Thresholds are determined using the Universal Outlier Detection criterion [1].' '' ...
        'Umin,Umax and Vmin,Vmax: "Umin,Umax" and "Vmin,Vmax"' 'Input the velocity thresholds on valid measurements.  Make sure that your flow does not exceed these limits for all frames.  Improper thresholding can be hard to detect by observing the vector field alone, use Eval matrix in output plt (Eval=100 for threshold limitation).' '' ...
        'Output basename' 'Input the basename used to create plt files (string).  Can be selected individually for each pass using the write output checkbox.  Outputs are numbered according to the number of the first image in the correlation pair.' '' ...
        'References:' ...
        '[1] J Westerweel, "Universal outlier detection for PIV data" Exp Fluids, 2005. 39: p. 1096-1100.' ''});
end
%Output
if A==7
    set(handles.edit1,'String',...
       {'Output basename' 'Input the basename used to create plt files (string).  Can be selected individually for each pass using the write output checkbox.  Outputs are numbered according to the number of the first image in the correlation pair.' ''});
end
%Interpolation
if A==8
    set(handles.edit1,'String',...
       {'General:' 'This section controls the velocity and image interpolation, which can only be activated for multigrid methods.' '' ...
        'Smoothing:' 'Input size for gaussian smoothing filter (positive real, vector units).  Smoothing is generally applied in order to remove high frequency noise that corrupts velocity interpolation.  Not recommended for general use.' '' ...
        'Velocity Interpolation Function:' 'Select desired method to be used in the velocity interpolation onto new grids.' '' ...
        'Image Interpolation Method:' 'Select the image interpolation function used to interpolate the deformed image back onto a rectilinear grid. [1]' '' ...
        'References:' ...
        '[2] F Scarano, Iterative image deformation methods in PIV" Meas Sci Tech, 2002. 13: p. R1-R19.' ''});
end
%Job list
if A==9
    set(handles.edit1,'String',...
       {'Job Listbox:' 'This lists the currently available PIV jobs.  To create, load, save, or delete jobs go to the file menu.  To run jobs, go to the execute menu.  To see the different settings of each job, select between the list of jobs.  Jobs cannot have names that begin with numbers or other invalid characters, as information is stored in a structure array.' ''});
end





