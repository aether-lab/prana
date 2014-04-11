function pranaprocessing(Data,I1,I2,maskname)
imClass = 'single';

SaveIMdeform = 1;

%% --- Read Formatted Parameters ---
%input/output directory
if ispc
    imbase=[Data.imdirec '\' Data.imbase];
    imbase2=[Data.imdirec2 '\' Data.imbase2];
    maskbase=[Data.maskdirec '\' Data.maskbase];
    pltdirec=[Data.outdirec '\'];
else
    imbase=[Data.imdirec '/' Data.imbase];
    imbase2=[Data.imdirec2 '/' Data.imbase2];
    maskbase=[Data.maskdirec '/' Data.maskbase];
    pltdirec=[Data.outdirec '/'];
end

if nargin<3
    I1 = str2double(Data.imfstart):str2double(Data.imfstep):str2double(Data.imfend);
    I2 = I1+str2double(Data.imcstep);
end

%processing mask
if strcmp(Data.masktype,'none')
    mask = 1+0*double(imread([imbase sprintf(['%0.' Data.imzeros 'i.' Data.imext],I1(1))]));
    maskname=[];
elseif strcmp(Data.masktype,'static')
    mask = double(imread(Data.staticmaskname));
    mask = flipud(mask);
    maskname=[];
elseif strcmp(Data.masktype,'dynamic')
    %     if nargin<4
    maskfend=str2double(Data.imfstart)+str2double(Data.imfstep)*length(str2double(Data.imfstart):str2double(Data.imfstep):str2double(Data.imfend))-1;
    maskname=str2double(Data.imfstart):str2double(Data.imfstep):maskfend;
    %     end
end

%method and passes
P = str2double(Data.passes);
Method = {'Multipass','Multigrid','Deform','Ensemble','EDeform','Multiframe','ForwardDeform'};
M = Method{str2double(Data.method)};
% Color channel
try
    if isfield(Data,'channel')
        channel = str2double(Data.channel);
        if isnan(channel);
            fprintf('Color channel returned nan.  Please check color designation\n and confirm that it takes a value of a string.\n Setting color channel to ''red.''\n')
            channel = 1;
            Data.channel = '1';
        end
    else
        channel = 1;
        Data.channel = '1';
        fprintf('Could not find color channel information, using ''red'' (first channel) as default\n')
    end
catch ME
    fprintf('Unknown issure with color channel, trying ''red'' (first channel) as default\n%s',ME.stack(1))
    channel = 1;
    Data.channel = '1';
end

%algorithm options
Velinterp = str2double(Data.velinterp);
Iminterp  = str2double(Data.iminterp);
Nmax      = str2double(Data.framestep);
ds        = str2double(Data.PIVerror);

%physical parameters
Mag  = str2double(Data.wrmag);
dt   = str2double(Data.wrsep);
Freq = str2double(Data.wrsamp);
%checking to makesure physical parameters are infact numbers.
if isnan(Mag) || isnan(dt) || isnan(Freq)
    if isnan(Mag)
        Mag = 1;
        fprintf('Magnifcation Value improporly set, changing the value to 1\n')
    end
    if isnan(dt)
        dt = 1;
        fprintf('Pulse Separation improporly set, changing the value to 1\n')
    end
    if isnan(Freq)
        Freq = 1;
        fprintf('Sampling Rate improporly set, changing the value to 1\n')
    end
end

% %initialization
% Effective size of window after Gaussian filtering
Wres = zeros(2, 2, P);

% Size of unfiltered window (pixels)
Wsize = zeros(P,2);

% Grid resolution
Gres = zeros(P,2);

% NOT SURE

% Not Sure
Gbuf            = zeros(P,2);
Corr            = cell(P,1);  %correlation type on each pass
D               = zeros(P,2);
Zeromean        = zeros(P,1);
Peaklocator     = zeros(P,1);
Velsmoothswitch = zeros(P,1);
Velsmoothfilt   = zeros(P,1);
Valswitch       = zeros(P,1);
UODswitch       = zeros(P,1);
Bootswitch      = zeros(P,1);
Threshswitch    = zeros(P,1);
Writeswitch     = zeros(P,1);
Peakswitch      = zeros(P,1);
UODwinsize      = zeros(P,2);
UODthresh       = zeros(P,1);
Bootper         = zeros(P,1);
Bootiter        = zeros(P,1);
Bootkmax        = zeros(P,1);
Uthresh         = zeros(P,2);
Vthresh         = zeros(P,2);
extrapeaks      = zeros(P,1);
PeakNum         = zeros(P,1);
PeakMag         = zeros(P,1);
PeakVel         = zeros(P,1);
wbase           = cell(0);
frac_filt       = zeros(P,1);
mindefloop      = zeros(P,1); % variables for the deformation convergences
maxdefloop      = zeros(P,1);
condefloop      = zeros(P,1);
saveplane       = zeros(P,1);
numDefPasses    = zeros(P,1); % This stores the number of interative deform passes that are performed

%read data info for each pass
for e=1:P
    
    %create structure for pass "e"
    eval(['A = Data.PIV' num2str(e) ';']);
    
    %store bulk window offset info
    if e==1
        BWO=[str2double(A.BWO(1:(strfind(A.BWO,',')-1))) str2double(A.BWO((strfind(A.BWO,',')+1):end))];
    end
    
    % Window Resolution (eventually move the str2double commands to the GUI callback)
    winres_com = regexp(A.winres,'[,;]');
    if length(winres_com) > 1%isfield(A,'winres1')
        %  Window resolutions for first image in correlation pair
        xwin_im1 = str2double(A.winres(1:winres_com(1)-1));
        ywin_im1 = str2double(A.winres(winres_com(1) + 1:winres_com(2)-1));
        
        %  Window resolutions for second image in correlation pair
        xwin_im2 = str2double(A.winres(winres_com(2)+1:winres_com(3)-1));
        ywin_im2 = str2double(A.winres(winres_com(3)+1:end));
    else
        xwin_im1 = str2double(A.winres(1:winres_com(1)-1));
        ywin_im1 = str2double(A.winres(winres_com(1)+1:end));
        xwin_im2 = xwin_im1;
        ywin_im2 = ywin_im1;
    end
    
    % Window resolution matrix
    Wres(:,:, e) = [xwin_im1 ywin_im1; xwin_im2 ywin_im2];
    
    % Window size and grid resolution
    Wsize(e,:) = [str2double(A.winsize(1:(strfind(A.winsize,',')-1))) str2double(A.winsize((strfind(A.winsize,',')+1):end))];
    Gres(e,:) = [str2double(A.gridres(1:(strfind(A.gridres,',')-1))) str2double(A.gridres((strfind(A.gridres,',')+1):end))];
    Gbuf(e,:) = [str2double(A.gridbuf(1:(strfind(A.gridbuf,',')-1))) str2double(A.gridbuf((strfind(A.gridbuf,',')+1):end))];
    Corr{e} = A.corr; %why do we subtract 1 from A.corr?  Just to make things more confusing? (SCC,RPC,GCC,FWC,SPC)
    D(e,:) = [str2double(A.RPCd(1:(strfind(A.RPCd,',')-1))) str2double(A.RPCd((strfind(A.RPCd,',')+1):end))];
    frac_filt(e) = str2double(A.frac_filt);
    Zeromean(e) = str2double(A.zeromean);
    Peaklocator(e) = str2double(A.peaklocator);
    Velsmoothswitch(e) = str2double(A.velsmooth);
    Velsmoothfilt(e) = str2double(A.velsmoothfilt);
    mindefloop(e) = str2double(A.deform_min);
    maxdefloop(e) = str2double(A.deform_max);
    condefloop(e) = str2double(A.deform_conv);
    
    if any(Wres(1,:,e)>Wsize(e,:)) || any(Wres(2,:,e)>Wsize(e,:))
        warning('warning:ResGraterThenSize','Pass %0.0f has a window resolution larger then the windown size!\n   [%0.0f,%0.0f;%0.0f %0.0f] > [%0.0f,%0.0f]\n',e,Wres(1,1,e),Wres(1,2,e),Wres(2,1,e),Wres(2,2,e),Wsize(e,1), Wsize(e,2))
    end
    
    %validation and thresholding
    Valswitch(e)=str2double(A.val);
    UODswitch(e)=str2double(A.uod);
    Bootswitch(e)=str2double(A.bootstrap);
    Threshswitch(e)=str2double(A.thresh);
    Writeswitch(e)=str2double(A.write);
    
    vpass=[0 strfind(A.uod_window,';') length(A.uod_window)+1];
    for q=1:(length(vpass)-1)
        B=A.uod_window((vpass(q)+1):(vpass(q+1)-1));
        UODwinsize(e,:,q)=[str2double(B(1:(strfind(B,',')-1))) str2double(B((strfind(B,',')+1):end))];
        UODthresh(e,q)=str2double(A.uod_thresh(1+2*(q-1)));
    end
    
    Bootper(e)=str2double(A.bootstrap_percentsampled);
    Bootiter(e)=str2double(A.bootstrap_iterations);
    Bootkmax(e)=str2double(A.bootstrap_passes);
    
    if str2double(A.thresh)==1
        Uthresh(e,:)=[str2double(A.valuthresh(1:(strfind(A.valuthresh,',')-1))) str2double(A.valuthresh((strfind(A.valuthresh,',')+1):end))];
        Vthresh(e,:)=[str2double(A.valvthresh(1:(strfind(A.valvthresh,',')-1))) str2double(A.valvthresh((strfind(A.valvthresh,',')+1):end))];
    else
        Uthresh(e,:)=[-inf,inf];
        Vthresh(e,:)=[-inf,inf];
    end
    
    extrapeaks(e)=str2double(A.valextrapeaks);
    
    %peak information
    Peakswitch(e)=str2double(A.savepeakinfo);
    PeakNum(e)=str2double(A.corrpeaknum);
    PeakMag(e)=str2double(A.savepeakmag);
    PeakVel(e)=str2double(A.savepeakvel);
    saveplane(e) = str2double(A.saveplane);
    
    %output directory
    wbase(e,:)={A.outbase};
    
end


maskSize=size(mask);
maskSize(3)=size(mask,3);
[XI,YI]=IMgrid(maskSize,[0 0]);
%cast the pixel coordinates to imClass for reduced memory
%index-based pixel coordinates are integers on center, but we need XI and YI to
%correspond to vector locations, which are centered on pixel edges for even sized
%windows.  This means the center of pixel index 1 is actually at coordinate
%0.5 from image edge, etc.
% We need to do this because the vector positions we will use
% for interpolation and deformation will be on the physical
% grid, but will will need to know where the pixels are located
% relative to them, not their indices.
%BUT they are pixel centered for ODD sized windows.  (FIX) but
%for now will ignore since even windows are more common.
XI = cast(XI,imClass) - 0.5;
YI = cast(YI,imClass) - 0.5;
if strcmpi(Data.input_vel_type,'static')  %'dynamic' is handled below in the processing loop
    Vel0 = load(Data.input_velocity);
    
    %velocity smoothing
    if Velsmoothswitch(1)==1
        [U,V]=VELfilt(Vel0.U(:,:,1),Vel0.V(:,:,1),UODwinsize(1,:,:),Velsmoothfilt(1));
        U = cast(U,imClass);
        V = cast(V,imClass);
        %when interpolating, assume saved coordinate system consistent with
        %vector locations (i.e. they are in vector grid positions, not
        %pixel-centered)
        
        %resample U, V and W from vector grid coordinates
        %V(X,Y,Z) onto UI,VI, and WI on pixel coordinates
        %[XI,YI] where XI and YI are a list of every
        %pixel centers in the image plane. Velinterp is the type of
        %interpolation to use.
        Vel0.U = VFinterp(Vel0.X,Vel0.Y,U,XI,YI,Velinterp);
        Vel0.V = VFinterp(Vel0.X,Vel0.Y,V,XI,YI,Velinterp);
    else
        %see above note about coordinate systems
        Vel0.U = cast(VFinterp(Vel0.X,Vel0.Y,Vel0.U(:,:,1),XI,YI,Velinterp),'single');
        Vel0.V = cast(VFinterp(Vel0.X,Vel0.Y,Vel0.V(:,:,1),XI,YI,Velinterp),'single');
    end
    VelInputFile = 1;
else
    Vel0.X = XI;
    Vel0.Y = YI;
    Vel0.U = BWO(1)*ones(size(XI),imClass);
    Vel0.V = BWO(2)*ones(size(XI),imClass);
    VelInputFile = 0;
end
wbase_org=wbase;

%% --- Evaluate Image Sequence ---
switch char(M)
    
    case {'Multipass','Multigrid','Deform','ForwardDeform'}
        %% --- Multipass, Multigrid, Deform
        frametime=zeros(length(I1),1);
        for q=1:length(I1)
            
            tf=tic;
            frametitle=['Frame' sprintf(['%0.' Data.imzeros 'i'],I1(q)) ' and Frame' sprintf(['%0.' Data.imzeros 'i'],I2(q))];
            
            %load image pair and flip coordinates
            if strcmpi(Data.imext,'mat') %read .mat file, image must be stored in variable 'I'
                loaddata=load([imbase sprintf(['%0.' Data.imzeros 'i.' Data.imext],I1(q))]);
                im1 = cast(loaddata.I,imClass);
                loaddata=load([imbase2 sprintf(['%0.' Data.imzeros 'i.' Data.imext],I2(q))]);
                im2 = cast(loaddata.I,imClass);
                loaddata =[];
            else
                im1 = cast(imread([imbase  sprintf(['%0.' Data.imzeros 'i.' Data.imext],I1(q))]),imClass);
                im2 = cast(imread([imbase2 sprintf(['%0.' Data.imzeros 'i.' Data.imext],I2(q))]),imClass);
            end
            
            % Specify which color channel(s) to consider
            % This was changed to greater then 2 because John had images
            % that were 4 channel with the last channel being a
            % transparency channel.
            if size(im1, 3) > 2
                %Extract only red channel
                if channel == 1;
                    im1 = im1(:,:,1);
                    im2 = im2(:,:,1);
                    %Extract only green channel
                elseif channel == 2;
                    im1 = im1(:,:,2);
                    im2 = im2(:,:,2);
                    %Extract only blue channel
                elseif channel == 3;
                    im1 = im1(:,:,3);
                    im2 = im2(:,:,3);
                    %Weighted average of channels (see rgb2gray for
                    %explanation of weighting factors)
                elseif channel == 4;
                    im1 = 0.2989 * im1(:, :, 1) + 0.5870 * im1(:, :, 2) + 0.1140 * im1(:, :, 3);
                    im2 = 0.2989 * im2(:, :, 1) + 0.5870 * im2(:, :, 2) + 0.1140 * im2(:, :, 3);
                    %Evenly weighted mean of channels
                elseif channel == 5;
                    im1 = (im1(:,:,1) + im1(:,:,2) + im1(:,:,3))/3;
                    im2 = (im2(:,:,1) + im2(:,:,2) + im2(:,:,3))/3;
                    %ensemble correlation of channels
                elseif channel == 6;
                    im1=im1(:,:,1:3);
                    im2=im2(:,:,1:3);
                end
                
            else
                %	Take only red channel
                im1 =im1(:,:,1);
                im2 =im2(:,:,1);
                channel = 1;
            end
            
            % Determine the number of channels in the image
            % to be deformed. This should be 3 for color
            % images or 1 for grayscale images.
            nChannels = size(im1, 3);
            
            % Flip images
            % flipud only works on 2D matices.  What about flipdim(im1,1) instead?
            im1 = im1(end:-1:1,:,:);%flipud(im1);
            im2 = im2(end:-1:1,:,:);%flipud(im2);
            
            % Determine size of images.
            imageSize = size(im1);
            
            % load dynamic mask and flip coordinates
            if strcmp(Data.masktype,'dynamic')
                mask = cast(imread([maskbase sprintf(['%0.' Data.maskzeros 'i.' Data.maskext],maskname(q))]),imClass);
                mask = flipud(mask);
            end
            
            % initialize grid and evaluation matrix for every pixel in image
            %cast the pixel coordinates to imClass for reduced memory
            %index-based pixel coordinates are integers on center, but we need XI and YI to
            %correspond to vector locations, which are centered on pixel edges for even sized
            %windows.  This means the center of pixel index 1 is actually at coordinate
            %0.5 from image edge, etc.
            % We need to do this because the vector positions we will use
            % for interpolation and deformation will be on the physical
            % grid, but will will need to know where the pixels are located
            % relative to them, not their indices.
            %BUT they are pixel centered for ODD sized windows.  (FIX) but
            %for now will ignore since even windows are more common.
            [XI,YI]=IMgrid(imageSize,[0 0]);
            XI = cast(XI,imClass) - 0.5;
            YI = cast(YI,imClass) - 0.5;
            
            %load dynamic input velocity, overwrite BWO or default
            if strcmpi(Data.input_vel_type,'dynamic')
                
                %for now, assume that input vector files are stored in .mat
                Vel0 = load(fullfile(Data.input_veldirec,[Data.input_velbase,sprintf(['%0.' Data.imzeros 'i.mat'],I1(q))]));
                
                %velocity smoothing, if first pass is smoothed, assumed
                %previous passs needed to be as well
                if Velsmoothswitch(1)==1
                    [U,V]=VELfilt(Vel0.U(:,:,1),Vel0.V(:,:,1),UODwinsize(1,:,:),Velsmoothfilt(1));
                    U = cast(U,imClass);
                    V = cast(V,imClass);
                    %when interpolating, assume saved coordinate system consistent with
                    %vector locations (i.e. they are in vector grid positions, not
                    %pixel-centered)
                    
                    %resample U, V and W from vector grid coordinates
                    %V(X,Y,Z) onto UI,VI, and WI on pixel coordinates
                    %[XI,YI] where XI and YI are a list of every
                    %pixel centers in the image plane. Velinterp is the type of
                    %interpolation to use.
                    Vel0.U = VFinterp(Vel0.X,Vel0.Y,U,XI,YI,Velinterp);
                    Vel0.V = VFinterp(Vel0.X,Vel0.Y,V,XI,YI,Velinterp);
                else
                    %see above note about coordinate systems
                    Vel0.U = cast(VFinterp(Vel0.X,Vel0.Y,Vel0.U(:,:,1),XI,YI,Velinterp),'single');
                    Vel0.V = cast(VFinterp(Vel0.X,Vel0.Y,Vel0.V(:,:,1),XI,YI,Velinterp),'single');
                end
                VelInputFile = 1;
            end
            
            UI = Vel0.U;
            VI = Vel0.V;
            %             UI = BWO(1)*ones(size(XI),imClass);
            %             VI = BWO(2)*ones(size(YI),imClass);
            
            % Preallocating variables
            corrtime=zeros(P,max(maxdefloop));
            valtime=zeros(P,max(maxdefloop));
            savetime=zeros(P,1);
            interptime=zeros(P,max(maxdefloop));
            deformtime=zeros(P,max(maxdefloop));
            defconvU = zeros(P,max(maxdefloop));
            defconvV = zeros(P,max(maxdefloop));
            
            e = 0; defloop = 1; % What is e?
            
            % before we start doing the new work, see if we need to deform
            % the input images to match the input velocity field
            if ( strcmpi(Data.input_vel_type,'dynamic') || strcmpi(Data.input_vel_type,'static'))  && ~isempty(regexpi(M,'Deform','once'))
                
                %Here, UI and VI are provided by the velocity input
                %file....
                %For deform, UI and VI will be used to deform the images, and
                %then downsampled in the next pass to the new grid
                %for addition to the next pass's displacement.
                %If not deform, then UI and VI just get downsampled
                %next pass and used for DWO.
                if strcmpi(M,'Deform')
                    t1=tic;
                    
                    % translate pixel locations, but
                    % since sincBlackmanInterp2 assumes
                    % coordinate system is pixel-centered, we need
                    % to convert back to index-coordinates for the
                    % deform.
                    %HOWEVER: we need to use these shifted
                    %coordinates later in vector coordinates, so
                    %move the -0.5 pixel correction to the call to
                    %the sincBlackmanInterp2 function
                    XD1 = XI-UI/2 ;
                    YD1 = YI-VI/2;
                    XD2 = XI+UI/2;
                    YD2 = YI+VI/2;
                    
                    % Preallocate memory for deformed images.
                    im1d = zeros(size(im1),imClass);
                    im2d = zeros(size(im2),imClass);
                    
                    % Deform images according to the interpolated velocity fields
                    for k = 1:nChannels % Loop over all of the color channels in the image
                        if Iminterp == 1 % Sinc interpolation (without blackman window)
                            im1d(:, :, k) = sincBlackmanInterp2(im1(:, :, k), XD1+0.5, YD1+0.5, 3, 'sinc');
                            im2d(:, :, k) = sincBlackmanInterp2(im2(:, :, k), XD2+0.5, YD2+0.5, 3, 'sinc');
                            
                        elseif Iminterp == 2 % Sinc interpolation with blackman filter
                            im1d(:, :, k) = sincBlackmanInterp2(im1(:, :, k), XD1+0.5, YD1+0.5, 3, 'blackman');
                            im2d(:, :, k) = sincBlackmanInterp2(im2(:, :, k), XD2+0.5, YD2+0.5, 3, 'blackman');
                        end
                    end
                    
                    %not sure this won't get overwritten below,
                    %but I think that normally the deformation is
                    %stored on the next pass?  Also, code can't get
                    %to deformation without finishing pass 1 (which
                    %resets defloop to 1 and stores the time at e+1)
                    %or incrementing defloop to 2, and storing time
                    %at e.  So it should always be safe to store at
                    %(e=1,defloop=1)
                    deformtime(1,defloop)=toc(t1);
                    
                elseif strcmpi(M,'ForwardDeform') %only the 2nd image is deformed
                    t1=tic;
                    
                    % translate pixel locations, but
                    % since sincBlackmanInterp2 assumes
                    % coordinate system is pixel-centered, we need
                    % to convert back to index-coordinates for the
                    % deform.
                    %HOWEVER: we need to use these shifted
                    %coordinates later in vector coordinates, so
                    %move the -0.5 pixel correction to the call to
                    %the sincBlackmanInterp2 function
                    XD1 = XI   ;
                    YD1 = YI   ;
                    XD2 = XI+UI;
                    YD2 = YI+VI;
                    
                    % Preallocate memory for deformed images.
                    im1d = im1;
                    im2d = zeros(size(im2),imClass);
                    
                    % Deform images according to the interpolated velocity fields
                    for k = 1:nChannels % Loop over all of the color channels in the image
                        if Iminterp == 1 % Sinc interpolation (without blackman window)
                            im2d(:, :, k) = sincBlackmanInterp2(im2(:, :, k), XD2+0.5, YD2+0.5, 3, 'sinc');
                        elseif Iminterp == 2 % Sinc interpolation with blackman filter
                            im2d(:, :, k) = sincBlackmanInterp2(im2(:, :, k), XD2+0.5, YD2+0.5, 3, 'blackman');
                        end
                    end
                    
                    %not sure this won't get overwritten below,
                    %but I think that normally the deformation is
                    %stored on the next pass?  Also, code can't get
                    %to deformation without finishing pass 1 (which
                    %resets defloop to 1 and stores the time at e+1)
                    %or incrementing defloop to 2, and storing time
                    %at e.  So it should always be safe to store at
                    %(e=1,defloop=1)
                    deformtime(1,defloop)=toc(t1);
                end
            end %pre-deformation for velocity input
            
            
            
            % This while statment is used to interatively move through the
            % deformations.  If the minimum number of loops hasn't been
            % reach it will keep iterating otherwise it should stop.
            while (e<P && defloop == 1) || (e<=P && defloop~=1)
                if defloop == 1
                    e=e+1;
                end
                
                t1=tic;
                [X,Y]=IMgrid(imageSize,Gres(e,:),Gbuf(e,:));
                S=size(X);X=X(:);Y=Y(:);
                
                if strcmp(M,'Multipass')
                    %UI and VI are already only defined on the grid points,
                    %which are the same at every iteration
                    Ub=UI(:);
                    Vb=VI(:);
                else
                    %resample UI(XI,YI) and VI at every pixel down to just
                    %the subset of interrogation points defined by Gres(e,:)
                    Ub = reshape(downsample(downsample( UI(Y(1):Y(end),X(1):X(end)),Gres(e,2))',Gres(e,1))',length(X),1);
                    Vb = reshape(downsample(downsample( VI(Y(1):Y(end),X(1):X(end)),Gres(e,2))',Gres(e,1))',length(X),1);
                end
                Eval=reshape(downsample(downsample( mask(Y(1):Y(end),X(1):X(end)),Gres(e,2))',Gres(e,1))',length(X),1);
                Eval(Eval==0)=-1;
                Eval(Eval>0)=0;
                
                %correlate image pair
                if (e~=1 || defloop~=1 || VelInputFile) && ~isempty(regexpi(M,'Deform','once'))          %then don't offset windows, images already deformed
                    %if Corr(e)<4
                    if ~strcmpi(Corr{e},'SPC')
                        [Xc,Yc,Uc,Vc,Cc,Dc,Cp]=PIVwindowed(im1d,im2d,Corr{e},Wsize(e,:),Wres(:, :, e),0,D(e,:),Zeromean(e),Peaklocator(e),Peakswitch(e) || (Valswitch(e) && extrapeaks(e)),frac_filt(e),saveplane(e),X(Eval>=0),Y(Eval>=0));
                        if strcmpi(M,'ForwardDeform')
                            % use coordinate system for deformed second
                            % frame to estimate deformation tensor and
                            % transform of subpixel corrector at t0 to
                            % deformed direction at t1.
                            
                            %                             %OPTION 1:
                            %                             % use the local deformation field in XD2 and
                            %                             % YD2 to estimate using central differences the
                            %                             % local strain rate tensor and use that to
                            %                             % transform Uc and Vc @ t0 to t1.
                            %                             %Question: what grid should the derivatives be
                            %                             % calculated over - at the pixel level, or the
                            %                             % vector level?
                            %
                            %                             % Need to figure out which XD2 and YD2, defined
                            %                             % over every pixel, correspond to the subset
                            %                             % embodied by Xc and Yc.  Don't forget to
                            %                             % make sure both use either pixel or grid
                            %                             % coordinate systems, not mixed!
                            %
                            %                             %loop over something similar to this for all
                            %                             %vectors positions Xc and Yc:
                            %
                            %                             %could substitute XD1 for XI here?
                            %                             Rest = [ (XD2(1,2)-XD2(1,1))/(XI(1,2)-XI(1,1)) , (XD2(2,1)-XD2(1,1))/(YI(2,1)-YI(1,1)) ; ...
                            %                                      (YD2(1,2)-YD2(1,1))/(XI(1,2)-XI(1,1)) , (YD2(2,1)-YD2(1,1))/(YI(2,1)-YI(1,1)) ];
                            %
                            %                             %this assumes Uc and Vc are column vectors
                            %                             UcVc = Rest*[Uc.';Vc.'];
                            %                             Uc = UcVc(1,:).';
                            %                             Vc = UcVc(2,:).';
                            
                            
                            %elseif strcmpi(M,'FowardDeformInterp')
                            % use coordinate system for deformed second
                            % frame to estimate deformation tensor and
                            % transform of subpixel corrector at t0 to
                            % deformed direction at t1.
                            
                            %OPTION 2:
                            % use interpolation to predict where a point
                            % shifted by (Uc,Vc) in (XD1,YD1) would fall in
                            % a deformed (XD2,YD2) system
                            
                            %first, figure out where vector centers are in (XD2,YD2)
                            XDc = interp2(XI,YI,XD2,Xc,Yc,'linear');
                            YDc = interp2(XI,YI,YD2,Xc,Yc,'linear');
                            %if shift takes us outside image domain, just
                            %return a NaN, we don't know where that point
                            %will go
                            if Peakswitch(e) || (Valswitch(e) && extrapeaks(e))
                                %there will be 3 velocity fields in Uc and Vc, so repmat vector origins
                                U2 = interp2(XI,YI,XD2,repmat(Xc,[1 3])+Uc,repmat(Yc,[1 3])+Vc,'linear',NaN) - repmat(XDc,[1,3]);
                                V2 = interp2(XI,YI,YD2,repmat(Xc,[1 3])+Uc,repmat(Yc,[1 3])+Vc,'linear',NaN) - repmat(YDc,[1,3]);
                            else
                                %only 1 velocity field is returned, use origins directly
                                U2 = interp2(XI,YI,YD2,Xc+Uc,Yc+Vc,'linear',NaN) - XDc;
                                V2 = interp2(XI,YI,YD2,Xc+Uc,Yc+Vc,'linear',NaN) - YDc;
                            end
                            
                            %if we have NaN values, we don't know how to
                            %correct the shift, so just use the raw value.
                            Uc(~isnan(U2)) = U2(~isnan(U2));
                            Vc(~isnan(V2)) = V2(~isnan(V2));
                            
                            clear U2 V2
                            
                            %Not sure what do with central difference
                            %deform correction yet - probably this section
                            %needs to be moved down to right before the
                            %deformation actually occurs?
                            %                         elseif strcmpi(M,'Deform')
                            %                             % use coordinate system for deformed 1st & 2nd
                            %                             % frames to estimate deformation tensor and
                            %                             % transform of subpixel corrector at t1/2 to
                            %                             % deformed direction at t0 and t1.
                            %
                            %                             %OPTION 2:
                            %                             % use interpolation to predict where a point
                            %                             % shifted by +/-(Uc/2,Vc/2) from (XI,YI) would fall in
                            %                             % a deformed (XD1,YD1) and (XD2,YD2) systems
                            %
                            %                             %first, figure out where vector centers are in (XD2,YD2)
                            %                             XDc1 = interp2(XI,YI,XD1,Xc,Yc,'linear');
                            %                             YDc1 = interp2(XI,YI,YD1,Xc,Yc,'linear');
                            %                             XDc2 = interp2(XI,YI,XD2,Xc,Yc,'linear');
                            %                             YDc2 = interp2(XI,YI,YD2,Xc,Yc,'linear');
                            %                             %if shift takes us outside image domain, just
                            %                             %return a NaN, we don't know where that point
                            %                             %will go
                            %                             if Peakswitch(e) || (Valswitch(e) && extrapeaks(e))
                            %                                 %there will be 3 velocity fields in Uc and Vc, so repmat vector origins
                            %                                 U1 = interp2(XI,YI,XD1,repmat(Xc,[1 3])+Uc,repmat(Yc,[1 3])+Vc,'linear',NaN) - repmat(XDc1,[1,3]);
                            %                                 V1 = interp2(XI,YI,YD1,repmat(Xc,[1 3])+Uc,repmat(Yc,[1 3])+Vc,'linear',NaN) - repmat(YDc1,[1,3]);
                            %                                 U2 = interp2(XI,YI,XD2,repmat(Xc,[1 3])+Uc,repmat(Yc,[1 3])+Vc,'linear',NaN) - repmat(XDc2,[1,3]);
                            %                                 V2 = interp2(XI,YI,YD2,repmat(Xc,[1 3])+Uc,repmat(Yc,[1 3])+Vc,'linear',NaN) - repmat(YDc2,[1,3]);
                            %                             else
                            %                                 %only 1 velocity field is returned, use origins directly
                            %                                 U1 = interp2(XI,YI,YD1,Xc+Uc,Yc+Vc,'linear',NaN) - XDc1;
                            %                                 V1 = interp2(XI,YI,YD1,Xc+Uc,Yc+Vc,'linear',NaN) - YDc1;
                            %                                 U2 = interp2(XI,YI,YD2,Xc+Uc,Yc+Vc,'linear',NaN) - XDc2;
                            %                                 V2 = interp2(XI,YI,YD2,Xc+Uc,Yc+Vc,'linear',NaN) - YDc2;
                            %                             end
                            %
                            %                             %if we have NaN values, we don't know how to
                            %                             %correct the shift, so just use the raw value.
                            %                             U1(isnan(U1)) = Uc(isnan(U1));
                            %                             V1(isnan(V1)) = Vc(isnan(V1));
                            %                             U2(isnan(U2)) = Uc(isnan(U2));
                            %                             V2(isnan(V2)) = Vc(isnan(V2));
                            %
                            %                             %clear U2 V2
                            
                        else
                            %other methods just add bulk offset predictor
                            %back to calculated subpixel corrector
                            % do nothing to Uc and Vc
                        end
                        
                        %reincorporate deformation as velocity for next pass
                        if Peakswitch(e) || (Valswitch(e) && extrapeaks(e))
                            %there will be 3 velocity fields in Uc an Vc, so repmat bulk offset
                            Uc = Uc + repmat(Ub(Eval>=0),[1 3]);
                            Vc = Vc + repmat(Vb(Eval>=0),[1 3]);
                        else
                            %only 1 velocity field is returned, use bulk offset directly
                            Uc = Uc + Ub(Eval>=0);
                            Vc = Vc + Vb(Eval>=0);
                        end
                    else
                        [Xc,Yc,Uc,Vc,Cc]=PIVphasecorr(im1d,im2d,Wsize(e,:),Wres(:, :, e),0,D(e,:),Zeromean(e),Peakswitch(e),X(Eval>=0),Y(Eval>=0));
                        %Sam deleted the Cc output from PIVPhaseCorr - why?  because we don't use it? But it's needed for Dc in next line?
                        %[Xc,Yc,Uc,Vc]=PIVphasecorr(im1d,im2d,Wsize(e,:),Wres(:, :, e),0,D(e),Zeromean(e),Peakswitch(e),X(Eval>=0),Y(Eval>=0));
                        Dc = zeros(size(Cc),imClass);
                        
                        Uc = Uc + Ub(Eval>=0);   %reincorporate deformation as velocity for next pass
                        Vc = Vc + Vb(Eval>=0);
                    end
                    
                else  %either first pass, or not deform
                    if ~strcmpi(Corr{e},'SPC')
                        if any(isnan(Ub(Eval>=0)))
                            keyboard
                        end
                        [Xc,Yc,Uc,Vc,Cc,Dc,Cp]=PIVwindowed(im1,im2,Corr{e},Wsize(e,:),Wres(:, :, e),0,D(e,:),Zeromean(e),Peaklocator(e),Peakswitch(e) || (Valswitch(e) && extrapeaks(e)),frac_filt(e),saveplane(e),X(Eval>=0),Y(Eval>=0),Ub(Eval>=0),Vb(Eval>=0));
                    else
                        [Xc,Yc,Uc,Vc,Cc]=PIVphasecorr(im1,im2,Wsize(e,:),Wres(:, :, e),0,D(e,:),Zeromean(e),Peakswitch(e),X(Eval>=0),Y(Eval>=0),Ub(Eval>=0),Vb(Eval>=0));
                        Dc = zeros(size(Cc),imClass);
                    end
                end
                
                %Where Eval<0, no correlation was performed and Uc, etc are
                %missing values.  Use Eval to fill in complete matrices U,V
                %over all grid points X,Y.
                if ~strcmpi(Corr{e},'SPC') %was not SPC=4
                    if Peakswitch(e) || (Valswitch(e) && extrapeaks(e))
                        U=zeros(size(X,1),3,imClass);
                        V=zeros(size(X,1),3,imClass);
                        U(repmat(Eval>=0,[1 3]))=Uc;V(repmat(Eval>=0,[1 3]))=Vc;
                        C=zeros(size(X,1),3,imClass);
                        Di=zeros(size(X,1),3,imClass);
                        C(repmat(Eval>=0,[1 3]))=Cc;
                        Di(repmat(Eval>=0,[1 3]))=Dc;
                    else
                        U=zeros(size(X),imClass);V=zeros(size(X),imClass);C=zeros(size(X),imClass);Di=zeros(size(X),imClass);
                        U(Eval>=0)=Uc;V(Eval>=0)=Vc;
                    end
                else %Corr was SPC=4
                    U=zeros(size(X),imClass);V=zeros(size(X),imClass);
                    U(Eval>=0)=Uc;V(Eval>=0)=Vc;
                    if Peakswitch(e)
                        C=zeros(size(X,1),3,imClass);
                        Di=zeros(size(X,1),3,imClass);
                        C(repmat(Eval>=0,[1 3]))=Cc;
                        Di(repmat(Eval>=0,[1 3]))=Dc;
                        
                    else
                        C=zeros(size(X),imClass);
                        Di=zeros(size(X),imClass);
                    end
                end
                
                corrtime(e,defloop)=toc(t1);
                
                %validation
                if Valswitch(e)
                    t1=tic;
                    
                    [Uval,Vval,Evalval,Cval,Dval]=VAL(X,Y,U,V,Eval,C,Di,Threshswitch(e),UODswitch(e),Bootswitch(e),extrapeaks(e),...
                        Uthresh(e,:),Vthresh(e,:),UODwinsize(e,:,:),UODthresh(e,UODthresh(e,:)~=0)',Bootper(e),Bootiter(e),Bootkmax(e));
                    
                    valtime(e,defloop)=toc(t1);
                else
                    Uval=U(:,1);Vval=V(:,1);Evalval=Eval(:,1);
                    if ~isempty(C)
                        Cval=C(:,1);
                        Dval=Di(:,1);
                    else
                        Cval=[];
                        Dval=[];
                    end
                end
                
                % --- Iterative Deformation Check ---
                % This block checks too see if the deformation has
                % converged or reach is max number of iterations.
                if ~isempty(regexpi(M,'Deform','once'))
                    if defloop == 1
                        Ud = Uval; Vd = Vval;
                    else
                        defconvU(e,defloop) = norm(Uval - Ud,2)./sqrt(numel(Ud));
                        defconvV(e,defloop) = norm(Vval - Vd,2)./sqrt(numel(Vd));
                        Ud = Uval; Vd = Vval;
                    end
                    if defloop == maxdefloop(e) || (defloop ~= 1 && defloop >= mindefloop(e) && defconvU(e,defloop) <= condefloop(e) && defconvV(e,defloop) <= condefloop(e))
                        if maxdefloop(e) ~= 1
                            % append the 'deform' and the pass number to
                            % the end of the file once the final number of
                            % iterations has been reached.
                            %wbase{e,:} = sprintf([wbase_org{e,:} 'deform' num2str(defloop) '_']);
                            wbase{e,:} = wbase_org{e,:};
                            numDefPasses(e) = defloop;
                        end
                        defloop = 1;
                    else
                        defloop = defloop+1;
                    end
                end
                
                
                %velocity smoothing
                % if Velsmoothswitch(e)==1
                %     U1=reshape(Uval(:,1),[S(1),S(2)]);
                %     V1=reshape(Vval(:,1),[S(1),S(2)]);
                %     [Us,Vs]=VELfilt(U1,V1,UODwinsize(e,:,:),Velsmoothfilt(e));
                %     Uval(:,1) = Us(:);
                %     Vval(:,1) = Vs(:);
                % end
                
                %write output
                if Writeswitch(e) && defloop == 1
                    t1=tic;
                    
                    %SPC only returns 1 peak right now?
                    if Peakswitch(e)
                        if PeakVel(e) && ~strcmpi(Corr{e},'SPC')
                            U=[Uval,U(:,1:PeakNum(e))];
                            V=[Vval,V(:,1:PeakNum(e))];
                        else
                            U=Uval; V=Vval;
                        end
                        if PeakMag(e)
                            C=[Cval,C(:,1:PeakNum(e))];
                            Di=[Dval,Di(:,1:PeakNum(e))];
                        else
                            C=Cval;
                            Di=Dval;
                        end
                    else
                        U=Uval; V=Vval; C=Cval; Di=Dval;
                    end
                    Eval=Evalval;
                    
                    %convert to physical units
                    Xval=X;Yval=Y;
                    X=X*Mag;Y=Y*Mag;
                    U=U*Mag/dt;V=V*Mag/dt;
                    
                    %convert to matrix if necessary
                    if size(X,2)==1
                        [X,Y,U,V,Eval,C,Di]=matrixform(X,Y,U,V,Eval,C,Di);
                    end
                    
                    %remove nans from data, replace with zeros
                    U(Eval<0|isinf(U))=0;V(Eval<0|isinf(V))=0;
                    
                    if str2double(Data.datout)
                        time=I1(q)/Freq;
                        if strcmpi(M,'Deform')
                            
                            % Here we are adding the information about the
                            % number of deformation passes into the title
                            % of the .dat file.
                            % First we convert the vector that contains the
                            % number of interations performed on each pass
                            % into a string
                            defPassString = num2str(numDefPasses(1:e)');
                            
                            % Then we find all of the characters that are
                            % not white space.
                            % tokens = regexp(defPassString,'([a-z_A-Z1-9])','tokenextents');
                            tokens = regexp(defPassString,'(\S)','tokenextents');
                            
                            % Next we reshape this cell into a matrix 2 by
                            % the number of locations.
                            tokens = cell2mat(tokens');
                            
                            % We will in the first part of our string with
                            % the first number of interations.
                            writeDefStr = defPassString(tokens(1,1):tokens(1,2));
                            
                            % Now we loop through the rest of the numbers
                            % adding a comma between numbers.
                            for sspaces = 2:size(tokens,1)
                                writeDefStr = [writeDefStr ',' defPassString(tokens(sspaces,1):tokens(sspaces,2))];
                            end
                            
                            % Now we append this information to the TITLE
                            % of the dat file (NOTE this doesn't change the
                            % file name only the title that is used in
                            % tecplot).
                            write_dat_val_C([pltdirec char(wbase(e,:)) sprintf(['%0.' Data.imzeros 'i.dat' ],I1(q))],X,Y,U,V,Eval,C,Di,e,time,char([wbase{e} ' Deform passes ' writeDefStr]));
                        else
                            write_dat_val_C([pltdirec char(wbase(e,:)) sprintf(['%0.' Data.imzeros 'i.dat' ],I1(q))],X,Y,U,V,Eval,C,Di,e,time,char(wbase(e,:)));
                        end
                    end
                    
                    if str2double(Data.multiplematout)
                        if ~isempty(regexpi(M,'Deform','once'))
                            save([pltdirec char(wbase(e,:)) sprintf(['%0.' Data.imzeros 'i.mat' ],I1(q))],'X','Y','U','V','Eval','C','Di','numDefPasses')
                        else
                            save([pltdirec char(wbase(e,:)) sprintf(['%0.' Data.imzeros 'i.mat' ],I1(q))],'X','Y','U','V','Eval','C','Di')
                        end
                    end
                    % This saves the correlation planes if that selection
                    % has been made in the job file.
                    if saveplane(e) && ~strcmpi(Corr{e},'SPC')
                        Xloc = Xc;Yloc=Yc;C_planes=Cp;%#ok
                        save(sprintf(['%s%scorrplanes_%0.' Data.imzeros 'i.mat' ],pltdirec,wbase{e,:},I1(q)),'Xloc','Yloc','C_planes')
                        clear Xloc Yloc C_planes
                    end
                        
                    X=Xval;Y=Yval;
                        
                      savetime(e,defloop)=toc(t1);
                end
                    
                    %reset any changes or scaling back to pixels
                    U=Uval; V=Vval;
                    
                    %If it isn't the lass pass, or the last pass hasn't
                    %converged yet for iterative deform, we need to prepare U
                    %and V from this pass for use as a predictor in the next
                    %pass
                    if e~=P || (~isempty(regexpi(M,'Deform','once')) && defloop ~=1)
                        %reshape from list of grid points to matrix
                        X=reshape(X,[S(1),S(2)]);
                        Y=reshape(Y,[S(1),S(2)]);
                        U=reshape(U(:,1),[S(1),S(2)]);
                        V=reshape(V(:,1),[S(1),S(2)]);
                        
                        %Multigrid and Deform need to interpolate this pass's displacements
                        %onto a different grid
                        if strcmpi(M,'Multigrid') || ~isempty(regexpi(M,'Deform','once'))
                            t1=tic;
                            
                            %velocity smoothing
                            if Velsmoothswitch(e)==1
                                [U,V]=VELfilt(U,V,UODwinsize(e,:,:),Velsmoothfilt(e));
                            end
                            
                            %velocity interpolation -
                            %resample U(X,Y) and V(X,Y) onto UI(XI,YI) and
                            %VI(XI,YI) where XI and YI are a list of every
                            %pixel in the image plane. Velinterp is the type of
                            %interpolation to use.
                            % We have previously made sure XI and YI correspond
                            % to the vector positions of the pixel centers by
                            % shifting by -0.5.
                            UI = VFinterp(X,Y,U,XI,YI,Velinterp);
                            VI = VFinterp(X,Y,V,XI,YI,Velinterp);
                            
                            if defloop == 1
                                interptime(e+1,defloop)=toc(t1);
                            else
                                interptime(e,defloop)=toc(t1);
                            end
                            
                            %For deform, UI and VI will be used to deform the images, and
                            %then downsampled in the next pass to the new grid
                            %for addition to the next pass's displacement.
                            %If not deform, then UI and VI just get downsampled
                            %next pass and used for DWO.
                            if strcmpi(M,'Deform')
                                t1=tic;
                                
                                % translate pixel locations, but
                                % since sincBlackmanInterp2 assumes
                                % coordinate system is pixel-centered, we need
                                % to convert back to index-coordinates for the
                                % deform.
                                %HOWEVER: we need to use these shifted
                                %coordinates later in vector coordinates, so
                                %move the -0.5 pixel correction to the call to
                                %the sincBlackmanInterp2 function
                                XD1 = XI-UI/2 ;
                                YD1 = YI-VI/2;
                                XD2 = XI+UI/2;
                                YD2 = YI+VI/2;
                                
                                % Preallocate memory for deformed images.
                                im1d = zeros(size(im1),imClass);
                                im2d = zeros(size(im2),imClass);
                                
                                % Deform images according to the interpolated velocity fields
                                for k = 1:nChannels % Loop over all of the color channels in the image
                                    if Iminterp == 1 % Sinc interpolation (without blackman window)
                                        im1d(:, :, k) = sincBlackmanInterp2(im1(:, :, k), XD1+0.5, YD1+0.5, 3, 'sinc');
                                        im2d(:, :, k) = sincBlackmanInterp2(im2(:, :, k), XD2+0.5, YD2+0.5, 3, 'sinc');
                                        
                                    elseif Iminterp == 2 % Sinc interpolation with blackman filter
                                        im1d(:, :, k) = sincBlackmanInterp2(im1(:, :, k), XD1+0.5, YD1+0.5, 3, 'blackman');
                                        im2d(:, :, k) = sincBlackmanInterp2(im2(:, :, k), XD2+0.5, YD2+0.5, 3, 'blackman');
                                    end
                                end
                                
                                % keyboard
                                % figure(1),imagesc(im1),colormap(gray),axis image xy,xlabel('im1')
                                % figure(2),imagesc(im2),colormap(gray),axis image xy,xlabel('im2')
                                % figure(3),imagesc(im1d),colormap(gray),axis image xy,xlabel('im1d')
                                % figure(4),imagesc(im2d),colormap(gray),axis image xy,xlabel('im2d')
                                % pause
                                % imwrite(uint8(im1d),[pltdirec char(wbase(e,:)) sprintf(['%0.' Data.imzeros 'ia.png' ],I1(q))]);
                                % imwrite(uint8(im2d),[pltdirec char(wbase(e,:)) sprintf(['%0.' Data.imzeros 'ib.png' ],I1(q))]);
                                
                                if defloop == 1
                                    deformtime(e+1,defloop)=toc(t1);
                                else
                                    deformtime(e,defloop)=toc(t1);
                                end
                            elseif strcmpi(M,'ForwardDeform') %only the 2nd image is deformed
                                t1=tic;
                                
                                % translate pixel locations, but
                                % since sincBlackmanInterp2 assumes
                                % coordinate system is pixel-centered, we need
                                % to convert back to index-coordinates for the
                                % deform.
                                %HOWEVER: we need to use these shifted
                                %coordinates later in vector coordinates, so
                                %move the -0.5 pixel correction to the call to
                                %the sincBlackmanInterp2 function
                                XD1 = XI   ;
                                YD1 = YI   ;
                                XD2 = XI+UI;
                                YD2 = YI+VI;
                                
                                % Preallocate memory for deformed images.
                                im1d = im1;
                                im2d = zeros(size(im2),imClass);
                                
                                % Deform images according to the interpolated velocity fields
                                for k = 1:nChannels % Loop over all of the color channels in the image
                                    if Iminterp == 1 % Sinc interpolation (without blackman window)
                                        im2d(:, :, k) = sincBlackmanInterp2(im2(:, :, k), XD2+0.5, YD2+0.5, 3, 'sinc');
                                    elseif Iminterp == 2 % Sinc interpolation with blackman filter
                                        im2d(:, :, k) = sincBlackmanInterp2(im2(:, :, k), XD2+0.5, YD2+0.5, 3, 'blackman');
                                    end
                                end
                                
                                % keyboard
                                % figure(1),imagesc(im1),colormap(gray),axis image xy,xlabel('im1')
                                % figure(2),imagesc(im2),colormap(gray),axis image xy,xlabel('im2')
                                % figure(3),imagesc(im1d),colormap(gray),axis image xy,xlabel('im1d')
                                % figure(4),imagesc(im2d),colormap(gray),axis image xy,xlabel('im2d')
                                % pause
                                % imwrite(uint8(im1d),[pltdirec char(wbase(e,:)) sprintf(['%0.' Data.imzeros 'ia.png' ],I1(q))]);
                                % imwrite(uint8(im2d),[pltdirec char(wbase(e,:)) sprintf(['%0.' Data.imzeros 'ib.png' ],I1(q))]);
                                
                                if defloop == 1
                                    deformtime(e+1,defloop)=toc(t1);
                                else
                                    deformtime(e,defloop)=toc(t1);
                                end
                                
                            end
                            
                            %save intermediate deformed images
                            if SaveIMdeform
                                %put the deformed images in a subdirectory of the vector directory to keep things tidy
                                saveIMDdir = fullfile(pltdirec,'imDeform');
                                [~,~,~] = mkdir(saveIMDdir); 
                                %build image names
                                saveIM1Dfile = [saveIMDdir, char(wbase(e,:)) ,'im1d_', sprintf(['%0.' Data.imzeros 'i.mat' ],I1(q))]; 
                                saveIM2Dfile = [saveIMDdir, char(wbase(e,:)) ,'im2d_', sprintf(['%0.' Data.imzeros 'i.mat' ],I1(q))]; %even though this is from I2, it corresponds to the vectors saved as I1, so use I1 in name
                                %check if images look like 8 bit or 12/16 bit
                                if max(im1d(:)) > 255 || max(im2d(:)) > 255
                                    saveclass = 'uint16';
                                else
                                    saveclass = 'uint8';
                                end
                                imwrite(cast(im1d,saveclass),saveIM1Dfile);
                                imwrite(cast(im2d,saveclass),saveIM2Dfile);
                             end
                            
                        else %must be Multipass - grid is same on every pass, no resampling needed
                            UI=U;VI=V;
                        end
                    end
                end
                
                eltime=toc(tf);
                %output text
                fprintf('\n----------------------------------------------------\n')
                fprintf(['Job: ',Data.batchname,'\n'])
                fprintf([frametitle ' Completed (' num2str(q) '/' num2str(length(I1)) ') at %s \n'], datestr(now));
                fprintf('----------------------------------------------------\n')
                for e=1:P
                    if e==1 && ~isempty(regexpi(M,'Deform','once')) && VelInputFile && mindefloop(e) == 1
                        %fprintf('velocity interpolation...        %0.2i:%0.2i.%0.0f\n',floor(sum(interptime(e,:))/60),floor(rem(sum(interptime(e,:)),60)),floor((rem(sum(interptime(e,:)),60)-floor(rem(sum(interptime(e,:)),60)))*10))
                        fprintf('image deformation...             %0.2i:%0.2i.%0.0f\n',floor(sum(deformtime(e,1))/60),floor(rem(sum(deformtime(e,1)),60)),floor((rem(sum(deformtime(e,1)),60)-floor(rem(sum(deformtime(e,1)),60)))*10))
                    end
                    
                    fprintf('correlation...                   %0.2i:%0.2i.%0.0f\n',floor(sum(corrtime(e,:))/60),floor(rem(sum(corrtime(e,:)),60)),floor((rem(sum(corrtime(e,:)),60)-floor(rem(sum(corrtime(e,:)),60)))*10))
                    if Valswitch(e)
                        fprintf('validation...                    %0.2i:%0.2i.%0.0f\n',floor(sum(valtime(e,:))/60),floor(rem(sum(valtime(e,:)),60)),floor((rem(sum(valtime(e,:)),60)-floor(rem(sum(valtime(e,:)),60)))*10))
                    end
                    if ~isempty(regexpi(M,'Deform','once')) && mindefloop(e) ~= 1
                        fprintf('velocity interpolation...        %0.2i:%0.2i.%0.0f\n',floor(sum(interptime(e,:))/60),floor(rem(sum(interptime(e,:)),60)),floor((rem(sum(interptime(e,:)),60)-floor(rem(sum(interptime(e,:)),60)))*10))
                        fprintf('image deformation...             %0.2i:%0.2i.%0.0f\n',floor(sum(deformtime(e,:))/60),floor(rem(sum(deformtime(e,:)),60)),floor((rem(sum(deformtime(e,:)),60)-floor(rem(sum(deformtime(e,:)),60)))*10))
                    end
                    if Writeswitch(e)
                        fprintf('save time...                     %0.2i:%0.2i.%0.0f\n',floor(savetime(e)/60),floor(rem(savetime(e),60)),floor((rem(savetime(e),60)-floor(rem(savetime(e),60)))*10))
                    end
                    if strcmp(M,'Multigrid') || (~isempty(regexpi(M,'Deform','once')) && mindefloop(e) == 1)
                        if e~=P
                            fprintf('velocity interpolation...        %0.2i:%0.2i.%0.0f\n',floor(sum(interptime(e,:))/60),floor(rem(sum(interptime(e,:)),60)),floor((rem(sum(interptime(e,:)),60)-floor(rem(sum(interptime(e,:)),60)))*10))
                            if ~isempty(regexpi(M,'Deform','once'))
                                fprintf('image deformation...             %0.2i:%0.2i.%0.0f\n',floor(sum(deformtime(e,:))/60),floor(rem(sum(deformtime(e,:)),60)),floor((rem(sum(deformtime(e,:)),60)-floor(rem(sum(deformtime(e,:)),60)))*10))
                            end
                        end
                    end
                end
                fprintf('total frame time...              %0.2i:%0.2i.%0.0f\n',floor(eltime/60),floor(rem(eltime,60)),floor((rem(eltime,60)-floor(rem(eltime,60)))*10))
                frametime(q)=eltime;
                comptime=mean(frametime(1:q))*(length(I1)-q);
                fprintf('estimated job completion time... %0.2i:%0.2i:%0.2i\n\n',floor(comptime/3600),floor(rem(comptime,3600)/60),floor(rem(comptime,60)))
            
        end
        
    case {'Ensemble','EDeform'}
        %% --- Ensemble and Ensemble Deform ---
        frametime=zeros(P,1);
        
        %initialize grid and evaluation matrix
        if strcmpi(Data.imext,'mat') %read .mat file, image must be stored in variable 'I'
            loaddata=load([imbase sprintf(['%0.' Data.imzeros 'i.' Data.imext],I1(1))]);
            im1 = cast(loaddata.I,imClass);
            loaddata =[];
        else
            im1=cast(imread([imbase sprintf(['%0.' Data.imzeros 'i.' Data.imext],I1(1))]),imClass);
        end
        imageSize=size(im1);
        imageSize(3)=size(im1,3);
        [XI,YI]=IMgrid(imageSize,[0 0]);
        %index-based pixel coordinates are integers on center, but we need XI and YI to
        %correspond to vector locations, which are centered on pixel edges for even sized
        %windows.  This means the center of pixel index 1 is actually at coordinate
        %0.5 from image edge, etc.
        % We need to do this because the vector positions we will use
        % for interpolation and deformation will be on the physical
        % grid, but will will need to know where the pixels are located
        % relative to them, not their indices.
        %BUT they are pixel centered for ODD sized windows.  (FIX) but
        %for now will ignore since even windows are more common.
        XI = XI - 0.5;
        YI = YI - 0.5;
        UI = Vel0.U;
        VI = Vel0.V;
        %UI = BWO(1)*ones(size(XI));
        %VI = BWO(2)*ones(size(XI));
        
        defconvU = zeros(P,max(maxdefloop));
        defconvV = zeros(P,max(maxdefloop));
        
        e = 0; defloop = 1;
        % This while statment is used to interatively move through the
        % deformations.  If the minimum number of loops hasn't been
        % reach it will keep iterating otherwise it should stop.
        while (e<P && defloop == 1) || (e<=P && defloop~=1)%for e=1:P
            if defloop == 1
                e=e+1;
                
                frametitle=['Frame' sprintf(['%0.' Data.imzeros 'i'],I1(1)) ' to Frame' sprintf(['%0.' Data.imzeros 'i'],I2(end))];
                fprintf('\n----------------------------------------------------\n')
                fprintf(['Job: ',Data.batchname,'\n'])
                fprintf([frametitle ' (Pass ' num2str(e) '/' num2str(P) ')\n'])
                fprintf('----------------------------------------------------\n')
                
            end
            tf=tic;
            
            [X,Y]=IMgrid(imageSize,Gres(e,:),Gbuf(e,:));
            S=size(X);X=X(:);Y=Y(:);
            Ub = reshape(downsample(downsample( UI(Y(1):Y(end),X(1):X(end)),Gres(e,2))',Gres(e,1))',length(X),1);
            Vb = reshape(downsample(downsample( VI(Y(1):Y(end),X(1):X(end)),Gres(e,2))',Gres(e,1))',length(X),1);
            Eval=reshape(downsample(downsample( mask(Y(1):Y(end),X(1):X(end)),Gres(e,2))',Gres(e,1))',length(X),1);
            Eval(Eval==0)=-1;
            Eval(Eval>0)=0;
            
            if Peakswitch(e) || (Valswitch(e) && extrapeaks(e))
                U=zeros(size(X,1),3,imClass);
                V=zeros(size(X,1),3,imClass);
                C=zeros(size(X,1),3,imClass);
                Di=zeros(size(X,1),3,imClass);
                DX=zeros(size(X,1),3,imClass);
                DY=zeros(size(X,1),3,imClass);
                ALPHA=zeros(size(X,1),3,imClass);
            else
                U=zeros(size(X),imClass);V=zeros(size(X),imClass);C=zeros(size(X),imClass);Di=zeros(size(X),imClass);
                DX=zeros(size(X),imClass);DY=zeros(size(X),imClass);ALPHA=zeros(size(X),imClass);
            end
            
            if str2double(Data.par) && matlabpool('size')>1
                
                spmd
                    verstr=version('-release');
                    if str2double(verstr(1:4))>=2010
                        I1dist=getLocalPart(codistributed(I1,codistributor('1d',2)));
                        I2dist=getLocalPart(codistributed(I2,codistributor('1d',2)));
                    else
                        I1dist=localPart(codistributed(I1,codistributor('1d',2),'convert'));
                        I2dist=localPart(codistributed(I2,codistributor('1d',2),'convert'));
                    end
                    
                    for q=1:length(I1dist)
                        
                        %load image pair and flip coordinates
                        if strcmpi(Data.imext,'mat') %read .mat file, image must be stored in variable 'I'
                            loaddata=load([imbase sprintf(['%0.' Data.imzeros 'i.' Data.imext],I1dist(q))]);
                            im1 = cast(loaddata.I,imClass);
                            loaddata=load([imbase2 sprintf(['%0.' Data.imzeros 'i.' Data.imext],I2dist(q))]);
                            im2 = cast(loaddata.I,imClass);
                            loaddata =[];
                        else
                            im1=cast(imread([imbase sprintf(['%0.' Data.imzeros 'i.' Data.imext],I1dist(q))]),imClass);
                            im2=cast(imread([imbase2 sprintf(['%0.' Data.imzeros 'i.' Data.imext],I2dist(q))]),imClass);
                        end
                        
                        if size(im1, 3) > 2
                            %Extract only red channel
                            if channel == 1;
                                im1 = im1(:,:,1);
                                im2 = im2(:,:,1);
                                %Extract only green channel
                            elseif channel == 2;
                                im1 = im1(:,:,2);
                                im2 = im2(:,:,2);
                                %Extract only blue channel
                            elseif channel == 3;
                                im1 = im1(:,:,3);
                                im2 = im2(:,:,3);
                                %Weighted average of channels (see rgb2gray for
                                %explanation of weighting factors)
                            elseif channel == 4;
                                im1 = 0.2989 * im1(:, :, 1) + 0.5870 * im1(:, :, 2) + 0.1140 * im1(:, :, 3);
                                im2 = 0.2989 * im2(:, :, 1) + 0.5870 * im2(:, :, 2) + 0.1140 * im2(:, :, 3);
                                %Evenly weighted mean of channels
                            elseif channel == 5;
                                im1 = (im1(:,:,1) + im1(:,:,2) + im1(:,:,3))/3;
                                im2 = (im2(:,:,1) + im2(:,:,2) + im2(:,:,3))/3;
                                %ensemble correlation of channels
                            elseif channel == 6;
                                im1=im1(:,:,1:3);
                                im2=im2(:,:,1:3);
                            end
                        else
                            %Take only red channel
                            im1 =im1(:,:,1);
                            im2 =im2(:,:,1);
                            channel = 1;
                        end
                        
                        % Determine the number of channels in the image
                        % to be deformed. This should be 3 for color
                        % images or 1 for grayscale images.
                        nChannels = size(im1, 3);
                        
                        %  Flip images
                        %flipud only works on 2D matices.
                        % im1=flipud(im1(:,:,1));
                        % im2=flipud(im2(:,:,1));
                        im1 = im1(end:-1:1,:,:);
                        im2 = im2(end:-1:1,:,:);
                        % L=size(im1);
                        
                        % The deformation for ensemble must be done before
                        % the correlation unlike in the instantaneous
                        % images where it is done after correlation
                        if strcmpi(M,'EDeform') && (e~=1 || defloop ~=1 || VelInputFile)
                            
                            t1=tic;
                            % translate pixel locations, but
                            % since sincBlackmanInterp2 assumes
                            % coordinate system is pixel-centered, we need
                            % to convert back to index-coordinates for the
                            % deform.
                            % Here this seems to be redundant ... Where is
                            % VFInterp called?
                            XD1 = XI-UI/2 + 0.5;
                            YD1 = YI-VI/2 + 0.5;
                            XD2 = XI+UI/2 + 0.5;
                            YD2 = YI+VI/2 + 0.5;
                            
                            % Preallocate memory for deformed images.
                            im1d = zeros(size(im1),imClass);
                            im2d = zeros(size(im2),imClass);
                            
                            % Deform images according to the interpolated
                            % velocity fields. Each color channel is
                            % deformed separately. They could have been
                            % done all together but I'm trying to save
                            % memory at the expense of some speed here.
                            for k = 1:nChannels % Loop over all of the color channels in the image
                                if Iminterp == 1 % Sinc interpolation (without blackman window)
                                    im1d(:, :, k) = sincBlackmanInterp2(im1(:, :, k), XD1, YD1, 3, 'sinc');
                                    im2d(:, :, k) = sincBlackmanInterp2(im2(:, :, k), XD2, YD2, 3, 'sinc');
                                    
                                elseif Iminterp == 2 % Sinc interpolation with blackman filter
                                    im1d(:, :, k) = sincBlackmanInterp2(im1(:, :, k), XD1, YD1, 3, 'blackman');
                                    im2d(:, :, k) = sincBlackmanInterp2(im2(:, :, k), XD2, YD2, 3, 'blackman');
                                end
                            end
                            
                            deformtime=toc(t1);
                        end
                        
                        t1=tic;
                        %correlate image pair and average correlations
                        %                      [Xc,Yc,CC]=PIVensemble(im1,im2,Corr(e),Wsize(e,:),Wres(e, :, :),0,D(e),Zeromean(e),X(Eval>=0),Y(Eval>=0),Ub(Eval>=0),Vb(Eval>=0));
                        if strcmpi(M,'EDeform') && (e~=1 || defloop ~=1  || VelInputFile)
                            [Xc,Yc,CC]=PIVensemble(im1d,im2d,Corr{e},Wsize(e,:),Wres(:, :, e),0,D(e,:),Zeromean(e),frac_filt(e),X(Eval>=0),Y(Eval>=0));
                        else
                            [Xc,Yc,CC]=PIVensemble(im1,im2,Corr{e},Wsize(e,:),Wres(:, :, e),0,D(e,:),Zeromean(e),frac_filt(e),X(Eval>=0),Y(Eval>=0),Ub(Eval>=0),Vb(Eval>=0));
                        end
                        
                        if ~strcmpi(Corr{e},'SPC')
                            if q==1
                                CCmdist=CC;
                                %cnvg_est = 0;
                                CC = []; %#ok% This clear is required for fine grids or big windows
                            else
                                % % cnvg_est = norm((CCmdist(:)*length(I1)/(q-1))-((CCmdist(:)*length(I1)+CC(:))/q),2);
                                % ave_pre = (CCmdist/(q-1));
                                CCmdist=CCmdist+CC;% Now includes the current frame
                                % ave_cur = ((CCmdist)/q);
                                % ave_cur(ave_cur==0)=nan; %This makes sure you don't divide by zeros
                                %cnvg_est = 0;%nanmean(mean(mean(abs(ave_pre-ave_cur),1),2)./nanmean(nanmean(abs(ave_cur),1),2));
                                CC = []; %#ok% This clear is required for fine grids or big windows
                            end
                        else %if Corr(e)==4 %SPC processor
                            error('SPC Ensemble does not work with parallel processing. Try running again on a single core.')
                        end
                        corrtime=toc(t1);
                        if strcmpi(M,'EDeform') && (e~=1 || defloop~=1 || VelInputFile)
                            fprintf('deformation %4.0f of %4.0f...      %0.2i:%0.2i.%0.0f\n',q,length(I1dist),floor(deformtime/60),floor(rem(deformtime,60)),floor((rem(deformtime,60)-floor(rem(deformtime,60)))*10))
                        end
                        %                         fprintf('correlation %4.0f of %4.0f...      %0.2i:%0.2i.%0.0f Ensemble %%change %0.2e\n',q,length(I1dist),floor(corrtime/60),floor(rem(corrtime,60)),rem(corrtime,60)-floor(rem(corrtime,60)),cnvg_est)
                        fprintf('correlation %4.0f of %4.0f...      %0.2i:%0.2i.%0.0f\n',q,length(I1dist),floor(corrtime/60),floor(rem(corrtime,60)),floor((rem(corrtime,60)-floor(rem(corrtime,60)))*10))
                    end
                end
                %                 if Corr(e)<4 %SCC or RPC processor
                CCm=zeros(size(CCmdist{1}),imClass);
                for i=1:length(CCmdist)
                    CCm=CCm+CCmdist{i}/length(I1);
                end
                %Xloc and Yloc need to be extracted from the distributed
                %variable in case we want to use them outside the spmd loop
                Xc = Xc{1};
                Yc = Yc{1};
                
                %                 elseif Corr(e)==2 %SPC processor
                %                     CCm=zeros(size(CCmdist{1},1),size(CCmdist{1},2),size(CCmdist{1},3),length(I1));
                %                     ind=1;
                %                     for i=1:length(CCmdist)
                %                         CCm(:,:,:,ind:ind+size(CCmdist{i},4)-1)=CCmdist{i};
                %                         ind=ind+size(CCmdist{i},4);
                %                     end
                %                 end
            else %Not running in parallel
                
                for q=1:length(I1)
                    
                    %load image pair and flip coordinates
                    if strcmpi(Data.imext,'mat') %read .mat file, image must be stored in variable 'I'
                        loaddata=load([imbase sprintf(['%0.' Data.imzeros 'i.' Data.imext],I1(q))]);
                        im1 = cast(loaddata.I,imClass);
                        loaddata=load([imbase2 sprintf(['%0.' Data.imzeros 'i.' Data.imext],I2(q))]);
                        im2 = cast(loaddata.I,imClass);
                        loaddata =[];
                    else
                        im1=cast(imread([imbase sprintf(['%0.' Data.imzeros 'i.' Data.imext],I1(q))]),imClass);
                        im2=cast(imread([imbase2 sprintf(['%0.' Data.imzeros 'i.' Data.imext],I2(q))]),imClass);
                    end
                    
                    if size(im1, 3) > 2
                        %Extract only red channel
                        if channel == 1;
                            im1 = im1(:,:,1);
                            im2 = im2(:,:,1);
                            %Extract only green channel
                        elseif channel == 2;
                            im1 = im1(:,:,2);
                            im2 = im2(:,:,2);
                            %Extract only blue channel
                        elseif channel == 3;
                            im1 = im1(:,:,3);
                            im2 = im2(:,:,3);
                            %Weighted average of channels (see rgb2gray for
                            %explanation of weighting factors)
                        elseif channel == 4;
                            im1 = 0.2989 * im1(:, :, 1) + 0.5870 * im1(:, :, 2) + 0.1140 * im1(:, :, 3);
                            im2 = 0.2989 * im2(:, :, 1) + 0.5870 * im2(:, :, 2) + 0.1140 * im2(:, :, 3);
                            %Evenly weighted mean of channels
                        elseif channel == 5;
                            im1 = (im1(:,:,1) + im1(:,:,2) + im1(:,:,3))/3;
                            im2 = (im2(:,:,1) + im2(:,:,2) + im2(:,:,3))/3;
                            %ensemble correlation of channels
                        elseif channel == 6;
                            im1=im1(:,:,1:3);
                            im2=im2(:,:,1:3);
                        end
                    else
                        %Take only red channel
                        im1 =im1(:,:,1);
                        im2 =im2(:,:,1);
                        channel = 1;
                    end
                    
                    % Determine the number of channels in the image
                    % to be deformed. This should be 3 for color
                    % images or 1 for grayscale images.
                    nChannels = size(im1, 3);
                    
                    %  Flip images.
                    im1 = im1(end:-1:1,:,:);
                    im2 = im2(end:-1:1,:,:);
                    
                    if strcmpi(M,'EDeform') && (e~=1 || defloop ~=1 || VelInputFile)
                        t1=tic;
                        
                        % translate pixel locations, but
                        % since sincBlackmanInterp2 assumes
                        % coordinate system is pixel-centered, we need
                        % to convert back to index-coordinates for the
                        % deform.
                        % Here this seems to be redundant ... Where is
                        % VFInterp called?
                        XD1 = XI-UI/2 + 0.5;
                        YD1 = YI-VI/2 + 0.5;
                        XD2 = XI+UI/2 + 0.5;
                        YD2 = YI+VI/2 + 0.5;
                        
                        % Preallocate memory for deformed images.
                        im1d = zeros(size(im1),imClass);
                        im2d = zeros(size(im2),imClass);
                        
                        % Deform images according to the interpolated velocity fields. Each color
                        % channel is deformed separately. They could have been done all together
                        % but I'm trying to save memory at the expense of some speed here.
                        for k = 1:nChannels % Loop over all of the color channels in the image
                            if Iminterp == 1 % Sinc interpolation (without blackman window)
                                im1d(:, :, k) = sincBlackmanInterp2(im1(:, :, k), XD1, YD1, 3, 'sinc');
                                im2d(:, :, k) = sincBlackmanInterp2(im2(:, :, k), XD2, YD2, 3, 'sinc');
                                
                            elseif Iminterp == 2 % Sinc interpolation with blackman filter
                                im1d(:, :, k) = sincBlackmanInterp2(im1(:, :, k), XD1, YD1, 3, 'blackman');
                                im2d(:, :, k) = sincBlackmanInterp2(im2(:, :, k), XD2, YD2, 3, 'blackman');
                            end
                        end
                        
                        deformtime=toc(t1);
                    end
                    
                    t1=tic;
                    %correlate image pair and average correlations
                    %                   [Xc,Yc,CC]=PIVensemble(im1,im2,Corr(e),Wsize(e,:),Wres(e, :, :),0,D(e),Zeromean(e),X(Eval>=0),Y(Eval>=0),Ub(Eval>=0),Vb(Eval>=0));
                    if strcmpi(M,'EDeform') && (e~=1 || defloop ~=1 || VelInputFile)
                        [Xc,Yc,CC]=PIVensemble(im1d,im2d,Corr{e},Wsize(e,:),Wres(:, :, e),0,D(e,:),Zeromean(e),frac_filt(e),X(Eval>=0),Y(Eval>=0));
                    else
                        [Xc,Yc,CC]=PIVensemble(im1,im2,Corr{e},Wsize(e,:),Wres(:, :, e),0,D(e,:),Zeromean(e),frac_filt(e),X(Eval>=0),Y(Eval>=0),Ub(Eval>=0),Vb(Eval>=0));
                    end
                    
                    if ~strcmpi(Corr{e},'SPC')   %SPC=4 %SCC or RPC
                        if q==1
                            CCm=CC/length(I1);
                            %cnvg_est = 0;
                            CC = []; %#ok% This clear is required for fine grids or big windows
                        else
                            % cnvg_est = norm((CCm(:)*length(I1)/(q-1))-((CCm(:)*length(I1)+CC(:))/q),2);
                            % ave_pre = (CCm*length(I1)/(q-1));
                            CCm=CCm+CC/length(I1);% Now adding the current pass
                            % ave_cur = ((CCm*length(I1))/q);
                            % ave_cur(ave_cur==0)=nan; %This makes sure you don't divide by zero.
                            %cnvg_est = 0;%nanmean(mean(mean(abs(ave_pre-ave_cur),1),2)./nanmean(nanmean(abs(ave_cur),1),2));
                            CC = []; %#ok% This clear is required for fine grids or big windows
                        end
                    else %if Corr(e)==4 %SPC processor, should this be just ELSE?
                        if q==1
                            CCm=CC;
                        else
                            CCm.U=[CCm.U,CC.U];
                            CCm.V=[CCm.V,CC.V];
                            CCm.C=[CCm.C,CC.C];
                        end
                        %cnvg_est = 0;
                    end
                    corrtime=toc(t1);
                    if strcmpi(M,'EDeform') && (e~=1 || defloop~=1 || VelInputFile)
                        fprintf('deformation %4.0f of %4.0f...      %0.2i:%0.2i.%0.0f\n',q,length(I1),floor(deformtime/60),floor(rem(deformtime,60)),floor((rem(deformtime,60)-floor(rem(deformtime,60)))*10))
                    end
                    %                     fprintf('correlation %4.0f of %4.0f...      %0.2i:%0.2i.%0.0f Ensemble %%change %0.2e\n',q,length(I1),floor(corrtime/60),floor(rem(corrtime,60)),rem(corrtime,60)-floor(rem(corrtime,60)),cnvg_est)
                    fprintf('correlation %4.0f of %4.0f...      %0.2i:%0.2i.%0.0f\n',q,length(I1),floor(corrtime/60),floor(rem(corrtime,60)),floor((rem(corrtime,60)-floor(rem(corrtime,60)))*10))
                end
            end
            
            Z=[size(CCm,1), size(CCm,2),length(X(Eval>=0))];
            cnorm=ones(Z(1),Z(2),imClass);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %{
            %fftshift indicies
            fftindy = [ceil(Wsize(e,2)/2)+1:Wsize(e,2) 1:ceil(Wsize(e,2)/2)];
            fftindx = [ceil(Wsize(e,1)/2)+1:Wsize(e,1) 1:ceil(Wsize(e,1)/2)];

            %window masking filter
            sfilt1 = windowmask([Wsize(e,1) Wsize(e,2)],[Wres(1, 1,e) Wres(1, 2,e)]);
            sfilt2 = windowmask([Wsize(e,1) Wsize(e,2)],[Wres(2, 1,e) Wres(2, 2,e)]);
            s1   = fftn(sfilt1,[Wsize(e,2) Wsize(e,1)]);
            s2   = fftn(sfilt2,[Wsize(e,2) Wsize(e,1)]);
            S21  = s2.*conj(s1);
            
            %Standard Fourier Based Cross-Correlation
            iS21 = ifftn(S21,'symmetric');
            iS21 = iS21(fftindy,fftindx);
            cnorm = 1./iS21;
            cnorm(isinf(cnorm)) = 0;
            %}
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            if Peakswitch(e) || (Valswitch(e) && extrapeaks(e))
                Uc=zeros(Z(3),3,imClass);
                Vc=zeros(Z(3),3,imClass);
                Cc=zeros(Z(3),3,imClass);
                Dc=zeros(Z(3),3,imClass);
                DXc=zeros(Z(3),3,imClass);
                DYc=zeros(Z(3),3,imClass);
                ALPHAc=zeros(Z(3),3,imClass);
                Ub=repmat(Ub,[1 3]);
                Vb=repmat(Vb,[1 3]);
                Eval=repmat(Eval,[1 3]);
            else
                Uc=zeros(Z(3),1,imClass);Vc=zeros(Z(3),1,imClass);Cc=zeros(Z(3),1,imClass);Dc=zeros(Z(3),1,imClass);
                DXc=zeros(Z(3),1,imClass);DYc=zeros(Z(3),1,imClass);ALPHAc=zeros(Z(3),1,imClass);
                
            end
            
            if ~strcmpi(Corr{e},'SPC')
                t1=tic;
                for s=1:Z(3) %Loop through grid points
                    %Find the subpixel fit of the average correlation matrix
                    [Uc(s,:),Vc(s,:),Cc(s,:),Dc(s,:),DXc(s,:),DYc(s,:),ALPHAc(s,:)]=subpixel(CCm(:,:,s),Z(2),Z(1),cnorm,Peaklocator(e),Peakswitch(e) || (Valswitch(e) && extrapeaks(e)),D(e,:));
                end
                peaktime=toc(t1);
                fprintf('peak fitting...                  %0.2i:%0.2i.%0.0f\n',floor(peaktime/60),floor(rem(peaktime,60)),floor((rem(peaktime,60)-floor(rem(peaktime,60)))*10))
            elseif strcmpi(Corr{e},'SPC') %SPC processor
                %RPC filter for weighting function
                wt = energyfilt(Z(2),Z(1),D(e,:),0);
                wtX=wt(Z(1)/2+1,:)';
                cutoff=2/pi/D(e,2);
                wtX(wtX<cutoff)=0;
                wtY=wt(:,Z(2)/2+1);
                cutoff=2/pi/D(e,1);
                wtY(wtY<cutoff)=0;
                lsqX=(0:Z(2)-1)-Z(2)/2;
                lsqY=(0:Z(1)-1)-Z(1)/2;
                lsqX=repmat(lsqX,[1 length(I1)]);
                lsqY=repmat(lsqY,[1 length(I1)]);
                Qp=squeeze(CCm.C(1,:,:)./CCm.C(2,:,:));
                Qp_norm=(Qp-repmat(min(Qp),[length(I1) 1]))./repmat(max(Qp)-min(Qp),[length(I1) 1]);
                
                for s=1:Z(3)
                    wtX_cum=reshape(repmat(wtX,[1 length(I1)]),[1 numel(wtX)*length(I1)]);
                    wtY_cum=reshape(repmat(wtY,[1 length(I1)]),[1 numel(wtY)*length(I1)]);
                    
                    for q=1:length(I1)
                        indX=(1:Z(2))+Z(2)*(q-1);
                        indY=(1:Z(1))+Z(1)*(q-1);
                        wtX_cum(indX)=wtX_cum(indX).*Qp_norm(q,s);
                        wtY_cum(indY)=wtX_cum(indY).*Qp_norm(q,s);
                    end
                    
                    %Perform the weighted lsq regression
                    Uc(s)= wlsq(CCm.U(1,:,s),lsqX,wtX_cum)*Z(2)/2/pi;
                    Vc(s)=-wlsq(CCm.V(1,:,s),lsqY,wtY_cum)*Z(1)/2/pi;
                    Cc(s)=max(max(CCm.C(:,:,s)));
                    Dc(s)=0;
                end
            end
            
            if strcmpi(M,'EDeform') && (e~=1 || defloop ~=1 || VelInputFile)
                U(Eval>=0)=Uc(:)+Ub(Eval>=0);
                V(Eval>=0)=Vc(:)+Vb(Eval>=0);
            else
                U(Eval>=0)=Uc(:)+round(Ub(Eval>=0));
                V(Eval>=0)=Vc(:)+round(Vb(Eval>=0));
            end
            if Peakswitch(e) || (Valswitch(e) && extrapeaks(e))%~isempty(Cc)
                C(Eval>=0)=Cc(:);
                Di(Eval>=0)=Dc(:);
                DX(Eval>=0)=DXc(:);
                DY(Eval>=0)=DYc(:);
                ALPHA(Eval>=0)=ALPHAc(:);
            end
            
            %validation
            if Valswitch(e)
                %keyboard
                t1=tic;
                
                [Uval,Vval,Evalval,Cval,Dval,DXval,DYval,ALPHAval]=VAL(X,Y,U,V,Eval,C,Di,Threshswitch(e),UODswitch(e),Bootswitch(e),extrapeaks(e),...
                    Uthresh(e,:),Vthresh(e,:),UODwinsize(e,:,:),UODthresh(e,UODthresh(e,:)~=0)',Bootper(e),Bootiter(e),Bootkmax(e),...
                    DX,DY,ALPHA);
                
                valtime=toc(t1);
                fprintf('validation...                    %0.2i:%0.2i.%0.0f\n',floor(valtime/60),floor(rem(valtime,60)),floor((rem(valtime,60)-floor(rem(valtime,60)))*10))
                
            else
                Uval=U(:,1);Vval=V(:,1);Evalval=Eval(:,1);
                if ~isempty(C)
                    Cval=C(:,1);
                    Dval=Di(:,1);
                else
                    Cval=[];
                    Dval=[];
                end
                DXval=DX(:,1);DYval=DY(:,1);ALPHAval=ALPHA(:,1);
            end
            
            % --- Iterative Deformation Check ---
            if strcmpi(M,'EDeform')
                if defloop == 1
                    Ud = Uval; Vd = Vval;
                else
                    defconvU(e,defloop) = norm(Uval - Ud,2)./sqrt(numel(Ud));
                    defconvV(e,defloop) = norm(Vval - Vd,2)./sqrt(numel(Ud));
                    Ud = Uval; Vd = Vval;
                end
                if defloop == maxdefloop(e) || (defloop ~= 1 && defloop >= mindefloop(e) && defconvU(e,defloop) <= condefloop(e) && defconvV(e,defloop) <= condefloop(e))
                    if maxdefloop(e) ~= 1
                        %wbase{e,:} = sprintf([wbase_org{e,:} 'deform' num2str(defloop) '_']);
                        wbase{e,:} = wbase_org{e,:};
                        numDefPasses(e) = defloop;
                    end
                    defloop = 1;
                else
                    defloop = defloop+1;
                end
            end
            
            %velocity smoothing
            % if Velsmoothswitch(e)==1
            %     U1=reshape(Uval(:,1),[S(1),S(2)]);
            %     V1=reshape(Vval(:,1),[S(1),S(2)]);
            %     [Us,Vs]=VELfilt(U1,V1,UODwinsize(e,:,:),Velsmoothfilt(e));
            %     Uval(:,1) = Us(:);
            %     Vval(:,1) = Vs(:);
            % end
            
            
            %write output
            if Writeswitch(e) && defloop == 1
                t1=tic;
                
                if Peakswitch(e)
                    if PeakVel(e) && ~strcmpi(Corr{e},'SPC')
                        U=[Uval,U(:,1:PeakNum(e))];
                        V=[Vval,V(:,1:PeakNum(e))];
                        Eval=[Evalval,Eval(:,1:PeakNum(e))];
                    else
                        U=Uval; V=Vval; Eval=Evalval;
                    end
                    if PeakMag(e)
                        C=[Cval,C(:,1:PeakNum(e))];
                        Di=[Dval,Di(:,1:PeakNum(e))];
                        DX=[DXval,DX(:,1:PeakNum(e))];
                        DY=[DYval,DY(:,1:PeakNum(e))];
                        ALPHA=[ALPHAval,ALPHA(:,1:PeakNum(e))];
                    else
                        C=Cval;
                        Di=Dval;
                        DX=DXval;
                        DY=DYval;
                        ALPHA=ALPHAval;
                    end
                else
                    U=Uval; V=Vval; Eval=Evalval; C=Cval; Di=Dval;
                    DX=DXval;DY=DYval;ALPHA=ALPHAval;
                end
                
                %convert to physical units
                Xval=X;Yval=Y;
                X=X*Mag;Y=Y*Mag;
                U=U*Mag/dt;V=V*Mag/dt;
                
                %convert to matrix if necessary
                if size(X,2)==1
                    [~,~,~,~,DX,DY,ALPHA]=matrixform(X,Y,U,V,DX,DY,ALPHA);
                    [X,Y,U,V,Eval,C,Di]=matrixform(X,Y,U,V,Eval,C,Di);
                end
                
                %remove nans from data, replace with zeros
                U(Eval<0 | isinf(U))=0;V(Eval<0 | isinf(V))=0;
                
                if str2double(Data.datout)
                    time=I1(1)/Freq;
                    if strcmpi(M,'EDeform')
                        
                        % Here we are adding the information about the
                        % number of deformation passes into the title
                        % of the .dat file.
                        % First we convert the vector that contains the
                        % number of interations performed on each pass
                        % into a string
                        defPassString = num2str(numDefPasses(1:e)');
                        
                        % Then we find all of the characters that are
                        % not white space.
                        % tokens = regexp(defPassString,'([a-z_A-Z1-9])','tokenextents');
                        tokens = regexp(defPassString,'(\S)','tokenextents');
                        
                        % Next we reshape this cell into a matrix 2 by
                        % the number of locations.
                        tokens = cell2mat(tokens');
                        
                        % We will in the first part of our string with
                        % the first number of interations.
                        writeDefStr = defPassString(tokens(1,1):tokens(1,2));
                        
                        % Now we loop through the rest of the numbers
                        % adding a comma between numbers.
                        for sspaces = 2:size(tokens,1)
                            writeDefStr = [writeDefStr ',' defPassString(tokens(sspaces,1):tokens(sspaces,2))];
                        end
                        
                        % Now we append this information to the TITLE
                        % of the dat file (NOTE this doesn't change the
                        % file name only the title that is used in
                        % tecplot).
                        write_dat_val_C([pltdirec char(wbase(e,:)) sprintf(['%0.' Data.imzeros 'i.dat' ],I1(1))],X,Y,U,V,Eval,C,Di,e,time,char([wbase{e} ' Deform passes ' writeDefStr]));
                    else
                        %write_dat_val_C([pltdirec char(wbase(e,:)) sprintf(['%0.' Data.imzeros 'i.dat' ],I1(1))],X,Y,U,V,Eval,C,e,0,frametitle);
                        write_dat_val_C([pltdirec char(wbase(e,:)) sprintf(['%0.' Data.imzeros 'i.dat' ],I1(1))],X,Y,U,V,Eval,C,Di,e,time,char(wbase(e,:)));
                    end
                end
                if str2double(Data.multiplematout)
                    if strcmpi(M,'EDeform')
                        save([pltdirec char(wbase(e,:)) sprintf(['%0.' Data.imzeros 'i.mat' ],I1(1))],'X','Y','U','V','Eval','C','Di','numDefPasses','DX','DY','ALPHA')
                    else
                        save([pltdirec char(wbase(e,:)) sprintf(['%0.' Data.imzeros 'i.mat' ],I1(1))],'X','Y','U','V','Eval','C','Di','DX','DY','ALPHA')
                    end
                end
                if saveplane(e) && ~strcmpi(Corr{e},'SPC')
                    Xloc = Xc;Yloc=Yc;%#ok
                    save(sprintf(['%s%scorrplanes_%0.' Data.imzeros 'i.mat' ],pltdirec,wbase{e,:},I1(1)),'Xloc','Yloc','CCm')
                    clear Xloc Yloc
                end
                X=Xval;Y=Yval;
                
                savetime=toc(t1);
                fprintf('save time...                     %0.2i:%0.2i.%0.0f\n',floor(savetime/60),floor(rem(savetime,60)),floor((rem(savetime,60)-floor(rem(savetime,60)))*10))
            end
            U=Uval; V=Vval;
            
            if e~=P || defloop ~= 1 %Not the last pass or not finished converging the final pass
                t1=tic;
                
                %reshape from list of grid points to matrix
                X=reshape(X,[S(1),S(2)]);
                Y=reshape(Y,[S(1),S(2)]);
                U=reshape(U(:,1),[S(1),S(2)]);
                V=reshape(V(:,1),[S(1),S(2)]);
                
                %velocity smoothing
                if Velsmoothswitch(e)==1
                    [U,V]=VELfilt(U,V,UODwinsize(e,:,:),Velsmoothfilt(e));
                end
                
                %velocity interpolation
                %resample U, V and W from vector grid coordinates
                %V(X,Y,Z) onto UI,VI, and WI on pixel coordinates
                %[XI,YI] where XI and YI are a list of every
                %pixel centers in the image plane. We adjusted XI and YI
                %previously.
                %Velinterp is the type of interpolation to use.
                UI = VFinterp(X,Y,U,XI,YI,Velinterp);
                VI = VFinterp(X,Y,V,XI,YI,Velinterp);
                
                interptime=toc(t1);
                fprintf('velocity interpolation...        %0.2i:%0.2i.%0.0f\n',floor(interptime/60),floor(rem(interptime,60)),floor((rem(interptime,60)-floor(rem(interptime,60)))*10))
            end
            
            if defloop == 1
                eltime=toc(tf);
                %output text
                fprintf('total pass time...               %0.2i:%0.2i.%0.0f\n',floor(eltime/60),floor(rem(eltime,60)),floor((rem(eltime,60)-floor(rem(eltime,60)))*10))
                frametime(e)=eltime;
                comptime=mean(frametime(1:e))*(P-e);
                fprintf('estimated job completion time... %0.2i:%0.2i:%0.2i\n',floor(comptime/3600),floor(rem(comptime,3600)/60),floor(rem(comptime,60)))
            end
            
        end
        
    case 'Multiframe'
        %% --- Multiframe ---
        I1_full=str2double(Data.imfstart):str2double(Data.imfstep):str2double(Data.imfend);
        time_full=str2double(Data.imfstart):(str2double(Data.imfend)+str2double(Data.imcstep));
        
        corrtime=zeros(P,1);
        valtime=zeros(P,1);
        savetime=zeros(P,1);
        interptime=zeros(P,1);
        
        %single-pulsed
        if round(1/Freq*10^6)==round(dt)
            time_full(2,:)=time_full(1,:);
        else
            %double-pulsed
            sample_t=1/Freq*10^6;
            for n=3:2:length(time_full)
                time_full(2,n)=floor(n/2)*sample_t/dt;
                time_full(2,n-1)=(floor((n-2)/2)*sample_t+dt)/dt;
            end
            time_full(2,end)=time_full(2,end-1)+1;
        end
        
        if I1(1)==I1_full(1)
            qstart=1;
        else
            qstart=Nmax+1;
        end
        if I1(end)==I1_full(end)
            qend=length(I1);
        else
            qend=length(I1)-Nmax;
        end
        frametime=nan(length(qstart:qend),1);
        
        for q=qstart:qend
            tf=tic;
            frametitle=['Frame' sprintf(['%0.' Data.imzeros 'i'],I1(q)) ' and Frame' sprintf(['%0.' Data.imzeros 'i'],I2(q))];
            
            %load dynamic mask and flip coordinates
            if strcmp(Data.masktype,'dynamic')
                mask = cast(imread([maskbase sprintf(['%0.' Data.maskzeros 'i.' Data.maskext],maskname(q))]),imClass);
                mask = flipud(mask);
            end
            
            %load image pairs, compute delta-t, and flip coordinates
            if q-Nmax<1 && q+Nmax>length(I1)
                N=min([q,length(I1)-q+1]);
            elseif q-Nmax<1
                N=q;
            elseif q+Nmax>length(I1)
                N=length(I1)-q+1;
            else
                N=Nmax+1;
            end
            im1=zeros(size(mask,1),size(mask,2),N,imClass); im2=im1;Dt=zeros(N,1,imClass);
            for n=1:N
                im1_temp=cast(imread([imbase sprintf(['%0.' Data.imzeros 'i.' Data.imext],I1(q)-(n-1))]),imClass);
                im2_temp=cast(imread([imbase2 sprintf(['%0.' Data.imzeros 'i.' Data.imext],I2(q)+(n-1))]),imClass);
                im1(:,:,n)=flipud(im1_temp(:,:,1));
                im2(:,:,n)=flipud(im2_temp(:,:,1));
                %                 if Zeromean(e)==1
                %                     im1(:,:,n)=im1(:,:,n)-mean(mean(im1(:,:,n)));
                %                     im2(:,:,n)=im2(:,:,n)-mean(mean(im2(:,:,n)));
                %                 end
                
                imind1= time_full(1,:)==I1(q)-(n-1);
                imind2= time_full(1,:)==I2(q)+(n-1);
                Dt(n)=time_full(2,imind2)-time_full(2,imind1);
            end
            imageSize=size(im1);
            
            % initialize grid and evaluation matrix for every pixel in image
            [XI,YI]=IMgrid(imageSize,[0 0]);
            %index-based pixel coordinates are integers on center, but we need XI and YI to
            %correspond to vector locations, which are centered on pixel edges for even sized
            %windows.  This means the center of pixel index 1 is actually at coordinate
            %0.5 from image edge, etc.
            % We need to do this because the vector positions we will use
            % for interpolation and deformation will be on the physical
            % grid, but will will need to know where the pixels are located
            % relative to them, not their indices.
            %BUT they are pixel centered for ODD sized windows.  (FIX) but
            %for now will ignore since even windows are more common.
            XI = XI - 0.5;
            YI = YI - 0.5;
            
            UI = Vel0.U;
            VI = Vel0.V;
            %             UI = BWO(1).*ones(size(XI),imClass);
            %             VI = BWO(2).*ones(size(YI),imClass);
            
            for e=1:P
                t1=tic;
                [X,Y]=IMgrid(imageSize,Gres(e,:),Gbuf(e,:));
                S=size(X);X=X(:);Y=Y(:);
                
                if ~strcmpi(Corr{e},'SPC')
                    U=zeros(size(X,1),3,N,imClass);
                    V=zeros(size(X,1),3,N,imClass);
                    C=zeros(size(X,1),3,N,imClass);
                    Di=zeros(size(X,1),3,N,imClass);
                    Cp=zeros(Wsize(e,1),Wsize(e,2),size(X,1),N,imClass);
                    Uval=zeros(size(X,1),3,imClass);
                    Vval=zeros(size(X,1),3,imClass);
                    Cval=zeros(size(X,1),3,imClass);
                    Dval=zeros(size(X,1),3,imClass);
                    Eval=repmat(reshape(downsample(downsample( mask(Y(1):Y(end),X(1):X(end)),Gres(e,2))',Gres(e,1))',length(X),1),[1 3]);
                    Eval(Eval==0)=-1;
                    Eval(Eval>0)=0;
                    Uc=zeros(sum(Eval(:,1)>=0),3,N,imClass);
                    Vc=zeros(sum(Eval(:,1)>=0),3,N,imClass);
                    Cc=zeros(sum(Eval(:,1)>=0),3,N,imClass);
                    Dc=zeros(sum(Eval(:,1)>=0),3,N,imClass);
                    
                    for t=1:N
                        Ub = reshape(downsample(downsample( UI(Y(1):Y(end),X(1):X(end)),Gres(e,2))',Gres(e,1))',length(X),1).*Dt(t);
                        Vb = reshape(downsample(downsample( VI(Y(1):Y(end),X(1):X(end)),Gres(e,2))',Gres(e,1))',length(X),1).*Dt(t);
                        
                        %correlate image pair
                        %                         [Xc,Yc,Uc(:,:,t),Vc(:,:,t),Cc(:,:,t)]=PIVwindowed(im1(:,:,t),im2(:,:,t),Corr(e),Wsize(e,:),Wres(e, :, :),0,D(e),Zeromean(e),Peaklocator(e),1,X(Eval(:,1)>=0),Y(Eval(:,1)>=0),Ub(Eval(:,1)>=0),Vb(Eval(:,1)>=0));
                        [Xc,Yc,Uc(:,:,t),Vc(:,:,t),Cc(:,:,t),Dc(:,:,t),Cp(:,:,:,t)]=PIVwindowed(im1(:,:,t),im2(:,:,t),Corr{e},Wsize(e,:),Wres(:, :, e),0,D(e,:),Zeromean(e),Peaklocator(e),1,frac_filt(e),saveplane(e),X(Eval(:,1)>=0),Y(Eval(:,1)>=0),Ub(Eval(:,1)>=0),Vb(Eval(:,1)>=0));
                    end
                    U(repmat(Eval>=0,[1 1 N]))=Uc;
                    V(repmat(Eval>=0,[1 1 N]))=Vc;
                    C(repmat(Eval>=0,[1 1 N]))=Cc;
                    Di(repmat(Eval>=0,[1 1 N]))=Dc;
                    
                    velmag=sqrt(U(:,1,:).^2+V(:,1,:).^2);
                    Qp=C(:,1,:)./C(:,2,:).*(1-ds./velmag);
                    %                     Qp=1-2.*exp(-0.5)./velmag.*(C(:,1,:)./C(:,2,:)-1).^(-1);
                    [Qmax,t_opt]=max(Qp,[],3); %#ok<ASGLU>
                    %                     for i=1:size(U,1)
                    %                         Uval(i,:)=U(i,:,t_opt(i));
                    %                         Vval(i,:)=V(i,:,t_opt(i));
                    %                         Cval(i,:)=C(i,:,t_opt(i));
                    %                         Dval(i,:)=Di(i,:,t_opt(i));
                    %                     end
                    %
                    %                     try
                    %                         U=Uval./repmat(Dt(t_opt)',[1 3]);
                    %                         V=Vval./repmat(Dt(t_opt)',[1 3]);
                    %                     catch
                    %                         U=Uval./repmat(Dt(t_opt),[1 3]);
                    %                         V=Vval./repmat(Dt(t_opt),[1 3]);
                    %                     end
                    for i=1:size(U,1)
                        U(i,:)=sum( (squeeze(U(i,1,:))./Dt).*squeeze(Qp(i,1,:)./sum(Qp(i,1,:))) );
                        V(i,:)=sum( (squeeze(V(i,1,:))./Dt).*squeeze(Qp(i,1,:)./sum(Qp(i,1,:))) );
                        Cval(i,:)=sum( (squeeze(C(i,1,:))    ).*squeeze(Qp(i,1,:)./sum(Qp(i,1,:))) );
                        Dval(i,:)=sum( (squeeze(Di(i,1,:))   ).*squeeze(Qp(i,1,:)./sum(Qp(i,1,:))) );
                    end
                    
                else
                    U=zeros(length(X),1);
                    V=zeros(length(X),1);
                    Eval=reshape(downsample(downsample( mask(Y(1):Y(end),X(1):X(end)),Gres(e,2))',Gres(e,1))',length(X),1);
                    Eval(Eval==0)=-1;
                    Eval(Eval>0)=0;
                    
                    Ub = reshape(downsample(downsample( UI(Y(1):Y(end),X(1):X(end)),Gres(e,2))',Gres(e,1))',length(X),1);
                    Vb = reshape(downsample(downsample( VI(Y(1):Y(end),X(1):X(end)),Gres(e,2))',Gres(e,1))',length(X),1);
                    %                     [Xc,Yc,Uc,Vc,Cc,t_optc]=PIVphasecorr(im1,im2,Wsize(e,:), Wres(e, :, :),0,D(e),Zeromean(e),Peakswitch(e),X(Eval>=0),Y(Eval>=0),Ub(Eval>=0),Vb(Eval>=0),Dt);
                    [Xc,Yc,Uc,Vc,Cc,t_optc]=PIVphasecorr(im1,im2,Wsize(e,:),Wres(:, :, e),0,D(e,:),Zeromean(e),Peakswitch(e),X(Eval>=0),Y(Eval>=0),Ub(Eval>=0),Vb(Eval>=0),Dt);
                    if Peakswitch(e)
                        C=zeros(length(X),3,imClass);
                        Di=zeros(length(X),3,imClass);
                        C(repmat(Eval,[1 3])>=0)=Cc;
                        t_opt=zeros(size(X),imClass);
                        t_opt(Eval>=0)=t_optc;
                    else
                        C=[];t_opt=[];Di=[];
                    end
                    U(Eval>=0)=Uc;V(Eval>=0)=Vc;
                end
                
                corrtime(e)=toc(t1);
                
                %validation
                if Valswitch(e)
                    t1=tic;
                    
                    [Uval,Vval,Evalval,Cval,Dval]=VAL(X,Y,U,V,Eval,C,Di,Threshswitch(e),UODswitch(e),Bootswitch(e),extrapeaks(e),...
                        Uthresh(e,:),Vthresh(e,:),UODwinsize(e,:,:),UODthresh(e,UODthresh(e,:)~=0)',Bootper(e),Bootiter(e),Bootkmax(e));
                    
                    valtime(e)=toc(t1);
                else
                    Uval=U(:,1);Vval=V(:,1);Evalval=Eval(:,1);
                    if ~isempty(C)
                        Cval=C(:,1);
                        Dval=Di(:,1);
                    else
                        Cval=[];
                        Dval=[];
                    end
                end
                
                %write output
                if Writeswitch(e)
                    t1=tic;
                    if Peakswitch(e)
                        if PeakVel(e) && ~strcmpi(Corr{e},'SPC')
                            U=[Uval(:,1),U(:,1:PeakNum(e))];
                            V=[Vval(:,1),V(:,1:PeakNum(e))];
                            Eval=[Evalval(:,1),Eval(:,1:PeakNum(e))];
                        else
                            U=Uval(:,1); V=Vval(:,1);Eval=Evalval(:,1);
                        end
                        if PeakMag(e)
                            C=[Cval(:,1),C(:,1:PeakNum(e))];
                            Di=[Dval(:,1),Di(:,1:PeakNum(e))];
                        else
                            C=[];
                            Di=[];
                        end
                    else
                        t_opt=[];
                    end
                    %convert to physical units
                    Xval=X;Yval=Y;
                    X=X*Mag;Y=Y*Mag;
                    U=U*Mag./dt;V=V*Mag./dt;
                    
                    %convert to matrix if necessary
                    if size(X,2)==1
                        [X,Y,U,V,Eval,C,Di]=matrixform(X,Y,U,V,Eval,C,Di);
                        if Peakswitch(e)
                            t_opt=reshape(t_opt,size(X,1),size(X,2));
                        end
                    end
                    
                    %remove nans from data, replace with zeros
                    U(Eval<0|isinf(U))=0;V(Eval<0|isinf(V))=0;
                    
                    if str2double(Data.datout)
                        %                         q_full=find(I1_full==I1(q),1,'first');
                        %                         time=(q_full-1)/Freq;
                        time=I1(q)/Freq;
                        %write_dat_val_C([pltdirec char(wbase(e,:)) sprintf(['%0.' Data.imzeros 'i.dat' ],I1(q))],X,Y,U,V,Eval,C,e,time,frametitle,t_opt);
                        write_dat_val_C([pltdirec char(wbase(e,:)) sprintf(['%0.' Data.imzeros 'i.dat' ],I1(q))],X,Y,U,V,Eval,C,Di,e,time,char(wbase(e,:)),t_opt);
                    end
                    if str2double(Data.multiplematout)
                        save([pltdirec char(wbase(e,:)) sprintf(['%0.' Data.imzeros 'i.mat' ],I1(q))],'X','Y','U','V','Eval','C','Di','t_opt')
                    end
                    if saveplane(e) && ~strcmpi(Corr{e},'SPC')
                        Xloc = Xc;Yloc=Yc;C_planes=Cp;%#ok
                        save(sprintf(['%s%scorrplanes_%0.' Data.imzeros 'i.mat' ],pltdirec,wbase{e,:},I1(q)),'Xloc','Yloc','C_planes')
                        clear Xloc Yloc C_planes
                    end
                    X=Xval;Y=Yval;
                    
                    savetime(e)=toc(t1);
                end
                U=Uval; V=Vval;
                
                if e~=P
                    %reshape from list of grid points to matrix
                    X=reshape(X,[S(1),S(2)]);
                    Y=reshape(Y,[S(1),S(2)]);
                    U=reshape(U(:,1),[S(1),S(2)]);
                    V=reshape(V(:,1),[S(1),S(2)]);
                    
                    t1=tic;
                    
                    %velocity smoothing
                    if Velsmoothswitch(e)==1
                        [U,V]=VELfilt(U,V,UODwinsize(e,:,:),Velsmoothfilt(e));
                    end
                    
                    %velocity interpolation -
                    %resample U(X,Y) and V(X,Y) onto UI(XI,YI) and
                    %VI(XI,YI) where XI and YI are a list of every
                    %pixel in the image plane. Velinterp is the type of
                    %interpolation to use.
                    % We have previously made sure XI and YI correspond
                    % to the vector positions of the pixel centers by
                    % shifting by -0.5.
                    UI = VFinterp(X,Y,U,XI,YI,Velinterp);
                    VI = VFinterp(X,Y,V,XI,YI,Velinterp);
                    
                    interptime(e)=toc(t1);
                end
                Uval=[];Vval=[];Cval=[];
            end
            
            eltime=toc(tf);
            %output text
            fprintf('\n----------------------------------------------------\n')
            fprintf(['Job: ',Data.batchname,'\n'])
            fprintf([frametitle ' Completed (' num2str(q+1-qstart) '/' num2str(length(qstart:qend)) ') at %s \n'], datestr(now));
            %             fprintf(1, 'Frame completed at %s \n', datestr(now)); % Print the date and time at which frame was completed
            fprintf('----------------------------------------------------\n')
            for e=1:P
                fprintf('correlation...                   %0.2i:%0.2i.%0.0f\n',floor(corrtime(e)/60),floor(rem(corrtime(e),60)),floor((rem(corrtime(e),60)-floor(rem(corrtime(e),60)))*10))
                if Valswitch(e)
                    fprintf('validation...                    %0.2i:%0.2i.%0.0f\n',floor(valtime(e)/60),floor(rem(valtime(e),60)),floor((rem(valtime(e),60)-floor(rem(valtime(e),60)))*10))
                end
                if Writeswitch(e)
                    fprintf('save time...                     %0.2i:%0.2i.%0.0f\n',floor(savetime(e)/60),floor(rem(savetime(e),60)),floor((rem(savetime(e),60)-floor(rem(savetime(e),60)))*10))
                end
                if e~=P
                    fprintf('velocity interpolation...        %0.2i:%0.2i.%0.0f\n',floor(interptime(e)/60),floor(rem(interptime(e),60)),floor((rem(interptime(e),60)-floor(rem(interptime(e),60)))*10))
                end
            end
            fprintf('total frame time...              %0.2i:%0.2i.%0.0f\n',floor(eltime/60),floor(rem(eltime,60)),floor((rem(eltime,60)-floor(rem(eltime,60)))*10))
            frametime(q+1-qstart)=eltime;
            comptime=nanmean(frametime(1:q+1-qstart))*(length(qstart:qend)-(q+1-qstart));
            fprintf('estimated job completion time... %0.2i:%0.2i:%0.2i\n\n',floor(comptime/3600),floor(rem(comptime,3600)/60),floor(rem(comptime,60)))
        end
end

%signal job complete
%beep,pause(0.2),beep

end


