function expsummary = write_expsummary(Data)

expsummary = [];

expsummary = [expsummary sprintf(['Summary for PIV Job: ',Data.batchname,'\n'])];
expsummary = [expsummary sprintf('\n--------------------Experiment Parameters--------------------\n')];
expsummary = [expsummary sprintf(['Date of Experiment:            ',Data.exp_date,'\n'])];
expsummary = [expsummary sprintf(['Laser Wavelength (um):         ',Data.exp_wavelength,'\n'])];
expsummary = [expsummary sprintf(['Laser Pulse Separation (us):   ',Data.wrsep,'\n'])];
expsummary = [expsummary sprintf(['Camera Frame Rate (Hz):        ',Data.wrsamp,'\n'])];
expsummary = [expsummary sprintf(['Pixel Size (um):               ',Data.exp_pixelsize,'\n'])];
expsummary = [expsummary sprintf(['Image Resolution (um/pix):     ',Data.wrmag,'\n'])];
expsummary = [expsummary sprintf(['Lens Focal Length (mm):        ',Data.exp_lensfocal,'\n'])];
if ~str2double(Data.exp_micro)
    expsummary = [expsummary sprintf(['Lens f#:                       ',Data.exp_lensfnum,'\n'])];
else
    expsummary = [expsummary sprintf(['Numerical Aperture:            ',Data.exp_NA,'\n'])];
    expsummary = [expsummary sprintf(['Index of Refraction:           ',Data.exp_n,'\n'])];
end
expsummary = [expsummary sprintf(['Particle Diameter (um):        ',Data.exp_partD,'\n'])];
expsummary = [expsummary sprintf(['Particle Density (kg/m^3):     ',Data.exp_partdensity,'\n'])];
expsummary = [expsummary sprintf(['Fluid Density (kg/m^3):        ',Data.exp_density,'\n'])];
expsummary = [expsummary sprintf(['Dynamic Viscosity (Pa*s):      ',Data.exp_viscosity,'\n'])];
expsummary = [expsummary sprintf(['Surface Tension (N/m):         ',Data.exp_surfacetension,'\n'])];
expsummary = [expsummary sprintf(['Characteristic Length (m):     ',Data.exp_L,'\n'])];
expsummary = [expsummary sprintf(['Characteristic Velocity (m/s): ',Data.exp_v0,'\n'])];
expsummary = [expsummary sprintf(['Reynolds Number:               ',Data.exp_Re,'\n'])];
expsummary = [expsummary sprintf(['Particle Stokes Number:        ',Data.exp_St,'\n'])];
expsummary = [expsummary sprintf(['Magnification Factor:          ',Data.exp_M,'\n'])];
expsummary = [expsummary sprintf(['Region of Interest (m):        ',Data.exp_ROI,'\n'])];
expsummary = [expsummary sprintf(['Diffraction Diameter (um):     ',Data.exp_diffractiondiameter,'\n'])];
expsummary = [expsummary sprintf(['Depth of Focus (um):           ',Data.exp_depthoffocus,'\n'])];

for i=1:size(Data.exp_notes,1)
    expsummary = [expsummary sprintf('Experimental Notes:\n')];
    expsummary = [expsummary sprintf([Data.exp_notes{i},'\n'])];
end

methods={'Multipass - DWO','Multigrid - DWO','Multigrid - Deform (DWO)',...
    'Multigrid - Ensemble (DWO)','Multigrid - Ensemble Deform (DWO)','Multigrid - Multiframe (DWO)'};
velinterp={'Nearest Neighbor','Bilinear','Cubic'};
iminterp={'Cardinal Function','Cardinal Function w/ Blackman Filter'};
expsummary = [expsummary sprintf('\n-----------------------PIV Processing------------------------\n')];
    expsummary = [expsummary sprintf(['Algorithm:                     ',methods{str2double(Data.method)},'\n'])];
if str2double(Data.method)~=1 && str2double(Data.method)~=6
    expsummary = [expsummary sprintf(['Velocity Interp Function:      ',velinterp{str2double(Data.velinterp)},'\n'])];
end
if any(str2double(Data.method)==[3 5])
    expsummary = [expsummary sprintf(['Image Interpolation Function:  ',iminterp{str2double(Data.iminterp)},'\n'])];
end
if str2double(Data.method)>=6
    expsummary = [expsummary sprintf(['PIV Error:                     ',Data.PIVerror,'\n'])];
    expsummary = [expsummary sprintf(['Maximum Framestep:             ',Data.framestep,'\n'])];
end

channel = {'Red (Grey Scale)','Green','Blue','Weighted Average','Mean','Color Ensemble'};
expsummary = [expsummary sprintf('\n----------------------Images and Masking---------------------\n')];
expsummary = [expsummary sprintf(['Image Directory:               ',Data.imdirec,'\n'])];
expsummary = [expsummary sprintf(['Image Basename:                ',Data.imbase,'\n'])];
expsummary = [expsummary sprintf(['Image Zeros:                   ',Data.imzeros,'\n'])];
expsummary = [expsummary sprintf(['Image Extension:               ',Data.imext,'\n'])];
expsummary = [expsummary sprintf(['Image Frame Start:             ',Data.imfstart,'\n'])];
expsummary = [expsummary sprintf(['Image Frame Step:              ',Data.imfstep,'\n'])];
expsummary = [expsummary sprintf(['Image Frame End:               ',Data.imfend,'\n'])];
expsummary = [expsummary sprintf(['Image Correlation Step:        ',Data.imcstep,'\n'])];
expsummary = [expsummary sprintf(['Image Color Channel:           ',channel{str2double(Data.channel)},'\n'])];
expsummary = [expsummary sprintf(['Output Directory:              ',Data.outdirec,'\n'])];
expsummary = [expsummary sprintf('Masking Type:                  ')];
if strcmp(Data.masktype,'static')
    if ispc
        slshind=strfind(Data.staticmaskname,'\');
    else
        slshind=strfind(Data.staticmaskname,'/');
    end
    expsummary = [expsummary sprintf(['Static\nMask File: ',Data.staticmaskname(slshind(end)+1:end),'\n'])];
elseif strcmp(Data.masktype,'dynamic')
    expsummary = [expsummary sprintf('Dynamic\n')];
    expsummary = [expsummary sprintf(['Mask Directory:               ',Data.maskdirec,'\n'])];
    expsummary = [expsummary sprintf(['Mask Basename:                ',Data.maskbase,'\n'])];
    expsummary = [expsummary sprintf(['Mask Zeros:                   ',Data.maskzeros,'\n'])];
    expsummary = [expsummary sprintf(['Mask Extension:               ',Data.maskext,'\n'])];
    expsummary = [expsummary sprintf(['Mask Frame Start:             ',Data.maskfstart,'\n'])];
    expsummary = [expsummary sprintf(['Mask Frame Step:              ',Data.maskfstep,'\n'])];
else
    expsummary = [expsummary sprintf('None\n')];
end        

for i=1:str2double(Data.passes)
    corr={'SCC','RPC','GCC','FWC','SPC'};
    peak={'Three-Point Gaussian','Four-Point Gaussian','Gaussian Least Squares'};
    y_n={'No','Yes'};
    A=eval(['Data.PIV' num2str(i)]);
    expsummary = [expsummary sprintf(['\n------------------------Pass ',num2str(i),' Setup-------------------------\n'])];
    if isfield(A,'winres1')
    expsummary = [expsummary sprintf(['Window Resolution (first image) (pix):        ',A.winres1,'\n'])];
    expsummary = [expsummary sprintf(['Window Resolution (second image) (pix):       ',A.winres2,'\n'])];
    else
    expsummary = [expsummary sprintf(['Window Resolution (first image) (pix):        ',A.winres,'\n'])];
    expsummary = [expsummary sprintf(['Window Resolution (second image) (pix):       ',A.winres,'\n'])];
    end        
    expsummary = [expsummary sprintf(['Window Size (pix):                            ',A.winsize,'\n'])];
    expsummary = [expsummary sprintf(['Grid Resolution (pix):                        ',A.gridres,'\n'])];
    expsummary = [expsummary sprintf(['Window Overlap Percentage:                    ',A.winoverlap,'\n'])];
    expsummary = [expsummary sprintf(['Grid Buffer (pix):                            ',A.gridbuf,'\n'])];
    if i==1
    expsummary = [expsummary sprintf(['Bulk Window Offset (pix):                     ',A.BWO,'\n'])];
    else
        expsummary = [expsummary sprintf(['Bulk Window Offset (pix):                     ','0,0','\n'])];
    end
    expsummary = [expsummary sprintf(['Correlation:                                  ',corr{str2double(A.corr)},'\n'])];
    if strcmpi(corr{str2double(A.corr)},'RPC')
        expsummary = [expsummary sprintf(['   RPC Diameter:                              ',A.RPCd,'\n'])];
    end
    expsummary = [expsummary sprintf(['Zero-Mean Image Windows:                      ',y_n{str2double(A.zeromean)+1},'\n'])];
    expsummary = [expsummary sprintf(['Subpixel Peak Location Method:                ',peak{str2double(A.peaklocator)},'\n'])];
    if str2double(Data.method)~=1 && str2double(Data.method)~=5
        expsummary = [expsummary sprintf(['Smoothing:                                    ',y_n{str2double(A.velsmooth)+1},'\n'])];
        if str2double(A.velsmooth)
            expsummary = [expsummary sprintf(['Smoothing Size:                               ',A.velsmoothfilt,'\n'])];
        end
    end
    expsummary = [expsummary sprintf('Validation Type(s):                           ')];
    if str2double(A.val)
        if str2double(A.thresh)
            expsummary = [expsummary sprintf('Thresholding ')];
        end
        if str2double(A.uod)
            expsummary = [expsummary sprintf('UOD ')];
        end
        if str2double(A.bootstrap)
            expsummary = [expsummary sprintf('Bootstrapping')];
        end
        expsummary = [expsummary sprintf('\n')];
        if str2double(A.thresh)
            expsummary = [expsummary sprintf(['Umin, Umax (pix):                             ',A.valuthresh,'\n'])];
            expsummary = [expsummary sprintf(['Vmin, Vmax (pix):                             ',A.valvthresh,'\n'])];
        end
        if str2double(A.uod)
            uod_type={'Mean','Median'};
            expsummary = [expsummary sprintf(['UOD Type:                                     ',uod_type{str2double(A.uod_type)},'\n'])];
            expsummary = [expsummary sprintf(['UOD Window Sizes:                             ',A.uod_window,'\n'])];
            expsummary = [expsummary sprintf(['UOD Thresholds:                               ',A.uod_thresh,'\n'])];
        end
        if str2double(A.bootstrap)
            expsummary = [expsummary sprintf(['Bootstrap Percent Sampled:                    ',A.bootstrap_percentsampled,'\n'])];
            expsummary = [expsummary sprintf(['Bootstrap Iterations:                         ',A.bootstrap_iterations,'\n'])];
            expsummary = [expsummary sprintf(['Bootstrap Passes:                             ',A.bootstrap_passes,'\n'])];
        end
        expsummary = [expsummary sprintf(['Try Additional Peaks:                         ',y_n{str2double(A.valextrapeaks)+1},'\n'])];
    else
        expsummary = [expsummary sprintf('None\n')];
    end
    expsummary = [expsummary sprintf(['Write Output:                                 ',y_n{str2double(A.write)+1},'\n'])];
    if str2double(A.write)
        expsummary = [expsummary sprintf(['Output Basename:                              ',A.outbase,'\n'])];
        expsummary = [expsummary sprintf(['Save Add. Peak Info:                          ',y_n{str2double(A.savepeakinfo)+1},'\n'])];
        if str2double(A.savepeakinfo)
            peaks={'1','1,2','1,2,3'};
            expsummary = [expsummary sprintf(['Save Data for Peaks:                          ',peaks{str2double(A.corrpeaknum)},'\n'])];
            expsummary = [expsummary sprintf(['Save Peak Magnitude:                          ',y_n{str2double(A.savepeakmag)+1},'\n'])];
            expsummary = [expsummary sprintf(['Save Resulting Vel.:                          ',y_n{str2double(A.savepeakvel)+1},'\n'])];
        end
        expsummary = [expsummary sprintf(['Save Correlation Planes:                      ',y_n{str2double(A.savepeakinfo)+1},'\n'])];
    end
end

% This string contains the date and time this function was called and
% appends it to the end of the expsummary.  This way when the jobs
% internals are changed but not the batch name a different text file will
% be created and not over written.
dateinfo = datestr(now);
dateinfo(12) = '-';
dateinfo([15 18]) = '.';
if ispc
    fname=[Data.outdirec,'\ExpSummary_',Data.batchname,'_',dateinfo,'.txt'];
else
    fname=[Data.outdirec,'/ExpSummary_',Data.batchname,'_',dateinfo,'.txt'];
end
fid=fopen(fname,'w');
if fid==-1
    try
        mkdir(Data.outdirec)
        fid=fopen(fname,'w');
        if fid==-1
            error(['error writing experiment summary ',fname])
        end
    catch ME
        error('Error writing experiment summary %s\n\n%s\n',fname,ME(1).message)
    end
end

fprintf(fid,expsummary);

fclose(fid);

fprintf(expsummary)

