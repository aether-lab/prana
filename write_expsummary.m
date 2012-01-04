function write_expsummary(Data)

if ispc
    fname=[Data.outdirec,'\ExpSummary_',Data.batchname,date,'.txt'];
else
    fname=[Data.outdirec,'/ExpSummary_',Data.batchname,date,'.txt'];
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
    fprintf(fid,'Experimental Notes:\n');
    fprintf(fid,[Data.exp_notes{i},'\n']);
end

methods={'Multipass - DWO','Multigrid - DWO','Multigrid - Deform (DWO)',...
    'Multigrid - Ensemble (DWO)','Multigrid - Multiframe (DWO)'};
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

channel = {'Red','Green','Blue','Weighted Average','Mean','Color Ensemble'};
fprintf(fid,'\n----------------------Images and Masking---------------------\n');
% fprintf(fid,['Image Directory:               ',Data.imdirec,'\n']);
fprintf(fid,['Image Basename:                ',Data.imbase,'\n']);
fprintf(fid,['Image Zeros:                   ',Data.imzeros,'\n']);
fprintf(fid,['Image Extension:               ',Data.imext,'\n']);
fprintf(fid,['Image Frame Start:             ',Data.imfstart,'\n']);
fprintf(fid,['Image Frame Step:              ',Data.imfstep,'\n']);
fprintf(fid,['Image Frame End:               ',Data.imfend,'\n']);
fprintf(fid,['Image Correlation Step:        ',Data.imcstep,'\n']);
fprintf(fid,['Image Color Channel:           ',channel{str2double(Data.channel)},'\n']);
fprintf(fid,'Masking Type:                  ');
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
    fprintf(fid,['Mask Basename:                ',Data.maskbase,'\n']);
    fprintf(fid,['Mask Zeros:                   ',Data.maskzeros,'\n']);
    fprintf(fid,['Mask Extension:               ',Data.maskext,'\n']);
    fprintf(fid,['Mask Frame Start:             ',Data.maskfstart,'\n']);
    fprintf(fid,['Mask Frame Step:              ',Data.maskfstep,'\n']);
else
    fprintf(fid,'None\n');
end

for i=1:str2double(Data.passes)
    corr={'SCC','RPC','SPC'};
    peak={'Three-Point Gaussian','Four-Point Gaussian','Gaussian Least Squares'};
    y_n={'No','Yes'};
    A=eval(['Data.PIV' num2str(i)]);
    fprintf(fid,['\n------------------------Pass ',num2str(i),' Setup-------------------------\n']);
    if isfield(A,'winres1')
    fprintf(fid,['Window Resolution (first image) (pix):        ',A.winres1,'\n']);
    fprintf(fid,['Window Resolution (second image) (pix):       ',A.winres2,'\n']);
    else
    fprintf(fid,['Window Resolution (first image) (pix):        ',A.winres,'\n']);
    fprintf(fid,['Window Resolution (second image) (pix):       ',A.winres,'\n']);
    end        
    fprintf(fid,['Window Size (pix):                            ',A.winsize,'\n']);
    fprintf(fid,['Grid Resolution (pix):                        ',A.gridres,'\n']);
    fprintf(fid,['Window Overlap Percentage:                    ',A.winoverlap,'\n']);
    fprintf(fid,['Grid Buffer (pix):                            ',A.gridbuf,'\n']);
    if i==1
    fprintf(fid,['Bulk Window Offset (pix):                     ',A.BWO,'\n']);
    else
        fprintf(fid,['Bulk Window Offset (pix):                     ','0,0','\n']);
    end
    fprintf(fid,['Correlation:                                  ',corr{str2double(A.corr)},'\n']);
    if strcmpi(corr{str2double(A.corr)},'RPC')
        fprintf(fid,['   RPC Diameter:                              ',A.RPCd,'\n']);
    end
    fprintf(fid,['Zero-Mean Image Windows:                      ',y_n{str2double(A.zeromean)+1},'\n']);
    fprintf(fid,['Subpixel Peak Location Method:                ',peak{str2double(A.peaklocator)},'\n']);
    if str2double(Data.method)~=1 && str2double(Data.method)~=5
        fprintf(fid,['Smoothing:                                    ',y_n{str2double(A.velsmooth)+1},'\n']);
        if str2double(A.velsmooth)
            fprintf(fid,['Smoothing Size:                               ',A.velsmoothfilt,'\n']);
        end
    end
    fprintf(fid,'Validation Type(s):                           ');
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
            fprintf(fid,['Umin, Umax (pix):                             ',A.valuthresh,'\n']);
            fprintf(fid,['Vmin, Vmax (pix):                             ',A.valvthresh,'\n']);
        end
        if str2double(A.uod)
            uod_type={'Mean','Median'};
            fprintf(fid,['UOD Type:                                     ',uod_type{str2double(A.uod_type)},'\n']);
            fprintf(fid,['UOD Window Sizes:                             ',A.uod_window,'\n']);
            fprintf(fid,['UOD Thresholds:                               ',A.uod_thresh,'\n']);
        end
        if str2double(A.bootstrap)
            fprintf(fid,['Bootstrap Percent Sampled:                    ',A.bootstrap_percentsampled,'\n']);
            fprintf(fid,['Bootstrap Iterations:                         ',A.bootstrap_iterations,'\n']);
            fprintf(fid,['Bootstrap Passes:                             ',A.bootstrap_passes,'\n']);
        end
        fprintf(fid,['Try Additional Peaks:                         ',y_n{str2double(A.valextrapeaks)+1},'\n']);
    else
        fprintf(fid,'None\n');
    end
    fprintf(fid,['Write Output:                                 ',y_n{str2double(A.write)+1},'\n']);
    if str2double(A.write)
        fprintf(fid,['Output Basename:                              ',A.outbase,'\n']);
        fprintf(fid,['Save Add. Peak Info:                          ',y_n{str2double(A.savepeakinfo)+1},'\n']);
        if str2double(A.savepeakinfo)
            peaks={'1','1,2','1,2,3'};
            fprintf(fid,['Save Data for Peaks:                          ',peaks{str2double(A.corrpeaknum)},'\n']);
            fprintf(fid,['Save Peak Magnitude:                          ',y_n{str2double(A.savepeakmag)+1},'\n']);
            fprintf(fid,['Save Resulting Vel.:                          ',y_n{str2double(A.savepeakvel)+1},'\n']);
        end
    end
end

fclose(fid);