%     This file is part of prana, an open-source GUI-driven program for
%     calculating velocity fields using PIV or PTV.
%     Copyright (C) 2012  Virginia Polytechnic Institute and State
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


% prana Runner
close all; clear all; clc; drawnow

data.imdirec = '~/Desktop'; % Image directory
data.maskdirec = ''; % Mask directory
data.imbase  = 'B004_'; %Image basename
data.imzeros = '1'; % Number of digits following basename
data.imext = 'tif'; % Image extension
data.imfstart = '1'; % Start frame
data.imfend = '2'; % End frame
data.imcstep = '1'; %Correlation step
data.imfstep = '2'; %Frame step
data.channel               = '1';       % Color Channel. 1 = Red, 2 = Green, 3 = blue, 4 = weighted average, 5 = mean.
data.par                        = '1';      % Use parallel processing? 1 for yes, 0 for no. 
data.parprocessors = '2';       % Number of processors to use in parallel
data.outdirec = '~/Desktop';    %Output directory


%% Experiment Info (these data do not affect processing)
data.clientversion          = '0.99'; 
data.version                = '4.3'; 
data.wrmag                  = '1';      % Image resolution (um/pix)
data.wrsamp                 = '1';      % Sampling rate (Hz)
data.wrsep                  = '1';      % Laser pulse separation (us)
data.batchname              = 'Proc1';      % Not sure about this one
data.datout                 = '1';      % Data format: write .DAT files?
data.multiplematout         = '0';      % Data format: write .MAT files?
data.exp_date               = '';       % Date of experiment    
data.exp_L                  = '';  % Characteristic length (m)
data.exp_v0                 = '';       % Characteristic velocity (m/s)
data.exp_notes              = {'Camera Description:'  ''  'Lens Description:'  ''  'Notes:'  ''};   %Camera description, lens description, notes
data.exp_density            = '1000';       % Fluid density (kg/m^3)
data.exp_viscosity          = '1.308e-3';       % Fluid dynamic viscosity (Pa * sec)
data.exp_surfacetension     = '0.07197';        % Fluid surface tension (N/m)
data.exp_partD              = '';       %Particle diameter (um)
data.exp_partdensity        = '';       % Particle density (kg/m^3)
data.exp_wavelength         = '.532';       %Laser wavelength (um)
data.exp_pixelsize          = '';       % Pixel size (um)
data.exp_lensfocal          = '';       % Lens focal length (mm)
data.exp_lensfnum           = '';   % Lens f#
data.exp_micro              = '0';     % Micro PIV?
data.exp_NA                 = '';       %Numerical aperature
data.exp_n                  = '';       % Refractive index
data.exp_Re                 = '';       % Reynolds number
data.exp_St                 = '';       % Particle stokes number
data.exp_M                  = '';       % Magnification
data.exp_ROI                = '';       % Region of interes (m)
data.exp_diffractiondiameter= '';   % Diffraction diameter (um)
data.exp_depthoffocus       = '';       %Depth of focus (um)
data.splash                 = '0';      % Display splash screen?

%% Mask Information
% Data masking type
data.masktype               = 'none';       % Type of mask. 'none', 'static', or 'dynamic'
data.staticmaskname         = '';          % Static mask name
data.maskbase               = 'maskfor_Img_';       % Basename for masks
data.maskzeros              = '6';      % Number of digits following mask basename
data.maskext                = 'tif';        %Mask image extension
data.maskfstep              = '1';      %Frame step for dynamic mask
data.maskfstart             = '1';      % Start frame for dynamic mask
data.passes                 = '2';      % Number of passes
data.method                 = '1';      % Processing method: 1 = Multipass DWO; 2 = Multigrid DWO; 3 = Multigrid Deform DWO; 4 = Ensemble DWO; 5 = Multiframe
data.velinterp              = '3';      % Velocity interpolation method: 1 = Nearest Neighbor; 2 = Bilinear; 3 = Bicubic
data.iminterp               = '1';      % Window deformation interpolation type. 1 = Cardinal function, 2 = Cardinal function w / Blackman filter
data.framestep              = '3';      % Multi-frame max step
data.PIVerror               = '0.1';      % Multi-frame error

%% data.PIV0
data.PIV0.outbase           = 'PIV_';       %Pass output basename
data.PIV0.winres            = '32,32; 32,32';      % Inhomogeneous window resolutions. Redundant field, needs to be removed.
data.PIV0.winres1          = '32, 32';      % Window resolution in first image
data.PIV0.winres2          = '32, 32';      % Window resolution in second image
data.PIV0.winsize           = '64,64';      % Size of interrogation region
data.PIV0.winauto           = '1';      % Auto window size? 
data.PIV0.gridres           = '8,8';      % Grid resolution 
data.PIV0. winoverlap       = '75,75';      % Window overlap
data.PIV0.gridtype          = '1';      % Grid type. Not sure what this means yet.
data.PIV0.gridbuf           = '8,8';      % Grid buffer
data.PIV0.BWO               = '0,0';      % Bulk window offset
data.PIV0.corr              = '2';      % Correlation type. '1' for SCC, '2' for RPC, '3' for SPC
data.PIV0.RPCd              = '2.8';      % RPC Filter Diameter
data.PIV0.zeromean          = '0';      % Perform zero-mean filtering?
data.PIV0.peaklocator       = '1';      % Peak locator type. '1' for three point gaussian, '2' for four point gaussian, '3' for Gaussian least squares regression 
data.PIV0.velsmooth         = '0';      % Perform velocity smoothing?
data.PIV0.velsmoothfilt     = '2';      % Velocity smoothing filter size
data.PIV0.val               = '0';      % Perform validation?
data.PIV0.uod               = '1';      % Perform universal outlier detection?
data.PIV0.bootstrap         = '0';      % Perform bootstrapping?
data.PIV0.thresh            = '0';      % Perform velocity thresholding?
data.PIV0.uod_type          = '2';      % UOD location parameter. '1' for mean, '2' for median.
data.PIV0.uod_window        = '3,3;3,3';      % UOD window size
data.PIV0.uod_thresh        = '3,2';      % UOD threshold
data.PIV0.bootstrap_percentsampled = '15';      % Bootstrapping percent sampled
data.PIV0.bootstrap_iterations = '700';      % Bootstrapping iterations per frame
data.PIV0.bootstrap_passes  = '12';      % Bootstraping number of passes
data.PIV0.valuthresh        = '-16,16';      % Velocity threshold (U-velocity)
data.PIV0.valvthresh        = '-16,16';      % Velocity threshold (V-velocity)
data.PIV0.valextrapeaks     = '0';      % Try additional peaks if validation fails?
data.PIV0.savepeakinfo      = '0';      % Save additional peak information?
data.PIV0.corrpeaknum       = '1';      % Peak info to save. '1' for peak 1, '2' for peaks 1 and 2, '3' for peaks 1, 2, and 3. 
data.PIV0.savepeakmag       = '0';      % Save peak magnitude?
data.PIV0.savepeakvel       = '0';      % Save peak velocity (location?)
data.PIV0.write             = '1';      % Write pass?


%% data.PIV1
data.PIV1.outbase           = 'PIV_';       %Pass output basename
data.PIV1.winres            = '32,32; 32,32';      % Inhomogeneous window resolutions. Redundant field, needs to be removed.
data.PIV1.winres1          = '32, 32';      % Window resolution in first image
data.PIV1.winres2          = '32, 32';      % Window resolution in second image
data.PIV1.winsize           = '64,64';      % Size of interrogation region
data.PIV1.winauto           = '1';      % Auto window size? 
data.PIV1.gridres           = '8,8';      % Grid resolution 
data.PIV1. winoverlap       = '75,75';      % Window overlap
data.PIV1.gridtype          = '1';      % Grid type. Not sure what this means yet.
data.PIV1.gridbuf           = '8,8';      % Grid buffer
data.PIV1.BWO               = '0,0';      % Bulk window offset
data.PIV1.corr              = '2';      % Correlation type. '1' for SCC, '2' for RPC, '3' for SPC
data.PIV1.RPCd              = '2.8';      % RPC Filter Diameter
data.PIV1.zeromean          = '0';      % Perform zero-mean filtering?
data.PIV1.peaklocator       = '1';      % Peak locator type. '1' for three point gaussian, '2' for four point gaussian, '3' for Gaussian least squares regression 
data.PIV1.velsmooth         = '0';      % Perform velocity smoothing?
data.PIV1.velsmoothfilt     = '2';      % Velocity smoothing filter size
data.PIV1.val               = '0';      % Perform validation?
data.PIV1.uod               = '1';      % Perform universal outlier detection?
data.PIV1.bootstrap         = '0';      % Perform bootstrapping?
data.PIV1.thresh            = '0';      % Perform velocity thresholding?
data.PIV1.uod_type          = '2';      % UOD location parameter. '1' for mean, '2' for median.
data.PIV1.uod_window        = '3,3;3,3';      % UOD window size
data.PIV1.uod_thresh        = '3,2';      % UOD threshold
data.PIV1.bootstrap_percentsampled = '15';      % Bootstrapping percent sampled
data.PIV1.bootstrap_iterations = '700';      % Bootstrapping iterations per frame
data.PIV1.bootstrap_passes  = '12';      % Bootstraping number of passes
data.PIV1.valuthresh        = '-16,16';      % Velocity threshold (U-velocity)
data.PIV1.valvthresh        = '-16,16';      % Velocity threshold (V-velocity)
data.PIV1.valextrapeaks     = '0';      % Try additional peaks if validation fails?
data.PIV1.savepeakinfo      = '0';      % Save additional peak information?
data.PIV1.corrpeaknum       = '1';      % Peak info to save. '1' for peak 1, '2' for peaks 1 and 2, '3' for peaks 1, 2, and 3. 
data.PIV1.savepeakmag       = '0';      % Save peak magnitude?
data.PIV1.savepeakvel       = '0';      % Save peak velocity (location?)
data.PIV1.write             = '1';      % Write pass?


%% data.PIV2
data.PIV2.outbase           = 'PIV_';       %Pass output basename
data.PIV2.winres            = '32,32; 32,32';      % Inhomogeneous window resolutions. Redundant field, needs to be removed.
data.PIV2.winres1          = '32, 32';      % Window resolution in first image
data.PIV2.winres2          = '32, 32';      % Window resolution in second image
data.PIV2.winsize           = '64,64';      % Size of interrogation region
data.PIV2.winauto           = '1';      % Auto window size? 
data.PIV2.gridres           = '8,8';      % Grid resolution 
data.PIV2. winoverlap       = '75,75';      % Window overlap
data.PIV2.gridtype          = '1';      % Grid type. Not sure what this means yet.
data.PIV2.gridbuf           = '8,8';      % Grid buffer
data.PIV2.BWO               = '0,0';      % Bulk window offset
data.PIV2.corr              = '2';      % Correlation type. '1' for SCC, '2' for RPC, '3' for SPC
data.PIV2.RPCd              = '2.8';      % RPC Filter Diameter
data.PIV2.zeromean          = '0';      % Perform zero-mean filtering?
data.PIV2.peaklocator       = '1';      % Peak locator type. '1' for three point gaussian, '2' for four point gaussian, '3' for Gaussian least squares regression 
data.PIV2.velsmooth         = '0';      % Perform velocity smoothing?
data.PIV2.velsmoothfilt     = '2';      % Velocity smoothing filter size
data.PIV2.val               = '0';      % Perform validation?
data.PIV2.uod               = '1';      % Perform universal outlier detection?
data.PIV2.bootstrap         = '0';      % Perform bootstrapping?
data.PIV2.thresh            = '0';      % Perform velocity thresholding?
data.PIV2.uod_type          = '2';      % UOD location parameter. '1' for mean, '2' for median.
data.PIV2.uod_window        = '3,3;3,3';      % UOD window size
data.PIV2.uod_thresh        = '3,2';      % UOD threshold
data.PIV2.bootstrap_percentsampled = '15';      % Bootstrapping percent sampled
data.PIV2.bootstrap_iterations = '700';      % Bootstrapping iterations per frame
data.PIV2.bootstrap_passes  = '12';      % Bootstraping number of passes
data.PIV2.valuthresh        = '-16,16';      % Velocity threshold (U-velocity)
data.PIV2.valvthresh        = '-16,16';      % Velocity threshold (V-velocity)
data.PIV2.valextrapeaks     = '0';      % Try additional peaks if validation fails?
data.PIV2.savepeakinfo      = '0';      % Save additional peak information?
data.PIV2.corrpeaknum       = '1';      % Peak info to save. '1' for peak 1, '2' for peaks 1 and 2, '3' for peaks 1, 2, and 3. 
data.PIV2.savepeakmag       = '0';      % Save peak magnitude?
data.PIV2.savepeakvel       = '0';      % Save peak velocity (location?)
data.PIV2.write             = '1';      % Write pass?


%% Run Prana
pranaPIVcode(data);
