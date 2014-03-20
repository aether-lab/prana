function [defaultdata] = buildDefaultPranaJob()

defaultdata.clientversion='2.0';
%defaultdata.version='2.0'; %gets set below in call to pranaPIVcode('version')
if ispc
    defaultdata.imdirec='C:\';
    defaultdata.outdirec='C:\';
    defaultdata.maskdirec='C:\';
else
    defaultdata.imdirec='/';
    defaultdata.outdirec='/';
    defaultdata.maskdirec='/';
end
defaultdata.imbase='Img_';
defaultdata.imzeros='6';
defaultdata.imext='tif';
defaultdata.imcstep='1';
defaultdata.imfstep='1';
defaultdata.imfstart='1';
defaultdata.imfend='1';

defaultdata.wrmag='1';
defaultdata.wrsamp='1';
defaultdata.wrsep='1';
defaultdata.batchname='Proc1';
defaultdata.datout='0';
defaultdata.multiplematout='1';

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

defaultdata.input_vel_type = 'none';
defaultdata.input_velocity = '';

defaultdata.masktype='none';
defaultdata.staticmaskname='';
defaultdata.maskbase='maskfor_Img_';
defaultdata.maskzeros='6';
defaultdata.maskext='tif';
defaultdata.maskfstep='1';
defaultdata.maskfstart='1';

defaultdata.runPIV='1';

defaultdata.PIV0.winres='32,32; 32,32';
%     defaultdata.PIV0.winres1='32,32';
%     defaultdata.PIV0.winres2='32,32';
defaultdata.PIV0.winsize='64,64';
defaultdata.PIV0.winauto='1';
defaultdata.PIV0.gridres='8,8';
defaultdata.PIV0.winoverlap='75,75';
defaultdata.PIV0.gridtype='1';
defaultdata.PIV0.gridbuf='8,8';
defaultdata.PIV0.BWO='0,0';
defaultdata.PIV0.corr='RPC';
defaultdata.PIV0.RPCd='2.8,2.8';
defaultdata.PIV0.frac_filt='1';
defaultdata.PIV0.zeromean='1';
defaultdata.PIV0.peaklocator='1';
defaultdata.PIV0.velsmooth='0';
defaultdata.PIV0.velsmoothfilt='2';
defaultdata.PIV0.deform_min ='1';
defaultdata.PIV0.deform_max ='1';
defaultdata.PIV0.deform_conv ='0.1';
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
defaultdata.PIV0.saveplane='0';
defaultdata.PIV0.outbase='PIV_';
defaultdata.PIV0.write='1';

defaultdata.PIV1=defaultdata.PIV0;
defaultdata.PIV2=defaultdata.PIV0;

defaultdata.passes='2';
defaultdata.method='2';
defaultdata.velinterp='3';
defaultdata.iminterp='2';
defaultdata.framestep='3';
defaultdata.PIVerror='0.1';
defaultdata.channel = '1';

% --- Tracking Default Values ---
% ID Default values
defaultdata.ID.runid        = '0';
defaultdata.ID.method       = '2';
defaultdata.ID.imthresh     = '10';
defaultdata.ID.savebase     = 'ID_';
% Sizing Default values
defaultdata.Size.runsize    = '0';
defaultdata.Size.method     = 'GEO';
defaultdata.Size.min_area   = '0';
defaultdata.Size.std        = '4';
defaultdata.Size.savebase   = 'SIZE_';
% Tracking Default values
defaultdata.Track.runtrack  = '0';
defaultdata.Track.method    = '2';
defaultdata.Track.prediction= '1';
defaultdata.Track.PIVweight = '0.5';
defaultdata.Track.radius    = '15';
defaultdata.Track.disweight = '1.0';
defaultdata.Track.sizeweight= '0.5';
defaultdata.Track.intensityweight = '0.5';
defaultdata.Track.estradius = '15';
defaultdata.Track.estweight = '.1';
defaultdata.Track.savebase  = 'Track_';
defaultdata.Track.trackdat  = '0';
defaultdata.Track.trackmat  = '1';
defaultdata.Track.vectors   = '3';
defaultdata.Track.iterations= '3';
% Tracking Validation Values
defaultdata.Track.valprops.run   = '1';
defaultdata.Track.valprops.valcoef = '0,0,0.2';
defaultdata.Track.valprops.valrad = '20,20,0';
defaultdata.Track.valprops.MAD_U = '1,0.75,0';
defaultdata.Track.valprops.MAD_V = '1,0.75,0';
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

%JJC: shouldn't this be in the catch statement as well?
if ~isfield(defaultdata,'outputpassbasename')
    defaultdata.outputpassbase = 'PIV_';
    defaultdata.PIV0.outbase = [defaultdata.outputpassbase 'pass0_'];
    defaultdata.PIV1.outbase = [defaultdata.outputpassbase 'pass1_'];
    defaultdata.PIV2.outbase = [defaultdata.outputpassbase 'pass2_'];
end

defaultdata.version=pranaPIVcode('version');  %why isn't this done in the catch statement above?
defaultdata.ptv_version=pranaPTVcode('version');
end