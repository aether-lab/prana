function [X2_est,Y2_est,Z2_est]=particle_estloc_MAIN_V2(completed_tracks,SIZE1,estlocprops)
%
% [X2_est,Y2_est,Z2_est]=particle_estloc_MAIN_V1(completed_tracks,SIZE1,...
%   estlocprops)
%
% PROGRAM DESCRIPTION
% This is a MAIN function for providing a location prediction for the
% motion of identfied particles between two images.  Currently the function
% utilizses three predictive methods: PIV-driven, PTV-driven, and a hybrid 
% PIV/PTV-driven mode.  It is important to note that the PIV mode is 
% designed to use information from the current image pair, while the PTV 
% mode must utilize information from a previous image pair (i.e. previous
% tracks).  The PTV location prediction will work well assuming that the
% data is sufficiently resolved in time, otherwise the PIV suggested.
%
% INPUT
%   completed_tracks (cell array - [X1 X2 Y1 Y2 Z1 Z2 d1 d2 I1 I2 p#1 p#2
%       match_probability])  particle tracks from the previous image pair
%   SIZE1 (structured array) - contains XYDiameter which is has information
%       on the location/size/intensity of particles in the 1st image
%       (X1 Y1 D1 I1 p#1 size-method)
%   estlocprops (structured array) - control parameters
%       .method - 'none' 'piv' 'ptv' 'piv-ptv'
%       .PIV_PTV_weight - weighting value for the combined method; (1-0)
%           where 1-PIV and 0-PTV; linear weighting scheme
%       .PIVprops - parameters related to PIV location prediction
%           .load_dir - full path of the PIV flowfields
%           .load_name - basename of the files
%           .slsh - '\' or '/'
%           .precision - number of zeros in the filename
%           .extension - usually '.plt'
%           .framelist - vector of the frame numbers [0 1 2....]
%           .frame1 - current frame number, set to [] initially
%       .PTVprops - parameters related to PTV location prediction
%           .predict_mode - 'static' or 'dynamic'; static mode used a fixed
%               search radius to id neighboring vectors and perform
%               weighted averaging; dynamic mode varies the
%               search/averaging radius until a set # of vectors is located
%           .r_weight - radius of the search/averaging window
%           .edgeval - value of the Gaussian weighting function at
%               r_weight; should vary BETWEEEN 0 and 1
%           .numvecs - min # of vectors for serach\averaging (dynamic only)
%           .max_iterations - prevents infinite loops (dynamic only)
%        .save_dir - saving filepath; make '0' if not saving
%        .slsh - '\' or '/' (PC or Linux)
%        .s_name - base saving name
%        .s_num - saving number
%       
% OUTPUT
%   X2_est, Y2_est, Z2_est - estimated location of the particles
%
%(v1) N.Cardwell - 11.18.2009
%(v2) N.Cardwell - 04.01.2010

%extract particle information from the SIZE1 structured array
X1=SIZE1.XYDiameter(:,1);  Y1=SIZE1.XYDiameter(:,2);  Z1=zeros(size(X1));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%MAIN LOGIC BLOCK FOR DETERMINING THE LOCATION ESTIMATION FOR IM2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nnz( strcmpi(estlocprops.method,{'none' 'piv' 'ptv' 'piv-ptv'}) ) == 0
    fprintf('ERROR - An incorrect estimation location method was provided');
    return
end

switch lower(estlocprops.method)
    case {'none'}
        X2_est=X1;   Y2_est=Y1;   Z2_est=Z1;
    
    case {'piv'}
        %call the PIVestloc subfunction
        [X2_est,Y2_est,Z2_est]=PIVestloc_V1(X1,Y1,Z1,estlocprops);
    
    case {'ptv'}
        %call the PTVestloc subfunction
        [X2_est,Y2_est,Z2_est]=PTVestloc_V1(X1,Y1,Z1,estlocprops,completed_tracks);
    
    case {'piv-ptv'}
        %call the both the PIVestloc and PTVestloc subfunctions
        [X2_estPIV,Y2_estPIV,Z2_estPIV]=PIVestloc_V1(X1,Y1,Z1,estlocprops);
        [X2_estPTV,Y2_estPTV,Z2_estPTV]=PTVestloc_V1(X1,Y1,Z1,estlocprops,completed_tracks);
        
        %apply the user-defined weighting scheme
        X2_est = X2_estPIV.*estlocprops.PIV_PTV_weight + ...
            X2_estPTV.*(1-estlocprops.PIV_PTV_weight);
        Y2_est = Y2_estPIV.*estlocprops.PIV_PTV_weight + ...
            Y2_estPTV.*(1-estlocprops.PIV_PTV_weight);
        Z2_est = Z2_estPIV.*estlocprops.PIV_PTV_weight + ...
            Z2_estPTV.*(1-estlocprops.PIV_PTV_weight);
end

%section for saving sized particle information
% if ~isempty(estlocprops.save_dir)
%     if exist(estlocprops.save_dir,'dir')~=7
%         fprintf('Making Save Directory %s \n',estlocprops.save_dir)
%         mkdir(estlocprops.save_dir)
%     end
%     sname = sprintf('%%s%%0%0.0fd',estlocprops.Data.imzeros);
%     save(fullfile(estlocprops.save_dir,sprintf(sname,estlocprops.s_name,estlocprops.s_num)),...
%         'estlocprops','X2_est','Y2_est','Z2_est');
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%plotting code - COMMENT OUT FOR NORMAL OPERATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% dummy=cell2mat(completed_tracks);
% figure;
% quiver(dummy(:,1),dummy(:,3),dummy(:,2)-dummy(:,1),dummy(:,4)-dummy(:,3),0,'Color','k');
% set(gca,'DataAspectRatio',[1 1 1]);
% hold on
% 
% %quiver(X1,Y1,X2_estPIV-X1,Y2_estPIV-Y1,0,'Color','r');
% %quiver(X1,Y1,X2_estPTV-X1,Y2_estPTV-Y1,0,'Color','b');
% quiver(X1,Y1,X2_est-X1,Y2_est-Y1,0,'Color','m');
% legend({'completed_tracks' 'PIV_est' 'PTV_est' 'Est_loc'});

end
