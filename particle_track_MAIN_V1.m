function [tracks]=particle_track_MAIN_V1(X2_est,Y2_est,Z2_est,SIZE1,SIZE2,trackprops,valprops)
%
% [tracks]=particle_track_MAIN_V1(X2_est,Y2_est,Z2_est,SIZE1,SIZE2,...
%   trackprops,valprops)
%
% PROGRAM DESCRIPTION
% This is a MAIN function for tracking particles between two consecutive
% images, where the particles have already been identified and sized using
% the particle_ID_MAIN and particle_size_MAIN functions.  The tracking
% fucntion also assumes that a displacement estimate has already been
% computed, typically using the particle_estloc_MAIN function.
%
% INPUT
%   X2_est, Y2_est, Z2_est - estimated location of the particles
%   SIZE1/SIZE2 - structured array containing the sizing results
%       .XYDiameter = 6-column matrix; 1st column is x-centorid location, 
%           2nd column is y-centroid location, 
%           3rd column is particle diameter, 
%           4th column is true max. intensity, 
%           5th column is particle id#, 
%           6th column is sizing method
%       .mapsizeinfo - (2 column array) defines the size (row col) of the 
%           each sized particle
%       .locxy - (2 column array) ROW/COL? location associated with the 
%           upper left pixel of the particle's rectangular projection
%   trackprops (structured array) tracking control parameters
%       .s_radius - radius to search for neighboring particles
%       .weights - [0-1 0-1 0-1] sets the relative emphasis for the pair
%               matching algorithm: inter-particle distance, size diff, &
%               max-intensity diff
%       .save_dir - full path of the save directory, make '0' for no saving
%       .slsh - '\' or '/'
%       .s_name - base save name
%       .s_num - intialize to [], controlled by the program
%       .plotfig - set to '1' to plot the results, otherwise set to '0'
%   valprops (structured array) validation control parameters
%       .num_pass - number of validation passes
%       .method - (vector of strings) 'none' 'coeff' 'mean' 'median'
%           'relaxation' sets the method to ID bad vectors and subsequently
%           revise the displacement estimate
%       .C_cutoff - max matching coeff to accept a track (between 0-1)
%       .s_radius - radius to search for neighboring tracks
%       .MAD_U/MAD_V - maximum allowable mean/median absolute deviation
%           before a track is considered 'bad'
%
%                               SAMPLE VALPROPS
%   valprops.numpass=4;
%   valprops.method={'median' 'median' 'median' 'coeff'};
%   valprops.C_cutoff=[1 1 1 0.2];
%   valprops.s_radius=[15 15 15 15];
%   valprops.MAD_U=[2 1.5 1 1];
%   valprops.MAD_V=[2 1.5 1 1];
%
% OUTPUT
%   tracks - main output array of the matched particle pairs:
%       [X1 X2 Y1 Y2 Z1 Z2 d1 d2 I1 I2 p#1 p#2 match_probability]
%       *where a lower match_coefficient is better, should vary
%        between 0 and 1*
%
%(v1) N.Cardwell - 11.18.2009

%set the locations for the particles in image 1 and 2
X1=SIZE1.XYDiameter(:,1);  X2=SIZE2.XYDiameter(:,1);
Y1=SIZE1.XYDiameter(:,2);  Y2=SIZE2.XYDiameter(:,2);
d1=SIZE1.XYDiameter(:,3);  d2=SIZE2.XYDiameter(:,3);
I1=SIZE1.XYDiameter(:,4);  I2=SIZE2.XYDiameter(:,4);
Z1=zeros(size(X1));        Z2=zeros(size(X2));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%MAIN TRACKING BLOCK
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%track the particles in the image pair using the 3D weighted
%nearest neighbor tracking method
[tracks]=weighted_nearest_neighbor3D(X1,X2,X2_est,Y1,Y2,Y2_est,...
    Z1,Z2,Z2_est,d1,d2,I1,I2,trackprops.weights,trackprops.s_radius);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Validation of the determined tracks
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%decision block - validate the measured tracks?
if valprops.run==1
for i=1:valprops.numpass
    %intialize temporary valprops array
    valprops_i.method=valprops.method{i};
    valprops_i.C_cutoff=valprops.C_cutoff(i);
    valprops_i.s_radius=valprops.s_radius(i);
    valprops_i.MAD_U=valprops.MAD_U(i);
    valprops_i.MAD_V=valprops.MAD_V(i);
    
    switch lower(valprops.method{i})
        case {'none'}
            %accepts the measured tracks with no post processing

        case {'coeff'}
            %retain tracks below the coefficient threshold - discard those above
            tracks=tracks((tracks(:,13) <= valprops_i.C_cutoff),:);

        case {'mean','median'}
            %call validation mean/median subfunction
            [MAD_ratio,MAD_ratio_hdr]=PTVval_meanandmedian(tracks,valprops_i);

            %refresh the location estimation with the validation results
            X2_est( MAD_ratio(:,5) ,1) = MAD_ratio(:,6);
            Y2_est( MAD_ratio(:,5) ,1) = MAD_ratio(:,7);

            %rerun tracking with new estimated locations
            [tracks]=weighted_nearest_neighbor3D(X1,X2,X2_est,Y1,Y2,Y2_est,...
                Z1,Z2,Z2_est,d1,d2,I1,I2,trackprops.weights,trackprops.s_radius);

        case {'relaxation'}
            %TO BE COMPLETED - see paper by Ohmi and Lee (2000) on the new
            %relaxation method for more information

    end
end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Saving of determined tracks
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%save tracking results and relevant info if called for by the user
if ~isempty(trackprops.save_dir)
    if exist(trackprops.save_dir,'dir')~=7
        fprintf('Making Save Directory %s \n',trackprops.save_dir)
        mkdir(trackprops.save_dir)
    end
    sname = sprintf('%%s%%0%0.0fd.mat',trackprops.Data.imzeros);
    save(fullfile(trackprops.save_dir,sprintf(sname,trackprops.s_name,trackprops.s_num)),...
        'tracks','trackprops','valprops');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Plotting of the determined tracks
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if trackprops.plotfig==1
    %plot the tracking results as a connected scatter plot of particle
    %positions in IM1 and IM2
    im_bounds=[1 800 1 1400 -1 1];
    %figure;  
    scatter(tracks(:,1),tracks(:,3),'.b');
    set(gca,'DataAspectRatio',[1 1 1]);
    set(gca,'YDir','reverse');
    hold on
    scatter(tracks(:,2),tracks(:,4),'+r');  hold on
    for j=1:size(tracks,1)
        line([tracks(j,1);tracks(j,2)],[tracks(j,3);tracks(j,4)]);
        hold on
    end
    axis(im_bounds(1:4))
    legend({'image1','image2'},'Location','NorthOutside','Orientation','horizontal')
end

end
