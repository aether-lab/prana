function [p_matrix,peaks,num_p]=particle_ID_MAIN_V1(im,particleIDprops,s_num)
%
% [p_matrix,peaks,num_p]=particle_ID_MAIN_V1(im,particleIDprops);
%
% This is a MAIN function for particle identification and segmentation from
% grayscale images.  It calls several differet sub-functions (contained
% below) depending the image parameters, specifically the particle image
% density.  For sparse-to-medium seeding densities, the 'blob' method is
% recommended.  For a densly seeded image, with many overlapped particle
% images, the 'dynamic' method is recommended.  The combined method is a
% beta-version and should be used with caution since it has not been
% rigorly tested
%
% INPUTS
%    im1,im2 - image to be processesed; should be uint8
%    particleIDprops - processing parameters related to particle ID
%        particleIDprops.method - selects the ID/segmentation method
%            'blob','dynamic','combined'
%        particleIDprops.v - base thresholding for the image
%        particleIDprops.contrast_ratio - ratio of max particle intensity to
%            possible pixel intensity, used to limit the intensity range a
%            particle can expand from its max value ('dynamic' only)
%            *typically set to 0, higher values can reduce processing time*
%        particleIDprops.save_dir - saveing filepath; make '0' if not saving
%        particleIDprops.slsh - '\' or '/' (PC or Linux)
%        particleIDprops.s_name - base saving name
%        particleIDprops.s_num - saving number
%
% OUTPUTS
%    p_matrix (2D array) - matrix of particle identification and extent
%    peaks (2D array) - matrix of identified image peaks (by image erosion)
%    num_p (num) - number of identified particles in p_matrix
%
%(v1) N. Cardwell - 10.29.2009
%(v1.1) N.Cardwell - 11.1.2009 (combined particleIDprops and saveprops into
%   a single input, particleIDprops

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% identify and segregate the particles %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%chose between different particle ID and segregation methods

switch lower(particleIDprops.method)
    case {'blob'}
        [p_matrix,peaks,num_p]=blob_segmentation(im,particleIDprops.v);
        
    case {'dynamic'}
        [p_matrix,peaks,num_p] = dynamic_threshold_segmentation_v3(im,...
            particleIDprops.v,particleIDprops.contrast_ratio);
        
    case {'combined'}
        [p_matrix,peaks,num_p]=combined_partID(im,particleIDprops.v);
    otherwise
        error('Unknown ID segmentation method\n')
end
%NOTE - currently the speed of each method is as follows:
%   'blob' (256x256) - 0.74sec per image pair - 1676 particles
%   'dynamic (256x256) - 9.07sec per image pair - 3592 particles
%   'combined (256x256) - 11.30sec per image pair - 4490 particles (hmm...)

if ~isempty(particleIDprops.save_dir)
    if exist(particleIDprops.save_dir,'dir')~=7
        fprintf('Making Save Directory %s \t',particleIDprops.save_dir)
        mkdir(particleIDprops.save_dir)
    end
    sname = sprintf('%%0%0.0fd',particleIDprops.Data.imzeros);
    save(fullfile(particleIDprops.save_dir,[particleIDprops.s_name sprintf(sname,s_num)]),...
        'particleIDprops','p_matrix','peaks','num_p');
end
end
