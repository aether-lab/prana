function [p_matrix,peaks,num_p]=combined_partID(im,v)
%
%This function uses a combination of the 'blob' and 'dynamic threshold
%segmentation' algorithms in an attempt to reduce the computational cost of
%the dynamic while still retaining its ability to accuartely segment 
%overlaped particles.  Special care must be taken in the intital
%thresholding of the image to partially segment the blobs without removing
%too many of the low intensity particle images
%
%(beta) N.Cardwell - 10.28.2009

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


%call the mapparticles subfunction ('blob' particle ID analysis)
[p_matrix,peaks,num_p]=mapparticles_v3(im,v);

%preallocate arrays
p_matrix_new=zeros(size(p_matrix));  num_p_new=0;  peaks_new=zeros(size(p_matrix));

%main program loop - to segment or not to segment, that is the question!
for i=1:num_p
    blob_i=double(p_matrix==i);
    
    %use 'regionprops to det. the extent of blob_i and segment the blob out
    STATS = regionprops(blob_i,'BoundingBox');
    B_Box=ceil(STATS.BoundingBox);
    if B_Box(3)==1;  width=0;  else width=B_Box(3)-1;  end
    if B_Box(4)==1;  height=0;  else height=B_Box(4)-1;  end
    part_i=im( B_Box(2):B_Box(2)+height , B_Box(1):B_Box(1)+width );
    part_i(part_i<v)=0;
    
    %determine if blob_i has multiple peaks, if so then pass it to the
    %dynamic function for additional segmentation, otherwise keep the
    %particle and move on
    BW_max=imregionalmax(part_i);
    if nnz(BW_max)~=1
        %crop out the blob and segment it
        im_crop=im(B_Box(2):B_Box(2)+height , B_Box(1):B_Box(1)+width);
        [p_matrix_crop,peaks_crop,num_p_crop]=dynamic_threshold_segmentation_v3(im_crop,0,0);
        
        %reset the 'particle number' of the segmented blod so it can be
        %directly inserted into the 'new' parameters
        for j=1:num_p_crop
            num_p_new=num_p_new+1;
            peaks_crop(peaks_crop==j) = num_p_new;
            p_matrix_crop(p_matrix_crop==j) = num_p_new;
        end
        
        %insert the segmented/cropped particle into the 'new' parameters
        peaks_new( B_Box(2):B_Box(2)+height ,...
            B_Box(1):B_Box(1)+width ) = peaks_crop;
        p_matrix_new( B_Box(2):B_Box(2)+height ,...
            B_Box(1):B_Box(1)+width ) = p_matrix_crop;
    else
        %assign the particle to the 'new' parameters
        num_p_new=num_p_new+1;
        [r_peak,c_peak]=find(BW_max==1);
        peaks_new(B_Box(2)+(r_peak-1) , B_Box(1)+(c_peak-1)) = num_p_new;
        p_matrix_new(B_Box(2):B_Box(2)+height , B_Box(1):B_Box(1)+width) = num_p_new;
    end
end

%reassign output variables
num_p=num_p_new;
peaks=peaks_new;
p_matrix=p_matrix_new;

end
