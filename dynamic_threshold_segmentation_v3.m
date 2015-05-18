function [p_matrix,peaks,num_p]=dynamic_threshold_segmentation_v3(im,v1,contrast_ratio)
%
%[p_matrix,peaks,num_p]=dynamic_threshold_segmentation_v3(im,v1,contrast_ratio)
%
%Uses an erosion/dilation process to identify peaks and then determine the
%extents of each particle.  Algorithm is very effective at separating
%overlapped particle images.  However, there is significant increase in
%processing time over the 'blob' method (apprximately 12X).
%
%INPUTS
%   im - image to be segmented (matrix-uint8)
%   v1 - intial threshold value (num)
%   contrast_ratio - may remove this...for now just set to zero
%
%OUTPUTS
%   p_matrix (2D array) - matrix of particle identification and extent
%   peaks (2D array) - matrix of identified image peaks (by image erosion)
%   num_p (num) - number of identified particles in p_matrix
%
%N.Cardwell (v3)   - 8.12.09
%N.Cardwell (v3.1) - 10.12.09 (replaced 'sub2ind' function call with a
%   direct calculation to increase speed)

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


%intially threshold the image using 'threshold'
im_thresh=im;  im_thresh(im_thresh<=v1)=0;

%use built in function to identify regional maxima (i.e peaks)
% tic;  fprintf('Locating image peaks...');
BW_max=imregionalmax(im_thresh);
[p_matrix,num_p]=bwlabel(BW_max,8);
% fprintf('DONE!----')
% fprintf(strcat('elapsed time=',num2str(toc),'seconds\n'));

%set the maximum intensity for each particle
Imax=zeros(num_p,1);
for i=1:num_p
    Imax(i,1)=max(im_thresh(p_matrix==i));
end;  clear i

%perform dilation on each particle until the contrast criterion is met for
%the boundary pixels of each particle; also has a check to make sure that
%the expanding pixels cannot grap "brighter" pixels
% tic;  fprintf('Expanding peaks........');
peaks=p_matrix;
p_matrix_temp=p_matrix;  s=size(p_matrix);
flags=ones(num_p,1);
%figure
while nnz(flags) > 0
%     imagesc(p_matrix_temp);  set(gca,'DataAspectRatio',[1 1 1]);
%     pause(0.5)
    for i=1:num_p

        %check to see if the particle has been flaged (ie. completly expanded)
        if flags(i,1)==1

            %initialize particle conditions
            l_row=size(p_matrix_temp,1);
            [r,c]=find(p_matrix_temp==i);  part_pixels=zeros(length(r),length(c));
            part_pixels(:,1)=r; part_pixels(:,2)=c;
            part_index=(part_pixels(:,2)-1).*l_row+part_pixels(:,1);
%             [part_pixels(:,1),part_pixels(:,2)]=find(p_matrix_temp==i);

%             part_index=find(p_matrix_temp==i);  part_pixels=zeros(length(part_index),2);
            

            %expand all particle pixels by one in each direction
            possible_pixels=[];            
            for j=1:size(part_pixels,1)
                p_pix_j=part_pixels(j,:);
                I_pix_j=im_thresh(p_pix_j(1),p_pix_j(2));
                if p_pix_j(2) > 1
%                     poss_pix=sub2ind(s,p_pix_j(1),p_pix_j(2)-1);
                    poss_pix=((p_pix_j(2)-1)-1)*l_row+p_pix_j(1);
                    if im_thresh(poss_pix) < I_pix_j
                        possible_pixels=[possible_pixels,poss_pix];
                    end
                end
                if p_pix_j(2) < s(2)
%                     poss_pix=sub2ind(s,p_pix_j(1),p_pix_j(2)+1);
                    poss_pix=((p_pix_j(2)-1)+1)*l_row+p_pix_j(1);
                    if im_thresh(poss_pix) < I_pix_j
                        possible_pixels=[possible_pixels,poss_pix];
                    end
                end
                if p_pix_j(1) > 1
%                     poss_pix=sub2ind(s,p_pix_j(1)-1,p_pix_j(2));
                    poss_pix=((p_pix_j(2)-1))*l_row+(p_pix_j(1)-1);
                    if im_thresh(poss_pix) < I_pix_j
                        possible_pixels=[possible_pixels,poss_pix];
                    end
                end
                if p_pix_j(1) < s(1)
%                     poss_pix=sub2ind(s,p_pix_j(1)+1,p_pix_j(2));
                    poss_pix=((p_pix_j(2)-1))*l_row+(p_pix_j(1)+1);
                    if im_thresh(poss_pix) < I_pix_j
                        possible_pixels=[possible_pixels,poss_pix];
                    end
                end
            end

            %remove all non-unique identifications in 'possible_pixels'
            possible_pixels=unique(possible_pixels);

            %check to see if the possible pixels are already part of the
            %particle or another particle (remove if true)
            check=p_matrix_temp(possible_pixels)~=0;
            possible_pixels=possible_pixels(check==0);

            %see if any of the border pixels satisy the contrast criterion
            %(if so then attach to the particle)
            check2=single(im_thresh(possible_pixels))./single(Imax(i,1)) > contrast_ratio;
            if nnz(check2)~=0
                %assign the border pixels to p_matrix final and border_pixels
                p_matrix_temp(possible_pixels(check2))=i;
            else
                flags(i,1)=0;
            end
        end
    end
end
% fprintf('DONE!----')
% fprintf(strcat('elapsed time=',num2str(toc),'seconds\n'));

p_matrix=p_matrix_temp;
end
