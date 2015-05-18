function [mapint,locxy,num_p,mapsize]=mapparticles_v3(im,p_matrix,num_p,p_area)
%
%[mapint,locxy,num_p,mapsize]=mapparticles_v3(im,p_matrix,num_p,p_area)
%
% This function segments individual sections of an image using the original
% image as well as a pre-constructed label matrix.  The function also
% employs a spatial filter (p_area) which allows the user to only return
% sections above a user-defined "size"
%
%INPUTS
%   im (2D array,uint8) - original image in intensity units
%   p_matrix (2D array) - matrix of particle identification and extent
%   num_p (num) - number of identified particles in p_matrix
%   p_area (num) - minimum area in pixels to identify a particle 
%
%OUTPUTS
%   mapint - (2D cell array) individual particle intensity profiles
%   locxy - (2 column array) ROW/COL? location associated with the upper 
%       left pixel of the particle's rectangular projection
%   mapsize - (2 column array) defines the size (row col) of the each
%       array in 'mapint'
%   num_p - number of identified particles
%
%N.Cardwell - 8.13.2009 (v3) - totally revamped the entire particle mapping
%   code; updated to utilize more built-in functions

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


%convert all inputs to single precision
im=single(im);  p_matrix=single(p_matrix);  
num_p=single(num_p);  p_area=single(p_area);

%check each particle to assure that its area is >= 'p_area'; remove if this
%condition is violated and reset the particle index number
% if p_area < 0
%     fprintf('\nWARNING!  Minimum particle identification area is < 0\n');
% elseif p_area > 0
%     STATS=regionprops(p_matrix,'Area');
%     idx=find([STATS.Area] >= p_area);
%     dummy=zeros(size(p_matrix));
%     for i=1:length(idx)
%         dummy=dummy+(p_matrix==idx(i)).*i;
%         %imagesc(dummy,[0 44]);  set(gca,'DataAspectRatio',[1 1 1]);
%         %pause(0.05);
%     end
%     num_p=length(idx);  p_matrix=dummy;
%     clear dummy idx STATS
% end

%determine the bounding box and location for each identified particle
STATS=regionprops(p_matrix,im,'BoundingBox','PixelIdxList','PixelList','Area');

%populate locxy (upper-left pixel of the bounding box)
locxy=zeros(length(STATS),4);
for i=1:length(STATS);  locxy(i,:)=STATS(i,1).BoundingBox(1:4);  end
locxy=ceil(locxy);

%populate mapint (particle confined by the bounding box)
% mapint=cell(1,length(STATS)); 
mapint=cell(1,1); keep=zeros(1,1); c=1;
for i=1:length(STATS)
    %check each particle to assure that its area is >= 'p_area'; remove if this
    %condition is violated and reset the particle index number; Also remove
    %particles that are the size of only 1 pixel.
    if length(STATS(i,1).PixelIdxList)>1 && STATS(i,1).Area >=p_area 
        mapint{1,c}=zeros(locxy(i,4),locxy(i,3));
        for j=1:length(STATS(i,1).PixelIdxList)
            pix_loc_row=STATS(i,1).PixelList(j,2)-locxy(i,2)+1;
            pix_loc_col=STATS(i,1).PixelList(j,1)-locxy(i,1)+1;
            mapint{1,c}(pix_loc_row,pix_loc_col)=im(STATS(i,1).PixelIdxList(j,1));
        end
        keep(c) = i;
        c=c+1;
    end
end

if keep(1) == 0
    mapsize = []; locxy = [];
    num_p = 0;
else
    %populate mapsize (size of each array in mapint) and also trim locxy
    mapsize=[locxy(keep,4),locxy(keep,3)];  locxy=locxy(keep,1:2);
    num_p = length(keep);
end
%code for plotting the results of this function 
%COMMENT OUT FOR NORMAL OPERATION
% figure
% for i=1:length(mapint)
%     subplot(1,2,1)
%     imagesc(im,[0 max(max(im))]);  set(gca,'DataAspectRatio',[1 1 1]); hold on
%     rectangle('Position',[locxy(i,:)-0.5,mapsize(i,2),mapsize(i,1)],'EdgeColor','w')
%     subplot(1,2,2)
%     imagesc(mapint{1,i},[0 max(max(im))]);  set(gca,'DataAspectRatio',[1 1 1]);
%     pause(0.5)
% end

end
