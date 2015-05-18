function [tracks]=weighted_nearest_neighbor3D(X1,X2,X2_est,Y1,Y2,Y2_est,Z1,Z2,Z2_est,d1,d2,I1,I2,weights,s_radius)
%
% This function tracks discrete locations (particles) between two successive
% realizations.  Pairing is based on a improved nearest-neighbor algorithm 
% which makes use of additional information about the particles to improve
% the pair-matching effecitiveness, specifically the diameter and peak
% intensity value.  Additionally, the program also makes use of an
% estimation of particle displacement if available
%
% INPUTS
%   (1-D number vectors of equal length)
%   X1,Y1,Z1 - spatial location of the particles in image 1
%   d1,I1 - diameter and maximum intensity of the particles in image 1
%   X2,Y2,Z2 - spatial location of the particles in image 2
%   d2,I2 - diameter and maximum intensity of the particles in image 2
%   X2_est,Y2_est,Z2_est - estimated particle location in image 2 
%
%   weights - [0-1 0-1 0-1] sets the relative emphasis for the pair
%       matching algorithm: inter-particle distance, size diff, & max-intensity diff
%   s_radius - (num, pixles) 3D search radius to identify particle pairs
%
% OUTPUTS
%   tracks - main output array of the matched particle pairs:
%       [X1 X2 Y1 Y2 Z1 Z2 d1 d2 I1 I2 p#1 p#2 match_probability]
%       *where a lower match_coefficient is better, should vary
%        between 0 and 1*
%
%(v6) N. Cardwell - 7.21.2009
%(v7) S. Raben    - 4.2011

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


%compute relevant parameters for each image
num_p1=length(X1);  num_p2=length(X2);
d_diff_max=max([d1;d2])-min([d1;d2]);
I_diff_max=max([I1;I2])-min([I1;I2]);
% If there is no difference between values (i.e. all the values are the
% same) then set variable to inf so that they are removed from the
% probablity calcuation.
if d_diff_max == 0;
    d_diff_max = inf;
end
if I_diff_max == 0
    I_diff_max = inf;
end

%adjust the particle locations in image 1 to match the estimated locations
%provided in X2_est, Y2_est, and Z2_est (save org. variables for later use)
X1_org=X1;  Y1_org=Y1;  Z1_org=Z1;
X1=X2_est;  Y1=Y2_est;  Z1=Z2_est;

%for each particle in image 1, create a "comparison" array for all particles
%in image 2 contained within the search radius and centered on the particle
%location in image 1
compare=cell(num_p1,1);
for i=1:num_p1
    dX=X2-X1(i,1);  dY=Y2-Y1(i,1);  dZ=Z2-Z1(i,1);
    distance=sqrt(dX.^2+dY.^2+dZ.^2);
    
    %build comparison array using particle whose linear distance is within
    %the user defined search radius (s_radius)
    p_index=find(distance<=s_radius);  compare_i=zeros(length(p_index),3);

    % This statement catches the possibility that no paritlces are with in
    % the search radius.
    if ~isempty(p_index)
        compare_i(:,1)=ones(size(p_index)).*i;  compare_i(:,2)=p_index;
        
        %compute the match probability for each possible pairing
        Prob_D=(distance(p_index)./s_radius).*weights(1);
        Prob_d=(abs(d2(p_index)-d1(i,1))./d_diff_max).*weights(2);
        Prob_I=(abs(I2(p_index)-I1(i,1))./I_diff_max).*weights(3);
        compare_i(:,3)=(Prob_D+Prob_d+Prob_I)./sum(weights);
        
        %populate main 'compare' array
        compare{i}=compare_i;
    end

    clear dX dY dZ distance p_index compare_i Prob_D Prob_d Prob_I
end
clear i
compare = cell2mat(compare);

% Check to make sure that the compare variable is not empty.  It will be
% empty if no particles where paired durring this processes.  The use
% should try increasing there search radius.
if ~isempty(compare)
    %sort the comparison array in ascending order based on the pairing
    %probability
    compare_sort=sortrows(compare,3);
    
    %determine the best particle matches and successively match particles
    %until the max number of pairs has been reached (or no more exist)
    p_pairs=[];
    max_matches=min(num_p1,num_p2);
    found_matches_im1=[];  found_matches_im2=[];
    c=0;  i=1;
    while (c<=max_matches) && (i<=size(compare_sort,1))
        test1=found_matches_im1==compare_sort(i,1);
        test2=found_matches_im2==compare_sort(i,2);
        if (any(test1) || any(test2))
            i=i+1;
        else
            found_matches_im1=vertcat(found_matches_im1,compare_sort(i,1));
            found_matches_im2=vertcat(found_matches_im2,compare_sort(i,2));
            p_pairs=vertcat(p_pairs,compare_sort(i,:));
            i=i+1;
            c=c+1;
        end
    end
    
    %populate the 'tracks' arrays with the following structure:
    %   tracks=[X1, X2, Y1, Y2, Z1, Z2, d1, d2, I1, I2, p#1 p#2 match_probability]
    tracks=[X1_org(p_pairs(:,1)), X2(p_pairs(:,2)), Y1_org(p_pairs(:,1)), Y2(p_pairs(:,2)), ...
        Z1_org(p_pairs(:,1)), Z2(p_pairs(:,2)), d1(p_pairs(:,1)), d2(p_pairs(:,2)), ...
        I1(p_pairs(:,1)), I2(p_pairs(:,2)), p_pairs(:,:)];
else
    tracks = [];
end
end
