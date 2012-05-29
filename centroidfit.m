function [x_centroid,y_centroid,diameter,I0] = centroidfit(mapint_i,locxy_i)

% 
% [x_centroid,y_centroid,diameter]=centroidfit(mapint_i,locxy_i)
% 
% given an input particle intensity profile, the function computes the x and
% y Intensity Weighted Centroid (IWC) location in index units relative to
% the upper left hand pixel of the entire input image
% 
% function also estimates a particle's diameter given an input particle 
% intensity profile(mapint_i), relative locator to the entire image(locxy_i), 
% and the previously derived particle centroid location(x_centroid,y_centroid).
% 
% The function extimates the diameter by meauring the total intensity of 
% the particle then calculating the portion of that associated with the 
% actual particle (86.4665%) from light scattering theory.  Starting with 
% the pixel closest to the particle center and working outward, the function 
% calculates a running intensity summation until 86.4665% of the total 
% intensity is reached.
% 
% mapint_i - input particle intensity profile
% locxy_i  - index location of the upper left pixel in the particles square
%           projection, used to orient the IWC to the entire image
% 
% N. Cardwell - 2.8.2008
% B. Drew - 7.31.2008

%The first section of the code calcualtes the centroid locations, the
%second calculates the particle diameter

%%%%%%%%%%%%%%%%%Section 1%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%determine extents of mapint_i (particle's square projection)
xintsize=size(mapint_i,2);
yintsize=size(mapint_i,1);

x_c=0;
y_c=0;

%calcualte the IWC x-index location
totalintensity=sum(sum(mapint_i));
for k=1:xintsize
    x_index=(k-0.5)+locxy_i(1,1);%+0.5
    x_c=x_c+sum(mapint_i(:,k))*x_index;
end
x_centroid=x_c/totalintensity+0.5;

%calcualte the IWC y-index location
for k=1:yintsize
    y_index=(k-0.5)+locxy_i(1,2);%+0.5
    y_c=y_c+sum(mapint_i(k,:))*y_index;
end
y_centroid=y_c/totalintensity+0.5;

%%%%%%%%%%%%%%%%%Section 2%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%preallocate distance array of pixel locations and distances from the centroid 
dist_cent=zeros(nnz(mapint_i),4);

%use the "regionprops" function to survey the particle intensity profile
bw_mapint_i=bwlabel(mapint_i);
[mapint_i_info]=regionprops(bw_mapint_i,'PixelIdxList','PixelList');

%populate the distance array
try
    % Check to see if mapint_i is all zeros
    if nnz(mapint_i) == 0
        diameter=NaN;
        I0=NaN;
        %keyboard
        return
    else
        %This loop will not run if size(mapint_i_info,1)>1 (ie - more than one
        %particle in mapint_i)
        for k=1:nnz(mapint_i);
            dist_cent(k,2)=mapint_i_info(1,1).PixelList(k,1)+(locxy_i(1,1)+0.5);%+1.5
            dist_cent(k,3)=mapint_i_info(1,1).PixelList(k,2)+(locxy_i(1,2)+0.5);%+1.5
            dist_cent(k,1)=(dist_cent(k,2)-x_centroid)^2+(dist_cent(k,3)-y_centroid)^2;
            dist_cent(k,4)=mapint_i_info(1,1).PixelIdxList(k);
        end
    end
catch
    %Return NaN's for the diameter and intensity and return to 'detectandsize'
    %    x_centroid=locxy_i(1);
    %    y_centroid=locxy_i(2);
    diameter=NaN;
    I0=NaN;
    %keyboard
    return
end

%sort the distance array from lowest to highest
sort_dist_cent=sortrows(dist_cent);

%calcuate required intensity amount, 86.4665% of total particle intensity
total_int=sum(sum(mapint_i));
req_int=0.864665*total_int;

%preallocate
diamc=0;
cum_int=0;
tag=1;
pp=0;

%sum closest pixel intensities until req_int is reached, also taken into
%account the "fraction pixel amount" needed to reach req_int(line 55)
while tag
    pp=pp+1;
    pixel_int=mapint_i(sort_dist_cent(pp,4));
    cum_int=cum_int + pixel_int;
    if cum_int<req_int
        diamc=diamc+1;
    else
        diamc=diamc+(req_int-(cum_int-pixel_int))/pixel_int;
        tag=0;
    end
end

%Uses the area of a circle to approximate the calculate the particle
%diameter.  If you particles are non-spherical, change this equation to
%ellipse/rectangle/whatever is most appropriate.  "regionprops" can calulate 
%other parameters to help with non=spherical particle diameter calculation, 
%such as orientation and the ratio of major and minor axis
diameter=2*sqrt(diamc/pi);

%The maximum intensity is currently set to the intensity of the brightest
%pixel. This could be modified to calculate the actual intensity if
%desired.
I0=max(max(mapint_i));
end
