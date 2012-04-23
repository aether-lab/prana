function [x_centroid,y_centroid,diameter,I0] = geometric_centroid(p_mat,im)
%
%
%
%

R = regionprops(p_mat,im,'Area','Centroid','MaxIntensity');
x_centroid = zeros(length(R),1);
y_centroid = zeros(length(R),1);
diameter = zeros(length(R),1);
I0 = zeros(length(R),1);

for i = 1:length(R)
    x_centroid(i) = R(i).Centroid(1);
    y_centroid(i) = R(i).Centroid(2);
    diameter(i)   = R(i).Area;
    I0(i)         = R(i).MaxIntensity;    
end
end
