function [Uf,Vf]=VELfilt(U,V,C)
% --- Velocity Smoothing Subfunction ---

%2D gaussian filtering
A=fspecial('gaussian',[7 7],C);
Uf=imfilter(U,A,'replicate');
Vf=imfilter(V,A,'replicate');