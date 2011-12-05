function [Uf,Vf]=VELfilt(U,V,H,C)
% --- Velocity Smoothing Subfunction ---
H = permute(H(:,:,end),[2 3 1]);
%2D gaussian filtering
A=fspecial('gaussian',[H(1) H(2)],C);
Uf=imfilter(U,A,'replicate');
Vf=imfilter(V,A,'replicate');

end