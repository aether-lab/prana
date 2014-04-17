function [Uf,Vf]=VELfilt(U,V,H,C)
H(H==0) = [];
% --- Velocity Smoothing Subfunction ---
H = permute(H(:,:,sum(H(1,1,:)~=0)),[2 3 1]);
%2D gaussian filtering
A=fspecial('gaussian',[H(1) H(2)],C);
%there seems to be a very strange bug in imfilter for when U is single and
%small in size (<18?).  Casts are very fast, and forcing it to double fixes 
%the problem.
Uf=cast(imfilter(double(U),A,'replicate'),class(U));
Vf=cast(imfilter(double(V),A,'replicate'),class(V));

end