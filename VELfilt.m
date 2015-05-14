function [Uf,Vf]=VELfilt(U,V,H,C)

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
