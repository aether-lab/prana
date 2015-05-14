function [Uf,Vf]=VELfilt(U,V,H,C)

%     This file is part of prana, an open-source GUI-driven program for
%     calculating velocity fields using PIV or PTV.
%
%     Copyright (C) 2012-2014  Virginia Polytechnic Institute and State
%     University
% 
%     Copyright 2014.  Los Alamos National Security, LLC. This material was
%     produced under U.S. Government contract DE-AC52-06NA25396 for Los 
%     Alamos National Laboratory (LANL), which is operated by Los Alamos 
%     National Security, LLC for the U.S. Department of Energy. The U.S. 
%     Government has rights to use, reproduce, and distribute this software.
%     NEITHER THE GOVERNMENT NOR LOS ALAMOS NATIONAL SECURITY, LLC MAKES ANY
%     WARRANTY, EXPRESS OR IMPLIED, OR ASSUMES ANY LIABILITY FOR THE USE OF
%     THIS SOFTWARE.  If software is modified to produce derivative works,
%     such modified software should be clearly marked, so as not to confuse
%     it with the version available from LANL.
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
