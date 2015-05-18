function F=poly_3xy_123z_2eqns(x,alldata)
% F=poly_3xy_123z_2eqns(x,alldata)
% this function solves for the xy object coordinates with input
% image coordiantes alldata.XYpoint.  the resulting x vector contains
% the (x y) object coords.  This is for S-PIV so the z coord. is 0.

% This function is called by reconstructvectorsfun.m

% Writen by M. Brady
% Edited and Commented by S. Raben

%     This file is part of prana, an open-source GUI-driven program for
%     calculating velocity fields using PIV or PTV.
%
%     Copyright (C) 2014  Virginia Polytechnic Institute and State
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

aX=alldata.aX;
aY=alldata.aY;
orderz=alldata.orderz;
XYpoint=alldata.XYpoint;

if orderz==1                % cubic xy, linear z
    polylist=[1 x(1) x(2) 0 x(1)^2 x(1)*x(2) x(2)^2 0  0 x(1)^3 x(1)^2*x(2) x(1)*x(2)^2 x(2)^3 0 0 0]';
    Fpoly=[aX*polylist;aY*polylist]-XYpoint;

elseif orderz==2            % cubic xy, quadratic z
    polylist=[1 x(1) x(2) 0 x(1)^2 x(1)*x(2) x(2)^2 0  0 0 x(1)^3 x(1)^2*x(2) x(1)*x(2)^2 x(2)^3 0 0 0 0 0]';
    Fpoly=[aX*polylist;aY*polylist]-XYpoint; 
    
else             % camera pinhole
    
end

F=Fpoly;
end
