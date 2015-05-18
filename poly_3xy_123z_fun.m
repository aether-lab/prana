function [Xgrid,Ygrid]=poly_3xy_123z_fun(xgrid,ygrid,orderz,aX,aY)
% [Xgrid Ygrid]=poly_3xy_123z_fun(xgrid,ygrid,orderz,aX,aY)
% 

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

x1=xgrid;
x2=ygrid;
[r,c]=size(xgrid);
x3=zeros(r,c);

if orderz==1                % cubic xy, linear z
    
    Xgrid=aX(1) + aX(2).*x1 + aX(3).*x2 + aX(4).*x3 + aX(5).*x1.^2 +...
        aX(6).*x1.*x2 + aX(7).*x2.^2 + aX(8).*x1.*x3 + aX(9).*x2.*x3 +...
        aX(10).*x1.^3 + aX(11).*x1.^2.*x2 + aX(12).*x1.*x2.^2 +...
        aX(13).*x2.^3 + aX(14).*x1.^2.*x3 + aX(15).*x1.*x2.*x3 +...
        aX(16).*x2.^2.*x3;
    
    Ygrid=aY(1) + aY(2).*x1 + aY(3).*x2 + aY(4).*x3 + aY(5).*x1.^2 +...
        aY(6).*x1.*x2 + aY(7).*x2.^2 + aY(8).*x1.*x3 + aY(9).*x2.*x3 +...
        aY(10).*x1.^3 + aY(11).*x1.^2.*x2 + aY(12).*x1.*x2.^2 +...
        aY(13).*x2.^3 + aY(14).*x1.^2.*x3 + aY(15).*x1.*x2.*x3 +...
        aY(16).*x2.^2.*x3;    

elseif orderz==2            % cubic xy, quadratic z
    
    Xgrid=aX(1) + aX(2).*x1 + aX(3).*x2 + aX(4).*x3 + aX(5).*x1.^2 +...
        aX(6).*x1.*x2 + aX(7).*x2.^2 + aX(8).*x1.*x3 + aX(9).*x2.*x3 +...
        aX(10).*x3.^2 + aX(11).*x1.^3 + aX(12).*x1.^2.*x2 + aX(13).*x1.*x2.^2 +...
        aX(14).*x2.^3 + aX(15).*x1.^2.*x3 + aX(16).*x1.*x2.*x3 +...
        aX(17).*x2.^2.*x3 + aX(18).*x1.*x3.^2 + aX(19).*x2.*x3.^2;   
    
    Ygrid=aY(1) + aY(2).*x1 + aY(3).*x2 + aY(4).*x3 + aY(5).*x1.^2 +...
        aY(6).*x1.*x2 + aY(7).*x2.^2 + aY(8).*x1.*x3 + aY(9).*x2.*x3 +...
        aY(10).*x3.^2 + aY(11).*x1.^3 + aY(12).*x1.^2.*x2 + aY(13).*x1.*x2.^2 +...
        aY(14).*x2.^3 + aY(15).*x1.^2.*x3 + aY(16).*x1.*x2.*x3 +...
        aY(17).*x2.^2.*x3 + aY(18).*x1.*x3.^2 + aY(19).*x2.*x3.^2;    
    
else             % pinhole
    
    
end

