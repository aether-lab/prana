function F=poly_3xy_123zmod(a,allx1data,allXdata,orderz)
% function F=poly_3xy_123z(a,data)

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

x1=allx1data(:,1);
x2=allx1data(:,2);
x3=allx1data(:,3);


if orderz==1                % cubic xy, linear z
    Fpoly=a(1) + a(2).*x1 + a(3).*x2 + a(4).*x3 + a(5).*x1.^2 +...
        a(6).*x1.*x2 + a(7).*x2.^2 + a(8).*x1.*x3 + a(9).*x2.*x3 +...
        a(10).*x1.^3 + a(11).*x1.^2.*x2 + a(12).*x1.*x2.^2 +...
        a(13).*x2.^3 + a(14).*x1.^2.*x3 + a(15).*x1.*x2.*x3 +...
        a(16).*x2.^2.*x3;

elseif orderz==2            % cubic xy, quadratic z
    Fpoly=a(1) + a(2).*x1 + a(3).*x2 + a(4).*x3 + a(5).*x1.^2 +...
        a(6).*x1.*x2 + a(7).*x2.^2 + a(8).*x1.*x3 + a(9).*x2.*x3 +...
        a(10).*x3.^2 + a(11).*x1.^3 + a(12).*x1.^2.*x2 + a(13).*x1.*x2.^2 +...
        a(14).*x2.^3 + a(15).*x1.^2.*x3 + a(16).*x1.*x2.*x3 +...
        a(17).*x2.^2.*x3 + a(18).*x1.*x3.^2 + a(19).*x2.*x3.^2;    
    
else             % cubic xyz
    
    
end

F=allXdata-Fpoly;

