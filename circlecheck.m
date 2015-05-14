function [locations] = circlecheck(locations,Icir)

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

[ymaxy,xmaxx] = size(Icir);
C_x = mean(locations(:,2));
C_y = mean(locations(:,1));
Cir = sqrt(((locations(:,1)-C_x).^2)+((locations(:,2)-C_y).^2));

while sum(abs(diff(Cir))) == 0
    if locations(1,1)-1 >= 1 && Icir(locations(1,1)-1,locations(1,2)) ~= 0
        locations(1,1) = locations(1,1)-1;
    elseif locations(2,1)+1 <= ymaxy  && Icir(locations(2,1)+1,locations(2,2)) ~= 0
        locations(2,1) = locations(2,1)+1;
    elseif locations(3,2)-1 >= 1 && Icir(locations(3,1),locations(3,2)-1) ~= 0
        locations(3,2) = locations(3,2)-1;
    elseif locations(4,2)+1 <= xmaxx && Icir(locations(4,1),locations(4,2)+1) ~= 0
        locations(4,2) = locations(4,2)+1;
    else
        %keyboard
    end
    
    C_x = mean(locations(:,2));
    C_y = mean(locations(:,1));
    Cir = sqrt(((locations(:,2)-C_x).^2)+((locations(:,1)-C_y).^2));
end
end
