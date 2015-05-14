function [Top,Bottom,Left,Right] = Edgesearch(Im,locmax)

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

[max_r max_c] = find(max(Im(:))==Im);
maxI          = max(Im(:));
[ymax xmax]   = size(Im);

if min(max_r)==1
    stop_ty = 1;
    Top = [round(mean(max_r)) round(mean(max_c))];
else
    stop_ty = 0;
end

if max(max_r)==ymax
    stop_by = 1;
    Bottom = [round(mean(max_r)) round(mean(max_c))];
else
    stop_by = 0;
end

if min(max_c)==1
    stop_lx = 1;
    Left = [round(mean(max_r)) round(mean(max_c))];
else
    stop_lx = 0;
end

if max(max_c)==xmax
    stop_rx = 1;
    Right = [round(mean(max_r)) round(mean(max_c))];
else
    stop_rx = 0;
end

%2 Logical loops to deal with saturated or equvalent pixels values surrounding
%max_int.  Incrementally moves outward until suitable values are found or
%the extents of I are reached
%stop_x=0;
%stop_y=0;
ii=1;
extent_left   = min(max_c);
extent_right  = max(max_c);
extent_top    = min(max_r);
extent_bottom = max(max_r);

%search in the x-dimension
while stop_lx==0 && extent_left>=1
    %locate x-index and intensity values for the nearest neighbors of max_int
    left_locxy  = [locmax(1) locmax(2)-ii];
    left_int    = Im(left_locxy(1),left_locxy(2));
    %logical check to assure that the three pixels chosen are not equal
    if left_int==maxI
        ii=ii+1;
    else
        stop_lx=1;%#ok
    end
    
    %make sure there are data points to the left & right
    Left  = left_locxy;%#ok
end
ii = 1;
while stop_rx==0 && extent_right<=xmax
    %locate x-index and intensity values for the nearest neighbors of max_int
    right_locxy = [locmax(1) locmax(2)+ii];
    right_int   = Im(right_locxy(1),right_locxy(2));
    %logical check to assure that the three pixels chosen are not equal
    if right_int==maxI
        ii=ii+1;
    else
        stop_rx=1;%#ok
    end
    
    %make sure there are data points to the left & right
    Right = right_locxy;%max(max_c)+i;%#ok
end
ii = 1;
%search in the y-dimension
while stop_ty==0 && extent_top>=1
    %locate y-index and intensity values for the nearest neighbors of max_int
    top_locxy    = [locmax(1)-ii locmax(2)];
    top_int      = Im(top_locxy(1),top_locxy(2));
    %logical check to assure that the three pixels chosen are not equal
    if top_int==maxI
        ii=ii+1;%#ok
    else
        stop_ty=1;%#ok
    end
    
    %make sure there are data points to the top & bottom
    Top    = top_locxy;%min(max_r)-j;
end
ii = 1;
while stop_by==0 && extent_bottom<=ymax
    %locate y-index and intensity values for the nearest neighbors of max_int
    bottom_locxy = [locmax(1)+ii locmax(2)];
    bottom_int   = Im(bottom_locxy(1),bottom_locxy(2));
    %logical check to assure that the three pixels chosen are not equal
    if bottom_int==maxI
        ii=ii+1;%#ok
    else
        stop_by=1;%#ok
    end
    
    %make sure there are data points to the top & bottom
    Bottom = bottom_locxy;%max(max_r)+j;
end
end
