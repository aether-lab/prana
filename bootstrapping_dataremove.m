function [M1] = bootstrapping_dataremove(DSIZE,ENUM,MASK)
% --- Bootstrapping Data Removal ---

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


Nx = DSIZE(1);
Ny = DSIZE(2);
Nt = 1;

M1   = zeros(DSIZE);
RMAT = rand(Nx,Ny,Nt);
EN   = 0;

while sum(M1(:))/(Nx*Ny) < ENUM && EN < 1
    M1 = RMAT<EN;
    M1(MASK<0) = 0; 
    EN = EN + 0.005;    
end

M1 = double(M1);

end
