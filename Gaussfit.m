function [x_c,y_c,D,P,E,Meth] = Gaussfit(intmap,method_in,sigma)
%
%[x_c,y_c,D,P,E] = Gaussfit(intmap,method,sigma)
%
%this function estimates a particle's diameter given an input particle
%intensity profile(intmap)
%
%The function calculated the particle centroid and diameter by assuming the
%light scattering profile (particle intensity profile) is Gaussian in
%shape.  Four points are chosen in both the x&y dimension which includes
%the maximum intensity pixel of the particle profile.
%
%intmap   - input particle intensity profile
%method   - switch: = 1 for standard 3-pt or = 2 for standard 4-pt or = 3 for continuous 4-pt
%sigma    - number of standard deviations in one diameter
%
%x_c      - X Centroid
%y_c      - Y Centroid
%D        - Diamter, for 3pt method it will return 2 Diameters (X and Y)
%P        - Peak Value (again for 3pt, 2 numbers are returned
%E        - Eccentricity of the 3pt fit that is always run first

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


%B.Drew  - 7.31.2008   (4-pt Gaussian method taken from the Master's
%                        Thesis of M. Brady)
%S.Raben - 9.20.2008


if strcmpi(method_in,'TPG')
    method = 1;
elseif strcmpi(method_in,'FPG')
    method = 3;
elseif strcmpi(method_in,'CFPG')
    method = 4;
else
% if method <1 || method > 4
    error('Unknown Method in Gaussfit Function')
end
Meth = method;
%Check to make sure intmap has at least four nonzero points
if method == 1
    if nnz(intmap)<3
        x_c = 1;
        y_c = 1;
        D   = NaN;
        P   = NaN;
        E   = NaN;
        fprintf('Not Enough Points for a Three Point Gauss Fit\n')
        return
    end
else 
    if nnz(intmap)<4
        x_c = 1;
        y_c = 1;
        D   = NaN;
        P   = NaN;
        E   = NaN;
        fprintf('Not Enough Points for a Four Point Gauss Fit\n')
        return
    end
end

%find all of the max intensity locations
[r c]=find(max(intmap(:)) == intmap);


[x_ct,y_ct,Dt,Pt,E] = threeptgaussfit(intmap,[r c],sigma);
Meth = method;
if sum(isnan(Dt)) > 0
    method = 0;
    Meth = method;
end

if method > 2
    [x_c,y_c,D,P,Meth] = fourptgaussfit(intmap,[r c],sigma,method);

    if sum(isnan([x_c,y_c,D,P]))>0
        x_c = x_ct;
        y_c = y_ct;
        D = Dt;
        P = Pt;
        Meth = 1;
    end
else
    x_c = x_ct;
    y_c = y_ct;
    D = Dt;
    P = Pt;
end
end
