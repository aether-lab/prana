function [R]=UOD_sub(W,p,q)
% --- Universal Outlier Detection Algorithm ---

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

%minimum variance assumption
e=0.1; 

%remove value from query point
x=W(p,q);
W(p,q)=nan;

%remove any erroneous points
P=W(:);
Ps = sort(P);
Psfull = Ps(~isnan(Ps));
N=length(Psfull);

if N<=floor(length(W)/3)
    %return negative threshold value if no valid vectors in block
    R = inf;
else
    %return the median deviation normalized to the MAD
    if mod(N,2)==0
        M = (Psfull(N/2)+Psfull(N/2+1))/2;
        MADfull = sort(abs(Psfull-M));
        Q = (MADfull(N/2)+MADfull(N/2+1))/2;
        R = abs(x-M)/(Q+e);
    else
        M = Psfull((N+1)/2);
        MADfull = sort(abs(Psfull-M));
        Q = MADfull((N+1)/2);
        R = abs(x-M)/(Q+e);
    end
end

end
