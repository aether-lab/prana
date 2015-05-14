function [u,v,eval,c,d]=vectorform(x,y,U,V,Eval,C,D)
% --- Matrix to Vector Subfunction ---

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


x=x(:);y=y(:);
%find unique x and y grid points
a=sort(unique(x));
b=sort(unique(y));
N=length(x(:));

imClass = 'double';

%initialize vectors
S=size(x(:));
u    = zeros(S,imClass);
v    = zeros(S,imClass);
eval = zeros(S,imClass);
if ~isempty(C)
    c = zeros(S,imClass);
    d = zeros(S,imClass);
else
    c = [];
    d = [];
end

%generate data vectors where data is available
for n=1:N
    I=find(b==y(n));
    J=find(a==x(n));
    u(n)    = U(I,J);
    v(n)    = V(I,J);
    eval(n) = Eval(I,J);
    if ~isempty(C)
        c(n)    = C(I,J);
        d(n)    = D(I,J);
    end
end

end
