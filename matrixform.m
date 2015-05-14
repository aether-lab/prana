function [X,Y,U,V,Eval,C,D]=matrixform(x,y,u,v,eval,c,d)
% --- Vector to Matrix Subfunction ---

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

imClass = 'double';

%find unique x and y grid points
a=sort(unique(x));
b=sort(unique(y));
N=length(x);

%initialize matrices
U=nan(length(b),length(a),size(u,2));
V=nan(length(b),length(a),size(v,2));
Eval=-1*ones(length(b),length(a),size(eval,2),imClass);

%generate grid matrix
[X,Y]=meshgrid(a,b);

%generate variable matrices (nans where no data available)
for i=1:size(U,3)
    for n=1:N
        I=find(b==y(n));
        J=find(a==x(n));
        U(I,J,i) = u(n,i);
        V(I,J,i) = v(n,i);
        Eval(I,J,i) = eval(n);
    end
end
if ~isempty(c)
    C=nan(length(b),length(a),size(c,2));
    for i=1:size(c,2)
        for n=1:N
            I= b==y(n);
            J= a==x(n);
            C(I,J,i)=c(n,i);
        end
    end
else
    C=[];
end
if ~isempty(d)
    D=nan(length(b),length(a),size(d,2));
    for i=1:size(d,2)
        for n=1:N
            I= b==y(n);
            J= a==x(n);
            D(I,J,i)=d(n,i);
        end
    end
else
    D=[];
end

end
