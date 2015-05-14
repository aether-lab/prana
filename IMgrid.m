function [X,Y]=IMgrid(L,S,G,imClass)
% --- Grid Generation Subfunction ---
if nargin<4
    imClass = 'double';
end

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


%grid buffer
if nargin<3
    G=[0 0 0 0];
end

S=[S(2) S(1)];
G=[G(2) G(1) L(1)-G(2)+1 L(2)-G(1)+1];

%form grid
if max(S)==0
    %pixel grid
    y=(1:L(1))';
    x=1:L(2);
else
    if G(1)==0
        %buffers 1/2 grid spacing
        y=(ceil((L(1)-(floor(L(1)/S(1))-2)*S(1))/2):S(1):(L(1)-S(1)))';
    else
        %predefined grid buffer
        y=(G(1):S(1):G(3))';
    end
    if G(2)==0
        %buffers 1/2 grid spacing
        x=ceil((L(2)-(floor(L(2)/S(2))-2)*S(2))/2):S(2):(L(2)-S(2));
    else
        %predefined grid buffer
        x=(G(2):S(2):G(4));
    end
end

%vector2matrix conversion
X=cast(x(ones(length(y),1),:),imClass);
Y=cast(y(:,ones(1,length(x))),imClass);


end
