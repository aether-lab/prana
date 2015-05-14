function [x_centroid,y_centroid,diameter,I0,Meth] = fourptgaussfit(I,locxy,sigma,method)
%
%[x_centroid,y_centroid,diameter]=fourptgausfit(mapint,locxy_i,method,sigma)
%
%this function estimates a particle's diameter given an input particle
%intensity profile(mapint), relative locator to the entire image(locxy_i),
%
%The function calculated the particle centroid and diameter by assuming the
%light scattering profile (particle intensity profile) is Gaussian in
%shape.  Four points are chosen in both the x&y dimension which includes
%the maximum intensity pixel of the particle profile.
%
%I        - input particle intensity profile
%locxy    - index location of the upper left pixel in the particles square
%           projection, used to orient the IWC to the entire image
%method   - switch: =2 for standard 4-pt or =3 for continuous 4-pt
%sigma    - number of standard deviations in one diameter
%
%B.Drew  - 7.31.2008   (4-pt Gaussian method taken from the Master's
%                        Thesis of M. Brady)
%S.Raben - 9.20.2008

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


% Find extents of I
[ymax_I xmax_I] = size(I);
Meth = method;
% Find the maxium value of I and all of its locations
[max_int,int_L] = max(I(:));%#ok
max_int_row = locxy(:,1);
max_int_col = locxy(:,2);
max_int_locxy(1) = round(median(locxy(:,1)));
max_int_locxy(2) = round(median(locxy(:,2)));

%Pick the four pixels to be used, starting with max_int
points=zeros(4,3);
points(1,:)=[max_int_locxy,max_int];

%         if numel(locxy(:,1)) > 1
%             max_try = max(I(I~=max_int));
%         else
%             max_try = max_int;
%         end
sat_int = max_int;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
I2 = zeros(size(I)+2);
I2(2:end-1,2:end-1) = I;
I = I2;
locxy = locxy + 1;
xmax_I = xmax_I+2;
ymax_I = ymax_I+2;
points(1,1:2) = points(1,1:2) + 1;
max_int_locxy(1) = max_int_locxy(1) + 1;
max_int_locxy(2) = max_int_locxy(2) + 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if min(locxy(:,1)) > 1 && max(locxy(:,1)) < ymax_I-2 && ...
        min(locxy(:,2)) > 1 && max(locxy(:,2)) < xmax_I-2
    %Try different combinations of points to find a suitable set
    %A suitable set of points is one where a circle cannot be drawn
    %which intersects all four, and no pixels are empty or
    %saturated
    
    [TE,BE,LE,RE] = Edgesearch(I,max_int_locxy);
    if numel(locxy(:,1))>1
        if I(TE(1)-1,TE(2))~=0
            [pts] = circlecheck([TE;[TE(1)-1 TE(2)];LE;RE],I);
        elseif I(BE(1)+1,BE(2))~=0
            [pts] = circlecheck([BE;[BE(1)+1 BE(2)];LE;RE],I);
        elseif I(LE(1),LE(2)-1)~=0
            [pts] = circlecheck([TE;[LE(1) LE(2)-1];LE;RE],I);
        elseif I(RE(1),RE(2)+1)~=0
            [pts] = circlecheck([TE;[RE(1) RE(2)+1];LE;RE],I);
        else
            %keyboard
        end
    else
        pts = [TE;locxy;BE;LE];
    end
    points(:,1:2) = pts(:,2:-1:1);
    
    for pp = 1:4
        points(pp,3)   = I(pts(pp,1),pts(pp,2));
    end
    
    if sum(diff(points(:,1))) == 0 || sum(diff(points(:,2))) == 0
        error('Points are in a Staight line')
    end
    
else
    %keyboard
    %If max_int is on the edge of the window, then try these point combinations
    if max_int_row==1 && max_int_col==1
        points(2,:)=[points(1,1)+1 , points(1,2)   , I(max_int_row,max_int_col+1)];
        points(3,:)=[points(1,1)   , points(1,2)+1 , I(max_int_row+1,max_int_col)];
        points(4,:)=[points(1,1)   , points(1,2)+2 , I(max_int_row+2,max_int_col)];
    elseif max_int_row==1 && max_int_col==xmax_I
        points(2,:)=[points(1,1)-1 , points(1,2)   , I(max_int_row,max_int_col-1)];
        points(3,:)=[points(1,1)   , points(1,2)+1 , I(max_int_row+1,max_int_col)];
        points(4,:)=[points(1,1)   , points(1,2)+2 , I(max_int_row+2,max_int_col)];
    elseif max_int_row==ymax_I && max_int_col==1
        points(2,:)=[points(1,1)+1 , points(1,2)   , I(max_int_row,max_int_col+1)];
        points(3,:)=[points(1,1)   , points(1,2)-1 , I(max_int_row-1,max_int_col)];
        points(4,:)=[points(1,1)   , points(1,2)-2 , I(max_int_row-2,max_int_col)];
    elseif max_int_row==ymax_I && max_int_col==xmax_I
        points(2,:)=[points(1,1)-1 , points(1,2)   , I(max_int_row,max_int_col-1)];
        points(3,:)=[points(1,1)   , points(1,2)-1 , I(max_int_row-1,max_int_col)];
        points(4,:)=[points(1,1)   , points(1,2)-2 , I(max_int_row-2,max_int_col)];
    elseif max_int_row==1
        points(2,:)=[points(1,1)+1 , points(1,2)   , I(max_int_row,max_int_col+1)];
        points(3,:)=[points(1,1)-1 , points(1,2)   , I(max_int_row,max_int_col-1)];
        points(4,:)=[points(1,1)   , points(1,2)+1 , I(max_int_row+1,max_int_col)];
    elseif max_int_row==ymax_I
        points(2,:)=[points(1,1)+1 , points(1,2)   , I(max_int_row,max_int_col+1)];
        points(3,:)=[points(1,1)-1 , points(1,2)   , I(max_int_row,max_int_col-1)];
        points(4,:)=[points(1,1)   , points(1,2)-1 , I(max_int_row-1,max_int_col)];
    elseif max_int_col==xmax_I
        points(2,:)=[points(1,1)-1 , points(1,2)   , I(max_int_row,max_int_col-1)];
        points(3,:)=[points(1,1)   , points(1,2)+1 , I(max_int_row+1,max_int_col)];
        points(4,:)=[points(1,1)   , points(1,2)-1 , I(max_int_row-1,max_int_col)];
    end
end%end

%If a suitable set of points could not be found:
if numel(max(points(:,3))==sat_int)>1 || min(points(:,3))==0
    x_centroid=locxy(1);
    y_centroid=locxy(2);
    %keyboard
    diameter=NaN;
    I0=NaN;
    return
end

%This code is specific to the standard four point gaussian estimator, but
%will be used to come up with guess values for the continuous 4-point
%method

[Isort,IsortI] = sort(points(:,3),'descend');
points = points(IsortI,:);

x1=points(1,1);
x2=points(2,1);
x3=points(3,1);
x4=points(4,1);
y1=points(1,2);
y2=points(2,2);
y3=points(3,2);
y4=points(4,2);
a1=points(1,3);
a2=points(2,3);
a3=points(3,3);
a4=points(4,3);

alpha(1) = (x4^2)*(y2 - y3) + (x3^2)*(y4 - y2) + ((x2^2) + (y2 - y3)*(y2 - y4))*(y3 - y4);
alpha(2) = (x4^2)*(y3 - y1) + (x3^2)*(y1 - y4) - ((x1^2) + (y1 - y3)*(y1 - y4))*(y3 - y4);
alpha(3) = (x4^2)*(y1 - y2) + (x2^2)*(y4 - y1) + ((x1^2) + (y1 - y2)*(y1 - y4))*(y2 - y4);
alpha(4) = (x3^2)*(y2 - y1) + (x2^2)*(y1 - y3) - ((x1^2) + (y1 - y2)*(y1 - y3))*(y2 - y3);

gamma(1) = (-x3^2)*x4 + (x2^2)*(x4 - x3) + x4*((y2^2) - (y3^2)) + x3*((x4^2) - (y2^2) + (y4^2)) + x2*(( x3^2) - (x4^2) + (y3^2) - (y4^2));
gamma(2) = ( x3^2)*x4 + (x1^2)*(x3 - x4) + x4*((y3^2) - (y1^2)) - x3*((x4^2) - (y1^2) + (y4^2)) + x1*((-x3^2) + (x4^2) - (y3^2) + (y4^2));
gamma(3) = (-x2^2)*x4 + (x1^2)*(x4 - x2) + x4*((y1^2) - (y2^2)) + x2*((x4^2) - (y1^2) + (y4^2)) + x1*(( x2^2) - (x4^2) + (y2^2) - (y4^2));
gamma(4) = ( x2^2)*x3 + (x1^2)*(x2 - x3) + x3*((y2^2) - (y1^2)) - x2*((x3^2) - (y1^2) + (y3^2)) + x1*((-x2^2) + (x3^2) - (y2^2) + (y3^2));

delta(1) = x4*(y2 - y3) + x2*(y3 - y4) + x3*(y4 - y2);
delta(2) = x4*(y3 - y1) + x3*(y1 - y4) + x1*(y4 - y3);
delta(3) = x4*(y1 - y2) + x1*(y2 - y4) + x2*(y4 - y1);
delta(4) = x3*(y2 - y1) + x2*(y1 - y3) + x1*(y3 - y2);

deno = 2*(log(a1)*delta(1) + log(a2)*delta(2) + log(a3)*delta(3) + log(a4)*delta(4));

x_centroid = (log(a1)*alpha(1) + log(a2)*alpha(2) + log(a3)*alpha(3) + log(a4)*alpha(4))/deno;

y_centroid = (log(a1)*gamma(1) + log(a2)*gamma(2) + log(a3)*gamma(3) + log(a4)*gamma(4))/deno;

if a2 ~= a1
    betas = abs((log(a2)-log(a1))/((x2-x_centroid)^2+(y2-y_centroid)^2-(x1-x_centroid)^2-(y1-y_centroid)^2));
elseif a3 ~= a1
    betas = abs((log(a3)-log(a1))/((x3-x_centroid)^2+(y3-y_centroid)^2-(x1-x_centroid)^2-(y1-y_centroid)^2));
elseif a4 ~= a1
    betas = abs((log(a4)-log(a1))/((x4-x_centroid)^2+(y4-y_centroid)^2-(x1-x_centroid)^2-(y1-y_centroid)^2));
else
    %keyboard
    I0 = NaN;
    diameter = NaN;
    return
end

I0=a1/exp(-betas*((x1-x_centroid)^2+(y1-y_centroid)^2));

if sum(isnan([x_centroid,y_centroid,betas])) > 0
    %keyboard
end
%Check solutions for errors - set guess values if the 4-point
%continuous method is being used
if (min([x_centroid,y_centroid,betas])<=0 || max([x_centroid,y_centroid,betas])>=10^4 ...
        || max(isnan([x_centroid,y_centroid,betas]))~=0) && method==4
    %keyboard
    x_centroid=points(1,1);
    y_centroid=points(1,2);
    betas=.25;
    I0=points(1,3);
elseif min([x_centroid,y_centroid,betas])<=0
    %keyboard
    x_centroid=locxy(1);
    y_centroid=locxy(2);
    diameter=NaN;
    I0=NaN;
    return
elseif x_centroid > xmax_I || y_centroid > ymax_I || betas >= 10^3%max([x_centroid,y_centroid,betas])>=10^3
    %keyboard
    x_centroid=locxy(1);
    y_centroid=locxy(2);
    diameter=NaN;
    I0=NaN;
    return
end

%This code is specific to the continuous four point gaussian estimator:
if method==4
    
    %Guess value for fsolve
    x0=[I0,betas,x_centroid,y_centroid];
    
    %Convert row/column measurements to center-of-pixel measurements
    points(:,1)=points(:,1)+0.5;
    points(:,2)=points(:,2)+0.5;
    options=optimset('MaxIter',400,'MaxFunEvals',2000,'LargeScale','off','Display','off');
    try
        %Call the fsolve function to find solutions to the equation in the
        %Powell_optimized_eq function
        [xvar,fval,exitflag]=fsolve(@Powell_optimized_eq,x0,options,points);%#ok
        I0=xvar(1);
        betas=xvar(2);
        x_centroid=xvar(3)-1;
        y_centroid=xvar(4)-1;
    catch %#ok
        %keyboard
        x_centroid=locxy(1);
        y_centroid=locxy(2);
        diameter=NaN;%#ok
        I0=NaN;
    end
end

%Calculate diameter
diameter=sqrt(sigma^2/(2*betas));
x_centroid = x_centroid - 1;
y_centroid = y_centroid - 1;

% Powel Optimized
    function F = Powell_optimized_eq(x,points)
        %This function is called by fsolve when the continuous four-pt method is
        %selected. X is a matrix containing the initial guesses [I0, betas, x_c,
        %y_c]. Points is a matrix containing the four selected points [x1 y1 I1;
        %...;x4 y4 I4]
        %
        %Adapted from M. Brady's 'fourptintgaussfit'
        %B.Drew - 7.18.2008
        
        xx1=points(1,1);
        xx2=points(2,1);
        xx3=points(3,1);
        xx4=points(4,1);
        yy1=points(1,2);
        yy2=points(2,2);
        yy3=points(3,2);
        yy4=points(4,2);
        aa1=points(1,3);
        aa2=points(2,3);
        aa3=points(3,3);
        aa4=points(4,3);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %For unknown reasons, the fsolve function tries negative values of x(2),
        %causing these equations to return an error. Putting the
        %abs() function in front of all the x(2)'s fixes this error, and fsolve's
        %final x(2) value ends up positive.
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        F = [pi*x(1)/4/x(2)*((erf(sqrt(abs(x(2)))*(xx1-x(3)))-erf(sqrt(abs(x(2)))*(xx1+1-x(3))))*(erf(sqrt(abs(x(2)))*(yy1-x(4)))-erf(sqrt(abs(x(2)))*(yy1+1-x(4)))))-aa1;
            pi*x(1)/4/x(2)*((erf(sqrt(abs(x(2)))*(xx2-x(3)))-erf(sqrt(abs(x(2)))*(xx2+1-x(3))))*(erf(sqrt(abs(x(2)))*(yy2-x(4)))-erf(sqrt(abs(x(2)))*(yy2+1-x(4)))))-aa2;
            pi*x(1)/4/x(2)*((erf(sqrt(abs(x(2)))*(xx3-x(3)))-erf(sqrt(abs(x(2)))*(xx3+1-x(3))))*(erf(sqrt(abs(x(2)))*(yy3-x(4)))-erf(sqrt(abs(x(2)))*(yy3+1-x(4)))))-aa3;
            pi*x(1)/4/x(2)*((erf(sqrt(abs(x(2)))*(xx4-x(3)))-erf(sqrt(abs(x(2)))*(xx4+1-x(3))))*(erf(sqrt(abs(x(2)))*(yy4-x(4)))-erf(sqrt(abs(x(2)))*(yy4+1-x(4)))))-aa4];
        
    end
end
