function [x_c,y_c,D,P,E,Meth] = Leastsqrfit(I_in,method_in,sigma_in)
% [x_c,y_c,D,P,E,Meth] = Leastsqrfit(I_in,method_in,sigma_in)
% This function is the front end call to the feature sizing algorithms that
% used a least squares Gaussian approach. The code takes an image of the
% feature along with the some relevant parameters and fits a non-linear
% least squares Gaussian to the surface.  There are two main options for
% the fit, either regular Gaussian or a continous Gaussian.  The continous
% Gaussian not only assumes a Gaussian shape to the object but also assumes
% that the surface is the results area integration (i.e. a pixel on a
% digital camera).  The continous method might not be useful for erregular
% features but has been shown to work well for digital particle images.

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


% Writen by: Sam Raben



% Based on the input information select the method to be used.
if strcmpi(method_in,'LSG')
    method = 3;
elseif strcmpi(method_in,'CLSG')
    method = 4;
else
    error('Unknow Least squares sizing method\n')    
end


% Options for the lsqnonlin solver
options=optimset('MaxIter',1200,'MaxFunEvals',5000,'TolX',5e-6,'TolFun',5e-6,...
    'Display','off','DiffMinChange',1e-7,'DiffMaxChange',1);%,'LevenbergMarquardt','off');%,'LargeScale','off');

% Find the center Max and all of the saturated points
[locxy_in(:,1) locxy_in(:,2)] = find(I_in == max(I_in(:)));
max_locxy_in(1) = round(median(locxy_in(:,1)));
max_locxy_in(2) = round(median(locxy_in(:,2)));

% If there are not enough points don't use the method
if nnz(I_in) - numel(find(I_in == max(I_in(:)))) + 1 < 5
    x_c = 1;
    y_c = 1;
    D   = NaN;
    P   = NaN;
    E   = NaN;
    Meth = 0;
    return
end

%Removes negitive values and puts in zeros
Ils = I_in;
Ils(I_in<0) = 0;

[x_c,y_c,D,P,E,Meth] = Leastsqrmethods(Ils,method,sigma_in,options,locxy_in,max_locxy_in);
end
