function [var1_est,var2_est,var3_est]=twoVar_spatial_weighting(r_weight,edgeval,loc,data)
%
% [var1_est,var2_est]=twoVar_spatial_weighting(r_weight,edgeval,loc,data)
%
% PROGRAM DESCRIPTION
% This function provides, for a given location, a weigthed estimate on 
% two variables based on the surrounding data, which may be random or 
% structured.  The weighting is based on a Gaussian distribution, of which 
% the user can vary.
%
% INPUTS
%   r_weight - radius of the spatial weighting function
%   edgeval - value of the spatial weighting function at r_weight
%       (between 0-1, but not equal to either 0 or 1); lower values weight the
%       center higher, higher values weight everything more evenly)
%   loc - location predicted point [x,y]
%   data - surrounding points and their values (two-variable)
%       [x,y,var1,var2]
%
% OUTPUT
%   var1_est - predicted value of var1 at 'loc'
%   var2_est - predicted value of var2 at 'loc'
%
%(v1) N.Cardwell - 11.16.2009
%(v2) S. Raben   - 05.25.2010

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


%compute the required standard deviation of the Gaussian curve so that the
%value ar r_weight equals edgeval
std_dev=sqrt(-(r_weight^2)/(2*log(edgeval)));

%compute the Gaussian weighting values at each point in data
gaus_weight_x2 = exp(- ((data(:,1)-loc(1)).^2) / (2*std_dev^2));
gaus_weight_y2 = exp(- ((data(:,2)-loc(2)).^2) / (2*std_dev^2));
gaus_weight_z2 = exp(- ((data(:,3)-loc(3)).^2) / (2*std_dev^2));

%predict the value of both variables at 'loc' (sum of the products over
%the sum of the weights)
var1_est=sum(gaus_weight_x2.*data(:,4))/sum(gaus_weight_x2);
var2_est=sum(gaus_weight_y2.*data(:,5))/sum(gaus_weight_y2);
var3_est=sum(gaus_weight_z2.*data(:,6))/sum(gaus_weight_z2);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%plotting code - COMMENT OUT FOR NORMAL OPERATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% figure;
% quiver(data(:,1),data(:,2),data(:,3),data(:,4),0);  hold on
% set(gca,'DataAspectRatio',[1 1 1]);
% axis([loc(1)-2*r_weight loc(1)+2*r_weight loc(2)-2*r_weight loc(2)+2*r_weight])
% scatter(loc(1),loc(2),'r');
% rectangle('Position',[loc(1)-r_weight,loc(2)-r_weight,r_weight*2,r_weight*2],...
%     'Curvature',[1,1],'LineStyle','--')
% scatter(loc(1)+var1_est,loc(2)+var2_est,'m','+');
% line([loc(1) loc(1)+var1_est],[loc(2) loc(2)+var2_est],'Color','m');
% 
% %compute the unweighted average and overlay on the plot
% avg_U=mean(data(:,3));  avg_V=mean(data(:,4));
% scatter(loc(1)+avg_U,loc(2)+avg_V,'g');
% line([loc(1) loc(1)+avg_U],[loc(2) loc(2)+avg_V],'Color','g');

% %set up the Gaussian weighting scheme
% mean_var=loc(1);
% x = (mean_var-r_weight : 2*r_weight/100 : mean_var+r_weight);
% std_dev=sqrt(-(r_weight^2)/(2*log(edgeval)));
% gaussian_weight_x=exp(- ((x-mean_var).^2) / (2*std_dev^2));
% plot(x,gaussian_weight_x.*r_weight,'b');

end
