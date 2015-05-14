function [XYDiameter,mapsizeinfo,locxy] = particle_size_MAIN_V1(im,p_matrix,num_p,sizeprops)
%
% [XYDiameter,mapsizeinfo,locxy] = particle_size_MAIN_V1(im,p_matrix,...
%     num_p,sizeprops);
%
% PROGRAM DESCRIPTION
% This is a MAIN function for the sizing of previously identified particles
% For a given image, this function reads in detected particles and then
% calculates the diameter and centroid location of each particle using 
% various fitting algorithms
%
% INPUTS
%   im - input image
%   p_matrix (2D array) - matrix of particle identification and extent
%   num_p (num) - number of identified particles in p_matrix
%   sizeprops - processing parameters related to particle sizing
%       sizeprops.method - choose sizing method
%           1-Intensity Weighted Centroid
%           2-Three Point Gaussian Method
%           3-Standard Four Point Gaussian
%           4-Continuous Four Point Gaussian
%           5-Least Squares Gaussian Fit
%           6-Continuous Least Squares Gaussian Fit
%       sizeprops.p_area - minimum area (in pixels) to identify a particle
%       sizeprops.sigma - number of standard deviations used to compute the diameter
%           if a gaussian sizing method (ie - anything but the intensity 
%           weighted centroid method) is selected
%       sizeprops.errors - switch: 
%           1-return NaN's if the selected method fails
%           0-return Intensity-Weighted Centroid results if method fails
%       sizeprops.save_dir - saveing filepath; make '0' if not saving
%       sizeprops.slsh - '\' or '/' (PC or Linux)
%       sizeprops.s_name - base saving name
%       sizeprops.s_num - saving number
%
% OUTPUTS
%   XYDiameter = 6-column matrix; 1st column is x-centorid location, 2nd
%       column is y-centroid location, 3rd column is particle diameter, 4th
%       column is true max. intensity, 5th column is particle id#, 6th
%       column is sizing method
%   mapsizeinfo - (2 column array) defines the size (row col) of the each
%       sized particle
%   locxy - (2 column array) ROW/COL? location associated with the upper 
%       left pixel of the particle's rectangular projection
%
%S. Raben - 9.15.08
%(v1) N. Cardwell - 10.29.2009 (now accepts generalized inputs from particleID 
%       methods; also revised header info and added save command)
%(v1.1) N. Cardwell - 11.1.2009 (combined particleIDprops and saveprops into
%       a single input, particleIDprops)
%(v1.2) N.Cardwell - 1.7.2010 (fixed the portions of the code related to
%       saving the IWC estimate whenever the user selected method fails)

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


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%identifies/locates all particles with (mapparticles)
%(adjacent pixels of greater than 0 intensity)
%   mapint - a cellular 2-D array of particle intensity profiles
%       locxy - a 2 column array of location associated with the upper
%       left pixel index of the particle intensity profile's
%       rectangular projection, the format of locxy is column then row 
%       (or x then y)
[mapint,locxy,num_p,mapsizeinfo]=mapparticles_v3(im,p_matrix,num_p,sizeprops.p_area);

%initialize particle properties array
particleprops=zeros(num_p,6);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Geometric Centroid (GEO)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%calculate the Geometeric Centorid (GEO) using (geometric_centroid)
%these results will be outputted if the other methods return NaN's and
%the user does not want errors in the output
%returns    *x_centroid - x (column) index location of the particles IWC
%           *y_centroid - y (row) index location of the particles IWC
%           *diameter   - calculated diameter of the particle
%           *I0         - intensity value of the brightest pixel
if strcmpi(sizeprops.method,'GEO')
    [particleprops(:,1),particleprops(:,2),particleprops(:,3),particleprops(:,4)]=...
        geometric_centroid(p_matrix,im,sizeprops.p_area);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Intensity Weigthed Centorid (IWC)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%calculate the Intensity Weigthed Centorid (IWC) using (centroidfit)
%these results will be outputted if the other methods return NaN's and
%the user does not want errors in the output
%returns    *x_centroid - x (column) index location of the particles IWC
%           *y_centroid - y (row) index location of the particles IWC
%           *diameter   - calculated diameter of the particle
%           *I0         - intensity value of the brightest pixel
c = 1;
if strcmpi(sizeprops.method,'IWC') % sizeprops.method<7
    for i=1:num_p
         [particleprops(c,1),particleprops(c,2),particleprops(c,3),particleprops(c,4)]=...
             centroidfit(mapint{i},locxy(i,:));
        particleprops(c,1)=particleprops(c,1) - 1;
        particleprops(c,2)=particleprops(c,2) - 1;
        particleprops(c,5)=i;
        particleprops(c,6)=1;
        c=c+1;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Three Point Gaussian (TPG)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%calculate the centroid location and particle diameter using the 3-point
%approximation method (threeptgausfit)
%returns    *x_centroid - x (column) index location of the particles IWC
%           *y_centroid - y (row) index location of the particles IWC
%           *diameter   - calculated diameter of the particle
%           *I0         - calculated maximum particle intensity
if strcmpi(sizeprops.method,'TPG')%sizeprops.method==2
    for i=1:num_p
        [x_cg,y_cg,D,I0,~,Meth]=Gaussfit(mapint{i},sizeprops.method,sizeprops.sigma);
        x_c = locxy(i,1)+x_cg-1;
        y_c = locxy(i,2)+y_cg-1;
        
        if nnz(isnan([x_c,y_c,D,I0])>1)
            if sizeprops.errors==0
                %retain the IWC estimate
            else
                particleprops(i,:)=[x_c,y_c,max(D),max(I0),i,Meth+1];
            end
        else
            particleprops(i,:)=[x_c,y_c,max(D),max(I0),i,Meth+1];
        end
    end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Four Point Gaussian (FPG and CFPG)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%calculate the centroid location and particle diameter using the 4-point
%or continuous 4-point approximation method (fourptgausfit)
%returns    *x_centroid - x (column) index location of the particles IWC
%           *y_centroid - y (row) index location of the particles IWC
%           *diameter   - calculated diameter of the particle
%           *I0         - calculated maximum particle intensity
elseif strcmpi(sizeprops.method,'FPG') || strcmpi(sizeprops.method,'CFPG') %sizeprops.method==3 || sizeprops.method==4
    for i=1:num_p
        [x_cg,y_cg,D,I0,~,Meth]=Gaussfit(mapint{i},sizeprops.method,sizeprops.sigma);
        x_c = locxy(i,1)+x_cg-1;
        y_c = locxy(i,2)+y_cg-1;

        if nnz(isnan([x_c,y_c,D,I0])>1)
            if sizeprops.errors==0
                %retain the IWC estimate
            else
                particleprops(i,:)=[x_c,y_c,max(D),max(I0),i,Meth];
            end
        else 
            particleprops(i,:)=[x_c,y_c,max(D),max(I0),i,Meth];
        end
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Least Squares Gaussian (LSG and CLSG)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%calculate the centroid location and particle diameter using the local
%least squares or continuous local least squares method (localleastsquares)
%returns    *x_centroid - x (column) index location of the particles IWC
%           *y_centroid - y (row) index location of the particles IWC
%           *diameter   - calculated diameter of the particle
%           *I0         - calculated maximum particle intensity
elseif strcmpi(sizeprops.method,'LSG') || strcmpi(sizeprops.method,'CLSG')%sizeprops.method==5 || sizeprops.method==6
    for i=1:num_p
        [x_cl,y_cl,D_l,I0_l,~,Meth] = Leastsqrfit(mapint{i},sizeprops.method,sizeprops.sigma);
        x_c = locxy(i,1)+x_cl-1;
        y_c = locxy(i,2)+y_cl-1;
        
        if nnz(isnan([x_c,y_c,D_l,I0_l]))>1
            if sizeprops.errors==0
                %retain the IWC estimate
            else
                particleprops(i,:)=[x_c,y_c,max(D_l),max(I0_l),i,Meth+2];
            end
        else
            particleprops(i,:)=[x_c,y_c,max(D_l),max(I0_l),i,Meth+2];
        end
    end
%     fprintf('\n\t%3.2f%% of particles sized with method %0.0f ',(sum(particleprops(:,6)==sizeprops.method)/num_p)*100,sizeprops.method)
end

XYDiameter=particleprops;

%section for saving sized particle information
if ~isempty(sizeprops.save_dir)
    if exist(sizeprops.save_dir,'dir')~=7
        fprintf('Making Save Directory %s \n',sizeprops.save_dir)
        mkdir(sizeprops.save_dir)
    end
    sname = sprintf('%%s%%0%0.0fd',sizeprops.Data.imzeros);
    save(fullfile(sizeprops.save_dir,sprintf(sname,sizeprops.s_name,sizeprops.s_num)),...
        'sizeprops','XYDiameter','mapsizeinfo','locxy');
end

% if ~isempty(particleIDprops.save_dir)
%     if exist(particleIDprops.save_dir,'dir')~=7
%         fprintf('Making Save Directory %s \t',particleIDprops.save_dir)
%         mkdir(particleIDprops.save_dir)
%     end
%     sname = sprintf('%%0%0.0fd',particleIDprops.Data.imzeros);
%     save(fullfile(particleIDprops.save_dir,[particleIDprops.s_name sprintf(sname,particleIDprops.s_num)]),...
%         'particleIDprops','p_matrix','peaks','num_p');
% end

end
