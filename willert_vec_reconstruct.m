function willert_vec_reconstruct(diroutlist,caldata,dewarp_grid,scaling,pulsesep)
% Performs geometric reconstruction on vector fields extracted using the
% Willert method (1997).  Method has also been modified to include angle
% calculations based on Scarano 2005 paper.
%
% Initialization of code requires a structured format containing
% predetermined factors from preceding stereo-piv steps
% to perform geometric reconstruction.
% (1) Code begins by loading in vector fields from each camera
% per time instance.
% (2) Mutual information region (overlap region between images) is
% calculated from de-warp image information.
% (3) Region is used to align the vector field data to ensure that fields
% overlap without production of addtional error.
% (4) Camera angles are calculated and geometric reconstruction is
% performed to obtain 2-D, 3 component flow field.
%
% Code written by Brett Meyers on June 16 2013
% Code modified by Brett Meyers on August 8 2013
%
% Reference Papers:
%
% Willert, C, "Stereoscopic digital particle image velocimetry for
% application in  tunnel flows", Meas. Sci. and Tech. 1997(8): pp 1465-1479
%
% Scarano, F, et al., "S-PIV comparative assessment: image dewarping +
% misalignment correction and pinhole + geometric back projection", Exp.
% Fluids 2005(39): 257-266

% Number of zeros in data vector name
% Original hard code time scalar (obsolete)
% Same as above

%     This file is part of prana, an open-source GUI-driven program for
%     calculating velocity fields using PIV or PTV.
%
%     Copyright (C) 2014  Virginia Polytechnic Institute and State
%     University
% 
%     Copyright 2014.  Los Alamos National Security, LLC. This material was
%     produced under U.S. Government contract DE-AC52-06NA25396 for Los 
%     Alamos National Laboratory (LANL), which is operated by Los Alamos 
%     National Security, LLC for the U.S. Department of Energy. The U.S. 
%     Government has rights to use, reproduce, and distribute this software.
%     NEITHER THE GOVERNMENT NOR LOS ALAMOS NATIONAL SECURITY, LLC MAKES ANY
%     WARRANTY, EXPRESS OR IMPLIED, OR ASSUMES ANY LIABILITY FOR THE USE OF
%     THIS SOFTWARE.  If software is modified to produce derivative works,
%     such modified software should be clearly marked, so as not to confuse
%     it with the version available from LANL.
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


%keyboard;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Dimension units determined from code during de-warp processing %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t = 1/pulsesep;
xscale=scaling.xscale;
yscale=scaling.yscale;

xgridc=dewarp_grid.xgrid;
ygridc=dewarp_grid.ygrid;
%zgrid = zeros(size(xgrid));

dir_struct1= dir(fullfile(diroutlist.willert2dcam1,['*.' 'mat']));
flname1={dir_struct1.name}';
dir_struct2= dir(fullfile(diroutlist.willert2dcam2,['*.' 'mat']));
flname2={dir_struct2.name}';

nof=length(flname1);

% Finals will not ALWAYS have "pass" in the name.  This could be a problem.
foutnamelist=regexp(flname1,'pass','split');

% Predefining variables helps things run faster.
vectorlist = cell(nof,1);

%keyboard;
for j=1:nof
    
    vectorlist{j}=[{fullfile(diroutlist.willert2dcam1,flname1{j})};{fullfile(diroutlist.willert2dcam2,flname2{j})}];
    
    vecfr1 = load(vectorlist{j}{1});
    if j==1 || (j>1 && ~strcmp(foutnamelist{j}{2}(1),foutnamelist{j-1}{2}(1)))
        %keyboard;
        % %this is the same error as in selfcalibration.m, fixed the same way        
        % X1 = xgridc(vecfr1.Y(:,1),vecfr1.X(1,:));
        % Y1 = ygridc(vecfr1.Y(:,1),vecfr1.X(1,:));
        X1 = xscale * (vecfr1.X+0.5 - 1) + min(xgridc(:));
        Y1 = yscale * (vecfr1.Y+0.5 - 1) + min(ygridc(:));
    end                   % Assign Camera 1 grids locs to local variable
    
    u1 = (xscale*t)*vecfr1.U(:,:,1);
    v1 = (yscale*t)*vecfr1.V(:,:,1);
    clear vecfr1;
    
    vecfr2 = load(vectorlist{j}{2});
    if j==1 || (j>1 && ~strcmp(foutnamelist{j}{2}(1),foutnamelist{j-1}{2}(1)))
        % X2 = xgridc(vecfr2.Y(:,1),vecfr2.X(1,:));
        % Y2 = ygridc(vecfr2.Y(:,1),vecfr2.X(1,:));
        X2 = xscale * (vecfr2.X+0.5 - 1) + min(xgridc(:));
        Y2 = yscale * (vecfr2.Y+0.5 - 1) + min(ygridc(:));
    end                   % Assign Camera 1 grids locs to local variable
    
    u2 = (xscale*t)*vecfr2.U(:,:,1);
    v2 = (yscale*t)*vecfr2.V(:,:,1);
    clear vecfr2;
    
    if j==1 || (j>1 && ~strcmp(foutnamelist{j}{2}(1),foutnamelist{j-1}{2}(1)))
        
        [rows,cols]=size(u1);
        
        xgrid=X1;
        ygrid=Y1;
        zgrid=zeros(size(X1));
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Compute gradients of calibration functions
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        aall=[caldata.aXcam1 caldata.aYcam1 caldata.aXcam2 caldata.aYcam2];
        dFdx1=zeros(rows,cols,4);       % the 3rd dimention corresponds to dFdx1 for (X1,Y1,X2,Y2)
        dFdx2=zeros(rows,cols,4);
        dFdx3=zeros(rows,cols,4);
        
        
        if caldata.modeltype==1
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Mapping the camera coord. to the World Coord. using 1sr order z
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            for gg=1:4
                a=aall(:,gg);
                dFdx1(:,:,gg) = a(2) + 2*a(5)*xgrid + a(6)*ygrid + a(8)*zgrid + 3*a(10)*xgrid.^2 + ...
                    2*a(11)*xgrid.*ygrid + a(12)*ygrid.^2 + 2*a(14)*xgrid.*zgrid + a(15)*ygrid.*zgrid;
                
                dFdx2(:,:,gg) = a(3) + a(6)*xgrid + 2*a(7)*ygrid + a(9)*zgrid + a(11)*xgrid.^2 + ...
                    2*a(12)*xgrid.*ygrid + 3*a(13)*ygrid.^2 + a(15)*xgrid.*zgrid + 2*a(16)*ygrid.*zgrid;
                
                dFdx3(:,:,gg) = a(4) + a(8)*xgrid + a(9)*ygrid + a(14)*xgrid.^2 + a(15)*xgrid.*ygrid + a(16)*ygrid.^2;
            end
        elseif caldata.modeltype==2
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Mapping the camera coord. to the World Coord. using 2nd order z
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            for gg=1:4
                a=aall(:,gg);
                dFdx1(:,:,gg) = a(2) + 2*a(5)*xgrid + a(6)*ygrid + a(8)*zgrid + 3*a(11)*xgrid.^2 + 2*a(12)*xgrid.*ygrid + ...
                    a(13)*ygrid.^2 + 2*a(15)*xgrid.*zgrid + a(16)*ygrid.*zgrid + a(18)*zgrid.^2;
                
                dFdx2(:,:,gg) = a(3) + a(6)*xgrid + 2*a(7)*ygrid + a(9)*zgrid + a(12)*xgrid.^2 + 2*a(13)*xgrid.*ygrid + ...
                    3*a(14)*ygrid.^2 + a(16)*xgrid.*zgrid + 2*a(17)*ygrid.*zgrid + a(19)*zgrid.^2;
                
                dFdx3(:,:,gg) = a(4) + a(8)*xgrid + a(9)*ygrid + 2*a(10)*zgrid + a(15)*xgrid.^2 + a(16)*xgrid.*ygrid + ...
                    a(17)*ygrid.^2 + 2*a(18)*xgrid.*zgrid + 2*a(19)*ygrid.*zgrid;
            end
        end
        
        alpha1 = ((dFdx3(:,:,2).*dFdx2(:,:,1))-(dFdx2(:,:,2).*dFdx3(:,:,1)))./(dFdx2(:,:,2).*dFdx1(:,:,1)-dFdx1(:,:,2).*dFdx2(:,:,1));
        beta1  = ((dFdx3(:,:,2).*dFdx1(:,:,1))-(dFdx1(:,:,2).*dFdx3(:,:,1)))./(dFdx1(:,:,2).*dFdx2(:,:,1)-dFdx2(:,:,2).*dFdx1(:,:,1));
        
        alpha2 = ((dFdx3(:,:,4).*dFdx2(:,:,3))-(dFdx2(:,:,4).*dFdx3(:,:,3)))./(dFdx2(:,:,4).*dFdx1(:,:,3)-dFdx1(:,:,4).*dFdx2(:,:,3));
        beta2  = ((dFdx3(:,:,4).*dFdx1(:,:,3))-(dFdx1(:,:,4).*dFdx3(:,:,3)))./(dFdx1(:,:,4).*dFdx2(:,:,3)-dFdx2(:,:,4).*dFdx1(:,:,3));
        
    end
    
    % Display camera angles for reference
%     
%     figure(100); subplot(2,2,1);
%     imagesc(atand(alpha1)); colorbar; %caxis([25 30]);
%     title('Camera 1 Angle \alpha1','FontSize',16)
%     subplot(2,2,2);
%     imagesc(atand(alpha2)); colorbar; %caxis([-30 -25]);
%     title('Camera 2 Angle \alpha2','FontSize',16)
%     subplot(2,2,3);
%     imagesc(atand(beta1)); colorbar; %caxis([-2 2]);(end-(zed+10):end-(zed+5))
%     title('Camera 1 Angle \beta1','FontSize',16)
%     subplot(2,2,4);
%     imagesc(atand(beta2)); colorbar; %caxis([-2 2]);
%     title('Camera 1 Angle \beta1','FontSize',16)
%     
%     keyboard;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Perform Geometric Reconstruction on Camera 1 & 2 PIV Fields   %
    % (The "Bread and Butter" of the Willert Method, the point where%
    % Stereo-PIV comes to be                                        %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %     keyboard
    if min(abs(atand(beta1(:))))< 8 && min(abs(atand(beta2(:))))< 8
        % This is implemented if there is no significant Y-axis tilt
        %         x = (X2.*alpha1-X1.*alpha2)./(alpha1-alpha2);
        x = (X1+X2)/2;
        %
        %         y = (Y1+Y2)/2+((X2-X1)/2).*((beta2-beta1)./(alpha1-alpha2));
        y = (Y1+Y2)/2;
        %
        u = (u2.*alpha1-u1.*alpha2)./(alpha1-alpha2);                       % U velocity reconstruction
        %
        v = (v1+v2)/2+((u2-u1)/2).*((beta2+beta1)./(alpha1-alpha2));        % V velocity reconstruction
        % yg
        w = (u2-u1)./(alpha1-alpha2);                                       % W velocity reconstruction
    elseif min(abs(atand(alpha1(:)))) < 8 && min(abs(atand(alpha2(:)))) < 8
        % This is implemented if there is no significant X-axis tilt
        %         x = (X1+X2)/2+((Y2-Y1)/2).*((alpha2-alpha1)./(beta1-beta2));
        x = (X1+X2)/2;
        %
        %         y = (Y2.*beta1-Y1.*beta2)./(beta1-beta2);
        y = (Y1+Y2)/2;
        %
        u = (u1+u2)/2+((v2-v1)/2).*((alpha2+alpha1)./(beta1-beta2));        % U velocity reconstruction
        %
        v = (v2.*beta1-v1.*beta2)./(beta1-beta2);                           % V velocity reconstruction
        %
        w = (v2-v1)./(beta1-beta2);                                         % W velocity reconstruction
    else
        %         keyboard
        %
        %         x = (X2.*alpha1-X1.*alpha2)./(alpha1-alpha2);
        x = (X1+X2)/2;
        %
        %         y = (Y2.*beta1-Y1.*beta2)./(beta1-beta2);
        y = (Y1+Y2)/2;
        %
        u = (u2.*alpha1-u1.*alpha2)./(alpha1-alpha2);                       % U velocity reconstruction
        %
        v = (v2.*beta1-v1.*beta2)./(beta1-beta2);                           % V velocity reconstruction
        %
        w = (u2-u1)./(alpha1-alpha2);                                       % W velocity reconstruction
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % output the plt file with all components %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    U=u/1000;                                                 % convert from mm/second (mm displacement) to m/sec
    V=v/1000;                                                 % pulsesep is in microsec
    W=w/1000;
    X= x/1000;
    Y= y/1000;
    Z= zeros(size(x));
    
    %foutname=regexp(flname1{j},'pass','split');
    foutname=foutnamelist{j}{2};
    
    stereo_output=fullfile(diroutlist.willert3cfields,['piv_2d3c_cam',num2str(caldata.camnumber(1)),'cam',num2str(caldata.camnumber(2)),'_pass_',foutname]);
    
    save(stereo_output,'X','Y','Z','U','V','W');
    
    fprintf(['stereo frame_pass_',foutname,'  done.\n']);
    %keyboard;
    clear X Y Z U V W;
end
