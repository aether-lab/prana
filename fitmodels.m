function [a_cam1, a_cam2, aXcam1, aYcam1, aXcam2, aYcam2, convergemessage] = fitmodels(allx1data,...
    allx2data,allX1data,allX2data,modeltype,optionsls)
% function [a_cam1 a_cam2 aXcam1 aYcam1 aXcam2 aYcam2 convergemessage]=fitcameramodels(allx1data,...
%     allx2data,allX1data,allX2data,modeltype,optionsls)
% This function ....
% Inputs:
%   allx1data:          Location of target markers (x,y,z) in world coord.
%                       for camera 1.
%   allx2data:          Location of target markers (x,y,z) in world coord.
%                       for camera 2.
%   allX1data:          Location of target markers (x,y) in pixel coord.
%                       for camera 1.
%   allX2data:          Location of target markers (x,y) in pixel coord.
%                       for camera 2.
%   Modeltype:          A switch to choose which camera model that will be
%                         used for the 3D fit.  1 = Cubic XY, Linear Z, 2 =
%                         Cubic XY, Quadratic Z, & 3 = Camera Pinhole Model
%                         (DLT).
%   optionsls:          These are the options for the least squares solver.
%                         This variable is determined by optimset
% 
%  Outputs:
%    a_cam1:            Polinomial coeff for camera 1 using the pinhole 
%                         camera model.
%    a_cam2:            Polinomial coeff for camera 2 using the pinhole 
%                         camera model.
%    a_Xcam1:           Polinomial coeff for the X component of camera 1
%                         using the cubic mapping function.
%    a_Ycam1:           Polinomial coeff for the Y component of camera 1
%                         using the cubic mapping function.
%    a_Xcam2:           Polinomial coeff for the X component of camera 2
%                         using the cubic mapping function.
%    a_Ycam2:           Polinomial coeff for the Y component of camera 2
%                         using the cubic mapping function.
%    convergemessage:   Provides a message about the convergence of the 
%                         non-linear least squares solver.  This variable 
%                         has 4 parts, (X and Y for both cameras).  The tag
%                         following the Yes or No is from the lsqnonlin
%                         function built into matlab.  The next to parts
%                         are RMS in pixels and the R value.

% Writen by M. Brady
% Edited and Commented by S. Raben

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


% save('allx1dataFLS.mat','allx1data');
% save('allx2dataFLS.mat','allx2data');
% save('allX1dataFLS.mat','allX1data');
% save('allX2dataFLS.mat','allX2data');

% save('allx1datanew.mat','allx1data');
% save('allx2datanew.mat','allx2data');
% save('allX1datanew.mat','allX1data');
% save('allX2datanew.mat','allX2data');

% save('allx1datagood.mat','allx1data');
% save('allx2datagood.mat','allx2data');
% save('allX1datagood.mat','allX1data');
% save('allX2datagood.mat','allX2data');

%keyboard;
%  allx1data(:,1)=(allx1data(:,1)-mean(allx1data(:,1)))./(1*max(abs(allx1data(:,1))));
%  allx1data(:,2)=(allx1data(:,2)-mean(allx1data(:,2)))./(1*max(abs(allx1data(:,2))));
%  allx1data(:,3)=(allx1data(:,3)-mean(allx1data(:,3)))./(1*max(abs(allx1data(:,3))));
%  
%  allx2data(:,1)=(allx2data(:,1)-mean(allx2data(:,1)))./(1*max(abs(allx2data(:,1))));
%  allx2data(:,2)=(allx2data(:,2)-mean(allx2data(:,2)))./(1*max(abs(allx2data(:,2))));
%  allx2data(:,3)=(allx2data(:,3)-mean(allx2data(:,3)))./(1*max(abs(allx2data(:,3))));
%  
%  allX1data(:,1)=(allX1data(:,1)-mean(allX1data(:,1)))./(1*max(abs(allX1data(:,1))));
%  allX1data(:,2)=(allX1data(:,2)-mean(allX1data(:,2)))./(1*max(abs(allX1data(:,2))));
%   
%  allX2data(:,1)=(allX2data(:,1)-mean(allX2data(:,1)))./(1*max(abs(allX2data(:,1))));
%  allX2data(:,2)=(allX2data(:,2)-mean(allX2data(:,2)))./(1*max(abs(allX2data(:,2))));

 
 
 
optionsls=optimset('MaxIter',30000,'MaxFunEvals',30000,'TolX',1e-11,'TolFun',1e-7,...
        'LargeScale','off','Display','off','Algorithm','levenberg-marquardt');


% Matlab perfers 'Levenberg-Marquardt' over 'levenbergmarquardt' but the
% optimset function hasn't been updated yet.  This should be fixed in later
% versions of matlab
% fprintf('LevenbergMarquart warning turned off in fitcameramodel.m\n')
% warning('off','optim:lsqncommon:AlgorithmConflict')

% Variable Initialization
a_cam1=[];
a_cam2=[];
aXcam1=[];
aXcam2=[];
aYcam1=[];
aYcam2=[];

% This is a switch that will allow for compulation of the camera pinhole
% coeff. (12) at the same time as the other coeff..  Otherwise it split it
% into 3 steps.
% 1 -> for simulatanious compulation 
% 0 -> for split
allsametime=0;

% Initalization of coeff. for the different models.
if modeltype==1
    a=[100; ones(15,1)];        % initial guesss for solver
elseif modeltype==2
    a=[100; ones(18,1)];        % initial guesss for solver
else    % camera pinhole model
    
end

% Defines the percision of the outputs in the fit table in the GUI.
printformat='%5.4f';        % for converge box message
printformat1='%7.6f';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% If either of the cubic XY camera models are going to be used.  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ((modeltype==1) || (modeltype==2))
    
    % Set model type
    alldata.orderz=modeltype;      
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % get coefficients for X, cam1
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    alldata.allxdata=allx1data;      
    alldata.allXdata=allX1data(:,1);
    % fit X for cam. 1 to all x data for cam1
    [aXcam1,resnorm,f_resid,exitflag]=lsqnonlin(@(a) poly_3xy_123z(a,alldata),...
        a,[],[],optionsls);   %#ok<ASGLU>
    rmsX1=sqrt(mean((f_resid).^2));%sqrt(resnorm/length(alldata.allXdata));
    corrcoef=1-sqrt(mean((f_resid./alldata.allXdata).^2));
    if ((exitflag ~= -1) && (exitflag ~= -2))
        msgX1=[' yes (' num2str(exitflag) ')     ' num2str(rmsX1,printformat) '       ' num2str(corrcoef,printformat1)];
    else
        msgX1=[' no (' num2str(exitflag) ')     ' num2str(rmsX1,printformat) '       ' num2str(corrcoef,printformat1)];
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % get coefficients for Y, cam1
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    alldata.allxdata=allx1data;      
    alldata.allXdata=allX1data(:,2);
    % fit Y for cam. 1 to all x data for cam1
    [aYcam1,resnorm,f_resid,exitflag]=lsqnonlin(@(a) poly_3xy_123z(a,alldata),...
        a,[],[],optionsls);     %#ok<ASGLU>
    rmsY1=sqrt(mean((f_resid).^2));%sqrt(resnorm/length(alldata.allXdata));
    corrcoef=1-sqrt(mean((f_resid./alldata.allXdata).^2));
    if ((exitflag ~= -1) && (exitflag ~= -2))
        msgY1=[' yes (' num2str(exitflag) ')     ' num2str(rmsY1,printformat) '       ' num2str(corrcoef,printformat1)];
    else
        msgY1=[' no (' num2str(exitflag) ')     ' num2str(rmsY1,printformat) '       ' num2str(corrcoef,printformat1)];
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % get coefficients for X, cam2
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    alldata.allxdata=allx2data;      
    alldata.allXdata=allX2data(:,1);
    % fit X for cam. 2 to all x data for cam1
    [aXcam2,resnorm,f_resid,exitflag]=lsqnonlin(@(a) poly_3xy_123z(a,alldata),...
        a,[],[],optionsls);   %#ok<ASGLU>  
    rmsX2=sqrt(mean((f_resid).^2));%sqrt(resnorm/length(alldata.allXdata));
    corrcoef=1-sqrt(mean((f_resid./alldata.allXdata).^2));
    if ((exitflag ~= -1) && (exitflag ~= -2))
        msgX2=[' yes (' num2str(exitflag) ')     ' num2str(rmsX2,printformat) '       ' num2str(corrcoef,printformat1)];
    else
        msgX2=[' no (' num2str(exitflag) ')     ' num2str(rmsX2,printformat) '       ' num2str(corrcoef,printformat1)];
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % get coefficients for Y, cam2
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    alldata.allxdata=allx2data;      
    alldata.allXdata=allX2data(:,2);
    % fit Y for cam. 2 to all x data for cam1
    [aYcam2,resnorm,f_resid,exitflag]=lsqnonlin(@(a) poly_3xy_123z(a,alldata),...
        a,[],[],optionsls);   %#ok<ASGLU>  
    rmsY2=sqrt(mean((f_resid).^2));%sqrt(resnorm/length(alldata.allXdata));
    corrcoef=1-sqrt(mean((f_resid./alldata.allXdata).^2));
    if ((exitflag ~= -1) && (exitflag ~= -2))
        msgY2=[' yes (' num2str(exitflag) ')     ' num2str(rmsY2,printformat) '       ' num2str(corrcoef,printformat1)];
    else
        msgY2=[' no (' num2str(exitflag) ')     ' num2str(rmsY2,printformat) '       ' num2str(corrcoef,printformat1)];
    end
    
    
    
else
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % If modeltype =3 (camera pinhole, direct linear transformation)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if allsametime         % this way calculates a1 thru a12 all at once
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Determination of Polynomial Coef. for camera 1 using the Pinhole
        % camera model
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        a=500*ones(11,1);
        alldata.allxdata=allx1data;
        alldata.allXdata=allX1data;
        [a_cam1,resnorm,f_resid,exitflag]=lsqnonlin(@(a) campinhole_LSfun(a,alldata),...
            a,[],[],optionsls);  % fit X for cam. 1 to all x data for cam1
        %[rows cols]=size(alldata.allxdata);
        %rmsX=sqrt(resnorm/rows/2);
        %corrcoef=1-sqrt(mean((f_resid./(reshape(alldata.allXdata,rows*2,1))).^2));
        %corrcoef=1-sqrt(mean((f_resid./([reshape(alldata.allXdata,rows*2,1)])).^2));
        Xd=alldata.allXdata(:,1);
        Yd=alldata.allXdata(:,2);
        xd=alldata.allxdata(:,1);
        yd=alldata.allxdata(:,2);
        zd=alldata.allxdata(:,3);
        ad=a_cam1;
        res1=(ad(1).*xd+ad(2).*yd+ad(3).*zd+ad(4))./(ad(9).*xd+ad(10).*yd+ad(11).*zd+1)-Xd;
        res2=(ad(5).*xd+ad(6).*yd+ad(7).*zd+ad(8))./(ad(9).*xd+ad(10).*yd+ad(11).*zd+1)-Yd;
        rmsX=sqrt(mean([res1;res2].^2));
        corrcoef=1-sqrt(mean([res1./Xd;res2./Yd].^2));
        if ((exitflag ~= -1) && (exitflag ~= -2))
            msgX1=[' yes (' num2str(exitflag) ')     ' num2str(rmsX,printformat) '       ' num2str(corrcoef,printformat1)];
            msgY1=[' yes (' num2str(exitflag) ')     ' num2str(rmsX,printformat) '       ' num2str(corrcoef,printformat1)];
        else
            msgX1=[' no (' num2str(exitflag) ')     ' num2str(rmsX,printformat) '       ' num2str(corrcoef,printformat1)];
            msgY1=[' no (' num2str(exitflag) ')     ' num2str(rmsX,printformat) '       ' num2str(corrcoef,printformat1)];
        end
        a_cam1=[a_cam1;1];      % add in a34 (or a12).  this was constrained to be 1 in the solver
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Determination of Polynomial Coef. for camera 2 using the Pinhole
        % camera model
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        a=500*ones(11,1);
        alldata.allxdata=allx2data;
        alldata.allXdata=allX2data;
        [a_cam2,resnorm,f_resid,exitflag]=lsqnonlin(@(a) campinhole_LSfun(a,alldata),...
            a,[],[],optionsls);  % fit X for cam. 1 to all x data for cam1
        %[rows cols]=size(alldata.allxdata);
        %rmsX=sqrt(resnorm/rows/2);
        %corrcoef=1-sqrt(mean((f_resid./(reshape(alldata.allXdata,rows*2,1))).^2));
        %corrcoef=1-sqrt(mean((f_resid./([reshape(alldata.allXdata,rows*2,1)])).^2));
        Xd=alldata.allXdata(:,1);
        Yd=alldata.allXdata(:,2);
        xd=alldata.allxdata(:,1);
        yd=alldata.allxdata(:,2);
        zd=alldata.allxdata(:,3);
        ad=a_cam2;
        res1=(ad(1).*xd+ad(2).*yd+ad(3).*zd+ad(4))./(ad(9).*xd+ad(10).*yd+ad(11).*zd+1)-Xd;
        res2=(ad(5).*xd+ad(6).*yd+ad(7).*zd+ad(8))./(ad(9).*xd+ad(10).*yd+ad(11).*zd+1)-Yd;
        rmsX=sqrt(mean([res1;res2].^2));
        corrcoef=1-sqrt(mean([res1./Xd;res2./Yd].^2));
        if ((exitflag ~= -1) && (exitflag ~= -2))
            msgX2=[' yes (' num2str(exitflag) ')     ' num2str(rmsX,printformat) '       ' num2str(corrcoef,printformat1)];
            msgY2=[' yes (' num2str(exitflag) ')     ' num2str(rmsX,printformat) '       ' num2str(corrcoef,printformat1)];
        else
            msgX2=[' no (' num2str(exitflag) ')     ' num2str(rmsX,printformat) '       ' num2str(corrcoef,printformat1)];
            msgY2=[' no (' num2str(exitflag) ')     ' num2str(rmsX,printformat) '       ' num2str(corrcoef,printformat1)];
        end
        a_cam2=[a_cam2;1];
        
        
    else
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % This way calculates a1 thru a12 on three separate times, they 
        % turn out to be the same
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
        
        a=ones(1,4);
        alldata.allxdata=allx1data;
        alldata.allXdata=allX1data(:,1);
        % fit X for cam. 1 to all x data for cam1
        [a1_cam1,resnorm1,f_resid1,~]=lsqnonlin(@(a) campinholemod_LSfun(a,alldata),...
            a,[],[],optionsls);  
        
        a=ones(1,4);
        alldata.allxdata=allx1data;
        alldata.allXdata=allX1data(:,2);
        % fit X for cam. 1 to all y data for cam1
        [a2_cam1,resnorm2,f_resid2,~]=lsqnonlin(@(a) campinholemod_LSfun(a,alldata),...
            a,[],[],optionsls);  
        
        a=ones(1,4);
        alldata.allxdata=allx1data;
        [rows ~]=size(alldata.allxdata);
        alldata.allXdata=ones(rows,1);
        % fit X for cam. 1 to all z data for cam1
        [a3_cam1,resnorm3,f_resid3,exitflag]=lsqnonlin(@(a) campinholemod_LSfun(a,alldata),...
            a,[],[],optionsls);  
        
        rmsX=sqrt((resnorm1+resnorm2+resnorm3)/3/rows);
        corrcoef=1-sqrt(mean([(f_resid1./(allX1data(:,1))).^2;(f_resid2./(allX1data(:,2))).^2;(f_resid3./(ones(rows,1))).^2]));
        
        msgX1=[' yes (' num2str(exitflag) ')     ' num2str(rmsX,printformat) '       ' num2str(corrcoef,printformat1)];
        msgY1=[' yes (' num2str(exitflag) ')     ' num2str(rmsX,printformat) '       ' num2str(corrcoef,printformat1)];
        
        %[f_resid1 f_resid2 f_resid3]
        
        a=ones(1,4);
        alldata.allxdata=allx2data;
        alldata.allXdata=allX2data(:,1);
        % fit X for cam. 2 to all x data for cam1
        [a1_cam2,resnorm1,f_resid1,~]=lsqnonlin(@(a) campinholemod_LSfun(a,alldata),...
            a,[],[],optionsls); 
        
        a=ones(1,4);
        alldata.allxdata=allx2data;
        alldata.allXdata=allX2data(:,2);
        % fit X for cam. 2 to all y data for cam1
        [a2_cam2,resnorm2,f_resid2,~]=lsqnonlin(@(a) campinholemod_LSfun(a,alldata),...
            a,[],[],optionsls);
        
        a=ones(1,4);
        alldata.allxdata=allx2data;
        [rows ~]=size(alldata.allxdata);
        alldata.allXdata=ones(rows,1);
        % fit X for cam. 2 to all z data for cam1
        [a3_cam2,resnorm3,f_resid3,exitflag]=lsqnonlin(@(a) campinholemod_LSfun(a,alldata),...
            a,[],[],optionsls);  
        
        rmsX=sqrt((resnorm1+resnorm2+resnorm3)/3/rows);
        corrcoef=1-sqrt(mean([(f_resid1./(allX2data(:,1))).^2;(f_resid2./(allX2data(:,2))).^2;(f_resid3./(ones(rows,1))).^2]));
        
        msgX2=[' yes (' num2str(exitflag) ')     ' num2str(rmsX,printformat) '       ' num2str(corrcoef,printformat1)];
        msgY2=[' yes (' num2str(exitflag) ')     ' num2str(rmsX,printformat) '       ' num2str(corrcoef,printformat1)];
        
        % [f_resid1 f_resid2 f_resid3]
        
        % Combine the camera coef for each camera
        a_cam1=[a1_cam1; a2_cam1; a3_cam1];
        a_cam2=[a1_cam2; a2_cam2; a3_cam2];
    end
end
convergemessage={msgX1;msgY1;msgX2;msgY2};

if modeltype==1 || modeltype==2
    
    [aXcam1 aYcam1 aXcam2 aYcam2]
    fprintf('\n Approximate Camera Angles: \n Alpha1  Alpha2  Beta1 Beta2\n');
    alpha1=atand((aYcam1(4)*aXcam1(3) - aYcam1(3)*aXcam1(4))/(aYcam1(3)*aXcam1(2) - aYcam1(2)*aXcam1(3)));
    alpha2=atand((aYcam2(4)*aXcam2(3) - aYcam2(3)*aXcam2(4))/(aYcam2(3)*aXcam2(2) - aYcam2(2)*aXcam2(3)));
    beta1=atand((aYcam1(4)*aXcam1(2) - aYcam1(2)*aXcam1(4))/(aYcam1(2)*aXcam1(3) - aYcam1(3)*aXcam1(2)));
    beta2=atand((aYcam2(4)*aXcam2(2) - aYcam2(2)*aXcam2(4))/(aYcam2(2)*aXcam2(3) - aYcam2(3)*aXcam2(2)));
    [alpha1 alpha2 beta1 beta2]
    
elseif modeltype==3
    [a_cam1;a_cam2]
    aXcam1=[a_cam1(1,4) a_cam1(1,1) a_cam1(1,2) a_cam1(1,3)]';
    aYcam1=[a_cam1(2,4) a_cam1(2,1) a_cam1(2,2) a_cam1(2,3)]';
    aXcam2=[a_cam2(1,4) a_cam2(1,1) a_cam2(1,2) a_cam2(1,3)]';
    aYcam2=[a_cam2(2,4) a_cam2(2,1) a_cam2(2,2) a_cam2(2,3)]';
    [aXcam1 aYcam1 aXcam2 aYcam2]
    
    fprintf('\n Approximate Camera Angles: \n Alpha1  Alpha2  Beta1 Beta2\n');
    alpha1=atand((aYcam1(4)*aXcam1(3) - aYcam1(3)*aXcam1(4))/(aYcam1(3)*aXcam1(2) - aYcam1(2)*aXcam1(3)));
    alpha2=atand((aYcam2(4)*aXcam2(3) - aYcam2(3)*aXcam2(4))/(aYcam2(3)*aXcam2(2) - aYcam2(2)*aXcam2(3)));
    beta1=atand((aYcam1(4)*aXcam1(2) - aYcam1(2)*aXcam1(4))/(aYcam1(2)*aXcam1(3) - aYcam1(3)*aXcam1(2)));
    beta2=atand((aYcam2(4)*aXcam2(2) - aYcam2(2)*aXcam2(4))/(aYcam2(2)*aXcam2(3) - aYcam2(3)*aXcam2(2)));
    [alpha1 alpha2 beta1 beta2]
    
    aXcam1=[];aYcam1=[];aXcam2=[];aYcam2=[];
    
end
% save('aXcam1.mat','aXcam1');
% save('aXcam2.mat','aXcam2');
% save('aYcam1.mat','aYcam1');
% save('aYcam2.mat','aYcam2');
% Turning the warning back on
%warning('on','optim:lsqncommon:AlgorithmConflict')
end
