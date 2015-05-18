function [X2_est,Y2_est,Z2_est]=PTVestloc_V1(X1,Y1,Z1,PTVprops,completed_tracks)
%
% [X2_est,Y2_est,Z2_est]=PTVestloc_V1(X1,Y1,Z1,PTVprops,completed_tracks)
%
% PROGRAM DESCRIPTION
% This function provides a diplacement estimation for particles in an image
% given previous tracks of the particles path.  A spatially weighted
% average (Gaussian) is used to provide the displacment prediction, based
% on the surrounding tracks.  Two main modes of operation are offered:
% static and dynamic.  The static used a fixed search raduis to search for
% neighboring particles.  The dynamic mode succesively increases the search
% radius until a user-defined number of vectors have been identified, or
% the max number of user-defined iterations has occured.  
% 
% The dynamic is suggested for sparse fields or when the spatial vector
% density varies thoroughout the image.  Otherwise the static is suggested.
%
% INPUTS
%   X1,Y1,Z1 - original particle locations
%   PTVprops (structured array)
%       .predict_mode - 'static' or 'dynamic'
%       .r_weight - radius (pixels) of the initial search window and the
%           gaussian weighting function
%       .edgeval - value of the GWF at r_weight (>0 & <1)
%       .numvecs - number of vectors to satisfy the dynamic mode
%       .max_iterations - prevents infinite runs in dynamic mode
%
% OUTPUTS
%   X2_est,Y2_est,Z2_est - estimated location of the particles
%
%(v1) N.Cardwell - 11.17.2009

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


if isempty(completed_tracks)==1
    %no completed tracks to use for PTV location prediction
    X2_est=X1;   Y2_est=Y1;   Z2_est=Z1;
else
    %completed tracks available for PTV location prediction
    switch lower(PTVprops.predict_mode)
        
        %static method uses a fixed radius to ID neighboring particles
        case {'static'}
            %intialize arrays
            X2_est=zeros(size(X1));  Y2_est=zeros(size(Y1));  Z2_est=zeros(size(Z1));
            dummy=cell2mat(completed_tracks);
            
            %perform static location prediction for each new particle
            for i=1:size(X1,1)
                loc=[X1(i,1),Y1(i,1),Z1(i,1)];

                %determine particles within the search radius 'r_weight'
                dX=dummy(:,1)-loc(1);  dY=dummy(:,3)-loc(2);  dZ=dummy(:,5)-loc(3);
                distance=sqrt(dX.^2+dY.^2+dZ.^2);
                check = (distance <= PTVprops.r_weight & distance ~= 0);

                %if no particles found within r_weight, initialize w/ the
                %origianl location
                if nnz(check)==0
                    X2_est(i,1)=loc(1);  Y2_est(i,1)=loc(2); Z2_est(i,1)=loc(3);

                %otherwise use the spatial weighting function to predict w/    
                else
                    data=zeros(nnz(check),4);
                    data(:,1)=dummy(check,1);  
                    data(:,2)=dummy(check,3);
                    data(:,3)=dummy(check,5);
                    data(:,4)=dummy(check,2)-dummy(check,1);
                    data(:,5)=dummy(check,4)-dummy(check,3);
                    data(:,6)=dummy(check,6)-dummy(check,5);
                    [U_est,V_est]=twoVar_spatial_weighting(...
                        PTVprops.r_weight,PTVprops.edgeval,loc,data);
                    X2_est(i,1)=X1(i,1)+U_est;  Y2_est(i,1)=Y1(i,1)+V_est;
                end
            end
%             figure; quiver(dummy(:,1),dummy(:,3),dummy(:,2)-dummy(:,1),dummy(:,4)-dummy(:,3),0,'Color','g');
%             hold on;  quiver(X1,Y1,X2_est-X1,Y2_est-Y1,0,'Color','r');
%             set(gca,'DataAspectRatio',[1 1 1]);

        %dynamic method varies to ID radius to get a specific # of part
        case {'dynamic'}
            %intialize arrays
            X2_est=zeros(size(X1));  Y2_est=zeros(size(Y1));  Z2_est=Z1;
            dummy=cell2mat(completed_tracks);

            %perform static location prediction for each new particle
            for i=1:size(X1,1)
                loc=[X1(i,1),Y1(i,1),Z1(i,1)];
                dX=dummy(:,1)-loc(1);  dY=dummy(:,3)-loc(2);  dZ=dummy(:,5)-loc(3);
                distance=sqrt(dX.^2+dY.^2+dZ.^2);
                
                %determine particles within the search radius 'r_weight',
                %increase the r_weight until numvecs is satified
                number_o_vecs=0;  count=1;  r_weight_temp=PTVprops.r_weight;
                while (number_o_vecs<PTVprops.numvecs && count<PTVprops.max_iterations)
                    check = (distance <= r_weight_temp & distance ~= 0);
                    number_o_vecs=nnz(check);
                    if number_o_vecs<PTVprops.numvecs
                        r_weight_temp=r_weight_temp+1;
                    end
                    count=count+1;
                end
                data=zeros(nnz(check),6);
                data(:,1)=dummy(check,1);  data(:,2)=dummy(check,3);  data(:,3)=dummy(check,5); 
                data(:,4)=dummy(check,2)-dummy(check,1);
                data(:,5)=dummy(check,4)-dummy(check,3);
                data(:,6)=dummy(check,6)-dummy(check,5);
                [U_est,V_est,W_est]=twoVar_spatial_weighting(...
                    r_weight_temp,PTVprops.edgeval,loc,data);
                X2_est(i,1)=X1(i,1)+U_est;  Y2_est(i,1)=Y1(i,1)+V_est; Z2_est(i,1)=Z1(i,1)+W_est;
            end
%             figure; quiver(dummy(:,1),dummy(:,3),dummy(:,2)-dummy(:,1),dummy(:,4)-dummy(:,3),0,'Color','g');
%             hold on;  quiver(X1,Y1,X2_est-X1,Y2_est-Y1,0,'Color','r');
%             set(gca,'DataAspectRatio',[1 1 1]); 
    end
end

end
