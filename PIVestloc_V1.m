function [X2_est,Y2_est,Z2_est]=PIVestloc_V1(X1,Y1,Z1,PIVprops)
%
% [X2_est,Y2_est,Z2_est]=PIVestloc_V1(X1,Y1,Z1,PIVprops);
%
% PROGRAM DESCRIPTION
% This function provides a diplacement estimation for particles in an image
% given a provided flowfield (typically from PIV).  Processing of PIV 
% flowfields is a separate step and must be completed BEFORE implementing 
% this function.  The user is cautioned to make sure that the supplied PIV 
% velocity field has the same orientation as the identified/sized/and 
% tracked particle fields.  This program assumes that the user will be
% supplying fields from the PIVadvance2 code (written by A.Eckstein) and 
% rotates/transforms the PIV field automatically.  This transformation can
% be edited in the 'PIV loading and transformation block' of the code if 
% desired (specifically line 41).
%
% INPUTS
%   X1,Y1,Z1 - original particle locations
%   PIVprops - loading parameters related to the PIV processed flowfield
%       PIVprops.load_dir - full path of the PIV-file directory
%       PIVprops.load_name - base name of the PIV-file
%       PIVprops.precision - base precision of the PIV-file
%       PIVprops.extension - base extension of the PIV-file (EX: '.plt')
%       PIVprops.frame1 - frame number of the PIV-file used to estimate the
%           particle displacment from SIZE1 to SIZE2
%       PIVprops.slsh - '\' or '/'
%
% OUTPUTS
%   X2_est,Y2_est,Z2_est - estimated location of the particles
%
%(v1) N.Cardwell - 11.1.2009

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%PIV loading and transformation block
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%load in the PIV processed data for the image pair 
%(rearrange coordiante system to match the output of the sizing program)
if any(strcmpi(PIVprops.PIVprops.extension,{'.plt' '.dat'}))
    [X,Y,U,V]=read_pltmod_NC(PIVprops.Data.PIV.Data.outbase,PIVprops.Data.Track.PIVprops.load_dir,...
        PIVprops.PIVprops.frame1,PIVprops.PIVprops.frame1,PIVprops.Data.imzeros);
    PIV1={X,Y,flipud(U),-1.*flipud(V)};  clear X Y U V
elseif strcmpi(PIVprops.PIVprops.extension,'.mat')
    lname = sprintf('%%s%%s%%0%0.0fd',PIVprops.Data.imzeros);
    PIV1t = load(sprintf(lname,PIVprops.Data.Track.PIVprops.load_dir,PIVprops.PIVprops.outbase,PIVprops.PIVprops.frame1,PIVprops.PIVprops.extension));
%     PIV1 = {PIV1t.flowvarunsteady(end:-1:1,:,1,1),PIV1t.flowvarunsteady(end:-1:1,:,1,2),...
%         PIV1t.flowvarunsteady(end:-1:1,:,1,3),PIV1t.flowvarunsteady(end:-1:1,:,1,4)};
    PIV1 = {PIV1t.X(end:-1:1,:),PIV1t.Y(end:-1:1,:),PIV1t.U(end:-1:1,:),PIV1t.V(end:-1:1,:)};

else
    error('Unknown PIV extension')
end

%interpolate the PIV flowfield at each particle location in im1
UI=interp2(PIV1{1},PIV1{2},PIV1{3},X1,Y1,'spline');
VI=interp2(PIV1{1},PIV1{2},PIV1{4},X1,Y1,'spline');

check=isnan(UI); %replace NaNs with the original particle location (occurs when interp2 fails)
X2_est=zeros(size(X1));  Y2_est=zeros(size(Y1));  Z2_est=zeros(size(Z1));
for j=1:length(UI)
    if check(j)==1
        X2_est(j)=X1(j);  Y2_est(j)=Y1(j);  Z2_est(j)=Z1(j)+0;
    else
        X2_est(j)=X1(j)+UI(j);  Y2_est(j)=Y1(j)+VI(j);  Z2_est(j)=Z1(j)+0;
    end
end

% %create overlaid quiver plot of the PIV and interpolated velocity field
% figure;  quiver(PIV1{1},PIV1{2},PIV1{3},PIV1{4},0);
% set(gca,'DataAspectRatio',[1 1 1]);
% set(gca,'YDir','reverse');
% hold on
% quiver(X1,Y1,UI,VI,0,'r');
% axis([1 64 1 64])
% 
% %create scatter plot of position and estimated position
% figure;  scatter(X1,Y1,'.r');
% set(gca,'DataAspectRatio',[1 1 1]);
% %set(gca,'YDir','reverse');
% hold on
% scatter(X2_est,Y2_est,'+b');
% scatter(X2,Y2,'xg');
% legend({'loc1','loc2'})
% axis([1 64 1 64])

end
