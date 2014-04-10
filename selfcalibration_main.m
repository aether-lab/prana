function [caldatamod]=selfcalibration_main(caldata,selfcaljob)

%THis function does self calibration and corrects for any disparity between
%two camera images and the existing calibration

%INPUTS

%caldata is a structure containing the world and image coordinates
%selected during calibration and also the calibration fit coefficients

%caldata.allx1data:- all world coordinates camera1 
%caldata.allx2data:- all world coordinates camera2 
%caldata.allX1data:- all image coordinates camera1 
%caldata.allX2data:- all image coordinates camera2 

%caldata.aXcam1
%caldata.aYcam1
%caldata.aXcam2
%caldata.aYcam2

%Selfcaljob is a structure containing the self calibration Job file

allx1data=caldata.allx1data;
allx2data=caldata.allx2data;
allX1data=caldata.allX1data;
allX2data=caldata.allX2data;

method=caldata.modeltype;
optionsls=caldata.optionsls;

aXcam1=caldata.aXcam1;
aYcam1=caldata.aYcam1;
aXcam2=caldata.aXcam2;
aYcam2=caldata.aYcam2;



%[aXcam1 aYcam1 aXcam2 aYcam2]

%keyboard;


% Doing dewarping and cross correlation to calculate the disparity map (Dux and Duy).
%imagelist=selfcaljob;
[outputdirlist,dewarp_grid,~]=imagedewarp(caldata,'Willert',selfcaljob);

%keyboard;

xgrid1=dewarp_grid.xgrid;
ygrid1=dewarp_grid.ygrid;

job1=selfcaljob;

job1.imdirec=outputdirlist.dewarpdir1;
job1.imdirec2=outputdirlist.dewarpdir2;
job1.imcstep='0';
job1.selfcal=1;
fprintf('Calculating Disparity.\n');
pranaPIVcode(job1);
%keyboard;
istring1=sprintf(['%%s%%s%%0%0.0fd.','mat'],str2double(job1.imzeros));

load(sprintf(istring1,job1.outdirec,[filesep,job1.outputpassbase,'pass',job1.passes,'_'],str2double(job1.imfstart)));

Dux=U(:,:,1);
Duy=V(:,:,1);%X and Y disparity matrices
X3=X;
Y3=Y;% correlation X and Y grid points
clear U V X Y;

%[X3,~,xgrid1,ygrid1,Dux,Duy,Imax1,Jmax1]=disparitycal(cameracal,selfcaljob);
% Dispx=Dux;
% Dispy=Duy;
figure('Name','Disparity Map');quiver(X3,Y3,Dux,Duy,1)
mdx=mean(mean(Dux));%mean x disparity
%figure(22);hist(Dux(:));
mdy=mean(mean(Duy));%mean y disparity
fprintf(['Average X Disparity in Pixels:',num2str(mdx),'\n','Average Y Disparity in Pixels:',num2str(mdy)])

%keyboard;
%calculating the length of the common grid in the object plane on which the
%images are dewarped
xmin=min(min(xgrid1));
xmax=max(max(xgrid1));
ymin=min(min(ygrid1));
ymax=max(max(ygrid1));

[Imax1,Jmax1]=size(xgrid1);
% Scaling factor = object plane grid length in x or y(in physical units) / pixel length of images
scalex=(xmax-xmin)/Jmax1; 
scaley=(ymax-ymin)/Imax1;


%figure(87);imagesc(Dux);title('X direction disparity');xlabel('x(pixels');ylabel('y(pixels)');%caxis on;
%figure(88);imagesc(Duy);title('Y direction disparity');xlabel('x(pixels');ylabel('y(pixels)');%caxis on;

%Making object grids of same size as the grid on which disparity is calculated
[a,b]=size(X3);
xg=xgrid1(Y3(:,1),X3(1,:));
yg=ygrid1(Y3(:,1),X3(1,:));
zgrid=zeros(size(xg));

%Just a check that if maximum absolute x disparity is less than 1 pixel then no need
%to go further.
% if max(max(abs(Dux)))<=1
%     fprintf('Disparity map converged\n');
%     return
% end


 
%scaling the disparity to physical dimensions 
Dux=Dux.*(scalex);
Duy=Duy.*(scaley);
%world grid points shifted by the amount of disparity to get the
%locations at which camera 2 local viewing angle are calculated.

xgrid=xg-Dux./2;
x2grid=xg+Dux./2;
ygrid=yg-Duy./2;
y2grid=yg+Duy./2; 

% ygrid=yg-Duy./2;
% y2grid=yg+Duy./2;
%Doing Triangulation

[z1grid]= geometricTriangulation(xgrid,x2grid,ygrid,y2grid,zgrid,Dux,Duy,aXcam1,aYcam1,aXcam2,aYcam2); %ouput is the projected z world points


x=reshape(xg,a*b,1);y=reshape(yg,a*b,1);z=reshape(z1grid,a*b,1);

%fitting a plane to the projected z points using linear regression

a1=length(x);a2=sum(x);a3=sum(y);a4=sum(x.*y);a5=sum(x.*z);a6=sum(y.*z);a7=sum(x.^2);a8=sum(y.^2);a9=sum(z);

AA=[a7 a4 a2; a4 a8 a3; a2 a3 a1]; bb=[a5;a6;a9];

coeff=AA\bb; % coefficients for the plane Z= Ax+By+C

%newz=coeff(1).*x + coeff(2).*y + coeff(3);
% figure(6);scatter3(x,y,newz);
%[XX,YY]=meshgrid(-6:12/15:6,-6:12/15:6);
%zzgrid=zeros(a,b);
%zprime=reshape(newz,a,b);
%z1grid=zprime;
%figure(7);surf(xgrid,ygrid,zprime);hold on;surf(xgrid,ygrid,zzgrid);title('corrected z plane');xlabel('x(mm)');ylabel('y(mm)');

alpha=atan(coeff(1));beta=atan(coeff(2));                                  % x-z and y-z plane slopes
Roty=[cos(alpha) 0 -sin(alpha);0 1 0;sin(alpha) 0 cos(alpha)];             %rotation about x and y axis
Rotx=[1 0 0;0 cos(beta) -sin(beta);0 sin(beta) cos(beta)];
tz=[0;0;coeff(3)];                                                         %translation along z
l1=length(allx1data(:,1));
l2=length(allx2data(:,1));
%Modifying the calibration plane world coordinates  for camera1 and camera2
%such that the new fitted z plane becomes z=0 or in other words the
%original calibration plane world z coordinates are calculated with respect
%to the new fitted plane.
ztrans1=Roty'*Rotx'*[allx1data(:,1)';allx1data(:,2)';allx1data(:,3)'] - [tz(1).*ones(1,l1);tz(2).*ones(1,l1);tz(3).*ones(1,l1)];
ztrans2=Roty'*Rotx'*[allx2data(:,1)';allx2data(:,2)';allx2data(:,3)'] - [tz(1).*ones(1,l2);tz(2).*ones(1,l2);tz(3).*ones(1,l2)];

%figure(8);scatter3(ztrans1(1,:),ztrans1(2,:),ztrans1(3,:));title('rotated calibration planes');xlabel('x(mm)');ylabel('y(mm)');


%projecting cal points to new z =0 plane

caldatamod.allx1data=ztrans1';
caldatamod.allx2data=ztrans2'; %outputs the modified planes
[~,~, aXcam1, aYcam1, aXcam2, aYcam2, ~]=fitmodels(caldatamod.allx1data,...
    caldatamod.allx2data,allX1data,allX2data,method,optionsls);
caldatamod.aXcam1=aXcam1;
caldatamod.aYcam1=aYcam1;
caldatamod.aXcam2=aXcam2;
caldatamod.aYcam2=aYcam2;
end


function [z1grid]= geometricTriangulation(xgrid,x2grid,ygrid,y2grid,zgrid,Dux,Duy,aXcam1,aYcam1,aXcam2,aYcam2)

%function for calculating local viewing angles and Triangulation
%input calibration matrices and common grid world coordinates

%outputs projected z coordinates.

%called in mainselfcal.m

%drawback does not take into account Duy for triangulation... should be
%modified later on;

%written by Sayantan Bhattacharya on 7/19/2013


[rows,cols]=size(zgrid);
%initializing local viewing angles 
alphatan=zeros(rows,cols,2);betatan=zeros(rows,cols,2);
%while (max(max(abs(dz))))>=th %&& (max(max(dz)))<=0.1

% Viewing angles for camera 1 calculated at gridpoints xgrid,ygrid;

 aall=[aXcam1 aYcam1 aXcam2 aYcam2];
    dFdx1=zeros(rows,cols,4);       % the 3rd dimention corresponds to dFdx1 for (X1,Y1,X2,Y2)
    dFdx2=zeros(rows,cols,4);
    dFdx3=zeros(rows,cols,4);

for gg=1:4
            a=aall(:,gg);
            dFdx1(:,:,gg) = a(2) + 2*a(5).*xgrid + a(6)*ygrid + a(8)*zgrid + 3*a(11)*xgrid.^2 + 2*a(12)*xgrid.*ygrid + ...
                a(13)*ygrid.^2 + 2*a(15)*xgrid.*zgrid + a(16)*ygrid.*zgrid + a(18)*zgrid.^2;

            dFdx2(:,:,gg) = a(3) + a(6)*xgrid + 2*a(7)*ygrid + a(9)*zgrid + a(12)*xgrid.^2 + 2*a(13)*xgrid.*ygrid + ...
                3*a(14)*ygrid.^2 + a(16)*xgrid.*zgrid + 2*a(17)*ygrid.*zgrid + a(19)*zgrid.^2;

            dFdx3(:,:,gg) = a(4) + a(8)*xgrid + a(9)*ygrid + 2*a(10)*zgrid + a(15)*xgrid.^2 + a(16)*xgrid.*ygrid + ...
                a(17)*ygrid.^2 + 2*a(18)*xgrid.*zgrid + 2*a(19)*ygrid.*zgrid;
end


%calculating the viewing angles using formula (7) and (8) from Giordano and Astarita's paper "Spatial resolution of the Stereo PIV technique" (2009)

alphatan(:,:,1)=(((dFdx3(:,:,2).*dFdx2(:,:,1)) - (dFdx2(:,:,2).*dFdx3(:,:,1)))./((dFdx2(:,:,2).*dFdx1(:,:,1)) - (dFdx1(:,:,2).*dFdx2(:,:,1))));
betatan(:,:,1)=(((dFdx3(:,:,2).*dFdx1(:,:,1)) - (dFdx1(:,:,2).*dFdx3(:,:,1)))./((dFdx1(:,:,2).*dFdx2(:,:,1)) - (dFdx2(:,:,2).*dFdx1(:,:,1))));

%aall=[aXcam1 aYcam1 aXcam2 aYcam2];
    dFdx1=zeros(rows,cols,4);       % the 3rd dimention corresponds to dFdx1 for (X1,Y1,X2,Y2)
    dFdx2=zeros(rows,cols,4);
    dFdx3=zeros(rows,cols,4);
%Xx1c1=aXcam1(2) + 2*aXcam1(5).*x1


%Viewing angles for camera 2 calculated at gridpoints x2grid=xgrid+Dux,y2grid=ygrid+Duy
for gg=1:4
            a=aall(:,gg);
            dFdx1(:,:,gg) = a(2) + 2*a(5).*x2grid + a(6)*y2grid + a(8)*zgrid + 3*a(11)*x2grid.^2 + 2*a(12)*x2grid.*y2grid + ...
                a(13)*y2grid.^2 + 2*a(15)*x2grid.*zgrid + a(16)*y2grid.*zgrid + a(18)*zgrid.^2;

            dFdx2(:,:,gg) = a(3) + a(6)*x2grid + 2*a(7)*y2grid + a(9)*zgrid + a(12)*x2grid.^2 + 2*a(13)*x2grid.*y2grid + ...
                3*a(14)*y2grid.^2 + a(16)*x2grid.*zgrid + 2*a(17)*y2grid.*zgrid + a(19)*zgrid.^2;

            dFdx3(:,:,gg) = a(4) + a(8)*x2grid + a(9)*y2grid + 2*a(10)*zgrid + a(15)*x2grid.^2 + a(16)*x2grid.*y2grid + ...
                a(17)*y2grid.^2 + 2*a(18)*x2grid.*zgrid + 2*a(19)*y2grid.*zgrid;
end

alphatan(:,:,2)=(((dFdx3(:,:,4).*dFdx2(:,:,3)) - (dFdx2(:,:,4).*dFdx3(:,:,3)))./((dFdx2(:,:,4).*dFdx1(:,:,3)) - (dFdx1(:,:,4).*dFdx2(:,:,3))));
betatan(:,:,2)=(((dFdx3(:,:,4).*dFdx1(:,:,3)) - (dFdx1(:,:,4).*dFdx3(:,:,3)))./((dFdx1(:,:,4).*dFdx2(:,:,3)) - (dFdx2(:,:,4).*dFdx1(:,:,3))));
%keyboard;

%Display camera angles for reference
    
%         figure(100); subplot(2,2,1);
%         imagesc(atand(alphatan(:,:,1))); colorbar; %caxis([25 30]);
%         title('Camera 1 Angle \alpha1','FontSize',16)
%         subplot(2,2,2);
%         imagesc(atand(alphatan(:,:,2))); colorbar; %caxis([-30 -25]);
%         title('Camera 2 Angle \alpha2','FontSize',16)
%         subplot(2,2,3);
%         imagesc(atand(betatan(:,:,1))); colorbar; %caxis([-2 2]);(end-(zed+10):end-(zed+5))
%         title('Camera 1 Angle \beta1','FontSize',16)
%         subplot(2,2,4);
%         imagesc(atand(betatan(:,:,2))); colorbar; %caxis([-2 2]);
%         title('Camera 1 Angle \beta1','FontSize',16)
%     
    %keyboard;
% Triangulation using formula (1) of that paper

if max(max(abs(alphatan(:,:,1)-abs(alphatan(:,:,2)))))>max(max(abs(betatan(:,:,1)-abs(betatan(:,:,2)))))
    dz1=(-sign(alphatan(1,1)).*Dux)./(abs(alphatan(:,:,1))+abs(alphatan(:,:,2)));
else
    %dz2=Duy./(abs(alphatan(:,:,1))+abs(alphatan(:,:,2)));%(abs(betatan(:,:,1))+abs(betatan(:,:,2)));
    dz1=(sign(betatan(1,1)).*Duy)./(abs(betatan(:,:,1))+abs(betatan(:,:,2)));
end
% if (mean(abs(Dux(:))))>(mean(abs(Duy(:))))
%     dz1=(-sign(alphatan(1,1)).*Dux)./(abs(alphatan(:,:,1))+abs(alphatan(:,:,2)));
% elseif (mean(abs(Duy(:))))>(mean(abs(Dux(:))))
%     %dz2=Duy./(abs(alphatan(:,:,1))+abs(alphatan(:,:,2)));%(abs(betatan(:,:,1))+abs(betatan(:,:,2)));
%     dz1=(sign(betatan(1,1)).*Duy)./(abs(betatan(:,:,1))+abs(betatan(:,:,2)));
% else
%     dz1=0;
%     fprintf('\n Disparity Map converged \n');
% end
%dz1
%if max(max(abs(Dux)))>max(max(abs(Duy)))
    
%max(max(abs(dz2)))

zgrid=zgrid-dz1;% the new z grid


% figure(21);surf(atand(alphatan(:,:,1)));
% figure(22);surf(atand(alphatan(:,:,2)));
% figure(23);surf(atand(betatan(:,:,1)));
% figure(24);surf(atand(betatan(:,:,2)));
% % % figure(13);surf(dz1);
% % % figure(14);surf(dz2);
% % % figure(15);surf(zgrid);
%  keyboard;




 



z1grid=zgrid;
%figure(3);surf(z1grid);title('change in z direction');xlabel('x(mm)');ylabel('y(mm)');
%figure(4);surf(z1grid);

end
