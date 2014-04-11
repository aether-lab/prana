function [scaling]=soloff_vec_reconstruction(diroutlist,caldata,pulsesep)

% This function performs the stereoreconstruction on the perviously
% processed vectors fields using the camera mapping information from the
% fitcameramodels.m function. The process is based off work by Soloff,
% Meas. Sci. Tech., 1997.
% Inputs:
t=1/pulsesep;
dir_struct1= dir(fullfile(diroutlist.soloff2dcam1,['*.' 'mat']));
flname1={dir_struct1.name}';
dir_struct2= dir(fullfile(diroutlist.soloff2dcam2,['*.' 'mat']));
flname2={dir_struct2.name}';
nof=length(flname1);

% Finals will not ALWAYS have "pass" in the name.  This could be a problem.
foutnamelist=regexp(flname1,'pass','split');
imagelist='';

% Predefining variables helps things run faster.
vectorlist = cell(nof,1);
%keyboard;
for j=1:nof
    
    vectorlist{j}=[{fullfile(diroutlist.soloff2dcam1,flname1{j})};{fullfile(diroutlist.soloff2dcam2,flname2{j})}];
    
    if j==1
        [~,dewarp_grid,scaling]=imagedewarp(caldata,'Soloff',imagelist,vectorlist{j});
        
        xgrid=dewarp_grid.xgrid;
        ygrid=dewarp_grid.ygrid;
        zgrid=zeros(size(xgrid));
        Xgrid1=dewarp_grid.Xgrid1;
        Ygrid1=dewarp_grid.Ygrid1;
        Xgrid2=dewarp_grid.Xgrid2;
        Ygrid2=dewarp_grid.Ygrid2;
    elseif j>1 && ~strcmp(foutnamelist{j}{2}(1),foutnamelist{j-1}{2}(1))
        [~,dewarp_grid,scaling]=imagedewarp(caldata,'Soloff',imagelist,vectorlist{j});
        
        xgrid=dewarp_grid.xgrid;
        ygrid=dewarp_grid.ygrid;
        zgrid=zeros(size(xgrid));
        Xgrid1=dewarp_grid.Xgrid1;
        Ygrid1=dewarp_grid.Ygrid1;
        Xgrid2=dewarp_grid.Xgrid2;
        Ygrid2=dewarp_grid.Ygrid2;
        
    end
    
    vecfr1 = load(vectorlist{j}{1});
    x1=vecfr1.X;
    y1=vecfr1.Y;
    u1=vecfr1.U(:,:,1);
    v1=vecfr1.V(:,:,1);
    clear vecfr1;
    vecfr2 = load(vectorlist{j}{2});
    x2=vecfr2.X;
    y2=vecfr2.Y;
    u2=vecfr2.U(:,:,1);
    v2=vecfr2.V(:,:,1);
    clear vecfr2;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Perform the interpolation in the image planes
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %     keyboard
    U1 = interp2(x1,y1,u1,Xgrid1,Ygrid1,'cubic',0);
    V1 = interp2(x1,y1,v1,Xgrid1,Ygrid1,'cubic',0);
    U2 = interp2(x2,y2,u2,Xgrid2,Ygrid2,'cubic',0);
    V2 = interp2(x2,y2,v2,Xgrid2,Ygrid2,'cubic',0);
    [rows,cols]=size(U1);
    
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
            dFdx1(:,:,gg) = a(2) + 2*a(5).*xgrid + a(6)*ygrid + a(8)*zgrid + 3*a(11)*xgrid.^2 + 2*a(12)*xgrid.*ygrid + ...
                a(13)*ygrid.^2 + 2*a(15)*xgrid.*zgrid + a(16)*ygrid.*zgrid + a(18)*zgrid.^2;
            
            dFdx2(:,:,gg) = a(3) + a(6)*xgrid + 2*a(7)*ygrid + a(9)*zgrid + a(12)*xgrid.^2 + 2*a(13)*xgrid.*ygrid + ...
                3*a(14)*ygrid.^2 + a(16)*xgrid.*zgrid + 2*a(17)*ygrid.*zgrid + a(19)*zgrid.^2;
            
            dFdx3(:,:,gg) = a(4) + a(8)*xgrid + a(9)*ygrid + 2*a(10)*zgrid + a(15)*xgrid.^2 + a(16)*xgrid.*ygrid + ...
                a(17)*ygrid.^2 + 2*a(18)*xgrid.*zgrid + 2*a(19)*ygrid.*zgrid;
        end
        
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Reconstruct the vectors according to Soloff, meas. sci. tech., 1997
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % solving using eqn 15 from Soloff's paper.
    u=zeros(rows,cols);
    v=zeros(rows,cols);
    w=zeros(rows,cols);
    for jj=1:rows
        for kk=1:cols
            d=[U1(jj,kk);V1(jj,kk);U2(jj,kk);V2(jj,kk)];
            C=[dFdx1(jj,kk,1) dFdx2(jj,kk,1) dFdx3(jj,kk,1);...
                dFdx1(jj,kk,2) dFdx2(jj,kk,2) dFdx3(jj,kk,2);...
                dFdx1(jj,kk,3) dFdx2(jj,kk,3) dFdx3(jj,kk,3);...
                dFdx1(jj,kk,4) dFdx2(jj,kk,4) dFdx3(jj,kk,4)];
            
            x=C\d;          % use lsqlin(C,d,...) for a constrained problem, this solves the linear system C*x=d
            
            u(jj,kk)=x(1);
            v(jj,kk)=x(2);
            w(jj,kk)=x(3);
        end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % output the plt file with all components
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    U=u.*(t/1000);% convert from mm/frame (mm displacement) to m/sec
    V=v.*(t/1000);%/data.pulsesep;%*1000;     % pulsesep is in microsec
    W=w.*(t/1000);%/data.pulsesep;%*1000;
    
    
    % Change from mm to m for SI output.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    X=xgrid./1000;
    Y=ygrid./1000;
    Z=zgrid./1000;
    
    %keyboard;
    
    %foutname=regexp(flname1{j},'pass','split');
    foutname=foutnamelist{j}{2};
    
    stereo_output=fullfile(diroutlist.soloff3cfields,['piv_2d3c_cam',num2str(caldata.camnumber(1)),'cam',num2str(caldata.camnumber(2)),'_pass_',foutname]);
    
    save(stereo_output,'X','Y','Z','U','V','W');
    
    fprintf(['stereo frame_pass_',foutname,'  done.\n']);
    %keyboard;
    clear X Y Z U V W;
end

end