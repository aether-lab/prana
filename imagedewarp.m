function [outputdirlist,dewarp_grid,scaling]=imagedewarp(caldata,dewarpmethod,imagelist,vectorlist)

%This code dewarps the images or vector grid depending on the
%reconstruction type and outputs the dewarped common grid coordinates and
%the modified magnification or scaling
%Inputs:-
%caldata=structure containing all the calibration information
%dewarpmethod='Soloff' or 'Willert
%imagelist=list containing directories of individual camera images (for Willert or selfcal)
%vectorlist=list containing directories of individual camera vector fields(Soloff)
%Outputs:-
%outputdirlist=dewarped image output directories
%dewarp_grid=dewarped grid coordinates
%scaling=magnification based on dewarped grid

%keyboard;
orderz=caldata.modeltype;
alldata.orderz=orderz;
optionsls=caldata.optionsls;

aXcam1=caldata.aXcam1;
aYcam1=caldata.aYcam1;
aXcam2=caldata.aXcam2;
aYcam2=caldata.aYcam2;

%keyboard;

if strcmp(dewarpmethod,'Willert')
    
    %keyboard;
    dir1=imagelist.imdirec;
    dir2=imagelist.imdirec2;
    base1=imagelist.imbase;
    base2=imagelist.imbase2;
    ext=imagelist.imext;
    zer=str2double(imagelist.imzeros);
    fstep=str2double(imagelist.imfstep);
    fstart=str2double(imagelist.imfstart);
    fend=str2double(imagelist.imfend);
    cstep=str2double(imagelist.imcstep);
    mkdir([dir1,'Dewarped Images',filesep]);
    mkdir([dir2,'Dewarped Images',filesep]);
    dirout1=([dir1,'Dewarped Images',filesep]);
    dirout2=([dir2,'Dewarped Images',filesep]);
    
    outputdirlist.dewarpdir1=dirout1;
    outputdirlist.dewarpdir2=dirout2;
    
    istring1=sprintf(['%%s%%s%%0%0.0fd.',ext],zer);
    
    IML=imread(sprintf(istring1,dir1,base1,1));
    IMR=imread(sprintf(istring1,dir2,base2,1));
    [Jmax1,Imax1] = size(IML);
    [Jmax2,Imax2] = size(IMR);
    
    X1points=[0.5 Jmax1-0.5 0.5 Jmax1-0.5];
    Y1points=[0.5 0.5 Imax1-0.5 Imax1-0.5];
    X2points=[0.5 Jmax2-0.5 0.5 Jmax2-0.5];
    Y2points=[0.5 0.5 Imax2-0.5 Imax2-0.5];
    
elseif strcmp(dewarpmethod,'Soloff')
    load(vectorlist{1});
    x1=X;
    y1=Y;
    clear X Y U V;
    load(vectorlist{2});
    x2=X;
    y2=Y;
    clear X Y U V;
    [Jmax1, Imax1]=size(x1);
    [Jmax2, Imax2]=size(x2);
    
    outputdirlist='';
    
    if caldata.targetside==1
            X1points=[min(min(x1)) max(max(x1)) min(min(x1)) max(max(x1))];
            Y1points=[min(min(y1)) min(min(y1)) max(max(y1)) max(max(y1))];
            X2points=[min(min(x2)) max(max(x2)) min(min(x2)) max(max(x2))];
            Y2points=[min(min(y2)) min(min(y2)) max(max(y2)) max(max(y2))];
    else    % different arrangement if camera 2 is on other side
            X1points=[min(min(x1)) max(max(x1)) min(min(x1)) max(max(x1))];
            Y1points=[min(min(y1)) min(min(y1)) max(max(y1)) max(max(y1))];
            X2points=[max(max(x2)) min(min(x2))  max(max(x2)) min(min(x2))];
            Y2points=[min(min(y2)) min(min(y2))  max(max(y2)) max(max(y2))];            
    end
    
end


    x0=[1 1];           % initial guess for solver
    % finds the over lap for the two cameras. This is only performed the
    % first time and then the information is used subsequently.
    %if c==1
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Determination of Area Overlap
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % find overlapping area using the corner points of the first loaded
        % images, the order for X1points and similar is [bl br tl tr] 
        
        fprintf('Computing overlapping area...\n');

        
%         X1points=[1 Jmax1 1 Jmax1];Y1points=[1 1 Imax1 Imax1];
%         X2points=[1 Jmax2 1 Jmax2];Y2points=[1 1 Imax2 Imax2];

%         X1points=[0 Jmax1 0 Jmax1];Y1points=[0 0 Imax1 Imax1];
%         X2points=[0 Jmax2 0 Jmax2];Y2points=[0 0 Imax2 Imax2];
        
        
      
          
        bottom1X=linspace(X1points(1),X1points(2),Jmax1)';
        top1X=linspace(X1points(3),X1points(4),Jmax1)';
        left1X=X1points(1)*ones(1,Imax1)';
        right1X=X1points(2)*ones(1,Imax1)';
        bottom1Y=Y1points(1)*ones(1,Jmax1)';
        top1Y=Y1points(3)*ones(1,Jmax1)';
        left1Y=linspace(Y1points(1),Y1points(3),Imax1)';
        right1Y=linspace(Y1points(2),Y1points(4),Imax1)';
        
        bottom2X=linspace(X2points(1),X2points(2),Jmax2)';
        top2X=linspace(X2points(3),X2points(4),Jmax2)';
        left2X=X2points(1)*ones(1,Imax2)';
        right2X=X2points(2)*ones(1,Imax2)';
        bottom2Y=Y2points(1)*ones(1,Jmax2)';
        top2Y=Y2points(3)*ones(1,Jmax2)';
        left2Y=linspace(Y2points(1),Y2points(3),Imax2)';
        right2Y=linspace(Y2points(2),Y2points(4),Imax2)';
        
        %keyboard;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % calculate x,y for the bottom,top,left,right vectors for X,Y from
        % camera 1
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        alldata.aX=aXcam1';
        alldata.aY=aYcam1';
        
        bottom1xy = zeros(2,Jmax1);
        top1xy    = zeros(2,Jmax1);
        points=[bottom1X bottom1Y];
        for k=1:Jmax1
            alldata.XYpoint=points(k,:)';
            % solve for x,y for camera 1
            [bottom1xy(:,k),~,~]=fsolve(@(x) poly_3xy_123z_2eqns(x,alldata),x0,optionsls);  
        end
        points=[top1X top1Y];
        for k=1:Jmax1
            alldata.XYpoint=points(k,:)';
            % solve for x,y for camera 1
            [top1xy(:,k),~,~]=fsolve(@(x) poly_3xy_123z_2eqns(x,alldata),x0,optionsls);
        end
        
        left1xy  = zeros(2,Imax1);
        right1xy = zeros(2,Imax1);
        points=[left1X left1Y];
        for k=1:Imax1
            alldata.XYpoint=points(k,:)';
             % solve for x,y for camera 1
            [left1xy(:,k),~,~]=fsolve(@(x) poly_3xy_123z_2eqns(x,alldata),x0,optionsls); 
        end
        points=[right1X right1Y];
        for k=1:Imax1
            alldata.XYpoint=points(k,:)';
            % solve for x,y for camera 1
            [right1xy(:,k),~,~]=fsolve(@(x) poly_3xy_123z_2eqns(x,alldata),x0,optionsls);  
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % calculate x,y for the bottom,top,left,right vectors for X,Y from
        % camera 2
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        alldata.aX=aXcam2';
        alldata.aY=aYcam2';
        
        bottom2xy = zeros(2,Jmax2);
        top2xy    = zeros(2,Jmax2);
        points=[bottom2X bottom2Y];
        for k=1:Jmax2
            alldata.XYpoint=points(k,:)';
            % solve for x,y for camera 2
            [bottom2xy(:,k),~,~]=fsolve(@(x) poly_3xy_123z_2eqns(x,alldata),x0,optionsls);  
        end
        points=[top2X top2Y];
        for k=1:Jmax2
            alldata.XYpoint=points(k,:)';
            % solve for x,y for camera 2
            [top2xy(:,k),~,~]=fsolve(@(x) poly_3xy_123z_2eqns(x,alldata),x0,optionsls);  
        end
        
        left2xy  = zeros(2,Imax2);
        right2xy = zeros(2,Imax2);
        points=[left2X left2Y];
        for k=1:Imax2
            alldata.XYpoint=points(k,:)';
            % solve for x,y for camera 2
            [left2xy(:,k),~,~]=fsolve(@(x) poly_3xy_123z_2eqns(x,alldata),x0,optionsls);  
        end
        points=[right2X right2Y];
        for k=1:Imax2
            alldata.XYpoint=points(k,:)';
            % solve for x,y for camera 2
            [right2xy(:,k),~,~]=fsolve(@(x) poly_3xy_123z_2eqns(x,alldata),x0,optionsls);  
        end
           
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % make object coordinate grid       
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        xlow=max([left1xy(1,:) left2xy(1,:)]);
        xhigh=min([right1xy(1,:) right2xy(1,:)]);
        ylow=max([bottom1xy(2,:) bottom2xy(2,:)]);
        yhigh=min([top1xy(2,:) top2xy(2,:)]);

        [xgrid,ygrid]=meshgrid(linspace(xlow,xhigh,Jmax1),linspace(ylow,yhigh,Imax1));
        %zgrid=zeros(size(xgrid));
        overplots=0;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Plot fig to check overlap
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if overplots == 1
            % use this to plot overlapping area in x,y
            figure(10);
            H(1) = plot(bottom1xy(1,:),bottom1xy(2,:),'+');hold on;
            plot(top1xy(1,:),top1xy(2,:),'+');hold on;
            plot(left1xy(1,:),left1xy(2,:),'+');hold on;
            plot(right1xy(1,:),right1xy(2,:),'+');hold on;
            H(2) = plot(bottom2xy(1,:),bottom2xy(2,:),'o');hold on;
            plot(top2xy(1,:),top2xy(2,:),'o');hold on;
            plot(left2xy(1,:),left2xy(2,:),'o');hold on;
            plot(right2xy(1,:),right2xy(2,:),'o');hold on;
            H(3) = plot([xlow xlow xhigh xhigh xlow],[ylow yhigh yhigh ylow ylow],'-k','LineWidth',2);hold on;
            H2 = plot(xgrid,ygrid,'.r','MarkerSize',4);xlabel('x (mm)');ylabel('y (mm)');
            H(4) = H2(1);
            title('Camera Overlap and new vector locations');
            Lstr = {'Camera 1 border','Camera 2 border','Overlap Border','Vector location'};
             legend(H,Lstr);
%             set(L,'Position',[0.4 0.4 0.2314 0.1869])
%             slashlocs = find(data.outputdirectory == '/');
%             set(gcf,'name',data.outputdirectory(slashlocs(end)+1:end))
            %         axis([-200 100 -150 250])
            drawnow
        end
        %keyboard;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Comput this grid in the IMAGE (object) plane to interpolate values
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
        [Xgrid1,Ygrid1]=poly_3xy_123z_fun(xgrid,ygrid,orderz,aXcam1,aYcam1);
        [Xgrid2,Ygrid2]=poly_3xy_123z_fun(xgrid,ygrid,orderz,aXcam2,aYcam2);
        
    dewarp_grid.Xgrid1=Xgrid1;
    dewarp_grid.Ygrid1=Ygrid1;
    dewarp_grid.Xgrid2=Xgrid2;
    dewarp_grid.Ygrid2=Ygrid2;
    dewarp_grid.xgrid=xgrid;
    dewarp_grid.ygrid=ygrid;
    if strcmp(dewarpmethod,'Willert')
        scaling.xscale =(max(xgrid(:))-min(xgrid(:)))/Imax1;
        scaling.yscale =(max(ygrid(:))-min(ygrid(:)))/Jmax1;
    elseif strcmp(dewarpmethod,'Soloff')
        scaling.xscale =(max(xgrid(:))-min(xgrid(:)))/(max(x1(:))-min(x1(:)));    
        scaling.yscale =(max(ygrid(:))-min(ygrid(:)))/(max(y1(:))-min(y1(:))); 
    end
    %keyboard;
    if strcmp(dewarpmethod,'Willert')
        [x1,y1] = meshgrid(0.5:1:Jmax1-0.5,0.5:1:Imax1-0.5);
        fprintf('Dewarping Images...\n');
                
        for k=fstart:fstep:fend+1
            %reading recorded images
            
            IMLi= im2double(imread(sprintf(istring1,dir1,base1,k)));
            IMRi= im2double(imread(sprintf(istring1,dir2,base2,k+cstep-1)));
            
            incl=class(IMLi);
            %flipping images
            IMLi=IMLi(end:-1:1,:);IMRi=IMRi(end:-1:1,:);
            %Interpolating on a common grid
            %     IMLo=double((sincBlackmanInterp2(IMLi, Xgrid1, Ygrid1, 3)));
            %     IMRo=double((sincBlackmanInterp2(IMRi, Xgrid2, Ygrid2, 3)));
            
            IMLo=double((interp2(x1,y1,IMLi, Xgrid1, Ygrid1, 'spline',0)));
            IMRo=double((interp2(x1,y1,IMRi, Xgrid2, Ygrid2, 'spline',0)));
            
            % figure(4);imagesc(IMLo);colormap('gray');%axis equal tight;
            % figure(5);imagesc(IMRo);colormap('gray');%axis equal tight;
            % figure(6);imagesc(IMLi);colormap('gray');%axis equal tight;
            % figure(7);imagesc(IMRi);colormap('gray');%axis equal tight;
            %  keyboard;
            %flipping them back
            IMLo=IMLo(end:-1:1,:);IMRo=IMRo(end:-1:1,:);
            IMLo=cast(IMLo,incl);IMRo=cast(IMRo,incl);
            %writing dewarped images to the output directory
            %writing the images because of two reasons:-
            %1)prana assumes all images to be in a sigle directory while processing
            %2)In prana processing just the image matrices cannot be given as input, have to give image location
            %     clear IM;
            %     IM(:,:,1)=IMLo;
            %     IM(:,:,2)=IMRo;
            %     IM(:,:,3)=IMRo;
            %     figure;imshow(IM);
            
            %keyboard;
            imwrite((IMLo),sprintf(istring1,dirout1,base1,k),'TIFF','WriteMode','overwrite','Compression','none');
            imwrite((IMRo),sprintf(istring1,dirout2,base2,k+cstep-1),'TIFF','WriteMode','overwrite','Compression','none');
            %keyboard;
            
        end
    end

end

function F=poly_3xy_123z_2eqns(x,alldata)
% F=poly_3xy_123z_2eqns(x,alldata)
% this function solves for the xy object coordinates with input
% image coordiantes alldata.XYpoint.  the resulting x vector contains
% the (x y) object coords.  This is for S-PIV so the z coord. is 0.

% This function is called by reconstructvectorsfun.m

% Writen by M. Brady
% Edited and Commented by S. Raben

aX=alldata.aX;
aY=alldata.aY;
orderz=alldata.orderz;
XYpoint=alldata.XYpoint;

if orderz==1                % cubic xy, linear z
    polylist=[1 x(1) x(2) 0 x(1)^2 x(1)*x(2) x(2)^2 0  0 x(1)^3 x(1)^2*x(2) x(1)*x(2)^2 x(2)^3 0 0 0]';
    Fpoly=[aX*polylist;aY*polylist]-XYpoint;

elseif orderz==2            % cubic xy, quadratic z
    polylist=[1 x(1) x(2) 0 x(1)^2 x(1)*x(2) x(2)^2 0  0 0 x(1)^3 x(1)^2*x(2) x(1)*x(2)^2 x(2)^3 0 0 0 0 0]';
    Fpoly=[aX*polylist;aY*polylist]-XYpoint; 
    
else             % camera pinhole
    
end

F=Fpoly;
end

function [Xgrid,Ygrid]=poly_3xy_123z_fun(xgrid,ygrid,orderz,aX,aY)
% [Xgrid Ygrid]=poly_3xy_123z_fun(xgrid,ygrid,orderz,aX,aY)
% 

% Writen by M. Brady
% Edited and Commented by S. Raben

x1=xgrid;
x2=ygrid;
[r,c]=size(xgrid);
x3=zeros(r,c);

if orderz==1                % cubic xy, linear z
    
    Xgrid=aX(1) + aX(2).*x1 + aX(3).*x2 + aX(4).*x3 + aX(5).*x1.^2 +...
        aX(6).*x1.*x2 + aX(7).*x2.^2 + aX(8).*x1.*x3 + aX(9).*x2.*x3 +...
        aX(10).*x1.^3 + aX(11).*x1.^2.*x2 + aX(12).*x1.*x2.^2 +...
        aX(13).*x2.^3 + aX(14).*x1.^2.*x3 + aX(15).*x1.*x2.*x3 +...
        aX(16).*x2.^2.*x3;
    
    Ygrid=aY(1) + aY(2).*x1 + aY(3).*x2 + aY(4).*x3 + aY(5).*x1.^2 +...
        aY(6).*x1.*x2 + aY(7).*x2.^2 + aY(8).*x1.*x3 + aY(9).*x2.*x3 +...
        aY(10).*x1.^3 + aY(11).*x1.^2.*x2 + aY(12).*x1.*x2.^2 +...
        aY(13).*x2.^3 + aY(14).*x1.^2.*x3 + aY(15).*x1.*x2.*x3 +...
        aY(16).*x2.^2.*x3;    

elseif orderz==2            % cubic xy, quadratic z
    
    Xgrid=aX(1) + aX(2).*x1 + aX(3).*x2 + aX(4).*x3 + aX(5).*x1.^2 +...
        aX(6).*x1.*x2 + aX(7).*x2.^2 + aX(8).*x1.*x3 + aX(9).*x2.*x3 +...
        aX(10).*x3.^2 + aX(11).*x1.^3 + aX(12).*x1.^2.*x2 + aX(13).*x1.*x2.^2 +...
        aX(14).*x2.^3 + aX(15).*x1.^2.*x3 + aX(16).*x1.*x2.*x3 +...
        aX(17).*x2.^2.*x3 + aX(18).*x1.*x3.^2 + aX(19).*x2.*x3.^2;   
    
    Ygrid=aY(1) + aY(2).*x1 + aY(3).*x2 + aY(4).*x3 + aY(5).*x1.^2 +...
        aY(6).*x1.*x2 + aY(7).*x2.^2 + aY(8).*x1.*x3 + aY(9).*x2.*x3 +...
        aY(10).*x3.^2 + aY(11).*x1.^3 + aY(12).*x1.^2.*x2 + aY(13).*x1.*x2.^2 +...
        aY(14).*x2.^3 + aY(15).*x1.^2.*x3 + aY(16).*x1.*x2.*x3 +...
        aY(17).*x2.^2.*x3 + aY(18).*x1.*x3.^2 + aY(19).*x2.*x3.^2;    
    
else             % pinhole
    
    
end
end

