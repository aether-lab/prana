function PIVadvance2code(Data)
%
% PIVadvance2code(Data)
%
% DPIV estimation and analysis tool based upon input structure "Data"
% generated using the PIVadvance2.fig GUI file

%% READ FORMATTED PARAMETERS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% clc;

%input/output directory
imbase=[Data.imdirec '/' Data.imbase];
pltdirec=[Data.outdirec '/'];

%image indices
I1 = str2double(Data.imfstart):str2double(Data.imfstep):str2double(Data.imfend);
I2 = I1+str2double(Data.imcstep);

%processing mask
if length(Data.maskname>0)
    mask = double(imread(Data.maskname));
else
    mask = 1+0*double(imread([imbase sprintf(['%0.' Data.imzeros 'i.' Data.imext],I1(1))]));
end
mask = flipud(mask);

%method and passes
P=str2double(Data.passes);
Method={'Multipass','Multigrid','Deform','Ensemble'};
M=Method(str2double(Data.method));

%interpolation and smoothing
Velsmoothswitch=str2double(Data.velsmooth);
Velsmoothfilt=str2double(Data.velsmoothfilt);
Velinterp=str2double(Data.velinterp);
Iminterp=str2double(Data.iminterp);

%physical parameters
Mag = str2double(Data.wrmag);
dt = str2double(Data.wrsep);
Freq = str2double(Data.wrsamp);

%initialization
Wres=zeros(P,2);
Wsize=zeros(P,2);
Gres=zeros(P,2);
Gbuf=zeros(P,2);
Corr=zeros(P);
D=zeros(P);
Valswitch=zeros(P);
Writeswitch=zeros(P);
Valsize=zeros(P,2,1);
Valthresh=zeros(P,1);
Uthresh=zeros(P,2);
Vthresh=zeros(P,2);
wbase=cell(0);

%read data info for each pass
for e=1:P
    
    %create structure for pass "e"
    eval(['A=Data.PIV' num2str(e) ';'])
    
    %store bulk window offset info
    if e==1
        BWO=[str2double(A.BWO(1:(strfind(A.BWO,',')-1))) str2double(A.BWO((strfind(A.BWO,',')+1):end))];
    end
    
    %window and grid resolution
    Wres(e,:)=[str2double(A.winres(1:(strfind(A.winres,',')-1))) str2double(A.winres((strfind(A.winres,',')+1):end))];
    Wsize(e,:)=[str2double(A.winsize(1:(strfind(A.winsize,',')-1))) str2double(A.winsize((strfind(A.winsize,',')+1):end))];
    Gres(e,:)=[str2double(A.gridres(1:(strfind(A.gridres,',')-1))) str2double(A.gridres((strfind(A.gridres,',')+1):end))];
    Gbuf(e,:)=[str2double(A.gridbuf(1:(strfind(A.gridbuf,',')-1))) str2double(A.gridbuf((strfind(A.gridbuf,',')+1):end))];
    
    %correlation
    Corr(e)=str2double(A.corr)-1;
    D(e)=str2double(A.RPCd);
    
    %validation and thresholding
    Valswitch(e)=str2double(A.val);
    Writeswitch(e)=str2double(A.write);
    
    vpass=[0 strfind(A.valsize,';') length(A.valsize)+1];
    for q=1:(length(vpass)-1)
        B=A.valsize((vpass(q)+1):(vpass(q+1)-1));
        Valsize(e,:,q)=[str2double(B(1:(strfind(B,',')-1))) str2double(B((strfind(B,',')+1):end))];
        Valthresh(e,q)=str2double(A.valthresh(1+2*(q-1)));
    end
    
    if str2double(A.thresh)==1
        Uthresh(e,:)=[str2double(A.valuthresh(1:(strfind(A.valuthresh,',')-1))) str2double(A.valuthresh((strfind(A.valuthresh,',')+1):end))];
        Vthresh(e,:)=[str2double(A.valvthresh(1:(strfind(A.valvthresh,',')-1))) str2double(A.valvthresh((strfind(A.valvthresh,',')+1):end))];
    else
        Uthresh(e,:)=[-inf,inf];
        Vthresh(e,:)=[-inf,inf];
    end
    
    %output directory
    wbase(e,:)={A.outbase};
    
end

%% IMAGE PREFILTER

%added to multipass, multigrid, deform, and ensemble for image loading
try
    IMmin=double(imread([Data.imdirec '/IMmin.tif']));
catch
    beep
    disp('Error reading image subtraction file...\nResuming without image prefilter...')
    IMmin=0*double(imread([imbase sprintf(['%0.' Data.imzeros 'i.' Data.imext],I1(1))]));
end

%% EVALUATE IMAGE SEQUENCE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

switch char(M)

    case 'Multipass'

        for q=1:length(I1)
            
            tf=cputime;
            
            %output text
            title=['Frame' sprintf(['%0.' Data.imzeros 'i'],I1(q)) ' and Frame' sprintf(['%0.' Data.imzeros 'i'],I2(q))];
            fprintf('\n----------------------------------------------------\n')
            fprintf(['Processing ' title ' (' num2str(q) '/' num2str(length(I1)) ')\n'])
            fprintf('----------------------------------------------------\n')

            %load image pair and flip coordinates
            im1=double(imread([imbase sprintf(['%0.' Data.imzeros 'i.' Data.imext],I1(q))]))-IMmin;
            im2=double(imread([imbase sprintf(['%0.' Data.imzeros 'i.' Data.imext],I2(q))]))-IMmin;
            im1=flipud(im1);
            im2=flipud(im2);
            L=size(im1);

            %initialize grid and evaluation matrix
            [X,Y]=IMgrid(L,Gres(e,:),Gbuf(e,:));
            X=X(:);Y=Y(:);
            Eval = reshape(downsample(downsample( mask(Y(1):Y(end),X(1):X(end)),Gres(e,2))',Gres(e,1))',length(X),1);
            Eval(Eval==0)=-1;
            Eval(Eval>0)=0;
            Ub = BWO(1)*ones(size(X));
            Vb = BWO(2)*ones(size(X));
            U  = zeros(size(X));
            V  = zeros(size(X));

            for e=1:P

                %correlate image pair
                [Xc,Yc,Uc,Vc]=PIVwindowed(im1,im2,Corr(e),Wsize(e,:),Wres(e,:),0,D(e),X(Eval>=0),Y(Eval>=0),Ub(Eval>=0),Vb(Eval>=0));
                U(Eval>=0)=Uc;
                V(Eval>=0)=Vc;

                %validation
                if Valswitch(e)==1
                    t=permute(Valsize(e,:,:),[2 3 1]);
                    t=t(:,t(1,:)~=0);
                    [U,V,Eval] = VAL(X,Y,U,V,t',Valthresh(e,Valthresh(e,:)~=0)',Uthresh(e,:),Vthresh(e,:),Eval);
                end

                %reset DWO for subsequent passes
                Ub=U;
                Vb=V;
                
                %write outut
                if Writeswitch(e)==1
                    time=(q-1)/Freq;
                    write_plt_val([pltdirec char(wbase(e,:)) sprintf(['%0.' Data.imzeros 'i.plt' ],I1(q))],X,Y,U,V,Eval,Mag,dt,time,title);
                end

            end
            
            eltime=cputime-tf;
            fprintf('total frame time...              %0.2i:%0.2i.%0.0f\n',floor(eltime/60),floor(rem(eltime,60)),rem(eltime,60)-floor(rem(eltime,60)))
            
        end


    case 'Multigrid'

        for q=1:length(I1)
            
            tf=cputime;
            
            %output text
            title=['Frame' sprintf(['%0.' Data.imzeros 'i'],I1(q)) ' and Frame' sprintf(['%0.' Data.imzeros 'i'],I2(q))];
            fprintf('\n----------------------------------------------------\n')
            fprintf(['Processing ' title ' (' num2str(q) '/' num2str(length(I1)) ')\n'])
            fprintf('----------------------------------------------------\n')
            
            %load image pair and flip coordinates
            im1=double(imread([imbase sprintf(['%0.' Data.imzeros 'i.' Data.imext],I1(q))]))-IMmin;
            im2=double(imread([imbase sprintf(['%0.' Data.imzeros 'i.' Data.imext],I2(q))]))-IMmin;
            im1=flipud(im1);
            im2=flipud(im2);
            L=size(im1);

            %initialize grid and evaluation matrix
            [XI,YI]=IMgrid(L,[0 0]);
            UI = BWO(1)*ones(size(XI));
            VI = BWO(2)*ones(size(XI));

            for e=1:P

                %find grid and evaluation matrix for each pass
                [X,Y]=IMgrid(L,Gres(e,:),Gbuf(e,:));
                S=size(X);X=X(:);Y=Y(:);
                Ub = reshape(downsample(downsample( UI(Y(1):Y(end),X(1):X(end)),Gres(e,2))',Gres(e,1))',length(X),1);
                Vb = reshape(downsample(downsample( VI(Y(1):Y(end),X(1):X(end)),Gres(e,2))',Gres(e,1))',length(X),1);
                Eval = reshape(downsample(downsample( mask(Y(1):Y(end),X(1):X(end)),Gres(e,2))',Gres(e,1))',length(X),1);
                Eval(Eval==0)=-1;
                Eval(Eval>0)=0;

                %correlate image pair
                U=zeros(size(X));
                V=zeros(size(X));
                [Xc,Yc,Uc,Vc]=PIVwindowed(im1,im2,Corr(e),Wsize(e,:),Wres(e,:),0,D(e),X(Eval>=0),Y(Eval>=0),Ub(Eval>=0),Vb(Eval>=0));
                U(Eval>=0)=Uc;
                V(Eval>=0)=Vc;

                %validation
                if Valswitch(e)==1
                    t=permute(Valsize(e,:,:),[2 3 1]);
                    t=t(:,t(1,:)~=0);
                    [U,V,Eval] = VAL(X,Y,U,V,t',Valthresh(e,Valthresh(e,:)~=0)',Uthresh(e,:),Vthresh(e,:),Eval);
                end
                
                %write output
                if Writeswitch(e)==1
                    time=(q-1)/Freq;
                    write_plt_val([pltdirec char(wbase(e,:)) sprintf(['%0.' Data.imzeros 'i.plt' ],I1(q))],X,Y,U,V,Eval,Mag,dt,time,title);
                end

                if e~=P

                    fprintf('interpolating velocity...        ')
                    t1=cputime;

                    %reshape from list of grid points to matrix
                    X=reshape(X,S(1),S(2));
                    Y=reshape(Y,S(1),S(2));
                    U=reshape(U,S(1),S(2));
                    V=reshape(V,S(1),S(2));

                    %velocity smoothing
                    if Velsmoothswitch==1
                        [U,V]=VELfilt(U,V,Velsmoothfilt);
                    end

                    %velocity interpolation
                    UI = VFinterp(X,Y,U,XI,YI,Velinterp);
                    VI = VFinterp(X,Y,V,XI,YI,Velinterp);

                    eltime=cputime-t1;
                    fprintf('%0.2i:%0.2i.%0.0f\n',floor(eltime/60),floor(rem(eltime,60)),rem(eltime,60)-floor(rem(eltime,60)))

                end

            end
            
            eltime=cputime-tf;
            fprintf('total frame time...              %0.2i:%0.2i.%0.0f\n',floor(eltime/60),floor(rem(eltime,60)),rem(eltime,60)-floor(rem(eltime,60)))
            
        end


    case 'Deform'

        for q=1:length(I1)
            
            tf=cputime;
            
            %output text
            title=['Frame' sprintf(['%0.' Data.imzeros 'i'],I1(q)) ' and Frame' sprintf(['%0.' Data.imzeros 'i'],I2(q))];
            fprintf('\n----------------------------------------------------\n')
            fprintf(['Processing ' title ' (' num2str(q) '/' num2str(length(I1)) ')\n'])
            fprintf('----------------------------------------------------\n')
            
            %load image pair and flip coordinates
            im1=double(imread([imbase sprintf(['%0.' Data.imzeros 'i.' Data.imext],I1(q))]))-IMmin;
            im2=double(imread([imbase sprintf(['%0.' Data.imzeros 'i.' Data.imext],I2(q))]))-IMmin;
            im1=flipud(im1);
            im2=flipud(im2);
            L=size(im1);

            %initialize deformed images
            im1d=im1;
            im2d=im2;

            %initialize grid and evaluation matrix
            [XI,YI]=IMgrid(L,[0 0]);
            UI = BWO(1)*ones(size(XI));
            VI = BWO(2)*ones(size(XI));

            for e=1:P

                %find grid and evaluation matrix for each pass
                [X,Y]=IMgrid(L,Gres(e,:),Gbuf(e,:));
                S=size(X);X=X(:);Y=Y(:);
                Ub = reshape(downsample(downsample( UI(Y(1):Y(end),X(1):X(end)),Gres(e,2))',Gres(e,1))',length(X),1);
                Vb = reshape(downsample(downsample( VI(Y(1):Y(end),X(1):X(end)),Gres(e,2))',Gres(e,1))',length(X),1);
                Eval = reshape(downsample(downsample( mask(Y(1):Y(end),X(1):X(end)),Gres(e,2))',Gres(e,1))',length(X),1);
                Eval(Eval==0)=-1;
                Eval(Eval>0)=0;

                %correlate image pair
                U=zeros(size(X));
                V=zeros(size(X));
                [Xc,Yc,Uc,Vc]=PIVwindowed(im1d,im2d,Corr(e),Wsize(e,:),Wres(e,:),0,D(e),X(Eval>=0),Y(Eval>=0));
                U(Eval>=0)=Uc;
                V(Eval>=0)=Vc;
                U=U+Ub;
                V=V+Vb;
                
                %validation
                if Valswitch(e)==1
                    t=permute(Valsize(e,:,:),[2 3 1]);
                    t=t(:,t(1,:)~=0);
                    [U,V,Eval] = VAL(X,Y,U,V,t',Valthresh(e,Valthresh(e,:)~=0)',Uthresh(e,:),Vthresh(e,:),Eval);
                end
                
                %write output
                if Writeswitch(e)==1
                    time=(q-1)/Freq;
                    write_plt_val([pltdirec char(wbase(e,:)) sprintf(['%0.' Data.imzeros 'i.plt' ],I1(q))],X,Y,U,V,Eval,Mag,dt,time,title);
                end

                if e~=P

                    fprintf('interpolating velocity...        ')
                    t1=cputime;

                    %reshape from list of grid points to matrix
                    X=reshape(X,S(1),S(2));
                    Y=reshape(Y,S(1),S(2));
                    U=reshape(U,S(1),S(2));
                    V=reshape(V,S(1),S(2));

                    %velocity smoothing
                    if Velsmoothswitch==1
                        [U,V]=VELfilt(U,V,Velsmoothfilt);
                    end

                    %velocity interpolation
                    UI = VFinterp(X,Y,U,XI,YI,Velinterp);
                    VI = VFinterp(X,Y,V,XI,YI,Velinterp);

                    %translate pixel locations
                    XD1 = XI+UI/2;
                    YD1 = YI+VI/2;
                    XD2 = XI-UI/2;
                    YD2 = YI-VI/2;

                    %output text
                    eltime=cputime-t1;
                    fprintf('%0.2i:%0.2i.%0.0f\n',floor(eltime/60),floor(rem(eltime,60)),rem(eltime,60)-floor(rem(eltime,60)))
                    fprintf('deforming images...              ')
                    t1=cputime;

                    %preallocate deformed images
                    im1d = zeros(L);
                    im2d = zeros(L);

                    %cardinal function interpolation
                    if Iminterp==1

                        for i=1:L(1)
                            for j=1:L(2)

                                %image 1 interpolation
                                nmin=max([1 (round(YD1(i,j))-3)]);
                                nmax=min([L(1) (round(YD1(i,j))+3)]);
                                mmin=max([1 (round(XD1(i,j))-3)]);
                                mmax=min([L(2) (round(XD1(i,j))+3)]);
                                for n=nmin:nmax
                                    for m=mmin:mmax
                                        wi = sin(pi*(m-XD1(i,j)))*sin(pi*(n-YD1(i,j)))/(pi^2*(m-XD1(i,j))*(n-YD1(i,j)));
                                        im1d(n,m)=im1d(n,m)+im1(i,j)*wi;
                                    end
                                end

                                %image 2 interpolation
                                nmin=max([1 (round(YD2(i,j))-3)]);
                                nmax=min([L(1) (round(YD2(i,j))+3)]);
                                mmin=max([1 (round(XD2(i,j))-3)]);
                                mmax=min([L(2) (round(XD2(i,j))+3)]);
                                for n=nmin:nmax
                                    for m=mmin:mmax
                                        wi = sin(pi*(m-XD2(i,j)))*sin(pi*(n-YD2(i,j)))/(pi^2*(m-XD2(i,j))*(n-YD2(i,j)));
                                        im2d(n,m)=im2d(n,m)+im2(i,j)*wi;
                                    end
                                end

                            end
                        end

                    %cardinal function interpolation with Blackman filter
                    elseif Iminterp==2

                        for i=1:L(1)
                            for j=1:L(2)

                                %image 1 interpolation
                                nmin=max([1 (round(YD1(i,j))-3)]);
                                nmax=min([L(1) (round(YD1(i,j))+3)]);
                                mmin=max([1 (round(XD1(i,j))-3)]);
                                mmax=min([L(2) (round(XD1(i,j))+3)]);
                                for n=nmin:nmax
                                    for m=mmin:mmax
                                        wi = sin(pi*(m-XD1(i,j)))*sin(pi*(n-YD1(i,j)))/(pi^2*(m-XD1(i,j))*(n-YD1(i,j)));
                                        bi = (0.42+0.5*cos(pi*(m-XD1(i,j))/3)+0.08*cos(2*pi*(m-XD1(i,j))/3))*(0.42+0.5*cos(pi*(n-YD1(i,j))/3)+0.08*cos(2*pi*(n-YD1(i,j))/3));
                                        im1d(n,m)=im1d(n,m)+im1(i,j)*wi*bi;
                                    end
                                end

                                %image 2 interpolation
                                nmin=max([1 (round(YD2(i,j))-3)]);
                                nmax=min([L(1) (round(YD2(i,j))+3)]);
                                mmin=max([1 (round(XD2(i,j))-3)]);
                                mmax=min([L(2) (round(XD2(i,j))+3)]);
                                for n=nmin:nmax
                                    for m=mmin:mmax
                                        wi = sin(pi*(m-XD2(i,j)))*sin(pi*(n-YD2(i,j)))/(pi^2*(m-XD2(i,j))*(n-YD2(i,j)));
                                        bi = (0.42+0.5*cos(pi*(m-XD2(i,j))/3)+0.08*cos(2*pi*(m-XD2(i,j))/3))*(0.42+0.5*cos(pi*(n-YD2(i,j))/3)+0.08*cos(2*pi*(n-YD2(i,j))/3));
                                        im2d(n,m)=im2d(n,m)+im2(i,j)*wi*bi;
                                    end
                                end

                            end
                        end

                    end

                    %clip lower values of deformed images
                    im1d(im1d<0)=0;
                    im2d(im2d<0)=0;

                    eltime=cputime-t1;
                    fprintf('%0.2i:%0.2i.%0.0f\n',floor(eltime/60),floor(rem(eltime,60)),rem(eltime,60)-floor(rem(eltime,60)))

                end

            end
            
            eltime=cputime-tf;
            fprintf('total frame time...              %0.2i:%0.2i.%0.0f\n',floor(eltime/60),floor(rem(eltime,60)),rem(eltime,60)-floor(rem(eltime,60)))
            
        end


    case 'Ensemble'
        
        %initialize grid and evaluation matrix
        im1=double(imread([imbase sprintf(['%0.' Data.imzeros 'i.' Data.imext],I1(1))]));
        L=size(im1);
        [XI,YI]=IMgrid(L,[0 0]);
        UI = BWO(1)*ones(size(XI));
        VI = BWO(2)*ones(size(XI));

        for e=1:P

            %find grid and evaluation matrix for each pass
            [X,Y]=IMgrid(L,Gres(e,:),Gbuf(e,:));
            S=size(X);X=X(:);Y=Y(:);
            Ub = reshape(downsample(downsample( UI(Y(1):Y(end),X(1):X(end)),Gres(e,2))',Gres(e,1))',length(X),1);
            Vb = reshape(downsample(downsample( VI(Y(1):Y(end),X(1):X(end)),Gres(e,2))',Gres(e,1))',length(X),1);
            Eval = reshape(downsample(downsample( mask(Y(1):Y(end),X(1):X(end)),Gres(e,2))',Gres(e,1))',length(X),1);
            Eval(Eval==0)=-1;
            Eval(Eval>0)=0;

            %initialize velocity outputs
            U=zeros(size(X));
            V=zeros(size(X));
            
            %output text
            title=['Frame' sprintf(['%0.' Data.imzeros 'i'],I1(1)) ' to Frame' sprintf(['%0.' Data.imzeros 'i'],I2(end))];
            fprintf('\n----------------------------------------------------\n')
            fprintf(['Ensemble Correlation ' title '\n'])
            fprintf('----------------------------------------------------\n')

            for q=1:length(I1)

                %load image pair and flip coordinates
                im1=double(imread([imbase sprintf(['%0.' Data.imzeros 'i.' Data.imext],I1(q))]))-IMmin;
                im2=double(imread([imbase sprintf(['%0.' Data.imzeros 'i.' Data.imext],I2(q))]))-IMmin;
                im1=flipud(im1);
                im2=flipud(im2);
                L=size(im1);
                
                %correlate image pair and average correlations
                [Xc,Yc,CC]=PIVensemble(im1,im2,Corr(e),Wsize(e,:),Wres(e,:),0,D(e),X(Eval>=0),Y(Eval>=0),Ub(Eval>=0),Vb(Eval>=0));
                if q==1
                    CCm=CC/length(I1);
                else
                    CCm=CCm+CC/length(I1);
                end
                fprintf(['(' sprintf(['%0.' num2str(length(num2str(length(I1)))) 'i' ],q) '/' num2str(length(I1)) ')\n'])

            end
            fprintf('----------------------------------------------------\n\n')

            %evaluate subpixel displacement of averaged correlation
            Z=size(CCm);
            Uc=zeros(Z(3),1);
            Vc=zeros(Z(3),1);
            ZZ=ones(Z(1),Z(2));
            for s=1:length(Xc)
                [Uc(s),Vc(s)]=subpixel(CCm(:,:,s),Z(2),Z(1),ZZ);
            end
            U(Eval>=0)=Uc+round(Ub(Eval>=0));
            V(Eval>=0)=Vc+round(Vb(Eval>=0));

            %validation
            if Valswitch(e)==1
                t=permute(Valsize(e,:,:),[2 3 1]);
                t=t(:,t(1,:)~=0);
                [U,V,Eval] = VAL(X,Y,U,V,t',Valthresh(e,Valthresh(e,:)~=0)',Uthresh(e,:),Vthresh(e,:),Eval);
            end

            %write output
            if Writeswitch(e)==1
                write_plt_val([pltdirec char(wbase(e,:)) sprintf(['%0.' Data.imzeros 'i_' ],I1(1)) sprintf(['%0.' Data.imzeros 'i.plt' ],I1(end))],X,Y,U,V,Eval,Mag,dt,0,title);
            end
            
            if e~=P

                fprintf('interpolating velocity...        ')
                t1=cputime;

                %reshape from list of grid points to matrix
                X=reshape(X,S(1),S(2));
                Y=reshape(Y,S(1),S(2));
                U=reshape(U,S(1),S(2));
                V=reshape(V,S(1),S(2));

                %velocity smoothing
                if Velsmoothswitch==1
                    [U,V]=VELfilt(U,V,Velsmoothfilt);
                end

                %velocity interpolation
                UI = VFinterp(X,Y,U,XI,YI,Velinterp);
                VI = VFinterp(X,Y,V,XI,YI,Velinterp);

                %output text
                eltime=cputime-t1;
                fprintf('%0.2i:%0.2i.%0.0f\n',floor(eltime/60),floor(rem(eltime,60)),rem(eltime,60)-floor(rem(eltime,60)))

            end

        end

end

%signal job complete
beep,pause(0.2),beep

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                           END MAIN FUNCTION                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




%% DPIV CORRELATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [X,Y,U,V]=PIVwindowed(im1,im2,corr,window,res,zpad,D,X,Y,Uin,Vin)

t1=cputime;

%convert input parameters
im1=double(im1);
im2=double(im2);
L = size(im1);

%convert to gridpoint list
X=X(:);
Y=Y(:);

%correlation and window mask types
ctype    = {'SCC','RPC'};
tcorr = char(ctype(corr+1)); 

%preallocate velocity fields and grid format
Nx = window(1);
Ny = window(2);
if nargin <=10
    Uin = zeros(length(X),1);
    Vin = zeros(length(X),1);
end
Uin = Uin(:);
Vin = Vin(:);
U = zeros(length(X),1);
V = zeros(length(X),1);

%sets up extended domain size
if zpad~=0
    Sy=2*Ny;
    Sx=2*Nx;
else
    Sy=Ny;
    Sx=Nx;
end

%window masking filter
sfilt = windowmask([Nx Ny],[res(1) res(2)]);

%correlation plane normalization function (always off)
cnorm = ones(Ny,Nx);

%RPC spectral energy filter
spectral = fftshift(energyfilt(Sx,Sy,D,0));

%fftshift indicies
fftindy = [Sy/2+1:Sy 1:Sy/2];
fftindx = [Sx/2+1:Sx 1:Sx/2];

fprintf('correlating...                   ')

switch upper(tcorr)

    %Standard Cross Correlation
    case 'SCC'

        for n=1:length(X)

            %apply the second order discrete window offset
            x1 = X(n) - floor(round(Uin(n))/2);
            x2 = X(n) +  ceil(round(Uin(n))/2);

            y1 = Y(n) - floor(round(Vin(n))/2);
            y2 = Y(n) +  ceil(round(Vin(n))/2);

            xmin1 = x1-Nx/2+1;
            xmax1 = x1+Nx/2;
            xmin2 = x2-Nx/2+1;
            xmax2 = x2+Nx/2;
            ymin1 = y1-Ny/2+1;
            ymax1 = y1+Ny/2;
            ymin2 = y2-Ny/2+1;
            ymax2 = y2+Ny/2;

            %find the image windows
            zone1 = im1( max([1 ymin1]):min([L(1) ymax1]),max([1 xmin1]):min([L(2) xmax1]) );
            zone2 = im2( max([1 ymin2]):min([L(1) ymax2]),max([1 xmin2]):min([L(2) xmax2]) );
            if size(zone1,1)~=Ny || size(zone1,2)~=Nx
                w1 = zeros(Ny,Nx);
                w1( 1+max([0 1-ymin1]):Ny-max([0 ymax1-L(1)]),1+max([0 1-xmin1]):Nx-max([0 xmax1-L(2)]) ) = zone1;
                zone1 = w1;
            end
            if size(zone2,1)~=Ny || size(zone2,2)~=Nx
                w2 = zeros(Ny,Nx);
                w2( 1+max([0 1-ymin2]):Ny-max([0 ymax2-L(1)]),1+max([0 1-xmin2]):Nx-max([0 xmax2-L(2)]) ) = zone2;
                zone2 = w2;
            end
            
            %apply the image spatial filter
            region1 = (zone1).*sfilt;
            region2 = (zone2).*sfilt;

            %FFTs and Cross-Correlation
            f1   = fftn(region1,[Sy Sx]);
            f2   = fftn(region2,[Sy Sx]);
            P21  = f2.*conj(f1);

            %Standard Fourier Based Cross-Correlation
            G = ifftn(P21,'symmetric');
            G = G(fftindy,fftindx);
            G = abs(G);
            
            %subpixel estimation
            [U(n),V(n)]=subpixel(G,Nx,Ny,cnorm);

        end

    %Robust Phase Correlation
    case 'RPC'
        
        for n=1:length(X)

            %apply the second order discrete window offset
            x1 = X(n) - floor(round(Uin(n))/2);
            x2 = X(n) +  ceil(round(Uin(n))/2);

            y1 = Y(n) - floor(round(Vin(n))/2);
            y2 = Y(n) +  ceil(round(Vin(n))/2);

            xmin1 = x1-Nx/2+1;
            xmax1 = x1+Nx/2;
            xmin2 = x2-Nx/2+1;
            xmax2 = x2+Nx/2;
            ymin1 = y1-Ny/2+1;
            ymax1 = y1+Ny/2;
            ymin2 = y2-Ny/2+1;
            ymax2 = y2+Ny/2;

            %find the image windows
            zone1 = im1( max([1 ymin1]):min([L(1) ymax1]),max([1 xmin1]):min([L(2) xmax1]) );
            zone2 = im2( max([1 ymin2]):min([L(1) ymax2]),max([1 xmin2]):min([L(2) xmax2]) );
            if size(zone1,1)~=Ny || size(zone1,2)~=Nx
                w1 = zeros(Ny,Nx);
                w1( 1+max([0 1-ymin1]):Ny-max([0 ymax1-L(1)]),1+max([0 1-xmin1]):Nx-max([0 xmax1-L(2)]) ) = zone1;
                zone1 = w1;
            end
            if size(zone2,1)~=Ny || size(zone2,2)~=Nx
                w2 = zeros(Ny,Nx);
                w2( 1+max([0 1-ymin2]):Ny-max([0 ymax2-L(1)]),1+max([0 1-xmin2]):Nx-max([0 xmax2-L(2)]) ) = zone2;
                zone2 = w2;
            end

            %apply the image spatial filter
            region1 = (zone1).*sfilt;
            region2 = (zone2).*sfilt;

            %FFTs and Cross-Correlation
            f1   = fftn(region1,[Sy Sx]);
            f2   = fftn(region2,[Sy Sx]);
            P21  = f2.*conj(f1);

            %Phase Correlation
            W = ones(Sy,Sx);
            Wden = sqrt(P21.*conj(P21));
            W(P21~=0) = Wden(P21~=0);
            R = P21./W;

            %Robust Phase Correlation with spectral energy filter
            G = ifftn(R.*spectral,'symmetric');
            G = G(fftindy,fftindx);
            G = abs(G);
            
            %subpixel estimation
            [U(n),V(n)]=subpixel(G,Nx,Ny,cnorm);

        end
end

%add DWO to estimation
U = round(Uin)+U;
V = round(Vin)+V;

eltime=cputime-t1;
fprintf('%0.2i:%0.2i.%0.0f\n',floor(eltime/60),floor(rem(eltime,60)),rem(eltime,60)-floor(rem(eltime,60)))


%% DPIV ENSEMBLE CORRELATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [X,Y,CC]=PIVensemble(im1,im2,corr,window,res,zpad,D,X,Y,Uin,Vin)

t1=cputime;

%convert input parameters
im1=double(im1);
im2=double(im2);
L = size(im1);

%convert to gridpoint list
X=X(:);
Y=Y(:);

%correlation and window mask types
ctype    = {'SCC','RPC'};
tcorr = char(ctype(corr+1)); 

%preallocate velocity fields and grid format
Nx = window(1);
Ny = window(2);
if nargin ==10
    Uin = zeros(length(X),1);
    Vin = zeros(length(X),1);
end
Uin = Uin(:);
Vin = Vin(:);

%sets up extended domain size
if zpad~=0
    Sy=2*Ny;
    Sx=2*Nx;
else
    Sy=Ny;
    Sx=Nx;
end

%window masking filter
sfilt = windowmask([Nx Ny],[res(1) res(2)]);

%RPC spectral energy filter
spectral = fftshift(energyfilt(Sx,Sy,D,0));

%fftshift indicies
fftindy = [Sy/2+1:Sy 1:Sy/2];
fftindx = [Sx/2+1:Sx 1:Sx/2];

fprintf('correlating...                   ')

%initialize correlation tensor
CC = zeros(Sy,Sx,length(X));

switch upper(tcorr)

    %Standard Cross Correlation
    case 'SCC'

        for n=1:length(X)

            %apply the second order discrete window offset
            x1 = X(n) - floor(round(Uin(n))/2);
            x2 = X(n) +  ceil(round(Uin(n))/2);

            y1 = Y(n) - floor(round(Vin(n))/2);
            y2 = Y(n) +  ceil(round(Vin(n))/2);

            xmin1 = x1-Nx/2+1;
            xmax1 = x1+Nx/2;
            xmin2 = x2-Nx/2+1;
            xmax2 = x2+Nx/2;
            ymin1 = y1-Ny/2+1;
            ymax1 = y1+Ny/2;
            ymin2 = y2-Ny/2+1;
            ymax2 = y2+Ny/2;

            %find the image windows
            zone1 = im1( max([1 ymin1]):min([L(1) ymax1]),max([1 xmin1]):min([L(2) xmax1]) );
            zone2 = im2( max([1 ymin2]):min([L(1) ymax2]),max([1 xmin2]):min([L(2) xmax2]) );
            if size(zone1,1)~=Ny || size(zone1,2)~=Nx
                w1 = zeros(Ny,Nx);
                w1( 1+max([0 1-ymin1]):Ny-max([0 ymax1-L(1)]),1+max([0 1-xmin1]):Nx-max([0 xmax1-L(2)]) ) = zone1;
                zone1 = w1;
            end
            if size(zone2,1)~=Ny || size(zone2,2)~=Nx
                w2 = zeros(Ny,Nx);
                w2( 1+max([0 1-ymin2]):Ny-max([0 ymax2-L(1)]),1+max([0 1-xmin2]):Nx-max([0 xmax2-L(2)]) ) = zone2;
                zone2 = w2;
            end
            
            %apply the image spatial filter
            region1 = (zone1).*sfilt;
            region2 = (zone2).*sfilt;

            %FFTs and Cross-Correlation
            f1   = fftn(region1-mean(region1(:)),[Sy Sx]);
            f2   = fftn(region2-mean(region2(:)),[Sy Sx]);
            P21  = f2.*conj(f1);

            %Standard Fourier Based Cross-Correlation
            G = ifftn(P21,'symmetric');
            G = G(fftindy,fftindx);
            G = abs(G);
            G = G/std(region1(:))/std(region2(:))/length(region1(:));
            
            %store correlation matrix
            CC(:,:,n) = G;

        end

    %Robust Phase Correlation
    case 'RPC'
        
        for n=1:length(X)

            %apply the second order discrete window offset
            x1 = X(n) - floor(round(Uin(n))/2);
            x2 = X(n) +  ceil(round(Uin(n))/2);

            y1 = Y(n) - floor(round(Vin(n))/2);
            y2 = Y(n) +  ceil(round(Vin(n))/2);

            xmin1 = x1-Nx/2+1;
            xmax1 = x1+Nx/2;
            xmin2 = x2-Nx/2+1;
            xmax2 = x2+Nx/2;
            ymin1 = y1-Ny/2+1;
            ymax1 = y1+Ny/2;
            ymin2 = y2-Ny/2+1;
            ymax2 = y2+Ny/2;

            %find the image windows
            zone1 = im1( max([1 ymin1]):min([L(1) ymax1]),max([1 xmin1]):min([L(2) xmax1]) );
            zone2 = im2( max([1 ymin2]):min([L(1) ymax2]),max([1 xmin2]):min([L(2) xmax2]) );
            if size(zone1,1)~=Ny || size(zone1,2)~=Nx
                w1 = zeros(Ny,Nx);
                w1( 1+max([0 1-ymin1]):Ny-max([0 ymax1-L(1)]),1+max([0 1-xmin1]):Nx-max([0 xmax1-L(2)]) ) = zone1;
                zone1 = w1;
            end
            if size(zone2,1)~=Ny || size(zone2,2)~=Nx
                w2 = zeros(Ny,Nx);
                w2( 1+max([0 1-ymin2]):Ny-max([0 ymax2-L(1)]),1+max([0 1-xmin2]):Nx-max([0 xmax2-L(2)]) ) = zone2;
                zone2 = w2;
            end

            %apply the image spatial filter
            region1 = zone1.*sfilt;
            region2 = zone2.*sfilt;

            %FFTs and Cross-Correlation
            f1   = fftn(region1,[Sy Sx]);
            f2   = fftn(region2,[Sy Sx]);
            P21  = f2.*conj(f1);

            %Phase Correlation
            W = ones(Sy,Sx);
            Wden = sqrt(P21.*conj(P21));
            W(P21~=0) = Wden(P21~=0);
            R = P21./W;

            %Robust Phase Correlation with spectral energy filter
            G = ifftn(R.*spectral,'symmetric');
            G = G(fftindy,fftindx);
            G = abs(G);
            
            %store correlation matrix
            CC(:,:,n) = G;

        end
end

eltime=cputime-t1;
fprintf('%0.2i:%0.2i.%0.0f ',floor(eltime/60),floor(rem(eltime,60)),rem(eltime,60)-floor(rem(eltime,60)))


%% GRID GENERATION SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [X,Y]=IMgrid(L,S,G)

%grid buffer
if nargin<3
    G=[0 0 0 0];
end

S=[S(2) S(1)];
G=[G(2) G(1) L(1)-G(2)+1 L(2)-G(1)+1];

%form grid
if max(S)==0
    %pixel grid
    y=(1:L(1))';
    x=1:L(2);
else
    if G(1)==0
        %buffers 1/2 grid spacing
        y=(ceil((L(1)-(floor(L(1)/S(1))-2)*S(1))/2):S(1):(L(1)-S(1)))';
    else
        %predefined grid buffer
        y=(G(1):S(1):G(3))';
    end
    if G(2)==0
        %buffers 1/2 grid spacing
        x=ceil((L(2)-(floor(L(2)/S(2))-2)*S(2))/2):S(2):(L(2)-S(2));
    else
        %predefined grid buffer
        x=(G(2):S(2):G(4));
    end
end

%vector2matrix conversion
X=x(ones(length(y),1),:);
Y=y(:,ones(1,length(x)));

% Y=Y+G(1);
% X=X+G(2);

%% VELOCITY INTERPOLATION SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ZI]=VFinterp(X,Y,Z,XI,YI,M)

%find grid sizes
Method={'nearest','linear','cubic'};
L=[max(max(YI)) max(max(XI))];
S=size(X);

%buffer matrix with nearest neighbor approximation for image boundaries
Xf = [1 X(1,:) L(2); ones(S(1),1) X L(2)*ones(S(1),1); 1 X(1,:) L(2);];
Yf = [ones(1,S(2)+2); Y(:,1) Y Y(:,1); L(1)*ones(1,S(2)+2)];
Zf = zeros(S+2);
Zf(2:end-1,2:end-1)=Z;
Zf(1,2:end-1)   = (Z(1,:)-Z(2,:))./(Y(1,:)-Y(2,:)).*(1-Y(2,:))+Z(1,:);
Zf(end,2:end-1) = (Z(end,:)-Z(end-1,:))./(Y(end,:)-Y(end-1,:)).*(L(1)-Y(end-1,:))+Z(end,:);
Zf(2:end-1,1)   = (Z(:,1)-Z(:,2))./(X(:,1)-X(:,2)).*(1-X(:,2))+Z(:,1);
Zf(2:end-1,end) = (Z(:,end)-Z(:,end-1))./(X(:,end)-X(:,end-1)).*(L(2)-X(:,end-1))+Z(:,end);
Zf(1,1)     = mean([Zf(2,1) Zf(1,2)]);
Zf(end,1)   = mean([Zf(end-1,1) Zf(end,2)]);
Zf(1,end)   = mean([Zf(2,end) Zf(1,end-1)]);
Zf(end,end) = mean([Zf(end-1,end) Zf(end,end-1)]);

%velocity interpolation
ZI=interp2(Xf,Yf,Zf,XI,YI,char(Method(M)));


%% VELOCITY SMOOTHING SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Uf,Vf]=VELfilt(U,V,C)

%2D gaussian filtering
A=fspecial('gaussian',[7 7],C);
Uf=imfilter(U,A,'replicate');
Vf=imfilter(V,A,'replicate');


%% GAUSSIAN WINDOW MASK SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [W]=windowmask(N,R)

% %generic indices
x  = -1:2/(N(1)-1):1;
y  = (-1:2/(N(2)-1):1)';
% 
% %gaussian window sizes
% px = (1.224*N(1)/R(1))^1.0172;
% py = (1.224*N(2)/R(2))^1.0172;
[px]=findwidth(R(1)/N(1));
[py]=findwidth(R(2)/N(2));
% 
% %generate 2D window
wx=exp(-px^2.*x.^2/2);
wy=exp(-py^2.*y.^2/2);

W  = wy*wx;


%% RPC SPECTRAL FILTER SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [W]=energyfilt(Nx,Ny,d,q)

%assume no aliasing
if nargin<4
    q = 0;
end

%initialize indices
[k1,k2]=meshgrid(-pi:2*pi/Ny:pi-2*pi/Ny,-pi:2*pi/Nx:pi-2*pi/Nx);

%particle-image spectrum
Ep = (pi*255*d^2/8)^2*exp(-d^2*k1.^2/16).*exp(-d^2*k2.^2/16);

%aliased particle-image spectrum
Ea = (pi*255*d^2/8)^2*exp(-d^2*(k1+2*pi).^2/16).*exp(-d^2*(k2+2*pi).^2/16)+...
     (pi*255*d^2/8)^2*exp(-d^2*(k1-2*pi).^2/16).*exp(-d^2*(k2+2*pi).^2/16)+...
     (pi*255*d^2/8)^2*exp(-d^2*(k1+2*pi).^2/16).*exp(-d^2*(k2-2*pi).^2/16)+...
     (pi*255*d^2/8)^2*exp(-d^2*(k1-2*pi).^2/16).*exp(-d^2*(k2-2*pi).^2/16)+...
     (pi*255*d^2/8)^2*exp(-d^2*(k1+0*pi).^2/16).*exp(-d^2*(k2+2*pi).^2/16)+...
     (pi*255*d^2/8)^2*exp(-d^2*(k1+0*pi).^2/16).*exp(-d^2*(k2-2*pi).^2/16)+...
     (pi*255*d^2/8)^2*exp(-d^2*(k1+2*pi).^2/16).*exp(-d^2*(k2+0*pi).^2/16)+...
     (pi*255*d^2/8)^2*exp(-d^2*(k1-2*pi).^2/16).*exp(-d^2*(k2+0*pi).^2/16);

%noise spectrum
En = pi/4*Nx*Ny;

%DPIV SNR spectral filter
W  = Ep./((1-q)*En+(q)*Ea);
W  = W'/max(max(W));


%% 3 POINT GAUSSIAN SUBPIXEL ESTIMATOR SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [u,v]=subpixel(G,ccsizex,ccsizey,W)

%intialize indices
cc_x = -ccsizex/2:ccsizex/2-1;
cc_y = -ccsizey/2:ccsizey/2-1;

%find maximum correlation value
[M,I] = max(G(:));

if M==0
    %if correlation empty
    u=0;
    v=0;
else
    %find x and y indices
    shift_locy = 1+mod(I-1,ccsizey);
    shift_locx = ceil(I/ccsizey);
    
    %find subpixel displacement in x
    if shift_locx == 1
        %boundary condition 1
        shift_err =  G( shift_locy , shift_locx+1 )/M;
    elseif shift_locx == ccsizex
        %boundary condition 2
        shift_err = -G( shift_locy , shift_locx-1 )/M;
    elseif G( shift_locy , shift_locx+1 ) == 0
        %endpoint discontinuity 1
        shift_err = -G( shift_locy , shift_locx-1 )/M;
    elseif G( shift_locy , shift_locx-1 ) == 0
        %endpoint discontinuity 2
        shift_err =  G( shift_locy , shift_locx+1 )/M;
    else
        %gaussian fit
        lCm1 = log(G( shift_locy , shift_locx-1 )*W( shift_locy , shift_locx-1 ));
        lC00 = log(G( shift_locy , shift_locx   )*W( shift_locy , shift_locx   ));
        lCp1 = log(G( shift_locy , shift_locx+1 )*W( shift_locy , shift_locx+1 ));
        if (2*(lCm1+lCp1-2*lC00)) == 0
            shift_err = 0;
        else
            shift_err = (lCm1-lCp1)/(2*(lCm1+lCp1-2*lC00));
        end
    end

    %add subpixel to discete pixel value
    u = cc_x(shift_locx) + shift_err;

    %find subpixel displacement in y
    if shift_locy == 1
        %boundary condition 1
        shift_err = -G( shift_locy+1 , shift_locx )/M;
    elseif shift_locy == ccsizey
        %boundary condition 2
        shift_err =  G( shift_locy-1 , shift_locx )/M;
    elseif G( shift_locy+1 , shift_locx ) == 0
        %endpoint discontinuity 1
        shift_err =  G( shift_locy-1 , shift_locx )/M;
    elseif G( shift_locy-1 , shift_locx ) == 0
        %endpoint discontinuity 2
        shift_err = -G( shift_locy+1 , shift_locx )/M;
    else
        %gaussian fit
        lCm1 = log(G( shift_locy-1 , shift_locx )*W( shift_locy-1 , shift_locx ));
        lC00 = log(G( shift_locy   , shift_locx )*W( shift_locy   , shift_locx ));
        lCp1 = log(G( shift_locy+1 , shift_locx )*W( shift_locy+1 , shift_locx ));
        if (2*(lCm1+lCp1-2*lC00)) == 0
            shift_err = 0;
        else
            shift_err = (lCm1-lCp1)/(2*(lCm1+lCp1-2*lC00));
        end
    end
    
    %add subpixel to discete pixel value
    v = (cc_y(shift_locy) + shift_err);

end


%% VALIDATION SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Uf,Vf,Eval] = VAL(xin,yin,uin,vin,t,tol,uthreshold,vthreshold,evalin)

%output text
fprintf('validating...                    ')
t1=cputime;

%preallocate evaluation matrix
if nargin<=8
    evalin = zeros(size(xin));
end

%neglect u and v threshold
if nargin<=6
    uthreshold = [-inf inf];
    vthreshold = [-inf inf];
end

%number of validation passes
pass = length(tol);

%vector size
Sin=size(xin);

%convert to matrix
if Sin(2)==1
    [X,Y,U,V,Eval]=matrixform(xin,yin,uin,vin,evalin);
    S=size(X);
else
    X=xin;
    Y=yin;
    U=uin;
    V=vin;
    Eval=evalin;
    S=Sin;
end

% [Evalprofile]=VALprofile2(U,V,2);
% Eval(Evalprofile>0)=100;
% U(Evalprofile>0)=nan;
% V(Evalprofile>0)=nan;

%thresholding
for i=1:S(1)
    for j=1:S(2)
        if Eval(i,j)==0
            %velocity threshold condition
            if U(i,j)<uthreshold(1) || U(i,j)>uthreshold(2) || V(i,j)<vthreshold(1) || V(i,j)>vthreshold(2)
                U(i,j)=nan;
                V(i,j)=nan;
                Eval(i,j)=100;
            end
        elseif Eval(i,j)==-1
            %boundary condition
            U(i,j)=nan;
            V(i,j)=nan;
        end
    end
end

%outlier searching
for k=1:pass
    
    q = (t(k,:)-1)/2;
    
    for i=1:S(1)
        for j=1:S(2)

            if Eval(i,j)==0
                
                %get statistical velocity evaluation block
                Imin = max([i-q(2) 1   ]);
                Imax = min([i+q(2) S(1)]);
                Jmin = max([j-q(1) 1   ]);
                Jmax = min([j+q(1) S(2)]);
                Iind = Imin:Imax;
                Jind = Jmin:Jmax;
                Ublock = U(Iind,Jind);
                Vblock = V(Iind,Jind);

                %universal outlier detection
                Ipos = find(Iind==i);
                Jpos = find(Jind==j);
                [Ru]=UOD(Ublock,Ipos,Jpos);
                [Rv]=UOD(Vblock,Ipos,Jpos);

                if Ru > tol(k) || Rv > tol(k)
                    %UOD threshold condition
                    U(i,j)=nan;
                    V(i,j)=nan;
                    Eval(i,j)=k;
                end

            end

        end
    end
end

%initialize output velocity
Uf=U;
Vf=V;

%replacement
for i=1:S(1)
    for j=1:S(2)
        
        if Eval(i,j) ~= 0

            %initialize replacement search size
            q=0;
            s=0;

            %get replacement block with at least 8 valid points
            while s==0
                q=q+1;
                Imin = max([i-q 1   ]);
                Imax = min([i+q S(1)]);
                Jmin = max([j-q 1   ]);
                Jmax = min([j+q S(2)]);
                Iind = Imin:Imax;
                Jind = Jmin:Jmax;
                Ublock = U(Iind,Jind);
                if length(Ublock(~isnan(Ublock)))>=8
                    Xblock = X(Iind,Jind)-X(i,j);
                    Yblock = Y(Iind,Jind)-Y(i,j);
                    Vblock = V(Iind,Jind);
                    s=1;
                end
            end
            
            %distance from erroneous vector
            Dblock = (Xblock.^2+Yblock.^2).^0.5;
            Dblock(isnan(Ublock))=nan;

            %validated vector
            Uf(i,j) = nansum(nansum(Dblock.*Ublock))/nansum(nansum(Dblock));
            Vf(i,j) = nansum(nansum(Dblock.*Vblock))/nansum(nansum(Dblock));
            
        end

    end
end

if Sin(2)==1
    %convert back to vector
    [Uf,Vf,Eval]=vectorform(xin,yin,Uf,Vf,Eval);
end

eltime=cputime-t1;
fprintf('%0.2i:%0.2i.%0.0f\n',floor(eltime/60),floor(rem(eltime,60)),rem(eltime,60)-floor(rem(eltime,60)))


%% UNVERSAL OUTLIER DETECTION SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [R]=UOD(W,p,q)

%minimum variance assumption
e=0.1; 

%remove value from query point
x=W(p,q);
W(p,q)=nan;

%remove any erroneous points
P=W(:);
Ps = sort(P);
Psfull = Ps(~isnan(Ps));
N=length(Psfull);

if N<=floor(length(W)/3)
    %return negative threshold value if no valid vectors in block
    R = inf;
else
    %return the median deviation normalized to the MAD
    if mod(N,2)==0
        M = (Psfull(N/2)+Psfull(N/2+1))/2;
        MADfull = sort(abs(Psfull-M));
        Q = (MADfull(N/2)+MADfull(N/2+1))/2;
        R = abs(x-M)/(Q+e);
    else
        M = Psfull((N+1)/2);
        MADfull = sort(abs(Psfull-M));
        Q = MADfull((N+1)/2);
        R = abs(x-M)/(Q+e);
    end
end


%% VECTOR TO MATRIX SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [X,Y,U,V,Eval]=matrixform(x,y,u,v,eval)

%find unique x and y grid points
a=sort(unique(x));
b=sort(unique(y));
N=length(x);

%initialize matrices
U=nan(length(b),length(a));
V=nan(length(b),length(a));
Eval=-1*ones(length(b),length(a));

%generate grid matrix
[X,Y]=meshgrid(a,b);

%generate variable matrices (nans where no data available)
for n=1:N
    I=find(b==y(n));
    J=find(a==x(n));
    U(I,J) = u(n);
    V(I,J) = v(n);
    Eval(I,J) = eval(n);
end


%% MATRIX TO VECTOR SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [u,v,eval]=vectorform(x,y,U,V,Eval)

%find unique x and y grid points
a=sort(unique(x));
b=sort(unique(y));
N=length(x);

%initialize vectors
S=size(x);
u    = zeros(S);
v    = zeros(S);
eval = zeros(S);

%generate data vectors where data is available
for n=1:N
    I=find(b==y(n));
    J=find(a==x(n));
    u(n)    = U(I,J);
    v(n)    = V(I,J);
    eval(n) = Eval(I,J);
end


%% PLT WRITER SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function []=write_plt_val(fname,X,Y,U,V,Eval,Mag,dt,T,title)

%output text
fprintf('saving...                        ')
t1=cputime;

%convert to physical units
X=X*Mag;
Y=Y*Mag;
U=U*Mag/dt;
V=V*Mag/dt;

%convert to matrix if necessary
if size(X,2)==1
    [X,Y,U,V,Eval]=matrixform(X,Y,U,V,Eval);
end

%remove nans from data, replace with zeros
U(Eval<0)=0;
V(Eval<0)=0;

%find I,J for plt
S = size(U);

%generate text file
fid = fopen(fname,'w');
if fid==-1
    error(['error creating file ',fname])
end

%header lines
fprintf(fid,['TITLE        = "' title '"\n']);
fprintf(fid,'VARIABLES    = "X", "Y", "U", "V", "Eval"\n');
fprintf(fid,'ZONE T="Time=%0.6f" I=%i J=%i C=BLACK STRANDID=1 SOLUTIONTIME = %0.6f\n',T,S(2),S(1),T);

%write data
for i=1:S(1)
    for j=1:S(2)
        if isnan(U(i,j)) || isnan(V(i,j))
            %second check to ensure no nans present
            fprintf(fid,'%14.6e %14.6e %14.6e %14.6e %14.6e\n',X(i,j),Y(i,j),0,0,-1);
        else
            %valid data points
            fprintf(fid,'%14.6e %14.6e %14.6e %14.6e %14.6e\n',X(i,j),Y(i,j),U(i,j),V(i,j),Eval(i,j));
        end
    end
end

fclose(fid);

%output text
eltime=cputime-t1;
fprintf('%0.2i:%0.2i.%0.0f\n',floor(eltime/60),floor(rem(eltime,60)),rem(eltime,60)-floor(rem(eltime,60)))


%% window size interpolation function
function [p]=findwidth(r)

R = [0.0000 0.0051 0.0052 0.0053 0.0055 0.0056 0.0057 0.0059 0.0060 ...
     0.0063 0.0064 0.0066 0.0067 0.0069 0.0070 0.0072 0.0074 0.0076 ...
     0.0079 0.0081 0.0083 0.0085 0.0087 0.0089 0.0091 0.0093 0.0095 ...
     0.0100 0.0102 0.0104 0.0107 0.0109 0.0112 0.0114 0.0117 0.0120 ...
     0.0125 0.0128 0.0131 0.0134 0.0137 0.0141 0.0144 0.0147 0.0151 ...
     0.0158 0.0161 0.0165 0.0169 0.0173 0.0177 0.0181 0.0185 0.0190 ...
     0.0199 0.0203 0.0208 0.0213 0.0218 0.0223 0.0228 0.0233 0.0239 ...
     0.0250 0.0256 0.0262 0.0268 0.0274 0.0281 0.0287 0.0294 0.0301 ...
     0.0315 0.0322 0.0330 0.0337 0.0345 0.0353 0.0361 0.0370 0.0378 ...
     0.0396 0.0406 0.0415 0.0425 0.0435 0.0445 0.0455 0.0466 0.0476 ...
     0.0499 0.0511 0.0522 0.0535 0.0547 0.0560 0.0573 0.0586 0.0600 ...
     0.0628 0.0643 0.0658 0.0673 0.0689 0.0705 0.0721 0.0738 0.0755 ...
     0.0791 0.0809 0.0828 0.0847 0.0867 0.0887 0.0908 0.0929 0.0951 ...
     0.0996 0.1019 0.1042 0.1067 0.1092 0.1117 0.1143 0.1170 0.1197 ...
     0.1253 0.1283 0.1312 0.1343 0.1374 0.1406 0.1439 0.1473 0.1507 ...
     0.1578 0.1615 0.1652 0.1691 0.1730 0.1770 0.1812 0.1854 0.1897 ...
     0.1986 0.2033 0.2080 0.2128 0.2178 0.2229 0.2281 0.2334 0.2388 ...
     0.2501 0.2559 0.2619 0.2680 0.2742 0.2806 0.2871 0.2938 0.3006 ...
     0.3148 0.3221 0.3296 0.3373 0.3451 0.3531 0.3613 0.3696 0.3781 ...
     0.3957 0.4048 0.4140 0.4233 0.4329 0.4425 0.4524 0.4623 0.4724 ...
     0.4930 0.5034 0.5139 0.5244 0.5351 0.5457 0.5564 0.5672 0.5779 ...
     0.5992 0.6099 0.6204 0.6309 0.6414 0.6517 0.6619 0.6720 0.6819 ...
     0.7014 0.7109 0.7203 0.7295 0.7385 0.7473 0.7559 0.7643 0.7726 ...
     0.7884 0.7960 0.8035 0.8107 0.8177 0.8245 0.8311 0.8376 0.8438 ...
     0.8556 0.8613 0.8667 0.8720 0.8771 0.8820 0.8867 0.8913 0.8957 ...
     0.9041 0.9080 0.9118 0.9155 0.9190 0.9224 0.9256 0.9288 0.9318 ...
     0.9374 0.9401 0.9426 0.9451 0.9474 0.9497 0.9519 0.9539 0.9559 ...
     0.9597 0.9614 0.9631 0.9647 0.9662 0.9677 0.9691 0.9705 0.9718 ...
     0.9742 0.9753 0.9764 0.9775 0.9785 0.9794 0.9803 0.9812 0.9820 ...
     0.9836 0.9843 0.9850 0.9857 0.9863 0.9869 0.9875 0.9881 0.9886 ...
     0.9896 0.9900 0.9905 0.9909 0.9913 0.9917 0.9921 0.9924 0.9928 ...
     0.9934 0.9937 0.9940 0.9943 0.9945 0.9948 0.9950 1.0000]';
 
P = [500.0000 245.4709 239.8833 234.4229 229.0868 223.8721 218.7762 213.7962 208.9296 ...
     199.5262 194.9845 190.5461 186.2087 181.9701 177.8279 173.7801 169.8244 165.9587 ...
     158.4893 154.8817 151.3561 147.9108 144.5440 141.2538 138.0384 134.8963 131.8257 ...
     125.8925 123.0269 120.2264 117.4898 114.8154 112.2018 109.6478 107.1519 104.7129 ...
     100.0000  97.7237  95.4993  93.3254  91.2011  89.1251  87.0964  85.1138  83.1764 ...
      79.4328  77.6247  75.8578  74.1310  72.4436  70.7946  69.1831  67.6083  66.0693 ...
      63.0957  61.6595  60.2560  58.8844  57.5440  56.2341  54.9541  53.7032  52.4807 ...
      50.1187  48.9779  47.8630  46.7735  45.7088  44.6684  43.6516  42.6580  41.6869 ...
      39.8107  38.9045  38.0189  37.1535  36.3078  35.4813  34.6737  33.8844  33.1131 ...
      31.6228  30.9030  30.1995  29.5121  28.8403  28.1838  27.5423  26.9153  26.3027 ...
      25.1189  24.5471  23.9883  23.4423  22.9087  22.3872  21.8776  21.3796  20.8930 ...
      19.9526  19.4984  19.0546  18.6209  18.1970  17.7828  17.3780  16.9824  16.5959 ...
      15.8489  15.4882  15.1356  14.7911  14.4544  14.1254  13.8038  13.4896  13.1826 ...
      12.5893  12.3027  12.0226  11.7490  11.4815  11.2202  10.9648  10.7152  10.4713 ...
      10.0000   9.7724   9.5499   9.3325   9.1201   8.9125   8.7096   8.5114   8.3176 ...
       7.9433   7.7625   7.5858   7.4131   7.2444   7.0795   6.9183   6.7608   6.6069 ...
       6.3096   6.1660   6.0256   5.8884   5.7544   5.6234   5.4954   5.3703   5.2481 ...
       5.0119   4.8978   4.7863   4.6774   4.5709   4.4668   4.3652   4.2658   4.1687 ...
       3.9811   3.8905   3.8019   3.7154   3.6308   3.5481   3.4674   3.3884   3.3113 ...
       3.1623   3.0903   3.0200   2.9512   2.8840   2.8184   2.7542   2.6915   2.6303 ...
       2.5119   2.4547   2.3988   2.3442   2.2909   2.2387   2.1878   2.1380   2.0893 ...
       1.9953   1.9498   1.9055   1.8621   1.8197   1.7783   1.7378   1.6982   1.6596 ...
       1.5849   1.5488   1.5136   1.4791   1.4454   1.4125   1.3804   1.3490   1.3183 ...
       1.2589   1.2303   1.2023   1.1749   1.1482   1.1220   1.0965   1.0715   1.0471 ...
       1.0000   0.9772   0.9550   0.9333   0.9120   0.8913   0.8710   0.8511   0.8318 ...
       0.7943   0.7762   0.7586   0.7413   0.7244   0.7079   0.6918   0.6761   0.6607 ...
       0.6310   0.6166   0.6026   0.5888   0.5754   0.5623   0.5495   0.5370   0.5248 ...
       0.5012   0.4898   0.4786   0.4677   0.4571   0.4467   0.4365   0.4266   0.4169 ...
       0.3981   0.3890   0.3802   0.3715   0.3631   0.3548   0.3467   0.3388   0.3311 ...
       0.3162   0.3090   0.3020   0.2951   0.2884   0.2818   0.2754   0.2692   0.2630 ...
       0.2512   0.2455   0.2399   0.2344   0.2291   0.2239   0.2188   0.2138   0.2089 ...
       0.1995   0.1950   0.1905   0.1862   0.1820   0.1778   0.1738   0.0000]';
 
p=interp1q(R,P,r);

%% axial profile validation thresholding

function [Eval]=VALprofile2(U,V,thresh)

um=median(U,2);
Um=um(:,ones(size(U,2),1));
Err=sqrt((U-Um).^2+V.^2);
Eval=zeros(size(U));
Eval(Err>thresh)=100;

