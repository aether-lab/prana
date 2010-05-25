function PIVadvance4code(Data)
%% --- Read Formatted Parameters ---

%input/output directory
if ispc
    imbase=[Data.imdirec '\' Data.imbase];
    maskbase=[Data.maskdirec '\' Data.maskbase];
    pltdirec=[Data.outdirec '\'];
else
    imbase=[Data.imdirec '/' Data.imbase];
    maskbase=[Data.maskdirec '/' Data.maskbase];
    pltdirec=[Data.outdirec '/'];
end

%image indices
I1 = str2double(Data.imfstart):str2double(Data.imfstep):str2double(Data.imfend);
I2 = I1+str2double(Data.imcstep);

%processing mask
if strcmp(Data.masktype,'none')
    mask = 1+0*double(imread([imbase sprintf(['%0.' Data.imzeros 'i.' Data.imext],I1(1))]));
elseif strcmp(Data.masktype,'static')
    mask = double(imread(Data.staticmaskname));
    mask = flipud(mask);
elseif strcmp(Data.masktype,'dynamic')
    maskfend=str2double(Data.maskfstart)+str2double(Data.maskfstep)*length(str2double(Data.imfstart):str2double(Data.imfstep):str2double(Data.imfend))-1;
    maskname=str2double(Data.maskfstart):str2double(Data.maskfstep):maskfend;
end

%method and passes
P=str2double(Data.passes);
Method={'Multipass','Multigrid','Deform','Ensemble','Multiframe'};
M=Method(str2double(Data.method));

%algorithm options
Velinterp=str2double(Data.velinterp);
Iminterp=str2double(Data.iminterp);
Nmax=str2double(Data.framestep);
ds=0.1;% PIV error (pix)

%physical parameters
Mag = str2double(Data.wrmag);
dt = str2double(Data.wrsep);
Freq = str2double(Data.wrsamp);

%initialization
Wres=zeros(P,2);
Wsize=zeros(P,2);
Gres=zeros(P,2);
Gbuf=zeros(P,2);
Corr=zeros(P,1);
D=zeros(P,1);
Velsmoothswitch=zeros(P,1);
Velsmoothfilt=zeros(P,1);
Valswitch=zeros(P,1);
UODswitch=zeros(P,1);
Bootswitch=zeros(P,1);
Threshswitch=zeros(P,1);
Writeswitch=zeros(P,1);
Peakswitch=zeros(P,1);
UODwinsize=zeros(P,2,1);
UODthresh=zeros(P,1);
Bootper=zeros(P,1);
Bootiter=zeros(P,1);
Bootkmax=zeros(P,1);
Uthresh=zeros(P,2);
Vthresh=zeros(P,2);
extrapeaks=zeros(P,1);
PeakNum=zeros(P,1);
PeakMag=zeros(P,1);
PeakVel=zeros(P,1);
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
    Corr(e)=str2double(A.corr)-1;
    D(e)=str2double(A.RPCd);
    Velsmoothswitch(e)=str2double(A.velsmooth);
    Velsmoothfilt(e)=str2double(A.velsmoothfilt);
    
    %validation and thresholding
    Valswitch(e)=str2double(A.val);
    UODswitch(e)=str2double(A.uod);
    Bootswitch(e)=str2double(A.bootstrap);
    Threshswitch(e)=str2double(A.thresh);
    Writeswitch(e)=str2double(A.write);

    vpass=[0 strfind(A.uod_window,';') length(A.uod_window)+1];
    for q=1:(length(vpass)-1)
        B=A.uod_window((vpass(q)+1):(vpass(q+1)-1));
        UODwinsize(e,:,q)=[str2double(B(1:(strfind(B,',')-1))) str2double(B((strfind(B,',')+1):end))];
        UODthresh(e,q)=str2double(A.uod_thresh(1+2*(q-1)));
    end
    
    Bootper(e)=str2double(A.bootstrap_percentsampled);
    Bootiter(e)=str2double(A.bootstrap_iterations);
    Bootkmax(e)=str2double(A.bootstrap_passes);
    
    if str2double(A.thresh)==1
        Uthresh(e,:)=[str2double(A.valuthresh(1:(strfind(A.valuthresh,',')-1))) str2double(A.valuthresh((strfind(A.valuthresh,',')+1):end))];
        Vthresh(e,:)=[str2double(A.valvthresh(1:(strfind(A.valvthresh,',')-1))) str2double(A.valvthresh((strfind(A.valvthresh,',')+1):end))];
    else
        Uthresh(e,:)=[-inf,inf];
        Vthresh(e,:)=[-inf,inf];
    end
    
    extrapeaks(e)=str2double(A.valextrapeaks);

    %peak information
    Peakswitch(e)=str2double(A.savepeakinfo);
    PeakNum(e)=str2double(A.corrpeaknum);
    PeakMag(e)=str2double(A.savepeakmag);
    PeakVel(e)=str2double(A.savepeakvel);
    
    %output directory
    wbase(e,:)={A.outbase};
    
end

%% --- Image Prefilter ---
%added to multipass, multigrid, deform, and ensemble for image loading
try
    if ispc
        IMmin=double(imread([Data.imdirec '\IMmin.tif']));
    else
        IMmin=double(imread([Data.imdirec '/IMmin.tif']));
    end
catch
    disp('Error reading image subtraction file: Resuming without image prefilter...')
    IMmin=0*double(imread([imbase sprintf(['%0.' Data.imzeros 'i.' Data.imext],I1(1))]));
end

%% --- Evaluate Image Sequence ---
switch char(M)

    case {'Multipass','Multigrid','Deform'}

        for q=1:length(I1)
            
            tf=cputime;
            
            %output text
            title=['Frame' sprintf(['%0.' Data.imzeros 'i'],I1(q)) ' and Frame' sprintf(['%0.' Data.imzeros 'i'],I2(q))];
            fprintf('\n----------------------------------------------------\n')
            fprintf(['Job: ',Data.batchname,'\n'])
            fprintf(['Processing ' title ' (' num2str(q) '/' num2str(length(I1)) ')\n'])
            fprintf('----------------------------------------------------\n')

            %load image pair and flip coordinates
            im1=double(imread([imbase sprintf(['%0.' Data.imzeros 'i.' Data.imext],I1(q))]))-IMmin;
            im2=double(imread([imbase sprintf(['%0.' Data.imzeros 'i.' Data.imext],I2(q))]))-IMmin;
            im1=flipud(im1);
            im2=flipud(im2);
            L=size(im1);
            
            %load dynamic mask and flip coordinates
            if strcmp(Data.masktype,'dynamic')
                mask = double(imread([maskbase sprintf(['%0.' Data.maskzeros 'i.' Data.maskext],maskname(q))]));
                mask = flipud(mask);
            end

            %initialize grid and evaluation matrix
            [XI,YI]=IMgrid(L,[0 0]);

            UI = BWO(1)*ones(size(XI));
            VI = BWO(2)*ones(size(YI));

            for e=1:P
                [X,Y]=IMgrid(L,Gres(e,:),Gbuf(e,:));
                S=size(X);X=X(:);Y=Y(:);
                
                Ub = reshape(downsample(downsample( UI(Y(1):Y(end),X(1):X(end)),Gres(e,2))',Gres(e,1))',length(X),1);
                Vb = reshape(downsample(downsample( VI(Y(1):Y(end),X(1):X(end)),Gres(e,2))',Gres(e,1))',length(X),1);
                maskds=downsample(downsample( mask(Y(1):Y(end),X(1):X(end)),Gres(e,2))',Gres(e,1))';
                Eval=reshape(maskds,length(X),1);
                Eval(Eval==0)=-1;
                Eval(Eval>0)=0;

                %correlate image pair
                [Xc,Yc,Uc,Vc,Cc]=PIVwindowed(im1,im2,Corr(e),Wsize(e,:),Wres(e,:),0,D(e),X(Eval>=0),Y(Eval>=0),Ub(Eval>=0),Vb(Eval>=0),Peakswitch(e) || (Valswitch(e) && extrapeaks(e)),ds);
                if Peakswitch(e) || (Valswitch(e) && extrapeaks(e))
                    U=zeros(size(X,1),3);
                    V=zeros(size(X,1),3);
                    C=zeros(size(X,1),3);
                    Eval=repmat(Eval,[1 3]);
                    C(Eval>=0)=Cc;
                else
                    U=zeros(size(X));V=zeros(size(X));C=[];
                end
                U(Eval>=0)=Uc;V(Eval>=0)=Vc;

                %validation
                if Valswitch(e)
                    %output text
                    fprintf('validating...                    ')
                    t1=cputime;
                    
                    [Uval,Vval,Evalval,Cval]=VAL(X,Y,U,V,Eval,C,Threshswitch(e),UODswitch(e),Bootswitch(e),extrapeaks(e),...
                        Uthresh(e,:),Vthresh(e,:),UODwinsize(e,:,:),UODthresh(e,UODthresh(e,:)~=0)',Bootper(e),Bootiter(e),Bootkmax(e));
                    
                    eltime=cputime-t1;
                    fprintf('%0.2i:%0.2i.%0.0f\n',floor(eltime/60),floor(rem(eltime,60)),rem(eltime,60)-floor(rem(eltime,60)))
                else
                    Uval=U(:,1);Vval=V(:,1);Evalval=Eval(:,1);
                    if ~isempty(C)
                        Cval=C(:,1);
                    else
                        Cval=[];
                    end
                end
                
                %write output
                if Writeswitch(e) 
                    if str2double(Data.datout) || str2double(Data.multiplematout)
                        fprintf('saving...                        ')
                        t1=cputime;
                    end         
                    if Peakswitch(e)
                        if PeakVel(e)
                            U=[Uval,U(:,1:PeakNum(e))];
                            V=[Vval,V(:,1:PeakNum(e))];
                            Eval=[Evalval,Eval(:,1:PeakNum(e))];
                        else
                            U=Uval; V=Vval;Eval=Evalval;
                        end
                        if PeakMag(e)
                            C=[Cval,C(:,1:PeakNum(e))];
                        else
                            C=Cval;
                        end
                    else
                        U=Uval; V=Vval; Eval=Evalval; C=Cval;
                    end

                    %convert to physical units
                    X=X*Mag;Y=Y*Mag;
                    U=U*Mag/dt;V=V*Mag/dt;

                    %convert to matrix if necessary
                    if size(X,2)==1
                        [X,Y,U,V,Eval,C]=matrixform(X,Y,U,V,Eval,C);
                    end

                    %remove nans from data, replace with zeros
                    U(Eval<0)=0;V(Eval<0)=0;
                    
                    if str2double(Data.datout)
                        time=(q-1)/Freq;
                        write_dat_val_C([pltdirec char(wbase(e,:)) sprintf(['%0.' Data.imzeros 'i.dat' ],I1(q))],X,Y,U,V,Eval,C,maskds,e,time,title);
                    end
                    
                    if str2double(Data.multiplematout)
                        save([pltdirec char(wbase(e,:)) sprintf(['%0.' Data.imzeros 'i.mat' ],I1(q))],'X','Y','U','V','Eval','C','maskds')
                    end

                    if str2double(Data.singlematout)
                        X_write{e}(:,:,q)=X;Y_write{e}(:,:,q)=Y;
                        U_write{e}(:,:,:,q)=U;V_write{e}(:,:,:,q)=V;
                        Eval_write{e}(:,:,:,q)=Eval;C_write{e}(:,:,:,q)=C;
                    end
                    
                    if str2double(Data.datout) || str2double(Data.multiplematout)
                        eltime=cputime-t1;
                        fprintf('%0.2i:%0.2i.%0.0f\n',floor(eltime/60),floor(rem(eltime,60)),rem(eltime,60)-floor(rem(eltime,60)))
                    end  
                end
                U=Uval; V=Vval;
        
                if e~=P
                    %reshape from list of grid points to matrix
                    X=reshape(X,[S(1),S(2)]);
                    Y=reshape(Y,[S(1),S(2)]);
                    U=reshape(U(:,1),[S(1),S(2)]);
                    V=reshape(V(:,1),[S(1),S(2)]);
                    
                    if strcmp(M,'Multigrid') || strcmp(M,'Deform')

                        fprintf('interpolating velocity...        ')
                        t1=cputime;

                        %velocity smoothing
                        if Velsmoothswitch(e)==1
                            [U,V]=VELfilt(U,V,Velsmoothfilt(e));
                        end

                        %velocity interpolation
                        UI = VFinterp(X,Y,U,XI,YI,Velinterp);
                        VI = VFinterp(X,Y,V,XI,YI,Velinterp);

                        eltime=cputime-t1;
                        fprintf('%0.2i:%0.2i.%0.0f\n',floor(eltime/60),floor(rem(eltime,60)),rem(eltime,60)-floor(rem(eltime,60)))
                        
                        if strcmp(M,'Deform')
                            fprintf('deforming images...              ')
                            t1=cputime;
                            
                            %translate pixel locations
                            XD1 = XI+UI/2;
                            YD1 = YI+VI/2;
                            XD2 = XI-UI/2;
                            YD2 = YI-VI/2;

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
                            im1d(im1d<0)=0; im1d(isnan(im1d))=0;
                            im2d(im2d<0)=0; im2d(isnan(im2d))=0;
                            
                            im1=im1d; im2=im2d;

                            eltime=cputime-t1;
                            fprintf('%0.2i:%0.2i.%0.0f\n',floor(eltime/60),floor(rem(eltime,60)),rem(eltime,60)-floor(rem(eltime,60)))
                        end
                    end
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
            [X,Y]=IMgrid(L,Gres(e,:),Gbuf(e,:));
            S=size(X);X=X(:);Y=Y(:);
            Ub = reshape(downsample(downsample( UI(Y(1):Y(end),X(1):X(end)),Gres(e,2))',Gres(e,1))',length(X),1);
            Vb = reshape(downsample(downsample( VI(Y(1):Y(end),X(1):X(end)),Gres(e,2))',Gres(e,1))',length(X),1);
            maskds=downsample(downsample( mask(Y(1):Y(end),X(1):X(end)),Gres(e,2))',Gres(e,1))';
            Eval=reshape(maskds,length(X),1);
            Eval(Eval==0)=-1;
            Eval(Eval>0)=0;
            
            if Peakswitch(e) || (Valswitch(e) && extrapeaks(e))
                U=zeros(size(X,1),3);
                V=zeros(size(X,1),3);
                C=zeros(size(X,1),3);
            else
                U=zeros(size(X));V=zeros(size(X));C=[];
            end
            
            %output text
            title=['Frame' sprintf(['%0.' Data.imzeros 'i'],I1(1)) ' to Frame' sprintf(['%0.' Data.imzeros 'i'],I2(end))];
            fprintf('\n----------------------------------------------------\n')
            fprintf(['Job: ',Data.batchname,'\n'])
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
            if Peakswitch(e) || (Valswitch(e) && extrapeaks(e))
                Uc=zeros(size(X,1),3);
                Vc=zeros(size(X,1),3);
                Cc=zeros(size(X,1),3);
                Ub=repmat(Ub,[1 3]);
                Vb=repmat(Vb,[1 3]);
                Eval=repmat(Eval,[1 3]);
            else
                Uc=zeros(size(X));Vc=zeros(size(X));Cc=[];
            end
            Z=size(CCm);
            ZZ=ones(Z(1),Z(2));
            for s=1:length(Xc)
                [Uc(s,:),Vc(s,:),Ctemp]=subpixel(CCm(:,:,s),Z(2),Z(1),ZZ,Peakswitch(e) || (Valswitch(e) && extrapeaks(e)));
                if ~isempty(Cc)
                    Cc(s,:)=Ctemp;
                end
            end

            U(Eval>=0)=Uc(Eval>=0)+round(Ub(Eval>=0));
            V(Eval>=0)=Vc(Eval>=0)+round(Vb(Eval>=0));
            if ~isempty(Cc)
                C(Eval>=0)=Cc(Eval>=0);
            end

            %validation
            if Valswitch(e)
                %output text
                fprintf('validating...                    ')
                t1=cputime;

                [Uval,Vval,Evalval,Cval]=VAL(X,Y,U,V,Eval,C,Threshswitch(e),UODswitch(e),Bootswitch(e),extrapeaks(e),...
                    Uthresh(e,:),Vthresh(e,:),UODwinsize(e,:,:),UODthresh(e,UODthresh(e,:)~=0)',Bootper(e),Bootiter(e),Bootkmax(e));

                eltime=cputime-t1;
                fprintf('%0.2i:%0.2i.%0.0f\n',floor(eltime/60),floor(rem(eltime,60)),rem(eltime,60)-floor(rem(eltime,60)))
            else
                Uval=U(:,1);Vval=V(:,1);Evalval=Eval(:,1);
                if ~isempty(C)
                    Cval=C(:,1);
                else
                    Cval=[];
                end
            end
                
            %write output
            if Writeswitch(e) 
                if str2double(Data.datout) || str2double(Data.multiplematout)
                    fprintf('saving...                        ')
                    t1=cputime;
                end

                if Peakswitch(e)
                    if PeakVel(e)
                        U=[Uval,U(:,1:PeakNum(e))];
                        V=[Vval,V(:,1:PeakNum(e))];
                        Eval=[Evalval,Eval(:,1:PeakNum(e))];
                    else
                        U=Uval; V=Vval; Eval=Evalval;
                    end
                    if PeakMag(e)
                        C=[Cval,C(:,1:PeakNum(e))];
                    else
                        C=Cval;
                    end
                else
                    U=Uval; V=Vval; Eval=Evalval; C=Cval;
                end

                %convert to physical units
                X=X*Mag;Y=Y*Mag;
                U=U*Mag/dt;V=V*Mag/dt;

                %convert to matrix if necessary
                if size(X,2)==1
                    [X,Y,U,V,Eval,C]=matrixform(X,Y,U,V,Eval,C);
                end

                %remove nans from data, replace with zeros
                U(Eval<0)=0;V(Eval<0)=0;

                if str2double(Data.datout)
                    time=(q-1)/Freq;
                    write_dat_val_C([pltdirec char(wbase(e,:)) sprintf(['%0.' Data.imzeros 'i.dat' ],I1(q))],X,Y,U,V,Eval,C,maskds,e,time,title);
                end
                if str2double(Data.multiplematout)
                    save([pltdirec char(wbase(e,:)) sprintf(['%0.' Data.imzeros 'i.mat' ],I1(q))],'X','Y','U','V','Eval','C','maskds')
                end

                if str2double(Data.singlematout)
                    X_write{e}(:,:,q)=X;Y_write{e}(:,:,q)=Y;
                    U_write{e}(:,:,:,q)=U;V_write{e}(:,:,:,q)=V;
                    Eval_write{e}(:,:,:,q)=Eval;C_write{e}(:,:,:,q)=C;
                    if strcmp(Data.masktype,'dynamic')
                        mask_write{e}(:,:,q)=maskds;
                    end
                end

                if str2double(Data.datout) || str2double(Data.multiplematout)
                    eltime=cputime-t1;
                    fprintf('%0.2i:%0.2i.%0.0f\n',floor(eltime/60),floor(rem(eltime,60)),rem(eltime,60)-floor(rem(eltime,60)))
                end  
            end
            U=Uval; V=Vval;
        
            if e~=P
                %reshape from list of grid points to matrix
                X=reshape(X,[S(1),S(2)]);
                Y=reshape(Y,[S(1),S(2)]);
                U=reshape(U(:,1),[S(1),S(2)]);
                V=reshape(V(:,1),[S(1),S(2)]);

                fprintf('interpolating velocity...        ')
                t1=cputime;

                %velocity smoothing
                if Velsmoothswitch(e)==1
                    [U,V]=VELfilt(U,V,Velsmoothfilt(e));
                end

                %velocity interpolation
                UI = VFinterp(X,Y,U,XI,YI,Velinterp);
                VI = VFinterp(X,Y,V,XI,YI,Velinterp);

                eltime=cputime-t1;
                fprintf('%0.2i:%0.2i.%0.0f\n',floor(eltime/60),floor(rem(eltime,60)),rem(eltime,60)-floor(rem(eltime,60)))
            end
        end
        
        
    case 'Multiframe'
        
        for q=1:length(I1)
            
            tf=cputime;
            
            %output text
            title=['Frame' sprintf(['%0.' Data.imzeros 'i'],I1(q)) ' and Frame' sprintf(['%0.' Data.imzeros 'i'],I2(q))];
            fprintf('\n----------------------------------------------------\n')
            fprintf(['Job: ',Data.batchname,'\n'])
            fprintf(['Processing ' title ' (' num2str(q) '/' num2str(length(I1)) ')\n'])
            fprintf('----------------------------------------------------\n')

            %load dynamic mask and flip coordinates
            if strcmp(Data.masktype,'dynamic')
                mask = double(imread([maskbase sprintf(['%0.' Data.maskzeros 'i.' Data.maskext],maskname(q))]));
                mask = flipud(mask);
            end
            
            %load image pairs and flip coordinates
            if q-Nmax<1 && q+Nmax>length(I1)
                j=min([q,length(I1)-q+1]);
            elseif q-Nmax<1
                j=q;
            elseif q+Nmax>length(I1)
                j=length(I1)-q+1;
            else
                j=Nmax+1;
            end
            im1=zeros(size(mask,1),size(mask,2),j); im2=im1;
            for i=1:j
                im1(:,:,i)=flipud(double(imread([imbase sprintf(['%0.' Data.imzeros 'i.' Data.imext],I1(q-i+1))]))-IMmin);
                im2(:,:,i)=flipud(double(imread([imbase sprintf(['%0.' Data.imzeros 'i.' Data.imext],I2(q+i-1))]))-IMmin);
            end
            L=size(im1);

            %initialize grid and evaluation matrix
            [XI,YI]=IMgrid(L,[0 0]);

            UI = zeros(size(XI));
            VI = zeros(size(YI));

            for e=1:P
                [X,Y]=IMgrid(L,Gres(e,:),Gbuf(e,:));
                X=X(:);Y=Y(:);
                
                Ub = reshape(downsample(downsample( UI(Y(1):Y(end),X(1):X(end)),Gres(e,2))',Gres(e,1))',length(X),1);
                Vb = reshape(downsample(downsample( VI(Y(1):Y(end),X(1):X(end)),Gres(e,2))',Gres(e,1))',length(X),1);
                maskds=downsample(downsample( mask(Y(1):Y(end),X(1):X(end)),Gres(e,2))',Gres(e,1))';
                Eval=reshape(maskds,length(X),1);
                Eval(Eval==0)=-1;
                Eval(Eval>0)=0;
                
                %correlate image pair
                [Xc,Yc,Uc,Vc,Cc]=PIVwindowed(im1,im2,Corr(e),Wsize(e,:),Wres(e,:),0,D(e),X(Eval>=0),Y(Eval>=0),Ub(Eval>=0),Vb(Eval>=0),Peakswitch(e) || (Valswitch(e) && extrapeaks(e)),ds);
                if Peakswitch(e) || (Valswitch(e) && extrapeaks(e))
                    U=zeros(size(X,1),3);
                    V=zeros(size(X,1),3);
                    C=zeros(size(X,1),3);
                    Eval=repmat(Eval,[1 3]);
                    C(Eval>=0)=Cc;
                else
                    U=zeros(size(X));V=zeros(size(X));C=[];
                end
                U(Eval>=0)=Uc;V(Eval>=0)=Vc;

                %validation
                if Valswitch(e)
                    %output text
                    fprintf('validating...                    ')
                    t1=cputime;
                    
                    [Uval,Vval,Evalval,Cval]=VAL(X,Y,U,V,Eval,C,Threshswitch(e),UODswitch(e),Bootswitch(e),extrapeaks(e),...
                        Uthresh(e,:),Vthresh(e,:),UODwinsize(e,:,:),UODthresh(e,UODthresh(e,:)~=0)',Bootper(e),Bootiter(e),Bootkmax(e));
                    
                    eltime=cputime-t1;
                    fprintf('%0.2i:%0.2i.%0.0f\n',floor(eltime/60),floor(rem(eltime,60)),rem(eltime,60)-floor(rem(eltime,60)))
                else
                    Uval=U(:,1);Vval=V(:,1);Evalval=Eval(:,1);
                    if ~isempty(C)
                        Cval=C(:,1);
                    else
                        Cval=[];
                    end
                end
                
                %write output
                if Writeswitch(e) 
                    if str2double(Data.datout) || str2double(Data.multiplematout)
                        fprintf('saving...                        ')
                        t1=cputime;
                    end
                                        
                    if Peakswitch(e)
                        if PeakVel(e)
                            U=[Uval,U(:,1:PeakNum(e))];
                            V=[Vval,V(:,1:PeakNum(e))];
                            Eval=[Evalval,Eval(:,1:PeakNum(e))];
                        else
                            U=Uval; V=Vval; Eval=Evalval;
                        end
                        if PeakMag(e)
                            C=[Cval,C(:,1:PeakNum(e))];
                        else
                            C=Cval;
                        end
                    else
                        U=Uval; V=Vval; Eval=Evalval; C=Cval;
                    end

                    %convert to physical units
                    X=X*Mag;Y=Y*Mag;
                    U=U*Mag/dt;V=V*Mag/dt;

                    %convert to matrix if necessary
                    if size(X,2)==1
                        [X,Y,U,V,Eval,C]=matrixform(X,Y,U,V,Eval,C);
                    end

                    %remove nans from data, replace with zeros
                    U(Eval<0)=0;V(Eval<0)=0;
                    
                    if str2double(Data.datout)
                        time=(q-1)/Freq;
                        write_dat_val_C([pltdirec char(wbase(e,:)) sprintf(['%0.' Data.imzeros 'i.dat' ],I1(q))],X,Y,U,V,Eval,C,maskds,e,time,title);
                    end
                    if str2double(Data.multiplematout)
                        save([pltdirec char(wbase(e,:)) sprintf(['%0.' Data.imzeros 'i.mat' ],I1(q))],'X','Y','U','V','Eval','C','maskds')
                    end

                    if str2double(Data.singlematout)
                        X_write{e}(:,:,q)=X;Y_write{e}(:,:,q)=Y;
                        U_write{e}(:,:,:,q)=U;V_write{e}(:,:,:,q)=V;
                        Eval_write{e}(:,:,:,q)=Eval;C_write{e}(:,:,:,q)=C;
                        if strcmp(Data.masktype,'dynamic')
                            mask_write{e}(:,:,q)=maskds;
                        end
                    end
                    
                    if str2double(Data.datout) || str2double(Data.multiplematout)
                        eltime=cputime-t1;
                        fprintf('%0.2i:%0.2i.%0.0f\n',floor(eltime/60),floor(rem(eltime,60)),rem(eltime,60)-floor(rem(eltime,60)))
                    end  
                end
            end
            
            eltime=cputime-tf;
            fprintf('total frame time...              %0.2i:%0.2i.%0.0f\n',floor(eltime/60),floor(rem(eltime,60)),rem(eltime,60)-floor(rem(eltime,60)))
            
        end
end

% Write Single .mat File
if str2double(Data.singlematout)
    fprintf('\n----------------------------------------------------\n')
    fprintf('saving .mat file(s)...           ')
    t1=cputime;
    for e=1:P
        X=X_write{e};Y=Y_write{e};
        U=U_write{e};V=V_write{e};
        Eval=Eval_write{e};C=C_write{e};
        save([pltdirec char(wbase(e,:)) 'all.mat'],'X','Y','U','V','Eval','C')
    end
    eltime=cputime-t1;
    fprintf('%0.2i:%0.2i.%0.0f\n',floor(eltime/60),floor(rem(eltime,60)),rem(eltime,60)-floor(rem(eltime,60)))
end

%signal job complete
beep,pause(0.2),beep

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                           END MAIN FUNCTION                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [X,Y,U,V,C]=PIVwindowed(im1,im2,corr,window,res,zpad,D,X,Y,Uin,Vin,Peakswitch,ds)
% --- DPIV Correlation ---
t1=cputime;

%convert input parameters
im1=double(im1);
im2=double(im2);
L = size(im1);
if size(im1,3)>1
    Peakreturn=1;
else
    Peakreturn=Peakswitch;
end

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

if Peakswitch
    Uin=repmat(Uin(:,1),[1 3]);
    Vin=repmat(Vin(:,1),[1 3]);
    U = zeros(length(X),3);
    V = zeros(length(X),3);
    C = zeros(length(X),3);
else
    U = zeros(length(X),1);
    V = zeros(length(X),1);
    C = [];
end

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

            for t=1:size(im1,3)
                dt=2*t-1;
                
                %find the image windows
                zone1 = im1( max([1 ymin1]):min([L(1) ymax1]),max([1 xmin1]):min([L(2) xmax1]), t);
                zone2 = im2( max([1 ymin2]):min([L(1) ymax2]),max([1 xmin2]):min([L(2) xmax2]), t);
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
                [Utemp(t,:),Vtemp(t,:),Ctemp(t,:)]=subpixel(G,Nx,Ny,cnorm,Peakreturn);
                Utemp(t,:)=Utemp(t,:)/dt;Vtemp(t,:)=Vtemp(t,:)/dt;
                if Peakreturn
                    velmag=sqrt(Utemp(t,1)^2+Vtemp(t,1)^2);
                    Qp(t)=Ctemp(t,1)/Ctemp(t,2)*(1-ds/velmag);
                else
                    Qp(t)=-1;
                end
            end
            [temp,t_opt]=max(Qp);
%             if t_opt~=1
%                 keyboard
%             end
            if Peakswitch
                U(n,:)=Utemp(t_opt,:);
                V(n,:)=Vtemp(t_opt,:);
                C(n,:)=Ctemp(t_opt,:);
            else
                U(n)=Utemp(t_opt,1);
                V(n)=Vtemp(t_opt,1);
            end
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
            
             for t=1:size(im1,3)
                dt=2*t-1;

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
                [Utemp(t,:),Vtemp(t,:),Ctemp(t,:)]=subpixel(G,Nx,Ny,cnorm,Peakreturn);
                Utemp(t,:)=Utemp(t,:)/dt;Vtemp(t,:)=Vtemp(t,:)/dt;
                if Peakreturn
                    velmag=sqrt(Utemp(t,1)^2+Vtemp(t,1)^2);
                    Qp(t)=Ctemp(t,1)/Ctemp(t,2)*(1-ds/velmag);
                else
                    Qp(t)=-1;
                end
             end
            [temp,t_opt]=max(Qp);
%             if t_opt~=1
%                 keyboard
%             end
            if Peakswitch
                U(n,:)=Utemp(t_opt,:);
                V(n,:)=Vtemp(t_opt,:);
                C(n,:)=Ctemp(t_opt,:);
            else
                U(n)=Utemp(t_opt,1);
                V(n)=Vtemp(t_opt,1);
            end
        end
end

%add DWO to estimation
U = round(Uin)+U;
V = round(Vin)+V;

eltime=cputime-t1;
fprintf('%0.2i:%0.2i.%0.0f\n',floor(eltime/60),floor(rem(eltime,60)),rem(eltime,60)-floor(rem(eltime,60)))

function [X,Y,CC]=PIVensemble(im1,im2,corr,window,res,zpad,D,X,Y,Uin,Vin)
% --- DPIV Ensemble Correlation ---

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

function [X,Y]=IMgrid(L,S,G)
% --- Grid Generation Subfunction ---

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

function [ZI]=VFinterp(X,Y,Z,XI,YI,M)
% --- Velocity Interpolation Subfunction

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

function [Uf,Vf]=VELfilt(U,V,C)
% --- Velocity Smoothing Subfunction ---

%2D gaussian filtering
A=fspecial('gaussian',[7 7],C);
Uf=imfilter(U,A,'replicate');
Vf=imfilter(V,A,'replicate');

function [W]=windowmask(N,R)
% --- Gaussian Window Mask Subfunction ---

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

function [W]=energyfilt(Nx,Ny,d,q)
% --- RPC Spectral Filter Subfunction ---

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

function [u,v,M]=subpixel(G,ccsizex,ccsizey,W,Peakswitch)
% --- 3 Point Gaussian subpixel Estimator Subfunction Plus Peak Return ---

%intialize indices
cc_x = -ccsizex/2:ccsizex/2-1;
cc_y = -ccsizey/2:ccsizey/2-1;

%find maximum correlation value
[M,I] = max(G(:));

%if correlation empty
if M==0
    if Peakswitch
        u=zeros(1,3);
        v=zeros(1,3);
        M=zeros(1,3);
    else
        u=0; v=0;
    end
else
    if Peakswitch
        A=imregionalmax(G);
        peakmat=G.*A;
        for i=2:3
            peakmat(peakmat==M(i-1))=0;
            [M(i),I(i)]=max(peakmat(:));
        end
        j=length(M);
    else
        j=1;    
    end

    for i=1:j
        %find x and y indices
        shift_locy = 1+mod(I(i)-1,ccsizey);
        shift_locx = ceil(I(i)/ccsizey);

        %find subpixel displacement in x
        if shift_locx == 1
            %boundary condition 1
            shift_err =  G( shift_locy , shift_locx+1 )/M(i);
        elseif shift_locx == ccsizex
            %boundary condition 2
            shift_err = -G( shift_locy , shift_locx-1 )/M(i);
        elseif G( shift_locy , shift_locx+1 ) == 0
            %endpoint discontinuity 1
            shift_err = -G( shift_locy , shift_locx-1 )/M(i);
        elseif G( shift_locy , shift_locx-1 ) == 0
            %endpoint discontinuity 2
            shift_err =  G( shift_locy , shift_locx+1 )/M(i);
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
        u(i) = cc_x(shift_locx) + shift_err;

        %find subpixel displacement in y
        if shift_locy == 1
            %boundary condition 1
            shift_err = -G( shift_locy+1 , shift_locx )/M(i);
        elseif shift_locy == ccsizey
            %boundary condition 2
            shift_err =  G( shift_locy-1 , shift_locx )/M(i);
        elseif G( shift_locy+1 , shift_locx ) == 0
            %endpoint discontinuity 1
            shift_err =  G( shift_locy-1 , shift_locx )/M(i);
        elseif G( shift_locy-1 , shift_locx ) == 0
            %endpoint discontinuity 2
            shift_err = -G( shift_locy+1 , shift_locx )/M(i);
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
        v(i) = (cc_y(shift_locy) + shift_err);
    end
end

function [Uval,Vval,Evalval,Cval]=VAL(X,Y,U,V,Eval,C,Threshswitch,UODswitch,Bootswitch,extrapeaks,Uthresh,Vthresh,UODwinsize,UODthresh,Bootper,Bootiter,Bootkmax)
% --- Validation Subfunction ---
if extrapeaks
    j=3;
else
    j=1;
end
Uval=U(:,1); Vval=V(:,1); Evalval=Eval(:,1);
if ~isempty(C)
    Cval=C(:,1);
end
for i=1:j
    if Threshswitch
        [Uval,Vval,Evalval] = Thresh(X,Y,Uval,Vval,Uthresh,Vthresh,Evalval);
    end
    if UODswitch
        t=permute(UODwinsize,[2 3 1]);
        t=t(:,t(1,:)~=0);
        [Uval,Vval,Evalval] = UOD(X,Y,Uval,Vval,t',UODthresh,Evalval);
    end
    if Bootswitch
        [Uval,Vval,Evalval] = bootstrapping(X,Y,Uval,Vval,Bootper,Bootiter,Bootkmax,Evalval);
    end
    disp(['';'Replaced ',num2str(sum(Evalval>0)),' Vectors'])
    if extrapeaks && i<3
        Uval(Evalval>0)=U(Evalval>0,i+1);
        Vval(Evalval>0)=V(Evalval>0,i+1);
        Evalval(Evalval>0)=Eval(Evalval>0,i+1);
        Cval(Evalval>0)=C(Evalval>0,i+1);
    end
end

function [Uf,Vf,Eval] = Thresh(xin,yin,uin,vin,uthreshold,vthreshold,evalin)
% --- Thresholding Validation Subfunction ---

%preallocate evaluation matrix
if nargin<=6
    evalin = zeros(size(xin));
end

%neglect u and v threshold
if nargin<=4
    uthreshold = [-inf inf];
    vthreshold = [-inf inf];
end

%vector size
Sin=size(xin);

%convert to matrix
if Sin(2)==1
    [X,Y,U,V,Eval]=matrixform(xin,yin,uin,vin,evalin,[]);
    S=size(X);
else
    U=uin;
    V=vin;
    Eval=evalin;
    S=Sin;
end

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

function [Uf,Vf,Eval] = UOD(xin,yin,uin,vin,t,tol,evalin)
% --- Universal Outlier Detection Validation Subfunction ---

%preallocate evaluation matrix
if nargin<=6
    evalin = zeros(size(xin));
end

%number of validation passes
pass = length(tol);

%vector size
Sin=size(xin);

%convert to matrix
if Sin(2)==1
    [X,Y,U,V,Eval]=matrixform(xin,yin,uin,vin,evalin,[]);
    S=size(X);
else
    X=xin;
    Y=yin;
    U=uin;
    V=vin;
    Eval=evalin;
    S=Sin;
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
                [Ru]=UOD_sub(Ublock,Ipos,Jpos);
                [Rv]=UOD_sub(Vblock,Ipos,Jpos);

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

function [R]=UOD_sub(W,p,q)
% --- Universal Outlier Detection Algorithm ---

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

function [U,V,Eval] = bootstrapping(xin,yin,uin,vin,per,iter,kmax,evalin)
% Bootstrapping Validation Subfunction 
%
% [U,V,Eval] = bootstraping(x,y,u,v,per,iter,kmax,Eval)
%
% per  = percent removed for each interpolation (0-1)
% iter = number of interpolations per frame (for histogram)
% kmax = number of passes 

if nargin<8
    evalin = zeros(size(xin));
end

%vector size
n_in=size(xin);

%convert to matrix
if n_in(2)==1
    [X,Y,U,V,Eval]=matrixform(xin,yin,uin,vin,evalin,[]);
    n=size(X);
else
    U=uin;
    V=vin;
    Eval=evalin;
    n=n_in;
end

M = zeros(n(1),n(2),iter);

tol = 0.3;
ktol = 1;

while tol > 0 && ktol <= kmax+1
    U = zeros(n(1),n(2),iter);
    V = zeros(n(1),n(2),iter);

    for i = 1:iter
        clear S m Up Vp Ui Vi
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Data Removal
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        [m]= bootstrapping_dataremove(size(U),per,sum(Eval,3));

        S(:,1) = X(m==1);
        S(:,2) = Y(m==1);

        Up = U(m==1);
        Vp = V(m==1);
        M(:,:,i) = m;
        
        Ui = gridfit(S(:,1),S(:,2),Up,X(1,:),Y(:,1));
        Vi = gridfit(S(:,1),S(:,2),Vp,X(1,:),Y(:,1));
        U(:,:,i) = Ui;
        V(:,:,i) = Vi;

    end

    PBad = 0;
    for j = 1:n(1)
        for k = 1:n(2)
            if sum(isnan(U(j,k,:))) == 0
                try              
                    [H.U,HX.U] = hist(U(j,k,:),iter/2);
                    [H.V,HX.V] = hist(V(j,k,:),iter/2);

                    modeU = HX.U(H.U==max(H.U));
                    modeV = HX.V(H.V==max(H.V));

                    tU = abs((modeU - U(j,k))/modeU);
                    tV = abs((modeV - V(j,k))/modeV);
                    if tU > tol || tV > tol && Eval(j,k) ~= -1
                        U(j,k) = modeU(1);
                        V(j,k) = modeV(1);
                        if Eval(j,k)<200
                            Eval(j,k) = 200;
                        else
                            Eval(j,k) = Eval(j,k)+1;
                        end
                        PBad = PBad+1;
                    end
                catch%#ok
                    Ems=lasterror;%#ok
                    fprintf('\n\n')
                    fprintf(Ems.message)
                end
            end
        end
    end
    ktol = ktol + 1;
    tol = tol-(tol/(kmax-1))*(ktol-1);
end

if n_in(2)==1
    %convert back to vector
    [U,V,Eval]=vectorform(xin,yin,U,V,Eval);
end

function [M1] = bootstrapping_dataremove(DSIZE,ENUM,MASK)
% --- Bootstrapping Data Removal ---

Nx = DSIZE(1);
Ny = DSIZE(2);
Nt = 1;

M1   = zeros(DSIZE);
RMAT = rand(Nx,Ny,Nt);
EN   = 0;

while sum(M1(:))/(Nx*Ny) < ENUM && EN < 1
    M1 = RMAT<EN;
    M1(MASK>0) = 0; 
    EN = EN + 0.005;    
end

M1 = double(M1);

function [zgrid,xgrid,ygrid] = gridfit(x,y,z,xnodes,ynodes,varargin)
% gridfit: estimates a surface on a 2d grid, based on scattered data
%          Replicates are allowed. All methods extrapolate to the grid
%          boundaries. Gridfit uses a modified ridge estimator to
%          generate the surface, where the bias is toward smoothness.
%
%          Gridfit is not an interpolant. Its goal is a smooth surface
%          that approximates your data, but allows you to control the
%          amount of smoothing.
%
% usage #1: zgrid = gridfit(x,y,z,xnodes,ynodes);
% usage #2: [zgrid,xgrid,ygrid] = gridfit(x,y,z,xnodes,ynodes);
% usage #3: zgrid = gridfit(x,y,z,xnodes,ynodes,prop,val,prop,val,...);
%
% Arguments: (input)
%  x,y,z - vectors of equal lengths, containing arbitrary scattered data
%          The only constraint on x and y is they cannot ALL fall on a
%          single line in the x-y plane. Replicate points will be treated
%          in a least squares sense.
%
%          ANY points containing a NaN are ignored in the estimation
%
%  xnodes - vector defining the nodes in the grid in the independent
%          variable (x). xnodes need not be equally spaced. xnodes
%          must completely span the data. If they do not, then the
%          'extend' property is applied, adjusting the first and last
%          nodes to be extended as necessary. See below for a complete
%          description of the 'extend' property.
%
%          If xnodes is a scalar integer, then it specifies the number
%          of equally spaced nodes between the min and max of the data.
%
%  ynodes - vector defining the nodes in the grid in the independent
%          variable (y). ynodes need not be equally spaced.
%
%          If ynodes is a scalar integer, then it specifies the number
%          of equally spaced nodes between the min and max of the data.

% set defaults
params.smoothness = 1;
params.interp = 'triangle';
params.regularizer = 'gradient';
params.solver = 'normal';
params.maxiter = [];
params.extend = 'warning';

% and check for any overrides
params = parse_pv_pairs(params,varargin);

% check the parameters for acceptability
% smoothness == 1 by default
if isempty(params.smoothness)
  params.smoothness = 1;
else
  if (params.smoothness<=0)
    error 'Smoothness must be real, finite, and positive.'
  end
end
% regularizer  - must be one of 4 options - the second and
% third are actually synonyms.
valid = {'springs', 'diffusion', 'laplacian', 'gradient'};
if isempty(params.regularizer)
  params.regularizer = 'diffusion';
end
ind = strmatch(lower(params.regularizer),valid);
if (length(ind)==1)
  params.regularizer = valid{ind};
else
  error(['Invalid regularization method: ',params.regularizer])
end

% interp must be one of:
%    'bilinear', 'nearest', or 'triangle'
% but accept any shortening thereof.
valid = {'bilinear', 'nearest', 'triangle'};
if isempty(params.interp)
  params.interp = 'triangle';
end
ind = strmatch(lower(params.interp),valid);
if (length(ind)==1)
  params.interp = valid{ind};
else
  error(['Invalid interpolation method: ',params.interp])
end

% solver must be one of:
%    'backslash', '\', 'symmlq', 'lsqr', or 'normal'
% but accept any shortening thereof.
valid = {'backslash', '\', 'symmlq', 'lsqr', 'normal'};
if isempty(params.solver)
  params.solver = '\';
end
ind = strmatch(lower(params.solver),valid);
if (length(ind)==1)
  params.solver = valid{ind};
else
  error(['Invalid solver option: ',params.solver])
end

% extend must be one of:
%    'never', 'warning', 'always'
% but accept any shortening thereof.
valid = {'never', 'warning', 'always'};
if isempty(params.extend)
  params.extend = 'warning';
end
ind = strmatch(lower(params.extend),valid);
if (length(ind)==1)
  params.extend = valid{ind};
else
  error(['Invalid extend option: ',params.extend])
end

% ensure all of x,y,z,xnodes,ynodes are column vectors,
% also drop any NaN data
x=x(:);
y=y(:);
z=z(:);
k = isnan(x) | isnan(y) | isnan(z);
if any(k)
  x(k)=[];
  y(k)=[];
  z(k)=[];
end
xmin = min(x);
xmax = max(x);
ymin = min(y);
ymax = max(y);

% did they supply a scalar for the nodes?
if length(xnodes)==1
  xnodes = linspace(xmin,xmax,xnodes)';
  xnodes(end) = xmax; % make sure it hits the max
end
if length(ynodes)==1
  ynodes = linspace(ymin,ymax,ynodes)';
  ynodes(end) = ymax; % make sure it hits the max
end

xnodes=xnodes(:);
ynodes=ynodes(:);
dx = diff(xnodes);
dy = diff(ynodes);
nx = length(xnodes);
ny = length(ynodes);
ngrid = nx*ny;

% default for maxiter?
if isempty(params.maxiter)
  params.maxiter = min(10000,nx*ny);
end

% check lengths of the data
n = length(x);
if (length(y)~=n)||(length(z)~=n)
  error 'Data vectors are incompatible in size.'
end
if n<3
  error 'Insufficient data for surface estimation.'
end

% verify the nodes are distinct
if any(diff(xnodes)<=0)||any(diff(ynodes)<=0)
  error 'xnodes and ynodes must be monotone increasing'
end

% do we need to tweak the first or last node in x or y?
if xmin<xnodes(1)
    xnodes(1) = xmin;
end
if xmax>xnodes(end)
    xnodes(end) = xmax;
end
if ymin<ynodes(1)
    ynodes(1) = ymin;
end
if ymax>ynodes(end)
    ynodes(end) = ymax;
end

% only generate xgrid and ygrid if requested.
if nargout>1
  [xgrid,ygrid]=meshgrid(xnodes,ynodes);
end

% determine which cell in the array each point lies in
[junk,indx] = histc(x,xnodes);
[junk,indy] = histc(y,ynodes);
% any point falling at the last node is taken to be
% inside the last cell in x or y.
k=(indx==nx);
indx(k)=indx(k)-1;
k=(indy==ny);
indy(k)=indy(k)-1;

% interpolation equations for each point
tx = min(1,max(0,(x - xnodes(indx))./dx(indx)));
ty = min(1,max(0,(y - ynodes(indy))./dy(indy)));
ind = indy + ny*(indx-1);
% Future enhancement: add cubic interpolant
switch params.interp
  case 'triangle'
    % linear interpolation inside each triangle
    k = (tx > ty);
    L = ones(n,1);
    L(k) = ny;
    
    t1 = min(tx,ty);
    t2 = max(tx,ty);
    A = sparse(repmat((1:n)',1,3),[ind,ind+ny+1,ind+L], ...
       [1-t2,t1,t2-t1],n,ngrid);
    
  case 'nearest'
    % nearest neighbor interpolation in a cell
    k = round(1-ty) + round(1-tx)*ny;
    A = sparse((1:n)',ind+k,ones(n,1),n,ngrid);
    
  case 'bilinear'
    % bilinear interpolation in a cell
    A = sparse(repmat((1:n)',1,4),[ind,ind+1,ind+ny,ind+ny+1], ...
       [(1-tx).*(1-ty), (1-tx).*ty, tx.*(1-ty), tx.*ty], ...
       n,ngrid);
    
end
rhs = z;

% Build regularizer. Add del^4 regularizer one day.
switch params.regularizer
  case 'springs'
    % zero "rest length" springs
    [i,j] = meshgrid(1:nx,1:(ny-1));
    ind = j(:) + ny*(i(:)-1);
    m = nx*(ny-1);
    stiffness = 1./dy;
    Areg = sparse(repmat((1:m)',1,2),[ind,ind+1], ...
       stiffness(j(:))*[-1 1],m,ngrid);
    
    [i,j] = meshgrid(1:(nx-1),1:ny);
    ind = j(:) + ny*(i(:)-1);
    m = (nx-1)*ny;
    stiffness = 1./dx;
    Areg = [Areg;sparse(repmat((1:m)',1,2),[ind,ind+ny], ...
       stiffness(i(:))*[-1 1],m,ngrid)];
    
    [i,j] = meshgrid(1:(nx-1),1:(ny-1));
    ind = j(:) + ny*(i(:)-1);
    m = (nx-1)*(ny-1);
    stiffness = 1./sqrt(dx(i(:)).^2 + dy(j(:)).^2);
    Areg = [Areg;sparse(repmat((1:m)',1,2),[ind,ind+ny+1], ...
       stiffness*[-1 1],m,ngrid)];
    
    Areg = [Areg;sparse(repmat((1:m)',1,2),[ind+1,ind+ny], ...
       stiffness*[-1 1],m,ngrid)];
    
  case {'diffusion' 'laplacian'}
    % thermal diffusion using Laplacian (del^2)
    [i,j] = meshgrid(1:nx,2:(ny-1));
    ind = j(:) + ny*(i(:)-1);
    dy1 = dy(j(:)-1);
    dy2 = dy(j(:));
    
    Areg = sparse(repmat(ind,1,3),[ind-1,ind,ind+1], ...
      [-2./(dy1.*(dy1+dy2)), 2./(dy1.*dy2), ...
       -2./(dy2.*(dy1+dy2))],ngrid,ngrid);
    
    [i,j] = meshgrid(2:(nx-1),1:ny);
    ind = j(:) + ny*(i(:)-1);
    dx1 = dx(i(:)-1);
    dx2 = dx(i(:));
    
    Areg = Areg + sparse(repmat(ind,1,3),[ind-ny,ind,ind+ny], ...
      [-2./(dx1.*(dx1+dx2)), 2./(dx1.*dx2), ...
       -2./(dx2.*(dx1+dx2))],ngrid,ngrid);
    
  case 'gradient'
    % Subtly different from the Laplacian. A point for future
    % enhancement is to do it better for the triangle interpolation
    % case.
    [i,j] = meshgrid(1:nx,2:(ny-1));
    ind = j(:) + ny*(i(:)-1);
    dy1 = dy(j(:)-1);
    dy2 = dy(j(:));

    Areg = sparse(repmat(ind,1,3),[ind-1,ind,ind+1], ...
      [-2./(dy1.*(dy1+dy2)), 2./(dy1.*dy2), ...
      -2./(dy2.*(dy1+dy2))],ngrid,ngrid);

    [i,j] = meshgrid(2:(nx-1),1:ny);
    ind = j(:) + ny*(i(:)-1);
    dx1 = dx(i(:)-1);
    dx2 = dx(i(:));

    Areg = [Areg;sparse(repmat(ind,1,3),[ind-ny,ind,ind+ny], ...
      [-2./(dx1.*(dx1+dx2)), 2./(dx1.*dx2), ...
      -2./(dx2.*(dx1+dx2))],ngrid,ngrid)];

end
nreg = size(Areg,1);

% Append the regularizer to the interpolation equations,
% scaling the problem first. Use the 1-norm for speed.
NA = norm(A,1);
NR = norm(Areg,1);
A = [A;Areg*(params.smoothness*NA/NR)];
rhs = [rhs;zeros(nreg,1)];

% solve the full system, with regularizer attached
switch params.solver
  case {'\' 'backslash'}
    % permute for minimum fill in for R (in the QR)
    p = colamd(A);
    zgrid=zeros(ny,nx);
    zgrid(p) = A(:,p)\rhs;
    
  case 'normal'
    % The normal equations, solved with \. Can be fast
    % for huge numbers of data points.
    
    % Permute for minimum fill-in for \ (in chol)
    APA = A'*A;
    p = symamd(APA);
    zgrid=zeros(ny,nx);
    zgrid(p) = APA(p,p)\(A(:,p)'*rhs);
    
  case 'symmlq'
    % iterative solver - symmlq - requires a symmetric matrix,
    % so use it to solve the normal equations. No preconditioner.
    tol = abs(max(z)-min(z))*1.e-13;
    [zgrid,flag] = symmlq(A'*A,A'*rhs,tol,params.maxiter);
    zgrid = reshape(zgrid,ny,nx);
    
    % display a warning if convergence problems
    switch flag
      case 0
        % no problems with convergence
      case 1
        % SYMMLQ iterated MAXIT times but did not converge.
        warning(['Symmlq performed ',num2str(params.maxiter), ...
          ' iterations but did not converge.'])
      case 3
        % SYMMLQ stagnated, successive iterates were the same
        warning 'Symmlq stagnated without apparent convergence.'
      otherwise
        warning(['One of the scalar quantities calculated in',...
          ' symmlq was too small or too large to continue computing.'])
    end
    
  case 'lsqr'
    % iterative solver - lsqr. No preconditioner here.
    tol = abs(max(z)-min(z))*1.e-13;
    [zgrid,flag] = lsqr(A,rhs,tol,params.maxiter);
    zgrid = reshape(zgrid,ny,nx);
    
    % display a warning if convergence problems
    switch flag
      case 0
        % no problems with convergence
      case 1
        % lsqr iterated MAXIT times but did not converge.
        warning(['Lsqr performed ',num2str(params.maxiter), ...
          ' iterations but did not converge.'])
      case 3
        % lsqr stagnated, successive iterates were the same
        warning 'Lsqr stagnated without apparent convergence.'
      case 4
        warning(['One of the scalar quantities calculated in',...
          ' LSQR was too small or too large to continue computing.'])
    end
    
end

function params=parse_pv_pairs(params,pv_pairs)
% parse_pv_pairs: parses sets of property value pairs
% usage: params=parse_pv_pairs(default_params,pv_pairs)
%
% arguments: (input)
%  default_params - structure, with one field for every potential
%             property/value pair. Each field will contain the default
%             value for that property. If no default is supplied for a
%             given property, then that field must be empty.
%
%  pv_array - cell array of property/value pairs.
%             Case is ignored when comparing properties to the list
%             of field names. Also, any unambiguous shortening of a
%             field/property name is allowed.
%
% arguments: (output)
%  params   - parameter struct that reflects any updated property/value
%             pairs in the pv_array.

npv = length(pv_pairs);
n = npv/2;

if n~=floor(n)
  error 'Property/value pairs must come in PAIRS.'
end
if n<=0
  % just return the defaults
  return
end

if ~isstruct(params)
  error 'No structure for defaults was supplied'
end

% there was at least one pv pair. process any supplied
propnames = fieldnames(params);
lpropnames = lower(propnames);
for i=1:n
  pi = lower(pv_pairs{2*i-1});
  vi = pv_pairs{2*i};
  
  ind = strmatch(pi,lpropnames,'exact');
  if isempty(ind)
    ind = strmatch(pi,lpropnames);
    if isempty(ind)
      error(['No matching property found for: ',pv_pairs{2*i-1}])
    elseif length(ind)>1
      error(['Ambiguous property name: ',pv_pairs{2*i-1}])
    end
  end
  pi = propnames{ind};
  
  % override the corresponding default in params
  params = setfield(params,pi,vi);
  
end

function [X,Y,U,V,Eval,C]=matrixform(x,y,u,v,eval,c)
% --- Vector to Matrix Subfunction ---

%find unique x and y grid points
a=sort(unique(x));
b=sort(unique(y));
N=length(x);

%initialize matrices
U=nan(length(b),length(a),size(u,2));
V=nan(length(b),length(a),size(v,2));
Eval=-1*ones(length(b),length(a),size(eval,2));

%generate grid matrix
[X,Y]=meshgrid(a,b);

%generate variable matrices (nans where no data available)
for i=1:size(U,3)
    for n=1:N
        I=find(b==y(n));
        J=find(a==x(n));
        U(I,J,i) = u(n,i);
        V(I,J,i) = v(n,i);
        Eval(I,J,i) = eval(n);
    end
end
if ~isempty(c)
    C=nan(length(b),length(a),size(c,2));
    for i=1:size(c,2)
        for n=1:N
            I= b==y(n);
            J= a==x(n);
            C(I,J,i)=c(n,i);
        end
    end
else
    C=[];
end

function [u,v,eval]=vectorform(x,y,U,V,Eval)
% --- Matrix to Vector Subfunction ---

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

function []=write_dat_val_C(fname,X,Y,U,V,Eval,C,mask,strand,T,title)
% --- .dat Writer Subfunction ---

%find I,J for plt
S = size(U);

%generate text file
fid = fopen(fname,'w');
if fid==-1
    error(['error creating file ',fname])
end

varlist='"X" "Y" "U" "V" "Eval"';
if ~isempty(C)
    varlist=[varlist,' "C"'];
    if size(U,3)>1
        for i=2:size(U,3)
            varlist=[varlist,' "U',num2str(i-1),'" "V',num2str(i-1),'" "C',num2str(i-1),'"'];
        end
    elseif size(C,3)>1
        for i=2:size(C,3)
            varlist=[varlist,' "C',num2str(i-1),'"'];
        end
    end
end
    

%header lines
fprintf(fid,['TITLE        = "' title '"\n']);
fprintf(fid,['VARIABLES    = ',varlist,'\n']);
fprintf(fid,'ZONE T="Time=%0.6f" I=%i J=%i C=BLACK STRANDID=%i SOLUTIONTIME = %0.6f\n',T,S(2),S(1),strand,T);
    

%write data
for i=1:S(1)
    for j=1:S(2)
        if isnan(U(i,j)) || isnan(V(i,j))
            %second check to ensure no nans present
            fprintf(fid,'%14.6e %14.6e %14.6e %14.6e %14.6e %14.6e',X(i,j),Y(i,j),0,0,-1);
        else
            %valid data points
            fprintf(fid,'%14.6e %14.6e %14.6e %14.6e %14.6e %14.6e',X(i,j),Y(i,j),U(i,j),V(i,j),Eval(i,j));
        end
        
        if ~isempty(C)
            if isnan(C(i,j,1))
                fprintf(fid,' %14.6e',0);
            else
                fprintf(fid,' %14.6e',C(i,j,1));
            end
            if size(U,3)>1
                for k=2:size(U,3)
                    if isnan(U(i,j,k)) || isnan(V(i,j,k)) || isnan(C(i,j,k))
                        fprintf(fid,' %14.6e %14.6e %14.6e',0,0,0);
                    else
                        fprintf(fid,' %14.6e %14.6e %14.6e',U(i,j,k),V(i,j,k),C(i,j,k));
                    end
                end
            elseif size(C,3)>1
                for k=2:size(C,3)
                    if isnan(C(i,j,k))
                        fprintf(fid,' %14.6e',0);
                    else
                        fprintf(fid,' %14.6e',C(i,j,k));
                    end
                end
            end
        end     
        fprintf(fid,'\n');
    end
end
    
fclose(fid);

function [p]=findwidth(r)
% --- Window Size Interpolation Function ---

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