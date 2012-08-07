function [X,Y,U,V,C,t_opt]=PIVphasecorr(im1,im2,window,res,zpad,D,Zeromean,Peakswitch,X,Y,Uin,Vin,dt)
% --- DPIV Correlation ---

%convert input parameters
im1=double(im1);
im2=double(im2);
L=size(im1);

if nargin<13
    dt=1;
end

%convert to gridpoint list
X=X(:);
Y=Y(:);

%preallocate velocity fields and grid format
Nx = window(1);
Ny = window(2);
if nargin <=11 || isempty(Uin) || isempty(Vin)
    Uin = zeros(length(X),1);
    Vin = zeros(length(X),1);
end

U = zeros(length(X),1);
V = zeros(length(X),1);
C = zeros(length(X),3);
t_opt=zeros(length(X),1);

%sets up extended domain size
if zpad~=0
    Sy=2*Ny;
    Sx=2*Nx;
else
    Sy=Ny;
    Sx=Nx;
end

%RPC cutoff weighting
% wt = energyfilt(Nx,Ny,D,0);
% wtX=wt(Nx/2+1,:);
% wtY=wt(:,Ny/2+1)';
% cutoff=2/pi/D;
% wtX(wtX<cutoff)=0;
% wtY(wtY<cutoff)=0;
wt = energyfilt(Sx,Sy,D,0);
wtX=wt(Sy/2+1,:);
wtY=wt(:,Sx/2+1,:)';
cutoff=2/pi/D(1);
wtX(wtX<cutoff)=0;
cutoff=2/pi/D(2);
wtY(wtY<cutoff)=0;

lsqX = cell(size(im1,3),1);
lsqY = cell(size(im1,3),1);
for i=1:size(im1,3)
    if sum(Uin)==0 || sum(Vin)==0
        DT=dt(i);
    else
        DT=1;
    end
    lsqX{i}=(0:Sx-1).*DT-Sx/2*DT;
    lsqY{i}=(0:Sy-1).*DT-Sy/2*DT;
end

%window masking filter
sfilt1 = windowmask([Sx Sy],[res(1, 1) res(1, 2)]);
sfilt2 = windowmask([Sx Sy],[res(2, 1) res(2, 2)]);

%fftshift indicies
fftindy = [Sy/2+1:Sy 1:Sy/2];
fftindx = [Sx/2+1:Sx 1:Sx/2];

for n=1:length(X)
    S=zeros(2,size(im1,3));um_cum=[];vm_cum=[];wtX_cum=[];wtY_cum=[];lsqX_cum=[];lsqY_cum=[];t_good=[];
    for t=1:size(im1,3)     
        %apply the second order discrete window offset
        x1 = X(n) - floor(round(Uin(n))*dt(t)/2);
        x2 = X(n) +  ceil(round(Uin(n))*dt(t)/2);

        y1 = Y(n) - floor(round(Vin(n))*dt(t)/2);
        y2 = Y(n) +  ceil(round(Vin(n))*dt(t)/2);

        xmin1 = x1-Nx/2+1;
        xmax1 = x1+Nx/2;
        xmin2 = x2-Nx/2+1;
        xmax2 = x2+Nx/2;
        ymin1 = y1-Ny/2+1;
        ymax1 = y1+Ny/2;
        ymin2 = y2-Ny/2+1;
        ymax2 = y2+Ny/2;
    
        %find the image windows
        zone1 = im1( max([1 ymin1]):min([L(1) ymax1]),max([1 xmin1]):min([L(2) xmax1]),t);
        zone2 = im2( max([1 ymin2]):min([L(1) ymax2]),max([1 xmin2]):min([L(2) xmax2]),t);
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
        
        if Zeromean==1
            zone1=zone1-mean(zone1(:));
            zone2=zone2-mean(zone2(:));
        end

        %apply the image spatial filter
        region1 = (zone1).*sfilt1;
        region2 = (zone2).*sfilt2;

        %FFTs
        f1   = fftn(region1,[Sy Sx]);
        f2   = fftn(region2,[Sy Sx]);
        P21  = f2.*conj(f1);
        W = ones(Ny,Nx);
        Wden = sqrt(P21.*conj(P21));
        W(P21~=0) = Wden(P21~=0);
        R = P21./W;
        R = R(fftindy,fftindx);
        
        %SVD-based Phase Correlation
        [u,s,v]=svd(R);
        v=unwrap(angle(v(:,1)));
        um=(v-v(Sx/2+1))';
        u=unwrap(angle(u(:,1)));
        vm=(u-u(Sy/2+1))';
                
        S(:,t)=[s(1,1) s(2,2)]';
        if (S(1,t)/S(2,t)>1.5 && t>1) || t==1
            if sum(abs(Uin))>0
                um=um./dt(t);
                vm=vm./dt(t);
            end
            um_cum=[um_cum,um];
            vm_cum=[vm_cum,vm];
            lsqX_cum=[lsqX_cum,lsqX{t}];
            lsqY_cum=[lsqY_cum,lsqY{t}];

            t_opt(n)=t;
            t_good=[t_good,t];
            if t>1
                velmag=sqrt((U(n).*dt(1:t)).^2+(V(n).*dt(1:t)).^2);
%                 Qp=1-exp(-2).*(S(1,:)./S(2,:)-1).^-1./velmag; %Uncertainty-based Quality
                Qp=S(1,:)./S(2,:).*(1-0.1./velmag); %Persoons Quality
                for T=1:length(t_good)
                    %Normalized by quality
                    wtX_cum((1:Sx)+Sx*(T-1))=wtX.*(Qp(t_good(T))-min(Qp))./(max(Qp)-min(Qp));
                    wtY_cum((1:Sy)+Sy*(T-1))=wtY.*(Qp(t_good(T))-min(Qp))./(max(Qp)-min(Qp));
                    %Un-normalized
%                     wtX_cum((1:Sx)+Sx*(T-1))=wtX;
%                     wtY_cum((1:Sy)+Sy*(T-1))=wtY;
                end
            else
                wtX_cum=wtX;
                wtY_cum=wtY;
            end
            
        end

        fit= wlsq(um_cum,lsqX_cum,wtX_cum);
%         um_refined=unwrap_refine(um_cum,lsqX_cum,fit);
%         fit= wlsq(um_refined,lsqX_cum,wtX_cum).*Sx./2./pi;
        U(n)=fit(1)*Sx/2/pi;
        
        fit=-wlsq(vm_cum,lsqY_cum,wtY_cum);
%         vm_refined=unwrap_refine(vm_cum,lsqY_cum,fit);
%         fit=-wlsq(vm_refined,lsqY_cum,wtY_cum).*Sy./2./pi;
        V(n)=fit(1)*Sy/2/pi;

        if t<size(im1,3)
            %Displacement cutoff
            if U(n)*dt(t+1)>res(1)/4 || V(n)*dt(t+1)>res(2)/4
                break
            end
        end
    end
    if Peakswitch
        C(n,:)=diag(s(1:3,1:3));
    end
end
%add DWO to estimation
U = round(Uin)+U;
V = round(Vin)+V;

end