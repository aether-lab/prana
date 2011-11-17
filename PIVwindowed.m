function [X,Y,U,V,C,Dia]=PIVwindowed(im1,im2,corr,window,res,zpad,D,Zeromean,Peaklocator,Peakswitch,X,Y,Uin,Vin)
% --- DPIV Correlation ---

%convert input parameters
im1=double(im1);
im2=double(im2);
L=size(im1);

%convert to gridpoint list
X=X(:);
Y=Y(:);

%correlation and window mask types
ctype    = {'SCC','RPC'};
tcorr = char(ctype(corr+1)); 

%preallocate velocity fields and grid format
Nx = window(1);
Ny = window(2);
if nargin <=13
    Uin = zeros(length(X),1);
    Vin = zeros(length(X),1);
end

if Peakswitch
    Uin=repmat(Uin(:,1),[1 3]);
    Vin=repmat(Vin(:,1),[1 3]);
    U = zeros(length(X),3);
    V = zeros(length(X),3);
    C = zeros(length(X),3);
    Dia = zeros(length(X),3);
else
    U = zeros(length(X),1);
    V = zeros(length(X),1);
    C = [];
    Dia = [];
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
sfilt1 = windowmask([Nx Ny],[res(1, 1) res(1, 2)]);
sfilt2 = windowmask([Nx Ny],[res(2, 1) res(2, 2)]);

%correlation plane normalization function (always off)
cnorm = ones(Ny,Nx);

%RPC spectral energy filter
spectral = fftshift(energyfilt(Sx,Sy,D,0));

%fftshift indicies
fftindy = [Sy/2+1:Sy 1:Sy/2];
fftindx = [Sx/2+1:Sx 1:Sx/2];

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
            zone1 = im1( max([1 ymin1]):min([L(1) ymax1]),max([1 xmin1]):min([L(2) xmax1]));
            zone2 = im2( max([1 ymin2]):min([L(1) ymax2]),max([1 xmin2]):min([L(2) xmax2]));
            
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
                zone1=zone1-mean(mean(zone1));
                zone2=zone2-mean(mean(zone2));
            end

            %apply the image spatial filter
            region1 = (zone1).*sfilt1;
            region2 = (zone2).*sfilt2;

            %FFTs and Cross-Correlation
            f1   = fftn(region1,[Sy Sx]);
            f2   = fftn(region2,[Sy Sx]);
            P21  = f2.*conj(f1);

            %Standard Fourier Based Cross-Correlation
            G = ifftn(P21,'symmetric');
            G = G(fftindy,fftindx);
            G = abs(G);
            
            %subpixel estimation
            [U(n,:),V(n,:),Ctemp,Dtemp]=subpixel(G,Nx,Ny,cnorm,Peaklocator,Peakswitch);
%             winmean=mean(mean(region1))*mean(mean(region2));
%             [U(n,:),V(n,:),Ctemp]=subpixel(G,Nx,Ny,cnorm,Peaklocator,Peakswitch,winmean);
            if Peakswitch
                C(n,:)=Ctemp;
                Dia(n,:)=Dtemp;
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

            %find the image windows
            zone1 = im1( max([1 ymin1]):min([L(1) ymax1]),max([1 xmin1]):min([L(2) xmax1]));
            zone2 = im2( max([1 ymin2]):min([L(1) ymax2]),max([1 xmin2]):min([L(2) xmax2]));
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
                zone1=zone1-mean(mean(zone1));
                zone2=zone2-mean(mean(zone2));
            end

            %apply the image spatial filter
            region1 = (zone1).*sfilt1;
            region2 = (zone2).*sfilt2;

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
            [U(n,:),V(n,:),Ctemp,Dtemp]=subpixel(G,Nx,Ny,cnorm,Peaklocator,Peakswitch);
%             winmean=mean(mean(region1))*mean(mean(region2));
%             [U(n,:),V(n,:),Ctemp]=subpixel(G,Nx,Ny,cnorm,Peaklocator,Peakswitch,winmean);
            if Peakswitch
                C(n,:)=Ctemp;
                Dia(n,:)=Dtemp;
            end
        end
end

%add DWO to estimation
U = round(Uin)+U;
V = round(Vin)+V;

end