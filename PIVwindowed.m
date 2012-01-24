function [X,Y,U,V,C,Dia]=PIVwindowed(im1,im2,corr,window,res,zpad,D,Zeromean,Peaklocator,Peakswitch,fracval,X,Y,Uin,Vin)
% --- DPIV Correlation ---

%convert input parameters
im1=double(im1);
im2=double(im2);
L=size(im1);

%convert to gridpoint list
X=X(:);
Y=Y(:);

%correlation and window mask types
%ctype    = {'SCC','RPC','GCC','FWC'};
ctype    = {'SCC','RPC','GCC','FWC','SPC','qRPC'};
tcorr = char(ctype(corr+1)); %and because corr = A.Corr-1, we have to fix it again here - why bother?

%preallocate velocity fields and grid format
Nx = window(1);
Ny = window(2);
if nargin <=14
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

if strcmpi(tcorr,'FWC')
    frac = fracval;
    spectral = ones(size(spectral));
elseif strcmpi(tcorr,'GCC')
    frac = 1;
    spectral = ones(size(spectral));
else
    frac = 1;
end

switch upper(tcorr)
    
    %Standard Cross Correlation
    case 'SCC'
        
        if size(im1,3) == 3
            Gens=zeros(Ny,Nx,3);
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
                
                for r=1:size(im1,3);
                    %find the image windows
                    zone1 = im1( max([1 ymin1]):min([L(1) ymax1]),max([1 xmin1]):min([L(2) xmax1]),r);
                    zone2 = im2( max([1 ymin2]):min([L(1) ymax2]),max([1 xmin2]):min([L(2) xmax2]),r);
                    
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
                    Gens(:,:,r) = abs(G);
                end
                G = mean(Gens,3);
                
                %subpixel estimation
                [U(n,:),V(n,:),Ctemp,Dtemp]=subpixel(G,Nx,Ny,cnorm,Peaklocator,Peakswitch);
                if Peakswitch
                    C(n,:)=Ctemp;
                    Dia(n,:)=Dtemp;
                end
                
            end
            
        else
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
                if Peakswitch
                    C(n,:)=Ctemp;
                    Dia(n,:)=Dtemp;
                end
            end
        end

    %Robust Phase Correlation
    case {'RPC','GCC','FWC'}
        
        if size(im1,3) == 3
            Gens=zeros(Ny,Nx,3);
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
                
                for r=1:size(im1,3);
                    %find the image windows
                    zone1 = im1( max([1 ymin1]):min([L(1) ymax1]),max([1 xmin1]):min([L(2) xmax1]),r);
                    zone2 = im2( max([1 ymin2]):min([L(1) ymax2]),max([1 xmin2]):min([L(2) xmax2]),r);
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
                    if frac ~= 1
                        R = P21./(W.^frac);
                    else
                        R = P21./W;
                    end
                    
                    %Robust Phase Correlation with spectral energy filter
                    G = ifftn(R.*spectral,'symmetric');
                    G = G(fftindy,fftindx);
                    Gens(:,:,r) = abs(G);
                end
                G=mean(Gens,3);
                
                %subpixel estimation
                [U(n,:),V(n,:),Ctemp,Dtemp]=subpixel(G,Nx,Ny,cnorm,Peaklocator,Peakswitch);
                if Peakswitch
                    C(n,:)=Ctemp;
                    Dia(n,:)=Dtemp;
                end
            end
            
        else
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
                if frac ~= 1
                    R = P21./(W.^frac);
                else
                    R = P21./W;
                end
                
                %Robust Phase Correlation with spectral energy filter
                G = ifftn(R.*spectral,'symmetric');
                G = G(fftindy,fftindx);
                G = abs(G);
                
                %subpixel estimation
                [U(n,:),V(n,:),Ctemp,Dtemp]=subpixel(G,Nx,Ny,cnorm,Peaklocator,Peakswitch);
                if Peakswitch
                    C(n,:)=Ctemp;
                    Dia(n,:)=Dtemp;
                end
            end
        end
    %Quaternion Robust Phase Correlation
    case {'QRPC'}
        
        %quaternion processing needs a 3 color matrix.
        %we could also throw an error, but this seems nicer
        if size(im1,3) ~=3
            error('not a color image')
            im1 = repmat(im1,[1,1,3]);
            im2 = repmat(im2,[1,1,3]);
        end
        
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
            
            
            zone1 = zeros(Ny,Nx,3);
            zone2 = zeros(Ny,Nx,3);
            
            
            %for r=1:3;
            %find the image windows
            w1 = im1( max([1 ymin1]):min([L(1) ymax1]),max([1 xmin1]):min([L(2) xmax1]),:);
            w2 = im2( max([1 ymin2]):min([L(1) ymax2]),max([1 xmin2]):min([L(2) xmax2]),:);
            
            %check to make sure enough data is present
            if size(w1,1)~=Ny || size(w1,2)~=Nx
                %w1 = zeros(Ny,Nx,size(im1,3));
                zone1( 1+max([0 1-ymin1]):Ny-max([0 ymax1-L(1)]),1+max([0 1-xmin1]):Nx-max([0 xmax1-L(2)]) ,:) = w1;
            else
                zone1 = w1;
            end
            if size(w2,1)~=Ny || size(w2,2)~=Nx
                %w2 = zeros(Ny,Nx);
                zone2( 1+max([0 1-ymin2]):Ny-max([0 ymax2-L(1)]),1+max([0 1-xmin2]):Nx-max([0 xmax2-L(2)]) ,:) = w2;
            else
                zone2 = w2;
            end
            
            %subtract mean from each region
            if Zeromean==1
                meanC1 = mean(mean(zone1,1),2);
                meanC2 = mean(mean(zone2,1),2);
                for r=1:3
                    zone1(:,:,r)=zone1(:,:,r)-meanC1(1,1,r);
                    zone2(:,:,r)=zone2(:,:,r)-meanC2(1,1,r);
                end
            end
                
            %apply the image spatial filter
            for r=1:3
                zone1(:,:,r) = (zone1(:,:,r)).*sfilt1;
                zone2(:,:,r) = (zone2(:,:,r)).*sfilt2;
            end
            
            
            %convert zones to quaternions
            region1 = convert(quaternion(zone1(:,:,1), zone1(:,:,2), zone1(:,:,3)), 'double');  %f or F below
            region2 = convert(quaternion(zone2(:,:,1), zone2(:,:,2), zone2(:,:,3)), 'double');  %g or G below
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %The following code is adpated from
            % $Id: vector_phase_correlation.m,v 1.4 2011/01/07 21:31:39 sangwine Exp $
            % Example script to calculate the vector phase correlation described in:
            %
            % Sangwine, S. J., Ell, T. A. and Moxey C. E., Vector Phase Correlation,
            % Electronics Letters, 37(25), 6 December 2001, 1513-5.
            % http://dx.doi.org/10.1049/el:20011035
            % GPL v3 or later
            % Copyright � : Steve Sangwine, May 2007
            
            
            
            mu = unit(quaternion(1,1,1)); % Transform axis, a unit pure quaternion.
            
            % Calculate the quantities needed to evaluate equation 4 in the paper. The
            % Matlab implementations of the FFTs differ in scale factor from the FFTs
            % in equations 2 and 3 of the paper (because they follow the Matlab
            % convention of applying the scale factor only to the inverse transform,
            % rather than splitting it equally between the forward and inverse) so we
            % have to make an adjustment for this:
            
            MN  = numel(region1); % This gives the number of pixels in each image, which
            % is the product MN in equations 2 and 3.
            
            %JJC: investigate optimizing qfft2 and iqfft2 using fftn
            FL  =  qfft2(region1, mu, 'L') ./ sqrt(MN); % This is F subscript L in the paper.
            FIL = iqfft2(region1, mu, 'L') .* sqrt(MN); % This is F subscript -L.
            GR  =  qfft2(region2, mu, 'R') ./ sqrt(MN); % This is G subscript R.
            
            % Now decompose GR into components parallel and perpendicular to the
            % transform axis, mu.
            
            GRpar  = (GR - mu .* GR .* mu) ./ 2; % Parallel/perpendicular
            GRperp =  GR - GRpar;                % decomposition.
            
            P21 = conj(FL) .* GRpar + conj(FIL) .* GRperp;
            
            % Now we have the hypercomplex cross-power spectrum P21. We compute the
            % hypercomplex phase correlation directly according to equation 6 in the
            % paper. Note that the unit function from the QTFM toolbox divides each
            % element of RR by its modulus. We don't have to write this out explicitly.
            
            %end adapted quaternion code
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
%                 %FFTs and Cross-Correlation
%                 f1   = fftn(region1,[Sy Sx]);
%                 f2   = fftn(region2,[Sy Sx]);
%                 P21  = f2.*conj(f1);

%                 %Phase Correlation
%                 W = ones(Sy,Sx);
%                 Wden = sqrt(P21.*conj(P21));
%                 W(P21~=0) = Wden(P21~=0);
%                 if frac ~= 1
%                     R = P21./(W.^frac);
%                 else
%                     R = P21./W;
%                 end
            
            %R = unit(P21);
            
            %Robust Phase Correlation with spectral energy filter
            %G = ifftn(R.*spectral,'symmetric');
            G = iqfft2(unit(P21).*spectral, mu, 'R') .* sqrt(MN); 
            
            G = abs(G); %convert from quaternions to magnitudes
            G = G(fftindy,fftindx); %should be 2-D now
            
            %subpixel estimation
            [U(n,:),V(n,:),Ctemp,Dtemp]=subpixel(G,Nx,Ny,cnorm,Peaklocator,Peakswitch);
            if Peakswitch
                C(n,:)=Ctemp;
                Dia(n,:)=Dtemp;
            end
        end
            
    otherwise
        %throw an error, we shouldn't be here (probably SPC?)
        keyboard
        error('invalid correlation type')
        
end

%add DWO to estimation
U = round(Uin)+U;
V = round(Vin)+V;
end