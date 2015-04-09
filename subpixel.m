function [u,v,M,D,major_axis_length, minor_axis_length, PEAK_ANGLE, PEAK_ECCENTRICITY] = subpixel(...
    SPATIAL_CORRELATION_PLANE,...
    CORRELATION_WIDTH, CORRELATION_HEIGHT, WEIGHTING_MATRIX, ...
    PEAK_FIT_METHOD, FIND_MULTIPLE_PEAKS, PARTICLE_DIAMETER_2D)

%intialize indices
cc_x = -floor(CORRELATION_WIDTH/2):ceil(CORRELATION_WIDTH/2)-1;
cc_y = -floor(CORRELATION_HEIGHT/2):ceil(CORRELATION_HEIGHT/2)-1;

% Set values for peak eccentricity and angle
% so that the function returns them properly
% even if the Guassian least-squares fit method isn't used.
PEAK_ECCENTRICITY = 0;
PEAK_ANGLE = 0;

%find maximum correlation value
[M,I] = max(SPATIAL_CORRELATION_PLANE(:));

% Use 4 standard deviations for the peak sizing (e^-2)
sigma = 4;

%if correlation empty
if M==0
    if FIND_MULTIPLE_PEAKS
        u=zeros(1,3);
        v=zeros(1,3);
        M=zeros(1,3);
        D=zeros(1,3);
        DX=zeros(1,3);
        DY=zeros(1,3);
        PEAK_ANGLE = zeros(1,3);
        
    else
        u=0; v=0; M=0; D=0; DX=0; DY=0; PEAK_ANGLE=0;
    end
else
    if FIND_MULTIPLE_PEAKS
        u=zeros(1,3);
        v=zeros(1,3);
        D=zeros(1,3);
        DX=zeros(1,3);
        DY=zeros(1,3);
        PEAK_ANGLE=zeros(1,3);

        %Locate peaks using imregionalmax
        A=imregionalmax(SPATIAL_CORRELATION_PLANE);
        peakmat=SPATIAL_CORRELATION_PLANE.*A;
        for i=2:3
            peakmat(peakmat==M(i-1))=0;
            [M(i),I(i)]=max(peakmat(:));
        end
        j=length(M);
    else
        u=zeros(1,1);
        v=zeros(1,1);
        D=zeros(1,1);
        DX=0; DY=0; PEAK_ANGLE=0;
        j=1;    
    end
    
    for i=1:j
        method=PEAK_FIT_METHOD;
        
        %find x and y indices
        shift_locy = 1+mod(I(i)-1,CORRELATION_HEIGHT);
        shift_locx = ceil(I(i)/CORRELATION_HEIGHT);

        shift_errx=[];
        shift_erry=[];
        %find subpixel displacement in x
        if CORRELATION_WIDTH == 1
            shift_errx = 1; method=1;
        elseif shift_locx == 1
            %boundary condition 1
            shift_errx =  SPATIAL_CORRELATION_PLANE( shift_locy , shift_locx+1 )/M(i); method=1;
        elseif shift_locx == CORRELATION_WIDTH
            %boundary condition 2
            shift_errx = -SPATIAL_CORRELATION_PLANE( shift_locy , shift_locx-1 )/M(i); method=1;
        elseif SPATIAL_CORRELATION_PLANE( shift_locy , shift_locx+1 ) == 0
            %endpoint discontinuity 1
            shift_errx = -SPATIAL_CORRELATION_PLANE( shift_locy , shift_locx-1 )/M(i); method=1;
        elseif SPATIAL_CORRELATION_PLANE( shift_locy , shift_locx-1 ) == 0
            %endpoint discontinuity 2
            shift_errx =  SPATIAL_CORRELATION_PLANE( shift_locy , shift_locx+1 )/M(i); method=1;
        end
        if CORRELATION_HEIGHT == 1
            shift_erry = 1; method=1;
        elseif shift_locy == 1
            %boundary condition 1
            shift_erry = -SPATIAL_CORRELATION_PLANE( shift_locy+1 , shift_locx )/M(i); method=1;
        elseif shift_locy == CORRELATION_HEIGHT
            %boundary condition 2
            shift_erry =  SPATIAL_CORRELATION_PLANE( shift_locy-1 , shift_locx )/M(i); method=1;
        elseif SPATIAL_CORRELATION_PLANE( shift_locy+1 , shift_locx ) == 0
            %endpoint discontinuity 1
            shift_erry =  SPATIAL_CORRELATION_PLANE( shift_locy-1 , shift_locx )/M(i); method=1;
        elseif SPATIAL_CORRELATION_PLANE( shift_locy-1 , shift_locx ) == 0
            %endpoint discontinuity 2
            shift_erry = -SPATIAL_CORRELATION_PLANE( shift_locy+1 , shift_locx )/M(i); method=1;
        end

        if method==2
            
            %%%%%%%%%%%%%%%%%%%%
            % 4-Point Gaussian %
            %%%%%%%%%%%%%%%%%%%%
            
            %Since the case where M is located at a border will default to
            %the 3-point gaussian and we don't have to deal with
            %saturation, just use 4 points in a tetris block formation:
            %
            %             *
            %            ***
            
            points=[shift_locy   shift_locx   SPATIAL_CORRELATION_PLANE(shift_locy  ,shift_locx  );...
                    shift_locy-1 shift_locx   SPATIAL_CORRELATION_PLANE(shift_locy-1,shift_locx  );...
                    shift_locy   shift_locx-1 SPATIAL_CORRELATION_PLANE(shift_locy  ,shift_locx-1);...
                    shift_locy   shift_locx+1 SPATIAL_CORRELATION_PLANE(shift_locy  ,shift_locx+1)];
                
            [~,IsortI] = sort(points(:,3),'descend');
            points = points(IsortI,:);

            x1=points(1,2); x2=points(2,2); x3=points(3,2); x4=points(4,2);
            y1=points(1,1); y2=points(2,1); y3=points(3,1); y4=points(4,1);
            a1=points(1,3); a2=points(2,3); a3=points(3,3); a4=points(4,3);

            peak_angle(1) = (x4^2)*(y2 - y3) + (x3^2)*(y4 - y2) + ((x2^2) + (y2 - y3)*(y2 - y4))*(y3 - y4);
            peak_angle(2) = (x4^2)*(y3 - y1) + (x3^2)*(y1 - y4) - ((x1^2) + (y1 - y3)*(y1 - y4))*(y3 - y4);
            peak_angle(3) = (x4^2)*(y1 - y2) + (x2^2)*(y4 - y1) + ((x1^2) + (y1 - y2)*(y1 - y4))*(y2 - y4);
            peak_angle(4) = (x3^2)*(y2 - y1) + (x2^2)*(y1 - y3) - ((x1^2) + (y1 - y2)*(y1 - y3))*(y2 - y3);

            gamma(1) = (-x3^2)*x4 + (x2^2)*(x4 - x3) + x4*((y2^2) - (y3^2)) + x3*((x4^2) - (y2^2) + (y4^2)) + x2*(( x3^2) - (x4^2) + (y3^2) - (y4^2));
            gamma(2) = ( x3^2)*x4 + (x1^2)*(x3 - x4) + x4*((y3^2) - (y1^2)) - x3*((x4^2) - (y1^2) + (y4^2)) + x1*((-x3^2) + (x4^2) - (y3^2) + (y4^2));
            gamma(3) = (-x2^2)*x4 + (x1^2)*(x4 - x2) + x4*((y1^2) - (y2^2)) + x2*((x4^2) - (y1^2) + (y4^2)) + x1*(( x2^2) - (x4^2) + (y2^2) - (y4^2));
            gamma(4) = ( x2^2)*x3 + (x1^2)*(x2 - x3) + x3*((y2^2) - (y1^2)) - x2*((x3^2) - (y1^2) + (y3^2)) + x1*((-x2^2) + (x3^2) - (y2^2) + (y3^2));

            delta(1) = x4*(y2 - y3) + x2*(y3 - y4) + x3*(y4 - y2);
            delta(2) = x4*(y3 - y1) + x3*(y1 - y4) + x1*(y4 - y3);
            delta(3) = x4*(y1 - y2) + x1*(y2 - y4) + x2*(y4 - y1);
            delta(4) = x3*(y2 - y1) + x2*(y1 - y3) + x1*(y3 - y2);

            deno = 2*(log(a1)*delta(1) + log(a2)*delta(2) + log(a3)*delta(3) + log(a4)*delta(4));

            x_centroid = (log(a1)*peak_angle(1) + log(a2)*peak_angle(2) + log(a3)*peak_angle(3) + log(a4)*peak_angle(4))/deno;
            y_centroid = (log(a1)*gamma(1) + log(a2)*gamma(2) + log(a3)*gamma(3) + log(a4)*gamma(4))/deno;
            shift_errx=x_centroid-shift_locx;
            shift_erry=y_centroid-shift_locy;
            
            betas = abs((log(a2)-log(a1))/((x2-x_centroid)^2+(y2-y_centroid)^2-(x1-x_centroid)^2-(y1-y_centroid)^2));
            D(i)=sqrt(sigma^2/(2*betas));
            
        elseif any(method==[3 4])
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%
            % Gaussian Least Squares %
            %%%%%%%%%%%%%%%%%%%%%%%%%%
            %convert the particle diameter to diameter of equivalent correlation peak
            D1 = sqrt(2) .* PARTICLE_DIAMETER_2D(1);
            D2 = sqrt(2) .* PARTICLE_DIAMETER_2D(2);
            goodSize = 0;  %gets set =1 after fit, but reset to 0 if betaX or betaY are bigger than 2*D1 or 2*D2
            
            %keep trying while method not 1 (G.3pt.fit), and the search diameter (2x expected diam.) is less than half the window size
            while ~goodSize && method~=1
                %Find a suitable window around the peak (+/- D1,D2)
                x_min=shift_locx-ceil(D1); x_max=shift_locx+ceil(D1);
                y_min=shift_locy-ceil(D2); y_max=shift_locy+ceil(D2);
                if x_min<1
                    x_min=1;
                end
                if x_max>CORRELATION_WIDTH
                    x_max=CORRELATION_WIDTH;
                end
                if y_min<1
                    y_min=1;
                end
                if y_max>CORRELATION_HEIGHT
                    y_max=CORRELATION_HEIGHT;
                end
                
                points = ...
                    double(SPATIAL_CORRELATION_PLANE(y_min:y_max,x_min:x_max).* ...
                    WEIGHTING_MATRIX(y_min:y_max,x_min:x_max));
                
                % Subtract the minimum value from the points matrix
                points_min_sub = points - min(points(:));
                
                % Normalize the points matrix so the max value is one
                points_norm = points_min_sub ./ max(points_min_sub(:));

                %Options for the lsqnonlin solver using Levenberg-Marquardt solver
                options=optimset('MaxIter',1200,'MaxFunEvals',5000,'TolX',1e-6,'TolFun',1e-6,...
                    'Display','off','DiffMinChange',1e-7,'DiffMaxChange',1,...
                    'Algorithm','levenberg-marquardt');
                
                %xvars is [M,betaX,betaY,CX,CY,alpha]
                
                % Set empty lower bounds (LB) and upper bounds (UB) 
                % for the least squares solver LSQNONLIN.
                LB = [];
                UB = [];

                %Initial values for the solver (have to convert D into Beta)
                %x0 must be double for lsqnonlin, so convert M
%                 x0 = [double(M(i)), ...
%                       0.5*(sigma/D1)^2, ...
%                       0.5*(sigma/D2)^2, ...
%                       shift_locx, ...
%                       shift_locy, ...
%                       0];
%                   
                  
              x0 = [1, ...
                      0.5*(sigma/D1)^2, ...
                      0.5*(sigma/D2)^2, ...
                      shift_locx, ...
                      shift_locy, ...
                      0];

                [xloc, yloc]=meshgrid(x_min:x_max,y_min:y_max);

                %Run solver; default to 3-point gauss if it fails
                try
                    %[xvars resnorm resid exitflag output]=lsqnonlin(@leastsquares2D,x0,LB,UB,options,points(:),[yloc(:),xloc(:)],method);
%                     [xvars]=lsqnonlin(@leastsquares2D,x0, LB, UB,options,points(:),[yloc(:),xloc(:)], method);
                    [xvars]=lsqnonlin(@leastsquares2D,x0, LB, UB,options,points_norm(:),[yloc(:),xloc(:)], method);
                    shift_errx=xvars(4)-shift_locx;
                    shift_erry=xvars(5)-shift_locy;
                    
                    %convert beta to diameter, diameter = 4*std.dev.
                    dA = sigma/sqrt(2*abs(xvars(2)));    %diameter of axis 1
                    dB = sigma/sqrt(2*abs(xvars(3)));    %diameter of axis 2
                    
                    % Find the equivalent diameter for a circle with 
                    % equal area and return that value
                    D(i) = sqrt(dA*dB);
                    peak_angle = mod(xvars(6),2*pi);
                                        
                    % Calculate the lengths of the major and minor axes 
                    % of the best-fit elliptical Gaussian 
                    major_axis_length = max(dA, dB);
                    minor_axis_length = min(dA, dB);
                    
                    % Calculate the eccentricity of the
                    % elliptical Gaussian peak.
                    PEAK_ECCENTRICITY = sqrt(1 - ...
                        minor_axis_length^2 / major_axis_length^2);
                    
                    % These are the lengths of the projections of the 
                    % elliptical Gaussian peak onto the horizontal and vertical
                    % axes.
                    dX = max( abs(dA*cos(peak_angle)), abs(dB*sin(peak_angle)) );
                    dY = max( abs(dA*sin(peak_angle)), abs(dB*cos(peak_angle)) );
                                        
                    DX(i) = dX;
                    DY(i) = dY;
                    PEAK_ANGLE(i) = peak_angle;
                    
                    %LSqF didn't fail...
                    goodSize = 1;
                    
 
                    %if D1 or D2 are already too big, just quit - it's
                    %the best we're going to do.
                    %Have to check in loop, if check at while statement,
                    %might never size it at all.
                    if 2*D1<CORRELATION_WIDTH/2 && 2*D2<CORRELATION_HEIGHT/2
                        goodSize = 1;
                    end

                catch err%#ok
                    %warning(err.message)
                    disp(err.message)
                   
                    method=1;
                end
            end %while trying to fit region
        end %if method==2,3,4
        if method==1

            %%%%%%%%%%%%%%%%%%%%
            % 3-Point Gaussian %
            %%%%%%%%%%%%%%%%%%%%
            
            if isempty(shift_errx)
                
                % Gaussian fit
                lCm1 = log(SPATIAL_CORRELATION_PLANE( shift_locy , shift_locx-1 )*WEIGHTING_MATRIX( shift_locy , shift_locx-1 ));
                lC00 = log(SPATIAL_CORRELATION_PLANE( shift_locy , shift_locx   )*WEIGHTING_MATRIX( shift_locy , shift_locx   ));
                lCp1 = log(SPATIAL_CORRELATION_PLANE( shift_locy , shift_locx+1 )*WEIGHTING_MATRIX( shift_locy , shift_locx+1 ));
                if (2*(lCm1+lCp1-2*lC00)) == 0
                    shift_errx = 0;
                    dX = nan;
                else
                    shift_errx = (lCm1-lCp1)/(2*(lCm1+lCp1-2*lC00));
                    betax = abs(lCm1-lC00)/((-1-shift_errx)^2-(shift_errx)^2);
                    dX = sigma./sqrt((2*betax));
                end
            else
                dX = nan;
            end
            
            if isempty(shift_erry)
                lCm1 = log(SPATIAL_CORRELATION_PLANE( shift_locy-1 , shift_locx )*WEIGHTING_MATRIX( shift_locy-1 , shift_locx ));
                lC00 = log(SPATIAL_CORRELATION_PLANE( shift_locy   , shift_locx )*WEIGHTING_MATRIX( shift_locy   , shift_locx ));
                lCp1 = log(SPATIAL_CORRELATION_PLANE( shift_locy+1 , shift_locx )*WEIGHTING_MATRIX( shift_locy+1 , shift_locx ));
                if (2*(lCm1+lCp1-2*lC00)) == 0
                    shift_erry = 0;
                    dY = nan;
                else
                    shift_erry = (lCm1-lCp1)/(2*(lCm1+lCp1-2*lC00));
                    betay = abs(lCm1-lC00)/((-1-shift_erry)^2-(shift_erry)^2);
                    dY = sigma./sqrt((2*betay));
                end
            else
                dY = nan;
            end
            
            D(i) = nanmean([dX dY]);
            
                    
            DX(i) = dX;
            DY(i) = dY;

            
        end
        u(i)=cc_x(shift_locx)+shift_errx;
        v(i)=cc_y(shift_locy)+shift_erry;
        
        if isinf(u(i)) || isinf(v(i))
            u(i)=0; v(i)=0;
        end
    end
end

end

function F = leastsquares2D(x,mapint_i,locxy_i,method)
% function F = leastsquares2D(x,mapint_i,locxy_i,method)
% This function is called by lsqnonlin if the least squares or continuous
% least squares method has been chosen. It solve (leastsqures) for a 
% Gaussian surface that best fist a list of sample points.  The code has
% been updated to now handle eliptical Gaussian shapes using a 
% trigonometric formulation for an arbitrary eliptical Gaussian function
% taken from ( Scharnowski (2012) Exp Fluids).  
% 
% Inputs: 
%  x:        Is a vectore containing an intial guess at the parameter values 
%            for the gaussian fit.  [Max Value, Beta in the X direction, 
%            Beta in the Y direction, Estimated Centroid for X, Estimated
%            Centroid for Y, Estimated Orientation Angle]
%  mapint_i: List of intensity values.
%  locxy_i:  Location of intensity samples for X and Y
%  method:   This switches between Standard Gaussian (3) and Continous
%            Gaussian (4).
%    
% Outputs:
% F:         Is the variable being minimized - the difference between the 
%            gaussian curve and the actual intensity values.
%
% Adapted from M. Brady's 'leastsquaresgaussfit' and 'mapintensity'
% Edited:
% B.Drew - 7.18.2008
% S. Raben - 7.24.2012


I0=x(1);
betasx=x(2);
betasy=x(3);
x_centroid=x(4);
y_centroid=x(5);
alpha = x(6);

if method==3

    
    xp = locxy_i(:,2);
    yp = locxy_i(:,1);

    
    % map an intensity profile of a gaussian function    
    gauss_int = I0   * exp(-abs(betasx).*(cos(alpha).*(xp-x_centroid) - ...
        sin(alpha)  .* (yp-y_centroid)).^2 - ...
        abs(betasy) .* (sin(alpha).*(xp-x_centroid) + ...
        cos(alpha)  .* (yp-y_centroid)).^2);
    
elseif method==4
    
    %Just like in the continuous four-point method, lsqnonlin tries negative
    %values for x(2) and x(3), which will return errors unless the abs() function is
    %used in front of all the x(2)'s.
    num1=(I0*pi)/4;
    num2=sqrt(abs(mean([betasx betasy])));
    
    S = size(mapint_i);
    gauss_int = zeros(S(1),S(2));
    xp = zeros(size(mapint_i));
    yp = zeros(size(mapint_i));
    for ii = 1:length(mapint_i)
        xp(ii) = locxy_i(ii,1)-0.5;
        yp(ii) = locxy_i(ii,2)-0.5;
        erfx1 = erf(num2*(xp(ii)-x_centroid));
        erfy1 = erf(num2*(yp(ii)-y_centroid));
        erfx2 = erf(num2*(xp(ii)+1-x_centroid));
        erfy2 = erf(num2*(yp(ii)+1-y_centroid));
        
        % map an intensity profile of a gaussian function:
        gauss_int(ii)=(num1/abs(betas))*(erfx1*(erfy1-erfy2)+erfx2*(-erfy1+erfy2));
        
    end
end
% compare the Gaussian curve to the actual pixel intensities
F = mapint_i-gauss_int;


end