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
    gauss_int = zeros(size(mapint_i));
    xp = zeros(size(mapint_i));
    yp = zeros(size(mapint_i));
    for ii = 1:length(mapint_i)
        xp(ii) = locxy_i(ii,2);
        yp(ii) = locxy_i(ii,1);
    end
    
    % map an intensity profile of a gaussian function:
    for rr = 1:size(xp,1)
%         gauss_int(rr)=I0*exp(-abs(betas)*(((xp(rr))-x_centroid)^2 + ...
%             ((yp(rr))-y_centroid)^2));
        gauss_int(rr)=I0*exp(-abs(betasx).*(cos(alpha).*(xp(rr)-x_centroid) - sin(alpha).*(yp(rr)-y_centroid)).^2 - ...
            abs(betasy).*(sin(alpha).*(xp(rr)-x_centroid) + cos(alpha).*(yp(rr)-y_centroid)).^2);
    end

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
F=mapint_i-gauss_int;

end