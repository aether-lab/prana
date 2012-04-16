function [x_c,y_c,D,P,E,Meth] = Leastsqrmethods(I,method,sigma,options,locxy,max_locxy)

% Find the maxium value of I and all of its locations
[max_int,int_L] = max(I(:));%#ok
max_row = locxy(:,1);
max_col = locxy(:,2);
E = 0;
Meth = method;
%Try making a 5x5 window around the maximum intensity
%pixel to reduce the chance of grouping particles together.
%If a 5x5 window isn't possible due to the particle's location on
%the image, then use the original window.

% Find extents of I
[ymax_I xmax_I] = size(I);

%This method cannot fit a Gaussian shape to a particle intensity profile
%where the highest intensity pixel is on the edge of I (assumes a
%Gaussian light scattering profile)
if min(max_row)==1
    stop_y = 1;
    edge_cky = 1;
elseif max(max_row)==ymax_I
    stop_y = 1;
    edge_cky = 1;
else
    stop_y = 0;
    edge_cky = 0;
end

if min(max_col)==1
    stop_x = 1;
    edge_ckx = 1;
elseif max(max_col)==xmax_I
    stop_x = 1;
    edge_ckx = 1;
else
    stop_x = 0;
    edge_ckx = 0;
end

%2 Logical loops to deal with saturated or equvalent pixels values surrounding
%max_int.  Incrementally moves outward until suitable values are found or
%the extents of I are reached
i=1;
j=1;
extent_left=min(max_col);
extent_right=max(max_col);
extent_top=min(max_row);
extent_bottom=max(max_row);
max_edge_left = [0 0];
max_edge_top = [0 0];
diameter_x = 0;
diameter_y = 0;
left_int = max_int;
left_locxy = max_locxy;
right_int = max_int;
right_locxy = max_locxy;
top_int = max_int;
top_locxy = max_locxy;
bottom_int = max_int;
bottom_locxy = max_locxy;

%search in the x-dimension
while stop_x==0 && extent_left>=1 && extent_right<=xmax_I
    %locate x-index and intensity values for the nearest neighbors of max_int
    left_locxy  = [max_locxy(1) max_locxy(2)-i];
    left_int    = I(left_locxy(1),left_locxy(2));
    right_locxy = [max_locxy(1) max_locxy(2)+i];
    right_int   = I(right_locxy(1),right_locxy(2));
    %logical check to assure that the three pixels chosen are not equal
    if left_int==max_int || right_int==max_int
        i=i+1;
        max_edge_left = left_locxy;
    else
        stop_x=1;
    end
    
    %make sure there are data points to the left & right
    extent_left=min(max_col)-i;
    extent_right=max(max_col)+i;
end

%search in the y-dimension
while stop_y==0 && extent_top>=1 && extent_bottom<=ymax_I
    %locate y-index and intensity values for the nearest neighbors of max_int
    top_locxy    = [max_locxy(1)-j max_locxy(2)];
    top_int      = I(top_locxy(1),top_locxy(2));
    bottom_locxy = [max_locxy(1)+j max_locxy(2)];
    bottom_int   = I(bottom_locxy(1),bottom_locxy(2));
    %logical check to assure that the three pixels chosen are not equal
    if top_int==max_int || bottom_int==max_int
        j=j+1;
        max_edge_top = top_locxy;
    else
        stop_y=1;
    end
    
    %make sure there are data points to the top & bottom
    extent_top=min(max_row)-j;
    extent_bottom=max(max_row)+j;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3pt Gauss Fit - Start
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Given that this method assumes a Gaussian shaped intensity profile, the
%function is designed to return NaN should any of the nearest neighbor
%pixels to max_int equal zero or if the particles have equal intensity
%values
if edge_ckx == 0;
    LR_denom = (2*(log(left_int)+log(right_int)-2*log(max_int)));
    if left_int==0 || right_int==0 || LR_denom==0
        x_centroid = max_locxy(2);
        %keyboard
        diameter_x=NaN;
        I0x=NaN;
    else
        %calculate the x&y centroid location using the 3-pt Gaussian Method
        x_centroid=max_locxy(2)+log(left_int/right_int)/LR_denom;
    end
else
    x_centroid = max_locxy(2);
    diameter_x=NaN;
    I0x=NaN;
end

if edge_cky == 0
    TB_denom = (2*(log(top_int)+log(bottom_int)-2*log(max_int)));
    if top_int==0 || bottom_int==0 || TB_denom==0
        y_centroid = max_locxy(1);
        %keyboard
        diameter_y=NaN;
        I0y=NaN;
    else
        %calculate the x&y centroid location using the 3-pt Gaussian Method
        y_centroid=max_locxy(1)+log(top_int/bottom_int)/TB_denom;
    end
else
    y_centroid = max_locxy(1);
    diameter_y=NaN;
    I0y=NaN;
end

%The diameter of the particle is directly related to the variance (where
%Beta=1/(variance)^2 and diameter=sigma*variance.  Logic statements deals with
%the possibility that one pixel adjacent to max_int is equal to max_int
%while the other adjacent pixel does not
if ~isnan(diameter_x)
    if nnz(max_edge_left) < 1
        betax=log(left_int/max_int)/((max_locxy(2)-x_centroid)^2-(left_locxy(2)-x_centroid)^2);
        I0x=left_int/exp(-betax*(left_locxy(2)-x_centroid)^2);
    else
        betax=log(left_int/max_int)/((max_edge_left(2)-x_centroid)^2-(left_locxy(2)-x_centroid)^2);
        I0x=left_int/exp(-betax*(left_locxy(2)-x_centroid)^2);
    end
    diameter_x=sigma./sqrt((2*betax));
end

if ~isnan(diameter_y)
    if nnz(max_edge_top) < 1
        betay=log(top_int/max_int)/((max_locxy(1)-y_centroid)^2-(top_locxy(1)-y_centroid)^2);
        I0y=top_int/exp(-betay*(top_locxy(1)-y_centroid)^2);
    else
        betay=log(top_int/max_int)/((max_edge_top(1)-y_centroid)^2-(top_locxy(1)-y_centroid)^2);
        I0y=top_int/exp(-betay*(top_locxy(1)-y_centroid)^2);
    end
    diameter_y=sigma./sqrt((2*betay));
end

%Check to make sure the diameter is not more than 1.5 times the size of the
%interrogation window
if diameter_x > 1.5*xmax_I
    diameter_x=NaN;
    I0x=NaN;
end
if diameter_y > 1.5*ymax_I
    diameter_y=NaN;
    I0y=NaN;
end

if isnan(diameter_x) && isnan(diameter_y)
    x_c = x_centroid;
    y_c = y_centroid;
    D   = NaN;
    P   = NaN;
    E   = NaN;
    Meth = 0;
    return
end

%Currently the diameter is approximated by a sphere of diameter equal to
%the max size in the x/y direction.  For non-spherical particles, this
%equation should be altered to an ellipse using diameterx and diametery as
%the major and minor axis
diameter=[diameter_x,diameter_y];
I0=[I0x,I0y];
E(1) = sqrt(max(diameter).^2 - min(diameter).^2)./max(diameter);
E(2) = sqrt(max(diameter).^2 - min(diameter).^2)./min(diameter);

%Validation to make sure the results are not negative
if min([x_centroid,y_centroid,diameter,I0])<=0
    x_c      = x_centroid; %x and y should not be negitive
    y_c      = y_centroid;
    D        = NaN;
    P        = NaN;
    Meth     = 0;
    return
end
if sum(isnan(diameter))>1 || sum(isnan(I0))>1
    x_c      = x_centroid; %x and y should not be negitive
    y_c      = y_centroid;
    D        = NaN;
    P        = NaN;
    Meth     = 0;
    return
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3pt Gauss Fit - End
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Initial guesses for lsqnonlin - [I0 beta x_c y_c]
if method < 3
    x0=[I0;(8./diameter.^2);repmat(x_centroid,[1 2]);repmat(y_centroid,[1 2])].';
else
    x0=[max(I0(diameter == max(diameter))),8/max(diameter)^2, x_centroid,y_centroid];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

try
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 1D Methods
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if method < 3
        %point Selection
        if left_locxy(2) > 2 && right_locxy(2) + 2 <= xmax_I
            if numel(max_row) > 1
                points_x = (I(max_locxy(1,1),[left_locxy(2)-2:left_locxy(2),right_locxy(2):right_locxy(2)+2])).';
                points_locx = ([left_locxy(2)-2:left_locxy(2),right_locxy(2):right_locxy(2)+2]).';
            else
                points_x = (I(max_locxy(1,1),max_locxy(1,2)-2:max_locxy(1,2)+2)).';
                points_locx = (max_locxy(1,2)-2:max_locxy(1,2)+2).';
                %If there are brighter particles in the window than the chosen
                %max_int, then use a 3x3 window
                if numel(find(points_x>max_int))>0
                    points_x = (I(max_locxy(1,1),max_locxy(1,2)-1:max_locxy(1,2)+1)).';
                    points_locx = (max_locxy(1,2)-1:max_locxy(1,2)-1).';
                    %keyboard
                end
            end
        else
            %keyboard
            points = I;
            points_locxy = max_locxy;
        end
        if top_locxy(1) > 2 && bottom_locxy(1) + 2 <= ymax_I
            if numel(max_row) > 1
                points_y = I([top_locxy(1,1)-2:top_locxy(1),bottom_locxy(1):bottom_locxy(1)+2],max_locxy(1,2));
                points_locy = ([top_locxy(1,1)-2:top_locxy(1),bottom_locxy(1):bottom_locxy(1)+2]).';
            else
                points_y = (I(max_locxy(1,1)-2:max_locxy(1,1)+2,max_locxy(1,2))).';
                points_locy = (max_locxy(1,1)-2:max_locxy(1,1)+2).';
                %If there are brighter particles in the window than the chosen
                %max_int, then use a 3x3 window
                if numel(find(points_y>max_int))>0
                    points_y = (I(max_locxy(1,1)-1:max_locxy(1,1)+1,max_locxy(1,2))).';
                    points_locy = (max_locxy(1,1)-1:max_locxy(1,1)+1).';
                    %keyboard
                end
            end
        else
            %keyboard
            points = I;
            points_locxy = max_locxy;
        end
        
        %Run a nonlinear least squares regression on the window -
        %method=1 will run the standard least squares equations, and
        %method=2 will use the continuous least squares equations
        xvars=lsqnonlin(@leastsquares1D,  x0(1,1:3)  ,[],[],options,points_x,points_locx,method);
        yvars=lsqnonlin(@leastsquares1D,x0(2,[1 2 4]),[],[],options,points_y,points_locy,method);
        x_c=xvars(3);
        y_c=yvars(3);
        
        %If a good solution could not be reached, try making the points
        %matrix smaller to buffer out any pixels that might be affected
        %by neighboring particles.
        if min([x_c,xvars(1),xvars(2)])<=0 || max(isnan([x_c,xvars(1),xvars(2)]))>0
            %keyboard
            points_x=points(2:size(points,1)-1,2:size(points,2)-1);
            points_locx=points_locxy+1;
            xvars=lsqnonlin(@leastsquares1D,x0,[],[],options,points_x,points_locx,method);
            x_c=xvars(3);
        end
        if min([y_c,yvars(1),yvars(2)])<=0 || max(isnan([y_c,yvars(1),yvars(2)]))>0
            %keyboard
            points_y=points(2:size(points_y,1)-1,2:size(points_y,2)-1);
            points_locy=points_locy+1;
            yvars=lsqnonlin(@leastsquares1D,x0,[],[],options,points_y,points_locy,method);
            y_c=yvars(3);
        end
        
        D(1)=sqrt(sigma^2/(2*xvars(2)));
        D(2)=sqrt(sigma^2/(2*yvars(2)));
        P(1)=xvars(1);
        P(2)=yvars(1);
        E(1) = sqrt(max(D).^2 - min(D).^2)./max(D);
        E(2) = sqrt(max(D).^2 - min(D).^2)./min(D);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % 2D Methods
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    else
        % Point Selection
        if left_locxy(2) > 2 && right_locxy(2) + 2 <= xmax_I && top_locxy(1) > 2 && bottom_locxy(1) + 2 <= ymax_I
            if numel(max_row) > 1
                [xloc yloc] = meshgrid(left_locxy(2)-2:right_locxy(2)+2,top_locxy(1)-2:bottom_locxy(1)+2);
                I2 = I(top_locxy(1)-2:bottom_locxy(1)+2,left_locxy(2)-2:right_locxy(2)+2);
                points = I2(I2 ~=max_int);
                points_locxy(:,1) = xloc(I2 ~=max_int);
                points_locxy(:,2) = yloc(I2 ~=max_int);
            else
                [xloc yloc] = meshgrid(1:xmax_I,1:ymax_I);
                points = reshape(I(max_locxy(1)-2:max_locxy(1)+2,max_locxy(2)-2:max_locxy(2)+2),[25 1]);
                points_locxy = [reshape(xloc(max_locxy(1)-2:max_locxy(1)+2,max_locxy(2)-2:max_locxy(2)+2),[25 1]),...
                    reshape(yloc(max_locxy(1)-2:max_locxy(1)+2,max_locxy(2)-2:max_locxy(2)+2),[25 1])];
                
                %If there are brighter particles in the window than the chosen
                %max_int, then use a 3x3 window
                if numel(find(points>max_int))>0
                    %keyboard
                    points = reshape(I(max_locxy(1)-1:max_locxy(1)+1,max_locxy(2)-1:max_locxy(2)+1),[25 1]);
                    points_locxy = [reshape(xloc(max_locxy(1)-1:max_locxy(1)+1,max_locxy(2)-1:max_locxy(2)+1),[25 1]),...
                        reshape(yloc(max_locxy(1)-1:max_locxy(1)+1,max_locxy(2)-1:max_locxy(2)+1),[25 1])];
                end
            end
        else
            %                     [xloc yloc] = meshgrid(1:xmax_I,1:ymax_I);
            %                     points = I(top_locxy(1):bottom_locxy(1),left_locxy(2):right_locxy(2));
            %                     points = reshape(points,[numel(points) 1]);
            %                     xloc = xloc(top_locxy(1):bottom_locxy(1),left_locxy(2):right_locxy(2));
            %                     yloc = yloc(top_locxy(1):bottom_locxy(1),left_locxy(2):right_locxy(2));
            %                     points_locxy = [reshape(xloc,[numel(xloc) 1]), reshape(yloc,[numel(yloc) 1])];
            [xloc yloc] = meshgrid(1:xmax_I,1:ymax_I);
            points = I(I~=0);
            xloc   = xloc(I~=0);
            yloc   = yloc(I~=0);
            points_locxy = [reshape(xloc,[numel(xloc) 1]), reshape(yloc,[numel(yloc) 1])];
        end
        % Need to have atleast 4 points to solve the nonlinear
        % leasquares problem.  Use the results from the 3pt Gauss
        % when this fails.
        if numel(points)<4
            x_c  = x_centroid;
            y_c  = y_centroid;
            D    = diameter;
            P    = I0;
            Meth = 0;
            return
        end
        
        %Run a nonlinear least squares regression on the window -
        %method=3 will run the standard least squares equations, and
        %method=4 will use the continuous least squares equations
        xvars=lsqnonlin(@leastsquares2D,x0,[],[],options,points,points_locxy,method);
        x_c=xvars(3);
        y_c=xvars(4);
        
        %If a good solution could not be reached, try making the points
        %matrix smaller to buffer out any pixels that might be affected
        %by neighboring particles.
        if min([x_c,y_c,xvars(1)])<=0 || max(isnan([x_c,y_c,xvars(1),xvars(2)]))>0
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %didn't check betas = xvars(2)
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            points=points(2:size(points,1)-1,2:size(points,2)-1);
            points_locxy=points_locxy+1;
            % Need to have atleast 4 points to solve the nonlinear
            % leasquares problem.  Use the results from the 3pt Gauss
            % when this fails.
            if numel(points)>=4
                xvars=lsqnonlin(@leastsquares2D,x0,[],[],options,points,points_locxy,method);
                x_c=xvars(3);
                y_c=xvars(4);
            else
                
                x_c  = x_centroid;
                y_c  = y_centroid;
                D    = diameter;
                P    = I0;
                Meth = 0;
                return
            end
        end
        
        D=sqrt(sigma^2/(2*abs(xvars(2))));
        P=xvars(1);
    end
    
    %Return Nan's if there was an error
catch errmess
    fprintf(errmess.message)
    keyboard
    x_c  = x_centroid;
    y_c  = y_centroid;
    D    = diameter;
    P    = I0;
    Meth = 0;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Sub Function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function F = leastsquares1D(x,mapint_i,locxy_i,method)
        %This function is called by lsqnonlin if the least squares or continuous
        %least squares method has been chosen. x contains initial guesses[I0, betas, x_c,
        %y_c]. mapint_i is a matrix containing pixel intensity values, and locxy_i
        %is a 1x2 vector containing the row/column coordinates of the top left
        %pixel in mapint_i
        %
        %F is the variable being minimized - the difference between the gaussian
        %curve and the actual intensity values.
        %
        %Adapted from M. Brady's 'leastsquaresgaussfit' and 'mapintensity'
        %B.Drew - 7.18.2008
        
        Io=x(1);
        betas=x(2);
        x_cent=x(3);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Just like in the continuous four-point method, lsqnonlin tries negative
        %values for x(2), which will return errors unless the abs() function is
        %used in front of all the x(2)'s.
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %num1=pi/4;
        num1=0.5*Io*sqrt(pi);
        num2=sqrt(abs(betas));
        S = size(mapint_i);
        
        if method==1
            gauss_int = zeros(S);
            mapint_mod = mapint_i(:);
            xp = zeros(S);
            for ii = 1:length(mapint_i)
                xp(ii) = locxy_i(ii);
                % map an intensity profile of a gaussian function:
                gauss_int(ii)=Io*exp(-abs(betas)*((xp(ii))-x_cent)^2);
            end
            
        elseif method==2
            gauss_int = zeros(S);
            mapint_mod = zeros(S);
            xp1 = zeros(S);
            erfx1 = zeros(S);
            erfx2 = zeros(S);
            for ii = 1:length(mapint_i)
                xp1(ii) = locxy_i(ii)-0.5;
                erfx1(ii) = erf(num2*(xp1(ii)-x_cent));
                erfx2(ii) = erf(num2*(xp1(ii)+1-x_cent));
                % map an intensity profile of a gaussian function:
                gauss_int(ii)=(num1/num2)*(erfx2(ii)-erfx1(ii));
                mapint_mod(ii)=mapint_i(ii);
            end
        end
        
        % compare the Gaussian curve to the actual pixel intensities
        F=mapint_mod-gauss_int;
    end

    function F = leastsquares2D(x,mapint_i,locxy_i,method)
        %This function is called by lsqnonlin if the least squares or continuous
        %least squares method has been chosen. x contains initial guesses[I0, betas, x_c,
        %y_c]. mapint_i is a matrix containing pixel intensity values, and locxy_i
        %is a 1x2 vector containing the row/column coordinates of the top left
        %pixel in mapint_i
        %
        %F is the variable being minimized - the difference between the gaussian
        %curve and the actual intensity values.
        %
        %Adapted from M. Brady's 'leastsquaresgaussfit' and 'mapintensity'
        %B.Drew - 7.18.2008
        
        Io=x(1);
        betas=x(2);
        x_cent=x(3);
        y_cent=x(4);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Just like in the continuous four-point method, lsqnonlin tries negative
        %values for x(2), which will return errors unless the abs() function is
        %used in front of all the x(2)'s.
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        num1=(Io*pi)/4;
        num2=sqrt(abs(betas));
        
        if method==3
            gauss_int = zeros(size(mapint_i));
            mapint_mod = zeros(size(mapint_i));
            xp = zeros(size(mapint_i));
            yp = zeros(size(mapint_i));
            for ii = 1:length(mapint_i)
                xp(ii) = locxy_i(ii,1);
                yp(ii) = locxy_i(ii,2);
            end
            
            % map an intensity profile of a gaussian function:
            for rr = 1:length(xp)
                gauss_int(rr)=Io*exp(-abs(betas)*(((xp(rr))-x_cent)^2 + ...
                    ((yp(rr))-y_cent)^2));
                mapint_mod(rr)=mapint_i(rr);
            end
            
        elseif method==4
            S = size(mapint_i);
            gauss_int = zeros(S(1),S(2));
            mapint_mod = zeros(S(1),S(2));
            xp = zeros(size(mapint_i));
            yp = zeros(size(mapint_i));
            for ii = 1:length(mapint_i)
                xp(ii) = locxy_i(ii,1)-0.5;
                yp(ii) = locxy_i(ii,2)-0.5;
                erfx1 = erf(num2*(xp(ii)-x_cent));
                erfy1 = erf(num2*(yp(ii)-y_cent));
                erfx2 = erf(num2*(xp(ii)+1-x_cent));
                erfy2 = erf(num2*(yp(ii)+1-y_cent));
                % map an intensity profile of a gaussian function:
                gauss_int(ii)=(num1/abs(betas))*(erfx1*(erfy1-erfy2)+erfx2*(-erfy1+erfy2));
                mapint_mod(ii)=mapint_i(ii);
            end
        end
        mapint_i=mapint_mod;
        % compare the Gaussian curve to the actual pixel intensities
        F=mapint_i-gauss_int;
    end
end