function [x_centroid,y_centroid,diameter,I0,eccentricity] = threeptgaussfit(I,locxy,sigma)
%
%[x_centroid,y_centroid,diameter]=threeptgausfit(I,locxy_i, sigma)
%
%this function estimates a particle's diameter given an input particle
%intensity profile(I), relative locator to the entire image(locxy_i),
%and previously derived particle centroid location(xcent,ycent)
%
%The function calculated the particle centroid and diameter by assuming the
%light scattering profile (particle intensity profile) is Gaussian in
%shape.  Three points are chosen in both the x&y dimension which includes
%the maximum intensity pixel of the particle profile.
%
%I        - input intensity profile
%locxy    - location on max(iums) intensities
%sigma    - number of standard deviations in one diameter
%
%N.Cardwell - 2.27.2008 (3-pt Gaussian method taken from the Master's Thesis of M. Brady)
%B.Drew     - 7.31.2008
%S.Raben    - 9.20.08

% Find extents of I
[ymax_I xmax_I] = size(I);

% Find the maxium value of I and all of its locations
[max_int,int_L] = max(I(:));%#ok
%        max_row = locxy(:,1);
%        max_col = locxy(:,2);
max_locxy(1) = round(median(locxy(:,1)));
max_locxy(2) = round(median(locxy(:,2)));
diameter_x = 0;
diameter_y = 0;

Better = 0;% if =0 then NO pixel selection just left 1 center and right 1.
%This method cannot fit a Gaussian shape to a particle intensity profile
%where the highest intensity pixel is on the edge of I (assumes a
%Gaussian light scattering profile)
if Better ~= 1
    LE = [max_locxy(1)   max_locxy(2)-1];
    RE = [max_locxy(1)   max_locxy(2)+1];
    TE = [max_locxy(1)-1 max_locxy(2)];
    BE = [max_locxy(1)+1 max_locxy(2)];
    
    %Given that this method assumes a Gaussian shaped intensity profile, the
    %function is designed to return NaN should any of the nearest neighbor
    %pixels to max_int equal zero or if the particles have equal intensity
    %values
    if LE(2) > 1 && RE(2) < xmax_I %edge_ckx == 0;
        left_int  = I(LE(1),LE(2));
        right_int = I(RE(1),RE(2));
        LR_denom = (2*(log(left_int)+log(right_int)-2*log(max_int)));
        if left_int==0 || right_int==0 || LR_denom==0
            x_centroid = max_locxy(2);
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
    
    if TE(1) > 1 && BE(1) < ymax_I %edge_cky == 0;
        top_int    = I(TE(1),TE(2));
        bottom_int = I(BE(1),BE(2));
        TB_denom = (2*(log(top_int)+log(bottom_int)-2*log(max_int)));
        if top_int==0 || bottom_int==0 || TB_denom==0
            y_centroid = max_locxy(1);
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
        betax=abs(log(left_int/max_int)/((LE(2)-x_centroid)^2-(max_locxy(2)-x_centroid)^2));
        I0x=left_int/exp(-betax*(LE(2)-x_centroid)^2);
        diameter_x=sigma./sqrt((2*betax));
    end
    
    if ~isnan(diameter_y)
        betay=abs(log(top_int/max_int)/((TE(1)-y_centroid)^2-(max_locxy(1)-y_centroid)^2));
        I0y=top_int/exp(-betay*(TE(1)-y_centroid)^2);
        diameter_y=sigma./sqrt((2*betay));
    end
    
    % General Solution to the 3pt Gausian Estimator
else
    
    [TE,BE,LE,RE] = Edgesearch(I,max_locxy);
    
    if numel(locxy(:,2)) > 1
        if TE(1) > 1 && BE(1) < ymax_I && I(TE(1),TE(2)) ~= 0 && I(BE(1),BE(2)) ~= 0
            if I(TE(1),TE(2)) > I(BE(1),BE(2))
                [TBLoc] = find(I(1:TE(1),TE(2)) < I(TE(1),TE(2)));%#ok matlab is just being stupid
                maxlocTB = [max(TBLoc) TE(2)];
            else
                [TBLoc] = find(I(BE(1):end,BE(2)) < I(BE(1),BE(2)));
                maxlocTB = [min(BE(1)-1+TBLoc) BE(2)];
            end
        else
            maxlocTB = max_locxy;
        end
        
        if LE(2) > 1 && RE(2) < xmax_I && I(LE(1),LE(2)) ~= 0 && I(RE(1),RE(2)) ~= 0
            if I(LE(1),LE(2)) > I(RE(1),RE(2))
                [LRLoc] = find(I(LE(1),1:LE(2)) < I(LE(1),LE(2)));%#ok matlab is just being stupid
                maxlocLR = [LE(1) max(LRLoc)];
            else
                [LRLoc] = find(I(RE(1),RE(2):end) < I(RE(1),RE(2)));
                maxlocLR = [RE(1) min(RE(2)-1+LRLoc)];
            end
        else
            maxlocLR = max_locxy;
        end
        
    else
        maxlocTB   = max_locxy;
        maxlocLR   = max_locxy;
    end
    try
        points_LR = [LE(1) LE(2) I(LE(1),LE(2)); RE(1) RE(2) I(RE(1),RE(2)); maxlocLR(1) maxlocLR(2) I(maxlocLR(1),maxlocLR(2))];
        points_TB = [TE(1) TE(2) I(TE(1),TE(2)); BE(1) BE(2) I(BE(1),BE(2)); maxlocTB(1) maxlocTB(2) I(maxlocTB(1),maxlocTB(2))];
    catch
        %keyboard
    end
    
    [LRD,LRI] = sort(points_LR(:,3),'ascend');
    [TBD,TBI] = sort(points_TB(:,3),'ascend');
    points_LR = points_LR(LRI,:);
    points_TB = points_TB(TBI,:);
    %Given that this method assumes a Gaussian shaped intensity profile, the
    %function is designed to return NaN should any of the nearest neighbor
    %pixels to max_int equal zero or if the particles have equal intensity
    %values
    if LE(2) ~= 1 && RE(2) ~= xmax_I %edge_ckx == 0;
        LR_denom = 2*((-points_LR(2,2) + points_LR(3,2))*log(points_LR(1,3))...
            + (points_LR(1,2) - points_LR(3,2))*log(points_LR(2,3)) + (-points_LR(1,2) + points_LR(2,2))*log(points_LR(3,3)));
        if points_LR(1,3)==0 || points_LR(3,3)==0 || LR_denom==0
            x_centroid = max_locxy(2);
            diameter_x=NaN;
            I0x=NaN;
        else
            %calculate the x&y centroid location using the 3-pt Gaussian Method
            %x_centroid=max_locxy(2)+log(left_int/right_int)/LR_denom;
            x_centroid = ((points_LR(1,2)^2 - points_LR(3,2)^2)*log(points_LR(2,3)/points_LR(1,3))...
                + (-points_LR(1,2)^2 + points_LR(2,2)^2)*log(points_LR(3,3)/points_LR(1,3)))/(LR_denom);
        end
    else
        x_centroid = max_locxy(2);
        diameter_x=NaN;
        I0x=NaN;
    end
    
    if TE(1) ~= 1 && BE(1) ~= ymax_I %edge_cky == 0;
        TB_denom = 2*((-points_TB(2,1) + points_TB(3,1))*log(points_TB(1,3))...
            + (points_TB(1,1) - points_TB(3,1))*log(points_TB(2,3)) + (-points_TB(1,1) + points_TB(2,1))*log(points_TB(3,3)));
        if points_TB(1,3)==0 || points_TB(3,3)==0 || TB_denom==0
            y_centroid = max_locxy(1);
            %keyboard
            diameter_y=NaN;
            I0y=NaN;
        else
            %calculate the x&y centroid location using the 3-pt Gaussian Method
            %y_centroid=max_locxy(1)+log(top_int/bottom_int)/TB_denom;
            y_centroid = ((points_TB(1,1)^2 - points_TB(3,1)^2)*log(points_TB(2,3)/points_TB(1,3))...
                + (-points_TB(1,1)^2 + points_TB(2,1)^2)*log(points_TB(3,3)/points_TB(1,3)))/(TB_denom);
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
        %                     betax=abs(log(left_int/max_int)/((max_locxy(2)-x_centroid)^2-(LE(2)-x_centroid)^2));
        %                     I0x=left_int/exp(-betax*(LE(2)-x_centroid)^2);
        sigma_x = sqrt(abs(((-points_LR(1,2)+points_LR(2,2))*(points_LR(1,2)-points_LR(3,2))...
            *(points_LR(2,2)-points_LR(3,2))/LR_denom)));
        I0x = points_LR(3,3)/exp((-(points_LR(3,2)-x_centroid)^2)/(2*sigma_x^2));
        diameter_x=sigma.*sigma_x;
    end
    
    if ~isnan(diameter_y)
        %                     betay=abs(log(top_int/max_int)/((max_locxy(1)-y_centroid)^2-(TE(1)-y_centroid)^2));
        %                     I0y=top_int/exp(-betay*(TE(1)-y_centroid)^2);
        sigma_y = sqrt(abs(((-points_TB(1,1)+points_TB(2,1))*(points_TB(1,1)-points_TB(3,1))...
            *(points_TB(2,1)-points_TB(3,1))/TB_denom)));
        I0y = points_TB(3,3)/exp((-(points_TB(3,1)-y_centroid)^2)/(2*sigma_y^2));
        diameter_y=sigma.*sigma_y;
    end
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

%Currently the diameter is approximated by a sphere of diameter equal to
%the max size in the x/y direction.  For non-spherical particles, this
%equation should be altered to an ellipse using diameterx and diametery as
%the major and minor axis
diameter=[diameter_x,diameter_y];
I0=[I0x,I0y];
eccentricity(1) = sqrt(max(diameter).^2 - min(diameter).^2)./max(diameter);
eccentricity(2) = sqrt(max(diameter).^2 - min(diameter).^2)./min(diameter);

%Validation to make sure the results are not negative
if min([x_centroid,y_centroid,diameter,I0])<=0
    x_centroid=max_locxy(2);
    y_centroid=max_locxy(1);
    diameter=NaN;
    I0=NaN;
    return
end
end
