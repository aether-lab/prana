function [u,v,M,D]=subpixel(G,ccsizex,ccsizey,W,Method,Peakswitch,d)
%intialize indices
cc_x = -floor(ccsizex/2):ceil(ccsizex/2)-1;
cc_y = -floor(ccsizey/2):ceil(ccsizey/2)-1;

%find maximum correlation value
[M,I] = max(G(:));

% Use 4 standard deviations for the peak sizing (e^-2)
sigma = 4;

%if correlation empty
if M==0
    if Peakswitch
        u=zeros(1,3);
        v=zeros(1,3);
        M=zeros(1,3);
        D=zeros(1,3);
    else
        u=0; v=0; M=0; D=0; 
    end
else
    if Peakswitch
        u=zeros(1,3);
        v=zeros(1,3);
        D=zeros(1,3);
        %Locate peaks using imregionalmax
        A=imregionalmax(G);
        peakmat=G.*A;
        for i=2:3
            peakmat(peakmat==M(i-1))=0;
            [M(i),I(i)]=max(peakmat(:));
        end
        j=length(M);
    else
        u=zeros(1,1);
        v=zeros(1,1);
        D=zeros(1,1);
        j=1;    
    end
    
    for i=1:j
        method=Method;
        
        %find x and y indices
        shift_locy = 1+mod(I(i)-1,ccsizey);
        shift_locx = ceil(I(i)/ccsizey);

        shift_errx=[];
        shift_erry=[];
        %find subpixel displacement in x
        if shift_locx == 1
            %boundary condition 1
            shift_errx =  G( shift_locy , shift_locx+1 )/M(i); method=1;
        elseif shift_locx == ccsizex
            %boundary condition 2
            shift_errx = -G( shift_locy , shift_locx-1 )/M(i); method=1;
        elseif G( shift_locy , shift_locx+1 ) == 0
            %endpoint discontinuity 1
            shift_errx = -G( shift_locy , shift_locx-1 )/M(i); method=1;
        elseif G( shift_locy , shift_locx-1 ) == 0
            %endpoint discontinuity 2
            shift_errx =  G( shift_locy , shift_locx+1 )/M(i); method=1;
        end
        if shift_locy == 1
            %boundary condition 1
            shift_erry = -G( shift_locy+1 , shift_locx )/M(i); method=1;
        elseif shift_locy == ccsizey
            %boundary condition 2
            shift_erry =  G( shift_locy-1 , shift_locx )/M(i); method=1;
        elseif G( shift_locy+1 , shift_locx ) == 0
            %endpoint discontinuity 1
            shift_erry =  G( shift_locy-1 , shift_locx )/M(i); method=1;
        elseif G( shift_locy-1 , shift_locx ) == 0
            %endpoint discontinuity 2
            shift_erry = -G( shift_locy+1 , shift_locx )/M(i); method=1;
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
            
            points=[shift_locy   shift_locx   G(shift_locy  ,shift_locx  );...
                    shift_locy-1 shift_locx   G(shift_locy-1,shift_locx  );...
                    shift_locy   shift_locx-1 G(shift_locy  ,shift_locx-1);...
                    shift_locy   shift_locx+1 G(shift_locy  ,shift_locx+1)];
                
            [~,IsortI] = sort(points(:,3),'descend');
            points = points(IsortI,:);

            x1=points(1,2); x2=points(2,2); x3=points(3,2); x4=points(4,2);
            y1=points(1,1); y2=points(2,1); y3=points(3,1); y4=points(4,1);
            a1=points(1,3); a2=points(2,3); a3=points(3,3); a4=points(4,3);

            alpha(1) = (x4^2)*(y2 - y3) + (x3^2)*(y4 - y2) + ((x2^2) + (y2 - y3)*(y2 - y4))*(y3 - y4);
            alpha(2) = (x4^2)*(y3 - y1) + (x3^2)*(y1 - y4) - ((x1^2) + (y1 - y3)*(y1 - y4))*(y3 - y4);
            alpha(3) = (x4^2)*(y1 - y2) + (x2^2)*(y4 - y1) + ((x1^2) + (y1 - y2)*(y1 - y4))*(y2 - y4);
            alpha(4) = (x3^2)*(y2 - y1) + (x2^2)*(y1 - y3) - ((x1^2) + (y1 - y2)*(y1 - y3))*(y2 - y3);

            gamma(1) = (-x3^2)*x4 + (x2^2)*(x4 - x3) + x4*((y2^2) - (y3^2)) + x3*((x4^2) - (y2^2) + (y4^2)) + x2*(( x3^2) - (x4^2) + (y3^2) - (y4^2));
            gamma(2) = ( x3^2)*x4 + (x1^2)*(x3 - x4) + x4*((y3^2) - (y1^2)) - x3*((x4^2) - (y1^2) + (y4^2)) + x1*((-x3^2) + (x4^2) - (y3^2) + (y4^2));
            gamma(3) = (-x2^2)*x4 + (x1^2)*(x4 - x2) + x4*((y1^2) - (y2^2)) + x2*((x4^2) - (y1^2) + (y4^2)) + x1*(( x2^2) - (x4^2) + (y2^2) - (y4^2));
            gamma(4) = ( x2^2)*x3 + (x1^2)*(x2 - x3) + x3*((y2^2) - (y1^2)) - x2*((x3^2) - (y1^2) + (y3^2)) + x1*((-x2^2) + (x3^2) - (y2^2) + (y3^2));

            delta(1) = x4*(y2 - y3) + x2*(y3 - y4) + x3*(y4 - y2);
            delta(2) = x4*(y3 - y1) + x3*(y1 - y4) + x1*(y4 - y3);
            delta(3) = x4*(y1 - y2) + x1*(y2 - y4) + x2*(y4 - y1);
            delta(4) = x3*(y2 - y1) + x2*(y1 - y3) + x1*(y3 - y2);

            deno = 2*(log(a1)*delta(1) + log(a2)*delta(2) + log(a3)*delta(3) + log(a4)*delta(4));

            x_centroid = (log(a1)*alpha(1) + log(a2)*alpha(2) + log(a3)*alpha(3) + log(a4)*alpha(4))/deno;
            y_centroid = (log(a1)*gamma(1) + log(a2)*gamma(2) + log(a3)*gamma(3) + log(a4)*gamma(4))/deno;
            shift_errx=x_centroid-shift_locx;
            shift_erry=y_centroid-shift_locy;
            
            betas = abs((log(a2)-log(a1))/((x2-x_centroid)^2+(y2-y_centroid)^2-(x1-x_centroid)^2-(y1-y_centroid)^2));
            D(i)=sqrt(sigma^2/(2*betas));
            
        elseif any(method==[3 4])
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%
            % Gaussian Least Squares %
            %%%%%%%%%%%%%%%%%%%%%%%%%%
            
            %Find a suitable window around the peak (5x5 preferred)
            x_min=shift_locx-ceil(sqrt(2).*d(1)/2); x_max=shift_locx+ceil(sqrt(2).*d(1)/2);
            y_min=shift_locy-ceil(sqrt(2).*d(2)/2); y_max=shift_locy+ceil(sqrt(2).*d(2)/2);
            if x_min<1
                x_min=1;
            end
            if x_max>ccsizex
                x_max=ccsizex;
            end
            if y_min<1
                y_min=1;
            end
            if y_max>ccsizey
                y_max=ccsizey;
            end
            points=G(y_min:y_max,x_min:x_max).*W(y_min:y_max,x_min:x_max);
            
            %Options for the lsqnonlin solver
            options=optimset('MaxIter',1200,'MaxFunEvals',5000,'TolX',5e-6,'TolFun',5e-6,...
                'LargeScale','off','Display','off','DiffMinChange',1e-7,'DiffMaxChange',1,...
                'Algorithm','levenberg-marquardt');
            
            %Initial values for the solver
            x0=[M(i) 1 1 shift_locx shift_locy 0];

            [xloc yloc]=meshgrid(x_min:x_max,y_min:y_max);

            %Run solver; default to 3-point gauss if it fails
            try
                xvars=lsqnonlin(@leastsquares2D,x0,[],[],options,points(:),[yloc(:),xloc(:)],method);
%                 shift_errx=xvars(3)-shift_locx;
%                 shift_erry=xvars(4)-shift_locy;
                shift_errx=xvars(4)-shift_locx;
                shift_erry=xvars(5)-shift_locy;
                D(i) = sqrt(sigma^2/(2*sqrt(xvars(2).^2 + xvars(3).^2)));
            catch %#ok
                method=1;
            end
        end
        if method==1

            %%%%%%%%%%%%%%%%%%%%
            % 3-Point Gaussian %
            %%%%%%%%%%%%%%%%%%%%
            
            if isempty(shift_errx)
                %gaussian fit
                lCm1 = log(G( shift_locy , shift_locx-1 )*W( shift_locy , shift_locx-1 ));
                lC00 = log(G( shift_locy , shift_locx   )*W( shift_locy , shift_locx   ));
                lCp1 = log(G( shift_locy , shift_locx+1 )*W( shift_locy , shift_locx+1 ));
                if (2*(lCm1+lCp1-2*lC00)) == 0
                    shift_errx = 0;
                    Dx = nan;
                else
                    shift_errx = (lCm1-lCp1)/(2*(lCm1+lCp1-2*lC00));
                    betax = abs(lCm1-lC00)/((-1-shift_errx)^2-(shift_errx)^2);
                    Dx = sigma./sqrt((2*betax));
                end
            else
                Dx = nan;
            end
            
            if isempty(shift_erry)
                lCm1 = log(G( shift_locy-1 , shift_locx )*W( shift_locy-1 , shift_locx ));
                lC00 = log(G( shift_locy   , shift_locx )*W( shift_locy   , shift_locx ));
                lCp1 = log(G( shift_locy+1 , shift_locx )*W( shift_locy+1 , shift_locx ));
                if (2*(lCm1+lCp1-2*lC00)) == 0
                    shift_erry = 0;
                    Dy = nan;
                else
                    shift_erry = (lCm1-lCp1)/(2*(lCm1+lCp1-2*lC00));
                    betay = abs(lCm1-lC00)/((-1-shift_erry)^2-(shift_erry)^2);
                    Dy = sigma./sqrt((2*betay));
                end
            else
                Dy = nan;
            end
            
            D(i) = nanmean([Dx Dy]);
        end
        
        u(i)=cc_x(shift_locx)+shift_errx;
        v(i)=cc_y(shift_locy)+shift_erry;
        
        if isinf(u(i)) || isinf(v(i))
            u(i)=0; v(i)=0;
        end
    end
end

end