function [Uval,Vval,Evalval,Cval]=VAL(X,Y,U,V,Eval,C,Threshswitch,UODswitch,Bootswitch,extrapeaks,Uthresh,Vthresh,UODwinsize,UODthresh,Bootper,Bootiter,Bootkmax)
% --- Validation Subfunction ---
if extrapeaks
    j=3;
else
    j=1;
end

[X,Y,U,V,Eval,C]=matrixform(X,Y,U,V,Eval,C);
Uval=U(:,:,1);Vval=V(:,:,1);Evalval=Eval(:,:,1);
if ~isempty(C)
    Cval=C(:,:,1);
else
    Cval=[];
end
S=size(X);

if Threshswitch || UODswitch
    for i=1:j
        %Thresholding
        if Threshswitch
            [Uval,Vval,Evalval] = Thresh(Uval,Vval,Uthresh,Vthresh,Evalval);
        end

        %Univeral Outlier Detection
        if UODswitch
            t=permute(UODwinsize,[2 3 1]);
            t=t(:,t(1,:)~=0);
            [Uval,Vval,Evalval] = UOD(Uval,Vval,t',UODthresh,Evalval);
        end
%         disp([num2str(sum(sum(Evalval>0))),' bad vectors'])
        %Try additional peaks where validation failed
        if i<j
            Utemp=U(:,:,i+1);Vtemp=V(:,:,i+1);Evaltemp=Eval(:,:,i+1);Ctemp=C(:,:,i+1);
            Uval(Evalval>0)=Utemp(Evalval>0);
            Vval(Evalval>0)=Vtemp(Evalval>0);
            Evalval(Evalval>0)=Evaltemp(Evalval>0);
            Cval(Evalval>0)=Ctemp(Evalval>0);
        end
    end
end

%replacement
for i=1:S(1)
    for j=1:S(2)
        if Evalval(i,j)~=0
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
                Ublock = Uval(Iind,Jind);
                if length(Ublock(~isnan(Ublock)))>=8
                    Xblock = X(Iind,Jind)-X(i,j);
                    Yblock = Y(Iind,Jind)-Y(i,j);
                    Vblock = Vval(Iind,Jind);
                    s=1;
                end
            end
            
            %distance from erroneous vector
            Dblock = (Xblock.^2+Yblock.^2).^0.5;
            Dblock(isnan(Ublock))=nan;

            %validated vector
            Uval(i,j) = nansum(nansum(Dblock.*Ublock))/nansum(nansum(Dblock));
            Vval(i,j) = nansum(nansum(Dblock.*Vblock))/nansum(nansum(Dblock));       
        end
    end
end

%Bootstrapping
if Bootswitch
    [Uval,Vval,Evalval] = bootstrapping(X,Y,Uval,Vval,Bootper,Bootiter,Bootkmax,Evalval);
end

%convert back to vector
[Uval,Vval,Evalval,Cval]=vectorform(X,Y,Uval,Vval,Evalval,Cval);

end