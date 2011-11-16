function [U,V,Eval] = UOD(U,V,t,tol,Eval)
% --- Universal Outlier Detection Validation Subfunction ---

%number of validation passes
pass = length(tol);

S=size(U);

%outlier searching
for k=1:pass
    
    q = (t(k,:)-1)/2;
    
    for i=1:S(1)
        for j=1:S(2)
            if Eval(i,j)==0           
                %get evaluation block with at least 8 valid points
                s=0;
                while s==0
                    Imin = max([i-q(2) 1   ]);
                    Imax = min([i+q(2) S(1)]);
                    Jmin = max([j-q(1) 1   ]);
                    Jmax = min([j+q(1) S(2)]);
                    Iind = Imin:Imax;
                    Jind = Jmin:Jmax;
                    Ublock = U(Iind,Jind);
                    if length(Ublock(~isnan(Ublock)))>=8
%                         Xblock = X(Iind,Jind)-X(i,j);
%                         Yblock = Y(Iind,Jind)-Y(i,j);
                        Vblock = V(Iind,Jind);
                        s=1;
                    else
                        q=q+1;
                    end
                end

%                 %distance from vector location
%                 Dblock = (Xblock.^2+Yblock.^2).^0.5;
%                 Dblock(isnan(Ublock))=nan;

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

end