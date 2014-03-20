function [U,V,Eval] = Thresh(U,V,uthreshold,vthreshold,Eval)
% --- Thresholding Validation Subfunction ---

%neglect u and v threshold
if nargin<=4
    uthreshold = [-inf inf];
    vthreshold = [-inf inf];
end

S=size(U);

%thresholding
for i=1:S(1)
    for j=1:S(2)
        if Eval(i,j)==0
            %velocity threshold condition
            if U(i,j)<uthreshold(1) || U(i,j)>uthreshold(2) || V(i,j)<vthreshold(1) || V(i,j)>vthreshold(2)
                U(i,j)=nan;
                V(i,j)=nan;
                Eval(i,j)=100;
            end
        elseif Eval(i,j)==-1
            %boundary condition
            U(i,j)=nan;
            V(i,j)=nan;
        end
    end
end

end