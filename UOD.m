function [U,V,Eval] = UOD(U,V,t,tol,Eval)
% --- Universal Outlier Detection Validation Subfunction ---

%     This file is part of prana, an open-source GUI-driven program for
%     calculating velocity fields using PIV or PTV.
%     Copyright (C) 2012  Virginia Polytechnic Institute and State
%     University
% 
%     prana is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     (at your option) any later version.
% 
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
% 
%     You should have received a copy of the GNU General Public License
%     along with this program.  If not, see <http://www.gnu.org/licenses/>.

%number of validation passes
pass = length(tol);

S=size(U);

%outlier searching
for k=1:pass
    
    q0 = floor((t(k,:)-1)/2);
    
    for i=1:S(1)
        for j=1:S(2)
            if Eval(i,j)==0           
                %get evaluation block with at least 8 valid points
                s=0;
                q = q0;
                while s==0
                    Imin = max([i-q(2) 1   ]);
                    Imax = min([i+q(2) S(1)]);
                    Jmin = max([j-q(1) 1   ]);
                    Jmax = min([j+q(1) S(2)]);
                    Iind = Imin:Imax;
                    Jind = Jmin:Jmax;
                    Ublock = U(Iind,Jind);
                    if length(Ublock(~isnan(Ublock)))>=8 || any(q >= 2*q0)
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
