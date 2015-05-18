function [U,V,Eval] = Thresh(U,V,uthreshold,vthreshold,Eval)
% --- Thresholding Validation Subfunction ---

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
