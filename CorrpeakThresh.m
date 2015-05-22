function [U,V,Eval] = CorrpeakThresh(U,V,C,PR,height_thresh,ratio_thresh,Eval)
% --- Correlation Peak Thresholding Validation Subfunction ---

%     This file is part of prana, an open-source GUI-driven program for
%     calculating velocity fields using PIV or PTV.
% 
%     Copyright 2015.  Los Alamos National Security, LLC. This material was
%     produced under U.S. Government contract DE-AC52-06NA25396 for Los 
%     Alamos National Laboratory (LANL), which is operated by Los Alamos 
%     National Security, LLC for the U.S. Department of Energy. The U.S. 
%     Government has rights to use, reproduce, and distribute this software.
%     NEITHER THE GOVERNMENT NOR LOS ALAMOS NATIONAL SECURITY, LLC MAKES ANY
%     WARRANTY, EXPRESS OR IMPLIED, OR ASSUMES ANY LIABILITY FOR THE USE OF
%     THIS SOFTWARE.  If software is modified to produce derivative works,
%     such modified software should be clearly marked, so as not to confuse
%     it with the version available from LANL.

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


%neglect thresholds if Eval isn't specified?
%copied from Thresh, but doesn't make complete sense
if nargin<=6
    height_thresh = 0;
    ratio_thresh  = 0;
end

S=size(U);

%thresholding
for i=1:S(1)
    for j=1:S(2)
        if Eval(i,j)==0
            %correlation threshold condition
            if C(i,j)<height_thresh
                U(i,j)=nan;
                V(i,j)=nan;
                Eval(i,j)=200;
            end
            if PR(i,j)<ratio_thresh
                U(i,j)=nan;
                V(i,j)=nan;
                Eval(i,j)=250;
            end
        elseif Eval(i,j)==-1
            %boundary condition
            U(i,j)=nan;
            V(i,j)=nan;
        end
    end
end

end
