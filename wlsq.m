function [xh]=wlsq(y,H,W)
% --- Weighted Least Squares Fit for Phase Correlation ---

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

tempmat=sortrows([y',H',W'],2);
y=tempmat(:,1);
H=tempmat(:,2);
W=diag(tempmat(:,3));

% xh=inv(H'*W*H)*H'*W*y;
xh=(H'*W*H)\(H'*W*y);

end
