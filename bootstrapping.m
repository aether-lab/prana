function [U,V,Eval] = bootstrapping(X,Y,u,v,per,iter,kmax,Eval)
% Bootstrapping Validation Subfunction 
%
% [U,V,Eval] = bootstraping(x,y,u,v,per,iter,kmax,Eval)
%
% per  = percent removed for each interpolation (0-1)
% iter = number of interpolations per frame (for histogram)
% kmax = number of passes 

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


n = size(X);

M = zeros(n(1),n(2),iter);

tol = 0.3;
ktol = 1;

while tol > 0 && ktol <= kmax+1
    U = zeros(n(1),n(2),iter);
    V = zeros(n(1),n(2),iter);

    for i = 1:iter
        clear S m Up Vp Ui Vi
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Data Removal
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        [m]= bootstrapping_dataremove(size(U),per,sum(Eval,3));

        S(:,1) = X(m==1);
        S(:,2) = Y(m==1);

        Up = u(m==1);
        Vp = v(m==1);
        M(:,:,i) = m;
        
        Ui = gridfit(S(:,1),S(:,2),Up,X(1,:),Y(:,1));
        Vi = gridfit(S(:,1),S(:,2),Vp,X(1,:),Y(:,1));
        U(:,:,i) = Ui;
        V(:,:,i) = Vi;

    end

    PBad = 0;
    for j = 1:n(1)
        for k = 1:n(2)
            if sum(isnan(U(j,k,:))) == 0
                try              
                    [H.U,HX.U] = hist(squeeze(U(j,k,:)),iter/2);
                    [H.V,HX.V] = hist(squeeze(V(j,k,:)),iter/2);

                    modeU = HX.U(H.U==max(H.U));
                    modeV = HX.V(H.V==max(H.V));

                    tU = abs((modeU - U(j,k))/modeU);
                    tV = abs((modeV - V(j,k))/modeV);
                    if tU > tol || tV > tol && Eval(j,k) ~= -1
                        U(j,k) = modeU(1);
                        V(j,k) = modeV(1);
                        Eval(j,k) = 200;
                        PBad = PBad+1;
                    end
                catch%#ok
                    Ems=lasterror;%#ok
                    fprintf('\n\n')
                    fprintf(Ems.message);
                end
            end
        end
    end
    ktol = ktol + 1;
    tol = tol-(tol/(kmax-1))*(ktol-1);
end

end
