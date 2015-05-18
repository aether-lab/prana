function []=write_dat_val_C(fname,X,Y,U,V,Eval,C,D,strand,T,frametitle,t_opt)
% --- .dat Writer Subfunction ---

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

if nargin<12
    t_opt=[];
end

%find I,J for plt
S = size(U);

%generate text file
fid = fopen(fname,'w');
if fid==-1
    error(['error creating file ',fname])
end

varlist='"X" "Y" "U" "V" "Eval"';

if ~isempty(C)
    varlist=[varlist,' "C" "D"'];
    if size(U,3)>1
        for i=2:size(U,3)
            varlist=[varlist,' "U',num2str(i-1),'" "V',num2str(i-1),'"'];
        end
    end
    if size(C,3)>1
        for i=2:size(C,3)
            varlist=[varlist,' "C',num2str(i-1),'"',' "D',num2str(i-1),'"'];
        end
    end
end
if ~isempty(t_opt)
    varlist=[varlist,' "t_opt"'];
end

%header lines
fprintf(fid,['TITLE        = "' frametitle '"\n']);
fprintf(fid,['VARIABLES    = ',varlist,'\n']);
fprintf(fid,'ZONE T="%s,Time=%0.6f" I=%i J=%i C=BLACK STRANDID=%i SOLUTIONTIME = %0.6f\n',frametitle,T,S(2),S(1),strand,T);
    

%write data
for i=1:S(1)
    for j=1:S(2)
        if isnan(U(i,j)) || isnan(V(i,j))
            %second check to ensure no nans present
            fprintf(fid,'%14.6e %14.6e %14.6e %14.6e %14.6e %14.6e',X(i,j),Y(i,j),0,0,-1);
        else
            %valid data points
            fprintf(fid,'%14.6e %14.6e %14.6e %14.6e %14.6e %14.6e',X(i,j),Y(i,j),U(i,j),V(i,j),Eval(i,j));
        end
        
        if ~isempty(C)
            if isnan(C(i,j,1)) || isnan(D(i,j,1))
                fprintf(fid,' %14.6e %14.6e',0,0);
            else
                fprintf(fid,' %14.6e %14.6e',C(i,j,1),D(i,j,1));
            end
            if size(U,3)>1
                for k=2:size(U,3)
                    if isnan(U(i,j,k)) || isnan(V(i,j,k))
                        fprintf(fid,' %14.6e %14.6e',0,0);
                    else
                        fprintf(fid,' %14.6e %14.6e',U(i,j,k),V(i,j,k));
                    end
                end
            end
            if size(C,3)>1
                for k=2:size(C,3)
                    if isnan(C(i,j,k)) || isnan(D(i,j,k))
                        fprintf(fid,' %14.6e %14.6e',0,0);
                    else
                        fprintf(fid,' %14.6e %14.6e',C(i,j,k),D(i,j,k));
                    end
                end
            end
        end
        
        if ~isempty(t_opt)
            if isnan(t_opt(i,j))
                fprintf(fid,' %14.6e',0);
            else
                fprintf(fid,' %14.6e',t_opt(i,j));
            end
        end
        
        fprintf(fid,'\n');
    end
end
    
fclose(fid);

end
