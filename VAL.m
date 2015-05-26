function [Uval,Vval,Evalval,Cval,Dval,DXval,DYval,ALPHAval]=VAL(X,Y,U,V,Eval,C,D,Valoptions,extrapeaks,DX,DY,ALPHA)
% --- Validation Subfunction ---

%     This file is part of prana, an open-source GUI-driven program for
%     calculating velocity fields using PIV or PTV.
%
%     Copyright (C) 2012-2014  Virginia Polytechnic Institute and State
%     University
% 
%     Copyright 2014-2015.  Los Alamos National Security, LLC. This material was
%     produced under U.S. Government contract DE-AC52-06NA25396 for Los 
%     Alamos National Laboratory (LANL), which is operated by Los Alamos 
%     National Security, LLC for the U.S. Department of Energy. The U.S. 
%     Government has rights to use, reproduce, and distribute this software.
%     NEITHER THE GOVERNMENT NOR LOS ALAMOS NATIONAL SECURITY, LLC MAKES ANY
%     WARRANTY, EXPRESS OR IMPLIED, OR ASSUMES ANY LIABILITY FOR THE USE OF
%     THIS SOFTWARE.  If software is modified to produce derivative works,
%     such modified software should be clearly marked, so as not to confuse
%     it with the version available from LANL.
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

%unpack the switches and options for the different validation types
Threshswitch    = Valoptions.Threshswitch;
Uthresh         = Valoptions.Uthresh;
Vthresh         = Valoptions.Vthresh;

UODswitch       = Valoptions.UODswitch;
UODwinsize      = Valoptions.UODwinsize';  %comes in arranged as [x/y, numpasses], UOD needs [numpasses, x/y]
UODthresh       = Valoptions.UODthresh';

Bootswitch      = Valoptions.Bootswitch;
Bootper         = Valoptions.Bootper;
Bootiter        = Valoptions.Bootiter;
Bootkmax        = Valoptions.Bootkmax;

Corrpeakswitch      = Valoptions.Corrpeakswitch;
Peakheight_thresh   = Valoptions.Peakheight_thresh;
Peakratio_thresh    = Valoptions.Peakratio_thresh;


imClass = 'double';

if exist('DX','var')
    [~,~,~,~,DX,DY,ALPHA]=matrixform(X,Y,U,V,DX,DY,ALPHA);
    [X,Y,U,V,Eval,C,D]=matrixform(X,Y,U,V,Eval,C,D);
else
    [X,Y,U,V,Eval,C,D]=matrixform(X,Y,U,V,Eval,C,D);
    DX    = zeros(size(D),imClass);
    DY    = zeros(size(D),imClass);
    ALPHA = zeros(size(D),imClass);
end

if extrapeaks
    numpeaks=3;
else
    numpeaks=1;
end

if Corrpeakswitch
    %create a peak ratio array, with values for the first two peaks, and
    %ones for the third peak (peak 3 will always fail ratio validation)
    PR = cat(3,C(:,:,1:2)./C(:,:,2:3),ones(size(C(:,:,1))));
else
    PR = ones(size(C(:,:,1)));
end  

Uval=U(:,:,1);Vval=V(:,:,1);Evalval=Eval(:,:,1);
if ~isempty(C)
    Cval  = C(:,:,1);
    Dval  = D(:,:,1);
    PRval = PR(:,:,1);
else
    Cval  = [];
    Dval  = [];
    PRval = [];
end

    DXval    = DX(:,:,1);
    DYval    = DY(:,:,1);
    ALPHAval = ALPHA(:,:,1);
    
    
    
S=size(X);

if Threshswitch || UODswitch
    for i=1:numpeaks
        %Thresholding
        if Threshswitch
            [Uval,Vval,Evalval] = Thresh(Uval,Vval,Uthresh,Vthresh,Evalval);
        end
        
        %Correlation peak thresholds
        if Corrpeakswitch && ~isempty(C)
            [Uval,Vval,Evalval] = CorrpeakThresh(Uval,Vval,Cval,PRval,Peakheight_thresh,Peakratio_thresh,Evalval);
        end

        %Univeral Outlier Detection
        if UODswitch
            %t=permute(UODwinsize,[2 3 1]);
            [Uval,Vval,Evalval] = UOD(Uval,Vval,UODwinsize,UODthresh,Evalval);
        end
%         disp([num2str(sum(sum(Evalval>0))),' bad vectors'])
        %Try additional peaks where validation failed
        if i<numpeaks %if there are extra peaks we haven't tested...
            %copy next peak data into temp arrays
            Utemp=U(:,:,i+1);Vtemp=V(:,:,i+1);Evaltemp=Eval(:,:,i+1);Ctemp=C(:,:,i+1);Dtemp=D(:,:,i+1);
            DXtemp=DX(:,:,i+1);DYtemp=DY(:,:,i+1);ALPHAtemp=ALPHA(:,:,i+1);
            PRtemp=PR(:,:,i+1);
            %copy next peak data onto sites that failed validation
            Uval(Evalval>0)=Utemp(Evalval>0);
            Vval(Evalval>0)=Vtemp(Evalval>0);
            Cval(Evalval>0)=Ctemp(Evalval>0);
            PRval(Evalval>0)=PRtemp(Evalval>0);
            Dval(Evalval>0)=Dtemp(Evalval>0); 
            Evalval(Evalval>0)=Evaltemp(Evalval>0);
            
            DXval(Evalval>0)=DXtemp(Evalval>0); 
            DYval(Evalval>0)=DYtemp(Evalval>0); 
            ALPHAval(Evalval>0)=ALPHAtemp(Evalval>0); 
        end
    end
end

%limit replacement search radius to largest region tested during UOD
maxSearch = floor( (max(UODwinsize(:))-1)/2 );

%replacement
for i=1:S(1)
    for numpeaks=1:S(2)
        if Evalval(i,numpeaks)>0
            %initialize replacement search size
            q=0;
            s=0;

            %get replacement block with at least 8 valid points, but stop if region grows bigger than UOD test region
            while s==0
                q=q+1;
                Imin = max([i-q 1   ]);
                Imax = min([i+q S(1)]);
                Jmin = max([numpeaks-q 1   ]);
                Jmax = min([numpeaks+q S(2)]);
                Iind = Imin:Imax;
                Jind = Jmin:Jmax;
                Ublock = Uval(Iind,Jind);
                if q >= maxSearch || length(Ublock(~isnan(Ublock)))>=8 
                    Xblock = X(Iind,Jind)-X(i,numpeaks);
                    Yblock = Y(Iind,Jind)-Y(i,numpeaks);
                    Vblock = Vval(Iind,Jind);
                    s=1;
                end
            end
            
            %distance from erroneous vector
            Dblock = (Xblock.^2+Yblock.^2).^-0.5;
            Dblock(isnan(Ublock))=nan;
            Dblock(isinf(Dblock))=nan;

            %validated vector
            Uval(i,numpeaks) = nansum(nansum(Dblock.*Ublock))/nansum(nansum(Dblock));
            Vval(i,numpeaks) = nansum(nansum(Dblock.*Vblock))/nansum(nansum(Dblock));       
        end
    end
end

%clean up any remaining NaNs
Uval(isnan(Uval)) = 0;
Vval(isnan(Vval)) = 0;

%Bootstrapping
if Bootswitch
    [Uval,Vval,Evalval] = bootstrapping(X,Y,Uval,Vval,Bootper,Bootiter,Bootkmax,Evalval);
end

%convert back to vector
[Uval,Vval,Evalval,Cval,Dval]=vectorform(X,Y,Uval,Vval,Evalval,Cval,Dval);
[DXval,DYval,ALPHAval]=vectorform(X,Y,DXval,DYval,ALPHAval,[],[]);

end
