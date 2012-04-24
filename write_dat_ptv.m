function write_dat_ptv(FILEBASE,XYZDATA,XYZNAME,VARDATA,VARNAME,TIME,TIMENAME,FORMAT,TIMEFORMAT)
%function write_plt(FILEBASE,XYZDATA,XYZNAME,VARDATA,VARNAME,TIME,TIMENAME,FORMAT,TIMEFORMAT)
% 
% write_plt('file',{X,Y},{'x','y'},{U},{'u'},T,'t','block,'volume')  
% write_plt('C:/dir/file',{X,Y,Z},{'x','y'},{U},{'u'},T,'t','point,'zone')   
%   
% FORMAT is either 'point' or 'block'.  'point' stores all variables point by
%  point, and is more readable, but 'block' groups data by variable, and is
%  more efficient to read and write.
%
% TIMEFORMAT can be either 'zone' or 'volume'.  In the case of 'zone'
%  output, each timestep, n,  will be stored in a separate zone with the
%  zone name "TIMENAME = TIME(n)".  This allows storage of either 2D or 3D
%  spatial information, with time being advanced by zone.
%
%  For 'volume', only 2D spatial information can be stored, the 3rd
%  dimension is time.  TIME is expanded to be a 3D matrix with time varying
%  with the k-index.  This is useful for plotting x-y format plots in
%  Tecplot a VARDATA vs TIME format.  In 'zone' format, Tecplot can't
%  follow a time history across multiple zones.  Using 'volume' is  
%  equivalent to 'zone' where TIME is substituted for Z in XYZDATA.
%
% TIME is always a 1D array of sequential time values.
%
% FILEBASE is the path and filename you want to use for you data.  ".dat"
%  will be appended.
%
% XYZDATA is {X, Y, Z} or {X,Y}.  Variables must be ordered so that X
%  varies with i, Y varies with j, and Z varies with k, i.e. X(i,j,k).
%
% VARDATA is {V1, V2, ..., VN} where V1, V2, etc are the data to be stored at
%  each of the XYZDATA points.  TIMEFORMAT is 'zone', each of VN must have
%  the same dimensions as X, Y, and Z.  If TIMEFORMAT is volume, each VN
%  must have i,j dimensions matching X and Y, and k dimensions matching the
%  length of TIME.
%
% XYZNAME and VARNAME are ordered lists of the names that the plt format
%  file will store for each matching variable.  They must have the same
%  number of entries as XYZDATA or VARDATA, respectively.

% v1 - 2007-11-02  - Charonko - 
%   need dimension checking and default handling, maybe some more comments?
%   Also, look into implementing variable sharing for zone output so we 
%    don't have to repeat all the values for the coordinates.  
%   Also, 'volume' may be able to be eliminated by just expanding TIME to 
%    fill Z, and creating a single 'zone' type file.
%   Another good idea is status messages or success/fail flags output for
%    file creation.
%   Change FILEBASE to require extension be included? (not assume .plt)?

% v2 - 2010-06-01 - Raben
%  Added three dimensional capabilities and if statments that allows for
%  checks before it crashes.

%write_plt(FILEBASE    ,{X(:,:,2),Y(:,:,2)},{'X','Y'},{Up          ,Vp          ,Wp          ,p_plot           ,I2,repmat(Z(:,:,2),[1,1,25])},{'U','V','W','P','Intensity','Z'},time,'T','block','volume')


% if nargin < 11
%     FORMAT = 'block';
% end

% NN = size(X);
% %keyboard
% if (size(Y)~=NN) | (size(U)~=NN) | (size(V)~=NN) | (size(W)~= NN) | (size(psi)~=NN) | (size(P)~=NN)
%     error('input fields must be the same size')
% end

NV = length(VARDATA);
NX = length(XYZDATA);
NT = length(TIME);

%NN = size(XYZDATA{1});

SX = size(XYZDATA{1});
SV = size(VARDATA{1});

% Checks to make sure that there the right amounts of names and vars
if length(XYZNAME) ~= length(XYZDATA)
    error('XYZDATA and XYZNAME are inconsistant')
end
if length(VARNAME) ~= length(VARDATA)
    error('VARDATA and VARNAME are inconsistant')
end

for pp = 1:NX
    if SX ~= size(XYZDATA{pp})
        error('XYZDATA size not Uniform')
    end
    if ~ischar(XYZNAME{pp})
        keyboard
        error('XYZNAME must be a string')
    end
    if any(ischar(XYZDATA{pp}))
        error('XYZData can not contain characters')
    end
end

for qq = 1:NV
    VARDATA{qq}(isnan(VARDATA{qq})) = 0;
    VARDATA{qq}(isinf(VARDATA{qq})) = 0;
    if SV ~= size(VARDATA{qq})
        error('VARDATA size not Uniform')
    end
    if ~ischar(VARNAME{qq})
        error('VARNAME must be a string')
    end
    if any(ischar(VARDATA{qq}))
        error('VARData can not contain characters')
    end    
end

if any(SX ~= SV(1:length(SX)))
    error('VARDATA is not the same shape as XYZDATA')
end

%This section doesn't allow for scattered three-D time dept data, i think
%maybe it should be removed
% if length(SX(SX~=1)) ~= length(XYZDATA)
%     error('XYZDATA incorrectly formated')
% end


NI = SX(1);
NJ = SX(2);
% if length(SV) == 2 && length(SX) == 2  && NT == NJ && NT > 1
%     NJ = 1;
%     NK = 1; 
%     SX(2) = [];
% %     NT = NT;
% elseif length(SV) == 2 && length(SX) == 2  && NT == 1
% %     NJ = NJ;
%     NK = 1;
% %     NT = NT;    
% else
if length(SV) == 3 && length(SX) == 2  && NT > 1
    NK = 1;
%     NT = NT;
elseif length(SV) == 3 && length(SX) == 2  && NT == 1
    NK = SX(3);
    NT = 1;
elseif length(SV) == 3 && length(SX) == 3
    NK = SX(3);
    NT = TIME;%SX(3);
elseif length(SV) == 4
    NK = SX(3);
    %NT = NT;
else
    NK = 1;
    NT = 1;
end

if ndims(VARDATA{1}) == 4
    TIMEFORMAT = 'zone';
    if length(TIME) ~= SV(4)
        error('Time Dimensions incorrect')
    end
end

filename = sprintf('%s.dat',FILEBASE);

fid = fopen(filename,'w');
%PARENTZONE=parentzone
%VARLOCATION=([varset]=varlocation, [varset]=varlocation),
%AUXDATA auxvar=�value�, VARSHARELIST=([varset]=zzz,[varset]=zzz)
%FACENEIGHBORMODE=faceneighbormode FACENEIGHBORCONNECTLIST=faceneighborconnections

% t0 = clock;
% fprintf('Writting %s...',filename)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Point
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmpi(FORMAT,'point')

    %%%%%%%%%%%%%%
    % Zone
    %%%%%%%%%%%%%%
    if strcmpi(TIMEFORMAT,'zone')
        fprintf(fid, 'TITLE="Matlab data - %s"\n',filename);                  %TITLE=�datasettitle�
        fprintf(fid,'VARIABLES =');
        for x=1:NX
            fprintf(fid,' "%s"',XYZNAME{x});
        end
        for v=1:length(VARNAME)
            fprintf(fid,' "%s"',VARNAME{v});
        end
        fprintf(fid,'\n');

        for k=1:NT
            
            if length(SX) == 1
                fprintf(fid, 'ZONE I=%g\n',NI);
            elseif length(SX) == 2
                fprintf(fid, 'ZONE I=%g J=%g\n',NI,NJ);
            else
                fprintf(fid, 'ZONE I=%g J=%g K=%g\n',NI,NJ,NK);
            end
            
            %DT=(datatypelist)C=color
            fprintf(fid,'ZONETYPE=ORDERED, DATAPACKING=%s\n',FORMAT);                           %ZONETYPE=ORDERED, DATAPACKING=datapacking,
            fprintf(fid,'T="%s=%f", SOLUTIONTIME=%f, STRANDID=1, C=BLACK\n',TIMENAME,TIME(k),TIME(k));   %T=�zonetitle�, SOLUTIONTIME=time, STRANDID=strandid
            

            for m=1:NK
                for j=1:NJ
                    for i=1:NI
                        for x=1:NX
                            fprintf(fid,' %14.7g',XYZDATA{x}(i,j,m));
                        end
                        for v=1:NV
                            if length(SX) == 2
                                fprintf(fid,' %14.7g',VARDATA{v}(i,j,k));
                            elseif length(SX) == 3
                                fprintf(fid,' %14.7g',VARDATA{v}(i,j,m,k));
                            end
                        end
                        fprintf(fid,'\n');
                    end
                end
            end
        end
        
    %%%%%%%%%%%%%%
    % Volume
    %%%%%%%%%%%%%%
    elseif strcmpi(TIMEFORMAT,'volume')
        fprintf(fid, 'TITLE="Matlab data - %s"\n',filename);                  %TITLE=�datasettitle�
        fprintf(fid,'VARIABLES =');
        for x=1:NX
            fprintf(fid,' "%s"',XYZNAME{x});
        end
        fprintf(fid,' "%s"',TIMENAME);
        for v=1:length(VARNAME)
            fprintf(fid,' "%s"',VARNAME{v});
        end
        fprintf(fid,'\n');

        fprintf(fid, 'ZONE I=%g J=%g K=%g\n',NI,NJ,NK);
        %DT=(datatypelist)C=color
        fprintf(fid,'ZONETYPE=ORDERED, DATAPACKING=%s\n',FORMAT);                           %ZONETYPE=ORDERED, DATAPACKING=datapacking,
        fprintf(fid,'T="%s=%f to %f", C=BLACK\n',TIMENAME,TIME(1),TIME(end));   %T=�zonetitle�, SOLUTIONTIME=time, STRANDID=strandid
        
        for k=1:NK
            for j=1:NJ
                for i=1:NI
                    for x=1:NX
                        fprintf(fid,' %14.7g',XYZDATA{x}(i,j));
                    end
                    fprintf(fid,' %14.7g',TIME(k));
                    for v=1:NV
                        fprintf(fid,' %14.7g',VARDATA{v}(i,j,k));
                    end
                    fprintf(fid,'\n');
                end
            end
        end
    else
        fclose(fid);
        error('unknown time organization format')
    end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Block
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
elseif strcmpi(FORMAT,'block')

    %%%%%%%%%%%%%%
    % Zone
    %%%%%%%%%%%%%%
    if strcmpi(TIMEFORMAT,'zone')
        fprintf(fid, 'TITLE="Matlab data - %s"\n',filename);                  %TITLE=�datasettitle�
        fprintf(fid,'VARIABLES =');
        for x=1:NX
            fprintf(fid,' "%s"',XYZNAME{x});
        end
        for v=1:length(VARNAME)
            fprintf(fid,' "%s"',VARNAME{v});
        end
        fprintf(fid,'\n');

        for m=1:NT
            
            if length(SX) == 1
                fprintf(fid, 'ZONE I=%g\n',NI);
            elseif length(SX) == 2
                fprintf(fid, 'ZONE I=%g J=%g\n',NI,NJ);
            else
                fprintf(fid, 'ZONE I=%g J=%g K=%g\n',NI,NJ,NK);
            end

            %DT=(datatypelist)C=color
            fprintf(fid,'ZONETYPE=ORDERED, DATAPACKING=%s\n',FORMAT);                     	%ZONETYPE=ORDERED, DATAPACKING=datapacking,
            fprintf(fid,'T="%s=%f", SOLUTIONTIME=%f, STRANDID=1, C=BLACK\n',TIMENAME,TIME(m),TIME(m));	%T=�zonetitle�, SOLUTIONTIME=time, STRANDID=strandid

            for x=1:NX
                fprintf(fid,'%12.5e %12.5e %12.5e %12.5e %12.5e %12.5e %12.5e\n',XYZDATA{x});
            end
            for v=1:NV
                if length(SX) == 2
                fprintf(fid,'%12.5e %12.5e %12.5e %12.5e %12.5e %12.5e %12.5e\n',VARDATA{v}(:,:,m));
                elseif length(SX) == 3
                    %This Needs Work but is heading in the right direction
                    fprintf(fid,'%12.5e %12.5e %12.5e %12.5e %12.5e %12.5e %12.5e\n',VARDATA{v}(:,:,:,m));                    
                end
            end

        end

    %%%%%%%%%%%%%%
    % Volume
    %%%%%%%%%%%%%%    
    elseif strcmpi(TIMEFORMAT,'volume')
        fprintf(fid, 'TITLE="Matlab data - %s"\n',filename);                  %TITLE=�datasettitle�
        fprintf(fid,'VARIABLES =');
        for x=1:NX
            fprintf(fid,' "%s"',XYZNAME{x});
        end
        fprintf(fid,' "%s"',TIMENAME);
        for v=1:length(VARNAME)
            fprintf(fid,' "%s"',VARNAME{v});
        end
        fprintf(fid,'\n');

        fprintf(fid, 'ZONE I=%g J=%g K=%g\n',NI,NJ,NK);
        %DT=(datatypelist)C=color
        fprintf(fid,'ZONETYPE=ORDERED, DATAPACKING=%s\n',FORMAT);                               %ZONETYPE=ORDERED, DATAPACKING=datapacking,
        fprintf(fid,'T="%s=%f to %f", C=BLACK\n',TIMENAME,TIME(1),TIME(end));	%T=�zonetitle�, SOLUTIONTIME=time, STRANDID=strandid

        for x=1:NX
            fprintf(fid,'%12.5e %12.5e %12.5e %12.5e %12.5e %12.5e %12.5e\n',repmat(XYZDATA{x},[1,1,NK]));
        end
        fprintf(fid,'%12.5e %12.5e %12.5e %12.5e %12.5e %12.5e %12.5e\n',repmat(permute(TIME(:),[3,2,1]),[NI,NJ,1]));
        for v=1:NV
            fprintf(fid,'%12.5e %12.5e %12.5e %12.5e %12.5e %12.5e %12.5e\n',VARDATA{v});
        end
        
    else
        fclose(fid);
        error('unknown time organization format')
    end
else
    fclose(fid);
    error('unknown tecplot datapacking format')
end

fclose(fid);

% eltime = etime(clock,t0);
% fprintf('Done\t')
% fprintf('%0.2i:%0.2i.%0.0f\n',floor(eltime/60),floor(rem(eltime,60)),rem(eltime,60)-floor(rem(eltime,60)))
