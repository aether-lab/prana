function [x,y,u,v,extravars,varlistnew]=read_pltmod_NC(testname,direc,startframe,endframe,numzeros,framestep,varlist)
%
%READ_PLT Read a Tecplot plt created in FlowIQ.
%   [X,Y,U,V]=READ_PLT(TESTNAME,DIREC,STARTFRAME,ENDFRAME,NUMZEROS,FRAMESTEP)
%   reads the series of files testnameXXXX.plt from directory DIREC, where
%   XXXX is some integer beginning with STARTFRAME and ending with ENDFRAME 
%   in increments of FRAMESTEP.  XXXX has a digit length equal to NUMZEROS.  
%
%   Velocity component data are returned in U and V matrices of unknown 
%   dimension MxNxT, where M is the size in the X direction, N is the size
%   in the Y direction, and T is the total number of frames read in. 
%
%   X and Y are coordinate matrices locating the positions of the vectors
%   stored in U and V.  They have dimension MxN.
%
%   NUMZEROS and FRAMESTEP are optional, and have default values of 4 and
%   1, respectively.

%%%
%   long term plan, allow variable output argument length
%   if length(vargout)=1
%   return contents in data.(variablelist(1))=x,etc
%   maybe allow for selecting fields from input list {'X','U','Correlation'}?
%   X,Y,U,V would be default, of course?
%%%

%allow for numzeros to default to 4, and framestep to 1
if nargin<7
    varlist = {'all'};
    if nargin<6
        framestep=1;
        if nargin==4
            numzeros=4;
        end
    end
end

extravars={};


try
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
numframes=endframe-startframe+1;

nameformat = sprintf('%%s%%0.%ui.dat',numzeros);
picfilename=sprintf(nameformat,testname,startframe);
fid = fopen(fullfile(direc,picfilename));
if fid == -1
    nameformat = sprintf('%%s%%0.%ui.plt',numzeros);
    picfilename=sprintf(nameformat,testname,startframe);
    fid = fopen(fullfile(direc,picfilename));
end


temp=fgetl(fid);%#ok                %TITLE="DPIV Data File"
temp=fgetl(fid);                    %VARIABLES="X" "Y" ....

%figure out what variables we have, and what their names are
[A,count,errmsg,nextindex]=sscanf(temp,'%*[^"]',1);
currindex = nextindex;
columns      = 0;
variablelist = cell(1);
var2         = cell(1);

while currindex<=length(temp)
    columns=columns+1;
    [variablelist{columns},count,errmsg,nextindex] = sscanf(temp(currindex:end),'%*["]%[^"]%*["]',1);    
    if isempty(variablelist{columns})
        currindex = currindex+1;
    else
        currindex = currindex+nextindex+1;
    end
end
r = 1;
for mm = 1:length(variablelist)
    if ~isempty(variablelist{mm})
        var2{r} = variablelist{mm};
        r = r+1;
    end
end
variablelist = var2;
columns = length(variablelist);

hh=0;varindexlist=[];varlistnew={};
if strcmp(varlist{1},'all')
    varlist=variablelist(5:end);
end
for i=1:length(varlist)
    for j=1:length(variablelist)
        if strcmp(varlist{i},variablelist{j})
            hh=hh+1;
            varindexlist(hh) = j;%#ok
            varlistnew{hh}=variablelist{j};%#ok     % in case a variable isnt typed in right
            break;
        end
    end
end

% if (nargout-4)~=length(varlist)
%     error('Number of outputs must match length of variable list')
% end

%determine size of arrays

%temp=fscanf(fid,'%*7c%u%*3c%u',2);   %Zone I=XXX J=XXX
%temp=fscanf(fid,'%*25c%u%*3c%u',2);   %Zone I=XXX J=XXX

tempstr=fgetl(fid);
 for rr=1:length(tempstr)-1
     if strcmp(tempstr(rr:rr+1),'I=')
         Iflag=rr;
     end
     if strcmp(tempstr(rr:rr+1),'J=')
         Jflag=rr;
     end
 end
Imax=sscanf(tempstr(Iflag+2:Jflag-1),'%u');
Jmax=sscanf(tempstr(Jflag+2:end),'%u');
        
fclose(fid);

%pre-allocate space for speed
u=zeros(Jmax,Imax,ceil(numframes/framestep));
v=zeros(Jmax,Imax,ceil(numframes/framestep));

numextravars = length(varindexlist);
for i = 1:numextravars
    vardata{i} = zeros(Jmax,Imax,ceil(numframes/framestep));%#ok
end

count=0;

for j=startframe:framestep:endframe
    count=count+1;

    picfilename=sprintf(nameformat,testname,j);
    fid = fopen(fullfile(direc,picfilename));
    
    temp=fgetl(fid);%#ok                    %TITLE="DPIV Data File"
    temp=fgetl(fid);%#ok                    %VARIABLES="X" "Y" ....
    temp=fgetl(fid);%#ok
%     temp=fscanf(fid,'%*7c%u%*3c%u',2);   %Zone I=XXX J=XXX
%     Imax=temp(1);
%     Jmax=temp(2);
    
    %read data for this frame
    variable=(fscanf(fid,'%g',[columns inf]))';
    numrows=size(variable,1);%#ok
    variabletemp=zeros(Jmax*Imax,columns);%#ok
%     if numrows~=Imax*Jmax               % if masking is used and not all points are printed in loaded .plt
%         tmpx=abs(diff(variable(:,1)));
%         tmpy=abs(diff(variable(:,2)));
%         indx=any(tmpx,2);
%         indy=any(tmpy,2);
%         dx=min(tmpx(indx));
%         dy=min(tmpy(indy));
%         rowindex=round((variable(:,2)-min(variable(:,2)))/dy+1);
%         colindex=round((variable(:,1)-min(variable(:,1)))/dx+1);
%         variabletemp((rowindex-1)*Jmax+colindex,:)=variable;
%         variable=variabletemp;
%     end
    
%     disp(['frame ' num2str(j) ' loaded']);
    fclose(fid);

        u(:,:,count) = reshape(variable(:,3),Imax,Jmax)';
        v(:,:,count) = reshape(variable(:,4),Imax,Jmax)';
        for i=1:numextravars
            vardata{i}(:,:,count) = reshape(variable(:,varindexlist(i)),Imax,Jmax)';%#ok
        end
    
end

    x(:,:) = reshape(variable(:,1),Imax,Jmax)';
    y(:,:) = reshape(variable(:,2),Imax,Jmax)';
    for i=1:numextravars
        extravars{i} = vardata{i};%#ok
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
catch %#ok
    if fid==-1
        error('%s does not exist',[direc picfilename])
    else
        rethrow(lasterror)%#ok
    end
end

end
