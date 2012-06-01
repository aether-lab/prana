function varargout=pranaPIVcode(Data)

if ischar(Data)
    if strcmpi(Data,'version')
        varargout{1} = '2.0.beta.r2012.05.30';
    else
        error('Error: String request must be ''version''.')
    end
else
    
    % Make sure the job file has all the required variables.  This will make
    % sure that jobs created from older versions have the necessary variables.
    [Data] = jobfile_validator(Data);
    
    % Write experimental summary for the given job.  This was previuosly in the
    % GUI but has been moved here so that it is always run even if prana is run
    % via a script.
    % This call also puts the out put on the screen for the user to evaluate.
    write_expsummary(Data);
    
    %% Set up a parallel job if needed
    if str2double(Data.par)
        fprintf('\n--- Initializing Processor Cores for Parallel Job ----\n')
        poolopen=1;
        
        % This section of code checks the job file for the number of
        % parallel processors requested for the job.  If the job requests
        % more then what is available in the machine it reduces the number
        % requested to 1 less then the machine maximum.  The reason for '1
        % less' is so that machine is not unintensionally over worked.  If
        % the users requests all of the machines processors this code will
        % not stop the user.  Only if a request is made OVER the avialable
        % number.
        compinfo=findResource('scheduler','configuration','local');
        avail_parprocs=num2str(compinfo.clustersize);
        if str2double(Data.parprocessors) > avail_parprocs
            fprintf('Job file requested %s processors while the machine only contains %s processors\n Updating job file to request %s processors',Data.parprocessors,num2str(avail_parprocs),num2str(avail_parprocs-1))
            Data.parprocessors = num2str(avail_parprocs-1);
        end
        
        %Don't open more processors than there are image pairs
        if length(str2double(Data.imfstart):str2double(Data.imfstep):str2double(Data.imfend)) < str2double(Data.parprocessors)
            Data.parprocessors=num2str(length(str2double(Data.imfstart):str2double(Data.imfstep):str2double(Data.imfend)));
        end
        
        try
            matlabpool('open','local',Data.parprocessors);
        catch
            try
                matlabpool close
                matlabpool('open','local',Data.parprocessors);
            catch
                beep
                disp('Error Running Job in Parallel - Defaulting to Single Processor')
                poolopen=0;
                fprintf('\n-------------- Processing Dataset (started at %s) ------------------\n', datestr(now));
                pranaprocessing(Data)
                fprintf('---------------- Job Completed at %s ---------------------\n', datestr(now));
            end
        end
        if poolopen
            I1=str2double(Data.imfstart):str2double(Data.imfstep):str2double(Data.imfend);
            I2=I1+str2double(Data.imcstep);
            if strcmp(Data.masktype,'dynamic')
                maskfend=str2double(Data.maskfstart)+str2double(Data.maskfstep)*length(str2double(Data.imfstart):str2double(Data.imfstep):str2double(Data.imfend))-1;
                maskname=str2double(Data.maskfstart):str2double(Data.maskfstep):maskfend;
            else
                maskname=nan(1,length(I1));
            end
            
            if any(str2double(Data.method)==[4 5])
                fprintf('\n-------------- Processing Dataset (started at %s) ------------------\n', datestr(now));
                pranaprocessing(Data)
                fprintf('---------------- Job Completed at %s ---------------------\n', datestr(now));
            else
                fprintf('\n-------------- Processing Dataset (started at %s) ------------------\n', datestr(now));
                spmd
                    verstr=version('-release');
                    if str2double(verstr(1:4))>=2010
                        I1dist=getLocalPart(codistributed(I1,codistributor('1d',2)));
                        I2dist=getLocalPart(codistributed(I2,codistributor('1d',2)));
                        masknamedist=getLocalPart(codistributed(maskname,codistributor('1d',2)));
                    else
                        I1dist=localPart(codistributed(I1,codistributor('1d',2),'convert'));
                        I2dist=localPart(codistributed(I2,codistributor('1d',2),'convert'));
                        masknamedist=localPart(codistributed(maskname,codistributor('1d',2),'convert'));
                    end
                    
                    if str2double(Data.method)==6
                        try
                            if labindex~=1
                                previous = labindex-1;
                            else
                                previous = numlabs;
                            end
                            if labindex~=numlabs
                                next = labindex+1;
                            else
                                next = 1;
                            end
                            
                            I1extra_end=labSendReceive(previous,next,I1dist(1:str2double(Data.framestep)));
                            I2extra_end=labSendReceive(previous,next,I2dist(1:str2double(Data.framestep)));
                            masknameextra_end=labSendReceive(previous,next,masknamedist(1:str2double(Data.framestep)));
                            
                            I1extra_beg=labSendReceive(next,previous,I1dist((end-str2double(Data.framestep)+1):end));
                            I2extra_beg=labSendReceive(next,previous,I2dist((end-str2double(Data.framestep)+1):end));
                            masknameextra_beg=labSendReceive(next,previous,masknamedist((end-str2double(Data.framestep)+1):end));
                            
                            if labindex<numlabs
                                I1dist = [I1dist,I1extra_end];
                                I2dist = [I2dist,I2extra_end];
                                masknamedist = [masknamedist,masknameextra_end];
                            end
                            if 1<labindex
                                I1dist = [I1extra_beg,I1dist];
                                I2dist = [I2extra_beg,I2dist];
                                masknamedist = [masknameextra_beg,masknamedist];
                            end
                        catch
                            beep
                            disp('Error Running Multiframe Job in Parallel (Not Enough Image Pairs) - Defaulting to Single Processor')
                            matlabpool close
                            poolopen=0;
                            pranaprocessing(Data)
                        end
                        
                    end
                    
                    pranaprocessing(Data,I1dist,I2dist,masknamedist);
                end
                fprintf('----------------- Job Completed at %s----------------------\n', datestr(now));
            end
            if poolopen
                matlabpool close
            end
        end
    else
        fprintf('\n-------------- Processing Dataset (started at %s)------------------\n', datestr(now));
        pranaprocessing(Data)
        fprintf('---------------- Job Completed at %s---------------------\n', datestr(now))
    end
end
end
