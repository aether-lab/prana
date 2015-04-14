function varargout=pranaPIVcode(Data)

if ischar(Data)
    if strcmpi(Data,'version')
        varargout{1} = '2.5.beta.r2014.04.14';
    else
        error('Error: String request must be ''version''.\m')
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
    
    % Determine whether a parallel job was specified
    run_parallel = str2double(Data.par);
    
    % Determine the current version of Matlab
    matlab_version_string = version('-release');
    
    % Convert the matlab version string to a number indicating
    % the release year
    matlab_version_year = str2double(matlab_version_string(1:4));
    
    %% Set up a parallel job if needed
    if run_parallel

        % Inform the user that a parallel pool is being initialized.
        fprintf('\n--- Initializing Processor Cores for Parallel Job ----\n')
        
        % Attempt to open a matlab pool
        %
        % Choose between matlabpool and parpool paradigms based on the 
        % release year of the currently instance of Matlab.
        if matlab_version_year < 2013
            
        else
            [Data, pool_opened] = open_prana_pool_2015(Data);
        end
        
        % If a pool wasn't successfully opened,
        % run as single core.
        if ~pool_opened
                    
            % Inform the user that processing has begun
            fprintf(['\n-------------- Processing Single-Core Dataset'...
                '(started at %s) ------------------\n'], datestr(now));

            % Run single-core processing.
            pranaprocessing(Data)
            
        end
        
        % Proceed with parallel processing if a parallel pool is open.
        if pool_opened
            I1=str2double(Data.imfstart):str2double(Data.imfstep):str2double(Data.imfend);
            I2=I1+str2double(Data.imcstep);
            if strcmp(Data.masktype,'dynamic')
                maskfend=str2double(Data.maskfstart)+str2double(Data.maskfstep)*length(str2double(Data.imfstart):str2double(Data.imfstep):str2double(Data.imfend))-1;
                maskname=str2double(Data.maskfstart):str2double(Data.maskfstep):maskfend;
            else
                maskname=nan(1,length(I1));
            end
            
            if any(str2double(Data.method)==[4 5])
                fprintf('\n-------------- Processing Parallel-Core Dataset (started at %s) ------------------\n', datestr(now));
                pranaprocessing(Data)
                fprintf('---------------- Job Completed at %s ---------------------\n', datestr(now));
            else
                fprintf('\n-------------- Processing Parallel-Core Dataset (started at %s) ------------------\n', datestr(now));
                spmd
                    
                    if matlab_version_year >= 2010
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
                            
                            % This closes the open matlab pool.
                            delete(gcp);
                            
                            % This runs the Prana job.
                            pranaprocessing(Data)
                        end
                        
                    end
                    
                    % Run Prana on parallel cores.
                    pranaprocessing(Data,I1dist,I2dist,masknamedist);
                end
                fprintf('----------------- Job Completed at %s----------------------\n', datestr(now));
            end
            
        end
    else
        
        % Inform the user that the processing is beginning.
        fprintf(['\n-------------- Processing Single-Core Dataset (started at %s)--'...
            '----------------\n'], datestr(now));
        pranaprocessing(Data)
        
        % Inform the user that processing has finished.
        fprintf(['---------------- Single-Core Job Completed at %s----------'...
            '-----------\n'], datestr(now))
        
    end
end


end
