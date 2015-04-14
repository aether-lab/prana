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
    
    %% Set up a parallel job if needed
    if run_parallel

        % Inform the user that a parallel pool is being initialized.
        fprintf('\n--- Initializing Processor Cores for Parallel Job ----\n')
        
        % Create a data structure containing the settings
        % regarding the matlab parallel pool profiles
        compinfo = parcluster('local');

        % Read the number of cores available to Matlab
        num_cores_available = compinfo.NumWorkers;
        
        % Read the number of cores requested
        num_cores_requested = round(str2double(Data.parprocessors));
        
        % If the number of cores requested exceeds the number available,
        % set the number requested to be one fewer than the number
        % available.
        if num_cores_requested > num_cores_available
            
            % Inform the user that the number of requested cores
            % exceeds the number of cores available to Matlab.
            fprintf(['Job file requested %d processors while the'...
                ' machine only contains %d processors\n Updating job'...
                 'file to request %d processors.\n'], ...
                 num_cores_requested, num_cores_available, ...
                 num_cores_available - 1);
            
            
            % Set the number of requested cores to one fewer
            % than the number available to Matlab         
            num_cores_requested = num_cores_available - 1;
                        
        end
        
        % Update the data structure with the new number
        % of requested cores.
        %
        % Doing this outside of the above if statement takes
        % care of non-integer numbers of requested cores.
        Data.parprocessors = num2str(num_cores_requested);
        
        % These lines deternube the number of image pairs that will be processed.
        %
        % This line reads the number of the first image to be processed.
        first_image_number = str2double(Data.imfstart);
        
        % This line reads the number of the last image to be processed.
        last_image_number = str2double(Data.imfend);
        
        % This line reads the frame step (number of images between 
        % subsequent pairs)
        frame_step = str2double(Data.imfstep);
        
        % This line creates a list of the numbers of the images to be
        % processed. These are the numbers of the first image in each pair.
        image_number_list = ...
            first_image_number : frame_step : last_image_number;
        
        % This line determines the number of image pairs that will be
        % processed.
        number_of_pairs = length(image_number_list);
        
        % This if-statement sets the number of cores requested to 
        % not exceed the number of image pairs to be processed.
        if number_of_pairs < num_cores_requested
            
            % This updates the "number of pairs" variable to equal
            % the number of frames.
            num_cores_requested = number_of_pairs;
            
            % This updates the data structure to reflect the new number
            % of cores.
            Data.parprocessors = num2str(num_cores_requested);
        end
        
        % This line gets parameters about any open pools.
        current_pool_info = gcp('nocreate');
        
        % This checks whether a pool is already open.
        pool_is_open = ~isempty(current_pool_info);
        
        % This checks the number of cores in an open pool
        if pool_is_open
            
            % This checks the number of cores open in the current pool
            num_open_cores = current_pool_info.NumWorkers;
            
        else
            % This sets the number of opened cores to zero.
            num_open_cores = 0;
        end
        
        
        % This if-statement opens a matlab pool only if either one doesn't
        % already exist or if a pool of the wrong size exists.
        if num_open_cores ~= num_cores_requested 
        
            % This try-catch statement attempts to open a parallel pool
            % with the number of processors requested.
            try

                % This opens a parallel pool with the specified number
                % of cores           
                parpool('local', num_cores_requested);

            catch

                % If the pool wasn't successfully opened, then set the number
                % of requested cores to zero and run the job.

                try

                    % This closes any existing matlab pools.
                    % gcp is a matlab command for "get current pool"
                   delete(gcp);

                    % This opens a parallel pool with the requested
                    % number of processors.
                    parpool('local', num_cores_requested);

                catch
                    beep
                    disp('Error Running Job in Parallel - Defaulting to Single Processor\n')

                    % This sets the boolean flag that indicates whether
                    % or not a parallel pool is open to false.
                    pool_is_open = 0;

                    % Inform the user that processing has begun
                    fprintf('\n-------------- Processing Dataset (started at %s) ------------------\n', datestr(now));

                    % Inform the user that processing has finished.
                    pranaprocessing(Data)
                    fprintf(['---------------- Job Completed at %s -------'...
                        '--------------\n'], datestr(now));
                end
            end

        end
        
        % Get info for the currently opened pool (if any)
        current_pool_info = gcp('nocreate');
        
        % Determine whether the pool is active.
        pool_is_open = ~isempty(current_pool_info);
        
        % Proceed with parallel processing if a parallel pool is open.
        if pool_is_open
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
                            
                            % This closes the open matlab pool.
                            delete(gcp);
                            
                            % This runs the Prana job.
                            pranaprocessing(Data)
                        end
                        
                    end
                    
                    pranaprocessing(Data,I1dist,I2dist,masknamedist);
                end
                fprintf('----------------- Job Completed at %s----------------------\n', datestr(now));
            end
            
        end
    else
        fprintf('\n-------------- Processing Dataset (started at %s)------------------\n', datestr(now));
        pranaprocessing(Data)
        fprintf('---------------- Job Completed at %s---------------------\n', datestr(now))
    end
end
end
