function [OUTPUT_DATA_STRUCTURE, POOL_IS_OPEN] = ...
    open_prana_pool_2015(INPUT_DATA_STRUCTURE)

% Copy the input data structure to the output data structure.
OUTPUT_DATA_STRUCTURE = INPUT_DATA_STRUCTURE;

% Create a data structure containing the settings
% regarding the matlab parallel pool profiles
compinfo = parcluster('local');

% Read the number of cores available to Matlab
num_cores_available = compinfo.NumWorkers;

% Read the number of cores requested
num_cores_requested = round(str2double(...
    INPUT_DATA_STRUCTURE.parprocessors));

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
OUTPUT_DATA_STRUCTURE.parprocessors = num2str(num_cores_requested);

% These lines deternube the number of image pairs that will be processed.
%
% This line reads the number of the first image to be processed.
first_image_number = str2double(INPUT_DATA_STRUCTURE.imfstart);

% This line reads the number of the last image to be processed.
last_image_number = str2double(INPUT_DATA_STRUCTURE.imfend);

% This line reads the frame step (number of images between 
% subsequent pairs)
frame_step = str2double(INPUT_DATA_STRUCTURE.imfstep);

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

    % This updates the output data structure to reflect the new number
    % of cores.
    OUTPUT_DATA_STRUCTURE.parprocessors = num2str(num_cores_requested);
end

% This line gets parameters about any open pools.
current_pool_info = gcp('nocreate');

% This checks whether a pool is already open.
POOL_IS_OPEN = ~isempty(current_pool_info);

% This checks the number of cores in an open pool
if POOL_IS_OPEN

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

            % Display warning if no pools were able to be created.
            disp(['Error Running Job in Parallel - ' ...
                'Defaulting to Single Processor\n'])

        end
    end

end

% Get info for the currently opened pool (if any)
current_pool_info = gcp('nocreate');

% Determine whether the pool is active.
POOL_IS_OPEN = ~isempty(current_pool_info);


end
        