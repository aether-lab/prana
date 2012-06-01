function varargout = pranaPTVcode(PTV_Data)

if ischar(PTV_Data)
    if strcmpi(PTV_Data,'version')
        varargout{1} = '2.0.beta.r2012.05.30';
    else
        error('Error: String request must be ''version''.')
    end
else
    
    % Make sure the job file has all the required variables.  This will make
    % sure that jobs created from older versions have the necessary variables.
    [PTV_Data] = jobfile_validator(PTV_Data);
    
    if ispc
        Data.slsh='\';
    else
        Data.slsh='/';
    end
    
    % --- Images ---
    Data.imdirec  = PTV_Data.imdirec;
    Data.imbase   = PTV_Data.imbase;
    Data.imext    = PTV_Data.imext;
    Data.imzeros  = str2double(PTV_Data.imzeros);
    Data.imcstep  = str2double(PTV_Data.imcstep);
    Data.imfstep  = str2double(PTV_Data.imfstep);
    Data.imfstart = str2double(PTV_Data.imfstart);
    Data.imfend   = str2double(PTV_Data.imfend);
    Data.channel  = str2double(PTV_Data.channel);
    if Data.channel==6
        fprintf(['Ensemble Color does not currently work with tracking code.\n'...
            'Resetting the color channel to 1 (grey scale)\n'])
        Data.channel=1;
    end
    
    % --- ID ---
    IDmethod = {'blob','dynamic','combined'};
    Data.ID.method         = IDmethod{str2double(PTV_Data.ID.method)};
    Data.ID.run            = str2double(PTV_Data.ID.runid);
    Data.ID.v              = str2double(PTV_Data.ID.imthresh);
    Data.ID.contrast_ratio = 0;%str2double(Data.ID.contrast_ratio);
    Data.ID.s_num          = 0;%str2double(Data.ID.s_num);
    Data.ID.s_name         = PTV_Data.ID.savebase;
    Data.ID.save_dir       = PTV_Data.ID.save_dir;
    
    % --- Sizing ---
    Data.Size.run      = str2double(PTV_Data.Size.runsize);
    Data.Size.thresh   = str2double(PTV_Data.ID.imthresh);
    Data.Size.method   = PTV_Data.Size.method;
    Data.Size.p_area   = str2double(PTV_Data.Size.min_area);
    Data.Size.sigma    = str2double(PTV_Data.Size.std);
    Data.Size.errors   = 1;%str2double(Data.Size.errors);
    Data.Size.s_name   = PTV_Data.Size.savebase;
    Data.Size.save_dir = PTV_Data.Size.save_dir;
    
    % --- PIV Info ---
    if str2double(PTV_Data.datout)
        Data.Track.PIVprops.extension='.dat';
    elseif str2double(PTV_Data.multiplematout)
        Data.Track.PIVprops.extension='.mat';
    end
    Data.Track.PIVprops.load_dir = PTV_Data.outdirec;
    Data.PIV.Data                = PTV_Data;
    eval(['Data.PIV.Data.outbase = PTV_Data.PIV' PTV_Data.passes '.outbase;'])
    
    % --- Tracking ---
    PTVmethod                 = {'none','ptv','piv','piv-ptv'};
    PTVpredict                = {'static','dynamic'};
    Data.Track.run            = str2double(PTV_Data.Track.runtrack);
    Data.Track.method         = PTVmethod{str2double(PTV_Data.Track.method)};
    Data.Track.predict_mode   = PTVpredict{str2double(PTV_Data.Track.prediction)};
    Data.Track.PIV_PTV_weight = str2double(PTV_Data.Track.PIVweight);
    Data.Track.plotfig        = 0;%str2double(Data.Track.plotfig);
    Data.Track.s_radius       = str2double(PTV_Data.Track.radius);
    Data.Track.r_weight       = str2double(PTV_Data.Track.estradius);
    Data.Track.edgeval        = str2double(PTV_Data.Track.estweight);
    Data.Track.numvecs        = str2double(PTV_Data.Track.vectors);
    Data.Track.max_iterations = str2double(PTV_Data.Track.iterations);
    Data.Track.s_name         = PTV_Data.Track.savebase;
    Data.Track.save_dir       = PTV_Data.Track.save_dir;
    Data.Track.savemat        = str2double(PTV_Data.Track.trackmat);
    Data.Track.savedat        = str2double(PTV_Data.Track.trackdat);
    
    cutoff_commas = strfind(PTV_Data.Track.valprops.valcoef,',');
    radius_commas = strfind(PTV_Data.Track.valprops.valrad,',');
    MAD_U_commas  = strfind(PTV_Data.Track.valprops.MAD_U,',');
    MAD_V_commas  = strfind(PTV_Data.Track.valprops.MAD_V,',');
    
    cutoff = zeros(numel(cutoff_commas)+1,1);
    radius = zeros(numel(radius_commas)+1,1);
    MAD_U  = zeros(numel(MAD_U_commas)+1,1);
    MAD_V  = zeros(numel(MAD_V_commas)+1,1);
    
    Data.Track.valprops.numpass = numel(cutoff_commas)+1;
    
    cutoff(1)=str2double(PTV_Data.Track.valprops.valcoef(1:cutoff_commas(1)-1));
    radius(1)=str2double(PTV_Data.Track.valprops.valrad(1:radius_commas(1)-1));
    MAD_U(1)=str2double(PTV_Data.Track.valprops.MAD_U(1:MAD_U_commas(1)-1));
    MAD_V(1)=str2double(PTV_Data.Track.valprops.MAD_V(1:MAD_V_commas(1)-1));
    method=cell(Data.Track.valprops.numpass,1);
    if radius(1)~=0
        method{1} = 'median';
    else
        method{1} = 'coeff';
    end
    for i=2:Data.Track.valprops.numpass
        if i ~= Data.Track.valprops.numpass
            cutoff(i)=str2double(PTV_Data.Track.valprops.valcoef(cutoff_commas(i-1)+1:cutoff_commas(i)-1));
            radius(i)=str2double(PTV_Data.Track.valprops.valrad(radius_commas(i-1)+1:radius_commas(i)-1));
            MAD_U(i)=str2double(PTV_Data.Track.valprops.MAD_U(MAD_U_commas(i-1)+1:MAD_U_commas(i)-1));
            MAD_V(i)=str2double(PTV_Data.Track.valprops.MAD_V(MAD_V_commas(i-1)+1:MAD_V_commas(i)-1));
        else
            cutoff(i)=str2double(PTV_Data.Track.valprops.valcoef(cutoff_commas(i-1)+1:end));
            radius(i)=str2double(PTV_Data.Track.valprops.valrad(radius_commas(i-1)+1:end));
            MAD_U(i)=str2double(PTV_Data.Track.valprops.MAD_U(MAD_U_commas(i-1)+1:end));
            MAD_V(i)=str2double(PTV_Data.Track.valprops.MAD_V(MAD_V_commas(i-1)+1:end));
        end
        if radius(i)~=0
            method{i} = 'median';
        else
            method{i} = 'coeff';
        end
    end
    Data.Track.weights=[str2double(PTV_Data.Track.disweight),str2double(PTV_Data.Track.sizeweight),str2double(PTV_Data.Track.intensityweight)];
    
    Data.Track.valprops.run=str2double(PTV_Data.Track.valprops.run);
    Data.Track.valprops.method = method;
    Data.Track.valprops.C_cutoff=cutoff;
    Data.Track.valprops.s_radius=radius;
    Data.Track.valprops.MAD_U=MAD_U;
    Data.Track.valprops.MAD_V=MAD_V;
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Particle ID Code
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if Data.ID.run
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Parallel Processing
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if str2double(PTV_Data.par)
            if length(str2double(PTV_Data.imfstart):str2double(PTV_Data.imfstep):str2double(PTV_Data.imfend)) < str2double(PTV_Data.parprocessors)
                PTV_Data.parprocessors=num2str(length(str2double(PTV_Data.imfstart):str2double(PTV_Data.imfstep):str2double(PTV_Data.imfend)));
            end
            fprintf('\n--- Initializing Processor Cores for Parallel Job ----\n')
            try
                matlabpool('open','local',PTV_Data.parprocessors);
            catch
                try
                    matlabpool close
                    matlabpool('open','local',PTV_Data.parprocessors);
                catch
                    beep
                    disp('Error Running Job in Parallel - Defaulting to Single Processor')
                    poolopen=0;
                end
            end
        end
        
        fprintf('\n-------------- Processing Dataset (started at %s) ------------------\n', datestr(now));
        
        % Write the Experimental Summary
        write_expsummary(PTV_Data);
        PTV_Data.ID.run='0';% Turn off ID so that it will write the next exp_summary in sizing folder
        
        I1 = Data.imfstart:Data.imfstep:Data.imfend;
        particleIDprops       = Data.ID;
        particleIDprops.Data  = Data;
        
        if str2double(PTV_Data.par) && matlabpool('size')>1
            spmd
                verstr=version('-release');
                if str2double(verstr(1:4))>=2010
                    I1dist=getLocalPart(codistributed(I1,codistributor('1d',2)));
                else
                    I1dist=localPart(codistributed(I1,codistributor('1d',2),'convert'));
                end
                for i= 1:length(I1dist)
                    
                    t0 = clock;
                    loadname = sprintf('%%s%%s%%s%%0%0.0fd.%s',Data.imzeros);
                    
                    if strcmpi(Data.imext,'tif')
                        IM = double(imread(sprintf(loadname,Data.imdirec,Data.slsh,Data.imbase,I1dist(i),'tif')));
                    elseif strcmpi(Data.imext,'mat')
                        IM = load(sprintf(loadname,Data.imdirec,Data.slsh,Data.imbase,I1dist(i),'mat'));
                    else
                        error('Unknown Extension type: .%s use, please use either ''.tif '' or ''.mat'' ',Data.imext)
                    end
                    
                    if size(IM, 3) > 2
                        %Extract only red channel
                        if Data.channel == 1;
                            IM = IM(:,:,1);
                            %Extract only green channel
                        elseif Data.channel == 2;
                            IM = IM(:,:,2);
                            %Extract only blue channel
                        elseif Data.channel == 3;
                            IM = IM(:,:,3);
                            %Weighted average of channels (see rgb2gray for
                            %explanation of weighting factors)
                        elseif Data.channel == 4;
                            IM = 0.2989 * IM(:, :, 1) + 0.5870 * IM(:, :, 2) + 0.1140 * IM(:, :, 3);
                            %Evenly weighted mean of channels
                        elseif Data.channel == 5;
                            IM = (IM(:,:,1) + IM(:,:,2) + IM(:,:,3))/3;
                        end
                    else
                        %	Take only red channel
                        IM =IM(:,:,1);
                    end
                    
                    s_num = I1dist(i);
                    
                    [p_matrix,peaks,num_p]=particle_ID_MAIN_V1(IM,particleIDprops,s_num);
                    %          figure; imagesc(p_matrix,[0 num_p]);  set(gca,'DataAspectRatio',[1 1 1])
                    
                    eltime=etime(clock,t0);
                    fprintf('%s...done!\t Time: %0.2i:%0.2i.%0.0f\n',sprintf(loadname,'processing image-',Data.slsh,Data.imbase,I1dist(i),Data.imext)...
                        ,floor(eltime/60),floor(rem(eltime,60)),rem(eltime,60)-floor(rem(eltime,60)))
                end
            end
            clear I1dist IM eltime i loadname num_p p_matrix peaks s_num t0 verstr
        else
            for i= 1:length(I1)
                
                t0 = clock;
                loadname = sprintf('%%s%%s%%s%%0%0.0fd.%s',Data.imzeros);
                fprintf(loadname,'processing image-',Data.slsh,Data.imbase,I1(i),Data.imext);
                
                if strcmpi(Data.imext,'tif')
                    IM = double(imread(sprintf(loadname,Data.imdirec,Data.slsh,Data.imbase,I1(i),'tif')));
                elseif strcmpi(Data.imext,'mat')
                    IM = load(sprintf(loadname,Data.imdirec,Data.slsh,Data.imbase,I1(i),'mat'));
                else
                    error('Unknown Extension type: .%s use, please use either ''.tif '' or ''.mat'' ',Data.imext)
                end
                if size(IM, 3) > 2
                    %Extract only red channel
                    if Data.channel == 1;
                        IM = IM(:,:,1);
                        %Extract only green channel
                    elseif Data.channel == 2;
                        IM = IM(:,:,2);
                        %Extract only blue channel
                    elseif Data.channel == 3;
                        IM = IM(:,:,3);
                        %Weighted average of channels (see rgb2gray for
                        %explanation of weighting factors)
                    elseif Data.channel == 4;
                        IM = 0.2989 * IM(:, :, 1) + 0.5870 * IM(:, :, 2) + 0.1140 * IM(:, :, 3);
                        %Evenly weighted mean of channels
                    elseif Data.channel == 5;
                        IM = (IM(:,:,1) + IM(:,:,2) + IM(:,:,3))/3;
                    end
                else
                    %	Take only red channel
                    IM =IM(:,:,1);
                end
                
                s_num = I1(i);
                
                [p_matrix,peaks,num_p]=particle_ID_MAIN_V1(IM,particleIDprops,s_num);
                %          figure; imagesc(p_matrix,[0 num_p]);  set(gca,'DataAspectRatio',[1 1 1])
                
                fprintf('...done!\t');
                eltime=etime(clock,t0);
                fprintf('Time: %0.2i:%0.2i.%0.0f\n',floor(eltime/60),floor(rem(eltime,60)),rem(eltime,60)-floor(rem(eltime,60)))
                
            end
        end
        
        save(fullfile(Data.ID.save_dir,'particle_ID_parameters.mat'));
        
        
        if str2double(PTV_Data.par)
            matlabpool('close')
        end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Particle Size Code
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if Data.Size.run
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Parallel Processing
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if str2double(PTV_Data.par)
            if length(str2double(PTV_Data.imfstart):str2double(PTV_Data.imfstep):str2double(PTV_Data.imfend)) < str2double(PTV_Data.parprocessors)
                PTV_Data.parprocessors=num2str(length(str2double(PTV_Data.imfstart):str2double(PTV_Data.imfstep):str2double(PTV_Data.imfend)));
            end
            fprintf('\n--- Initializing Processor Cores for Parallel Job ----\n')
            try
                matlabpool('open','local',PTV_Data.parprocessors);
            catch
                try
                    matlabpool close
                    matlabpool('open','local',PTV_Data.parprocessors);
                catch
                    beep
                    disp('Error Running Job in Parallel - Defaulting to Single Processor')
                    poolopen=0;
                end
            end
        end
        fprintf('\n-------------- Processing Dataset (started at %s) ------------------\n', datestr(now));
        
        % Write the Experimental Summary
        write_expsummary(PTV_Data);
        PTV_Data.Size.run='0';% Turn off ID so that it will write the next exp_summary in sizing folder
        
        I1 = Data.imfstart:Data.imfstep:Data.imfend;
        
        if str2double(PTV_Data.par) && matlabpool('size')>1
            spmd
                verstr=version('-release');
                if str2double(verstr(1:4))>=2010
                    I1dist=getLocalPart(codistributed(I1,codistributor('1d',2)));
                else
                    I1dist=localPart(codistributed(I1,codistributor('1d',2),'convert'));
                end
                for i= 1:length(I1dist)
                    t0 = clock;
                    loadname = sprintf('%%s%%s%%s%%0%0.0fd.%%s',Data.imzeros);
                    
                    if strcmpi(Data.imext,'tif')
                        IM = imread(sprintf(loadname,Data.imdirec,Data.slsh,Data.imbase,I1dist(i),'tif'));
                    elseif strcmpi(Data.imext,'mat')
                        IM = load(sprintf(loadname,Data.imdirec,Data.slsh,Data.imbase,I1dist(i),'mat'));
                    else
                        error('Unknown Extension type: .%s use, please use either ''.tif '' or ''.mat'' ',Data.imext)
                    end
                    
                    if size(IM, 3) > 2
                        %Extract only red channel
                        if Data.channel == 1;
                            IM = IM(:,:,1);
                            %Extract only green channel
                        elseif Data.channel == 2;
                            IM = IM(:,:,2);
                            %Extract only blue channel
                        elseif Data.channel == 3;
                            IM = IM(:,:,3);
                            %Weighted average of channels (see rgb2gray for
                            %explanation of weighting factors)
                        elseif Data.channel == 4;
                            IM = 0.2989 * IM(:, :, 1) + 0.5870 * IM(:, :, 2) + 0.1140 * IM(:, :, 3);
                            %Evenly weighted mean of channels
                        elseif Data.channel == 5;
                            IM = (IM(:,:,1) + IM(:,:,2) + IM(:,:,3))/3;
                        end
                    else
                        %	Take only red channel
                        IM =IM(:,:,1);
                    end
                    
                    im1=IM(:,:);
                    im1(im1<=Data.Size.thresh) = 0;
                    
                    %load in the identified particles
                    sname = sprintf('%%s%%0%0.0fd',Data.imzeros);
                    
                    ID_info = load(fullfile(Data.ID.save_dir,sprintf(sname,Data.ID.s_name,I1dist(i),'.mat')));
                    
                    %size the particles
                    sizeprops       = Data.Size;
                    sizeprops.Data  = Data;
                    sizeprops.s_num =I1dist(i);
                    
                    [SIZE1.XYDiameter,SIZE1.mapsizeinfo,SIZE1.locxy]=particle_size_MAIN_V1(im1,ID_info.p_matrix,...
                        ID_info.num_p,sizeprops);
                    
                    eltime=etime(clock,t0);
                    fprintf('%s...done!\t Time: %0.2i:%0.2i.%0.0f\n',sprintf(loadname,'Sizing Frame-',Data.slsh,Data.imbase,I1dist(i),Data.imext),...
                        floor(eltime/60),floor(rem(eltime,60)),rem(eltime,60)-floor(rem(eltime,60)))
                    
                end
            end
            clear I1dist ID_info IM SIZE1 eltime i im1 loadname sizeprops sname t0 verstr
        else
            for i = 1:length(I1)
                
                t0 = clock;
                loadname = sprintf('%%s%%s%%s%%0%0.0fd.%%s',Data.imzeros);
                fprintf(loadname,'Sizing Frame-',Data.slsh,Data.imbase,I1(i),Data.imext)
                
                if strcmpi(Data.imext,'tif')
                    IM = imread(sprintf(loadname,Data.imdirec,Data.slsh,Data.imbase,I1(i),'tif'));
                elseif strcmpi(Data.imext,'mat')
                    IM = load(sprintf(loadname,Data.imdirec,Data.slsh,Data.imbase,I1(i),'mat'));
                else
                    error('Unknown Extension type: .%s use, please use either ''.tif '' or ''.mat'' ',Data.imext)
                end
                
                if size(IM, 3) > 2
                    %Extract only red channel
                    if Data.channel == 1;
                        IM = IM(:,:,1);
                        %Extract only green channel
                    elseif Data.channel == 2;
                        IM = IM(:,:,2);
                        %Extract only blue channel
                    elseif Data.channel == 3;
                        IM = IM(:,:,3);
                        %Weighted average of channels (see rgb2gray for
                        %explanation of weighting factors)
                    elseif Data.channel == 4;
                        IM = 0.2989 * IM(:, :, 1) + 0.5870 * IM(:, :, 2) + 0.1140 * IM(:, :, 3);
                        %Evenly weighted mean of channels
                    elseif Data.channel == 5;
                        IM = (IM(:,:,1) + IM(:,:,2) + IM(:,:,3))/3;
                    end
                else
                    %	Take only red channel
                    IM =IM(:,:,1);
                end
                
                im1=IM(:,:);
                im1(im1<=Data.Size.thresh) = 0;
                
                %load in the identified particles
                sname = sprintf('%%s%%0%0.0fd',Data.imzeros);
                
                ID_info = load(fullfile(Data.ID.save_dir,sprintf(sname,Data.ID.s_name,I1(i),'.mat')));
                
                %size the particles
                sizeprops       = Data.Size;
                sizeprops.Data  = Data;
                sizeprops.s_num = I1(i);
                
                [SIZE1.XYDiameter,SIZE1.mapsizeinfo,SIZE1.locxy]=particle_size_MAIN_V1(im1,ID_info.p_matrix,...
                    ID_info.num_p,sizeprops);
                
                fprintf('...done!\t');
                eltime=etime(clock,t0);
                fprintf('Time: %0.2i:%0.2i.%0.0f\n',floor(eltime/60),floor(rem(eltime,60)),rem(eltime,60)-floor(rem(eltime,60)))
                
            end
        end
        
        save(fullfile(Data.Size.save_dir,'particle_SIZE_parameters.mat'));
        
        if str2double(PTV_Data.par)
            matlabpool close
        end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Particle Tracking Code
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if Data.Track.run
        im_list = Data.imfstart:Data.imfstep:Data.imfend;
        completed_tracks = cell(length(im_list),1);
        
        % Write the Experimental Summary
        write_expsummary(PTV_Data);
        
        Data.Track.Data  = Data;
        
        fprintf('\n-------------- Processing Dataset (started at %s) ------------------\n', datestr(now));
        
        for k=1:length(im_list)-1
            t0 = clock;
            
            fprintf('Tracking Frame %s %6.0f to %6.0f\t',Data.imbase,im_list(k),im_list(k)+Data.imcstep)
            
            sname = sprintf('%%s%%0%0.0fd',Data.imzeros);
            SIZE1=load(fullfile(Data.Size.save_dir,sprintf(sname,Data.Size.s_name,im_list(k),'.mat')));
            SIZE2=load(fullfile(Data.Size.save_dir,sprintf(sname,Data.Size.s_name,im_list(k)+Data.imcstep,'.mat')));
            
            %remove all NaNs from the particle arrays (indicated a failed
            %sizing method)
            check1=isnan(SIZE1.XYDiameter(:,3))==0;
            check2=isnan(SIZE2.XYDiameter(:,3))==0;
            SIZE1.XYDiameter=SIZE1.XYDiameter(check1,:);
            SIZE2.XYDiameter=SIZE2.XYDiameter(check2,:);
            
            if ~isempty(SIZE1.XYDiameter) && ~isempty(SIZE2.XYDiameter)
                
                estlocprops      = Data.Track;
                estlocprops.Data = Data;
                
                %compute location prediction for the particles in im1
                if nnz( strcmpi(Data.Track.method,{'piv' 'piv-ptv'}) )
                    estlocprops.PIVprops.frame1 = im_list(k);
                    if (k==1 && strcmpi(Data.Track.method,{'piv-ptv'}))
                        org_weight = Data.Track.PIV_PTV_weight;
                        estlocprops.PIV_PTV_weight=1;
                    elseif (k>1 && strcmpi(Data.Track.method,{'piv-ptv'}))
                        estlocprops.PIV_PTV_weight=org_weight;
                    end
                    
                elseif nnz( k==1 && strcmpi(estlocprops.method,{'ptv'})  )
                    %THIS SECTION OF CODE PROVIDES INTIALIZATION FOR THE PTV LOCATION
                    %ESTIMATION - USER SHOULD SUPPLY THIS IF POSSIBLE TO INCREASE
                    %TRACKING RELIABILITY
                end
                
                estlocprops.s_num=im_list(k);
                %     [X2_est,Y2_est,Z2_est]=particle_estloc_MAIN_V1(completed_tracks{k},SIZE1,estlocprops);
                [X2_est,Y2_est,Z2_est]=particle_estloc_MAIN_V2(completed_tracks{k},SIZE1,estlocprops);
                
                %track images (run for all four pair-matching methods)
                Data.Track.s_num = im_list(k);
                
                [new_tracks]=particle_track_MAIN_V1(X2_est,Y2_est,Z2_est,SIZE1,SIZE2,Data.Track,estlocprops.valprops);
                
                
                %         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %         %% Plot the tracking results as a connected scatter plot of particle
                %         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %         if k == 1%length(im_list)-1
                %             % %     figure(1)
                %             % %     quiver(SIZE1.XYDiameter(:,1),SIZE1.XYDiameter(:,2),...
                %             % %         X2_est-SIZE1.XYDiameter(:,1),Y2_est-SIZE1.XYDiameter(:,2),0,'Color','BLACK');
                %             % %     axis image
                %             %
                %             %positions in IM1 and IM2
                %             figure(2);
                %             scatter(new_tracks(:,1),new_tracks(:,3),'ob','sizedata',10);
                %             hold on
                %             scatter(new_tracks(:,2),new_tracks(:,4),'+r','sizedata',10);  hold on
                %             scatter(X2_est,Y2_est,'.BLACK');  hold on
                %             axis tight
                %             for j=1:size(new_tracks,1)
                %                 line([new_tracks(j,1);new_tracks(j,2)],[new_tracks(j,3);new_tracks(j,4)]);
                %                 hold on
                %             end
                %             axis image
                %             legend({'image1','image2'},'Location','SouthOutside','Orientation','horizontal')
                %             title(sprintf('Track frame %0.0f to %0.0f',im_list(k),im_list(k+1)))
                %             saveas(gcf,sprintf('%s%sTrack_Image_Lag_Tracks_%06d.tif',save_dir,slsh,k))
                %
                %             close(figure(3))
                %             figure(3);
                %             X=new_tracks(:,1); Y=new_tracks(:,3);
                %             U=new_tracks(:,2)-new_tracks(:,1); V=new_tracks(:,4)-new_tracks(:,3);
                %             quiver(X,Y,U,V,1.5);
                %             if k > im_list(1)
                %                 temp = cell2mat(completed_tracks{k});
                %                 X=temp(:,1); Y=temp(:,3);
                %                 U=temp(:,2)-temp(:,1); V=temp(:,4)-temp(:,3);
                %                 hold on
                %                 quiver(X,Y,U,V,1.5,'RED');
                %             end
                %             axis image
                %             title(sprintf('Track frame %0.0f to %0.0f',im_list(k),im_list(k+1)))
                %             saveas(gcf,sprintf('%s%sTrack_Image_Vec_Refine_%06d.tif',save_dir,slsh,k))
                %             drawnow
                %             keyboard
                %         end
                %         %%
                %         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                %restructure and reassign new_tracks to completed_tracks
                new_tracks=num2cell(new_tracks,2);
                completed_tracks{k+1}=new_tracks;
                
                fprintf('...done!\t');
                eltime=etime(clock,t0);
                fprintf('Time: %0.2i:%0.2i.%0.0f\n',floor(eltime/60),floor(rem(eltime,60)),rem(eltime,60)-floor(rem(eltime,60)))
                
            else
                fprintf('...SKIPPED Due to lack of particles!\n');
            end
            
            %save processing parameters
            save(sprintf('%s%s%s',Data.Track.save_dir,'particle_TRACK_Val_PIV_parameters.mat'));
        end
        fprintf('---------------- Job Completed at %s ---------------------\n', datestr(now));
    end
end
end
