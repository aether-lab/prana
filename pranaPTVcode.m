function pranaPTVcode(PTV_Data)

if ispc
    Data.slsh='\';
else
    Data.slsh='/';
end
% --- Images ---
Data.imdirec=PTV_Data.imdirec;
Data.imbase=PTV_Data.imbase;
Data.imext=PTV_Data.imext;
Data.imzeros=str2double(PTV_Data.imzeros);
Data.imcstep=str2double(PTV_Data.imcstep);
Data.imfstep=str2double(PTV_Data.imfstep);
Data.imfstart=str2double(PTV_Data.imfstart);
Data.imfend=str2double(PTV_Data.imfend);
Data.channel=str2double(PTV_Data.channel);
if Data.channel==6
    fprintf(['Ensemble Color does not currently work with tracking code.\n'...
        'Resetting the color channel to 1 (grey scale)\n'])
    Data.channel=1;
end

% --- ID ---
IDmethod = {'blob','dynamic','combined'};
Data.ID.method=IDmethod{str2double(PTV_Data.ID.method)};
Data.ID.run=str2double(PTV_Data.ID.runid);
Data.ID.v=str2double(PTV_Data.ID.imthresh);
Data.ID.contrast_ratio=0;%str2double(Data.ID.contrast_ratio);
Data.ID.s_num=0;%str2double(Data.ID.s_num);
Data.ID.s_name=PTV_Data.ID.savebase;
Data.ID.save_dir=PTV_Data.ID.save_dir;

% --- Sizing ---
Data.Size.run=str2double(PTV_Data.Size.runsize);
Data.Size.thresh=str2double(PTV_Data.ID.imthresh);
Data.Size.method=str2double(PTV_Data.Size.method);
Data.Size.p_area=str2double(PTV_Data.Size.min_area);
Data.Size.sigma=str2double(PTV_Data.Size.std);
Data.Size.errors=0;%str2double(Data.Size.errors);
Data.Size.s_name=PTV_Data.Size.savebase;
Data.Size.save_dir=PTV_Data.Size.save_dir;

% --- PIV Info ---
if str2double(PTV_Data.datout)
     Data.Track.PIVprops.extension='.dat';
elseif str2double(PTV_Data.multiplematout)
    Data.Track.PIVprops.extension='.mat';
end
Data.Track.PIVprops.load_dir       = PTV_Data.outdirec;
Data.PIV.Data=PTV_Data;
eval(['Data.PIV.Data.outbase = PTV_Data.PIV' PTV_Data.passes '.outbase;'])

% --- Tracking ---
Data.Track.run=str2double(PTV_Data.Track.runtrack);
PTVmethod = {'ptv','piv','piv-ptv'};
PTVpredict = {'static','dynamic'};
Data.Track.method = PTVmethod{str2double(PTV_Data.Track.method)};
Data.Track.predict_mode = PTVpredict{str2double(PTV_Data.Track.prediction)};
Data.Track.PIV_PTV_weight=str2double(PTV_Data.Track.PIVweight);
Data.Track.plotfig=0;%str2double(Data.Track.plotfig);
Data.Track.s_radius=str2double(PTV_Data.Track.radius);
Data.Track.r_weight=str2double(PTV_Data.Track.estradius);
Data.Track.edgeval=str2double(PTV_Data.Track.estweight);
Data.Track.numvecs=str2double(PTV_Data.Track.vectors);
Data.Track.max_iterations=str2double(PTV_Data.Track.iterations);
Data.Track.s_name=PTV_Data.Track.savebase;
Data.Track.save_dir=PTV_Data.Track.save_dir;

cutoff_commas=strfind(PTV_Data.Track.valprops.valcoef,',');
radius_commas=strfind(PTV_Data.Track.valprops.valrad,',');
MAD_U_commas=strfind(PTV_Data.Track.valprops.MAD_U,',');
MAD_V_commas=strfind(PTV_Data.Track.valprops.MAD_V,',');

cutoff = zeros(numel(cutoff_commas)+1,1);
radius = zeros(numel(radius_commas)+1,1);
MAD_U = zeros(numel(MAD_U_commas)+1,1);
MAD_V = zeros(numel(MAD_V_commas)+1,1);

Data.Track.valprops.numpass=numel(cutoff_commas)+1;

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
%% Parallel Processing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
    fprintf('\n-------------- Processing Dataset (started at %s) ------------------\n', datestr(now));
else
    fprintf('\n-------------- Processing Dataset (started at %s) ------------------\n', datestr(now));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Particle ID Code
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if Data.ID.run
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
                        %ensemble correlation of channels
%                     elseif channel == 6;
%                         im1=im1(:,:,1:3);
%                         im2=im2(:,:,1:3);
                    end
                else
                    %	Take only red channel
                    IM =IM(:,:,1);
%                     Data.channel = 1;
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
                    %ensemble correlation of channels
%                 elseif channel == 6;
%                     im1=im1(:,:,1:3);
%                     im2=im2(:,:,1:3);
                end
            else
                %	Take only red channel
                IM =IM(:,:,1);
                Data.channel = 1;
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
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if str2double(PTV_Data.par)
    matlabpool('close')
end

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
    fprintf('\n-------------- Processing Dataset (started at %s) ------------------\n', datestr(now));
else
    fprintf('\n-------------- Processing Dataset (started at %s) ------------------\n', datestr(now));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Particle Size Code
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if Data.Size.run
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
                        %ensemble correlation of channels
%                     elseif channel == 6;
%                         im1=im1(:,:,1:3);
%                         im2=im2(:,:,1:3);
                    end
                else
                    %	Take only red channel
                    IM =IM(:,:,1);
%                     Data.channel = 1;
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
                %ensemble correlation of channels
%             elseif channel == 6;
%                 im1=im1(:,:,1:3);
%                 im2=im2(:,:,1:3);
            end
        else
            %	Take only red channel
            IM =IM(:,:,1);
            Data.channel = 1;
        end

        im1=IM(:,:);  
        im1(im1<=Data.Size.thresh) = 0;

        %load in the identified particles
        sname = sprintf('%%s%%0%0.0fd',Data.imzeros);

        ID_info = load(fullfile(Data.ID.save_dir,sprintf(sname,Data.ID.s_name,I1(i),'.mat')));

        %size the particles
        sizeprops       = Data.Size;
        sizeprops.Data  = Data;
        sizeprops.s_num =I1(i);
        
        [SIZE1.XYDiameter,SIZE1.mapsizeinfo,SIZE1.locxy]=particle_size_MAIN_V1(im1,ID_info.p_matrix,...
            ID_info.num_p,sizeprops);

        fprintf('...done!\t');
        eltime=etime(clock,t0);
        fprintf('Time: %0.2i:%0.2i.%0.0f\n',floor(eltime/60),floor(rem(eltime,60)),rem(eltime,60)-floor(rem(eltime,60)))

    end
    end

    save(fullfile(Data.Size.save_dir,'particle_SIZE_parameters.mat'));
end

if str2double(PTV_Data.par)
    matlabpool close
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Particle Tracking Code
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if Data.Track.run
    im_list = Data.imfstart:Data.imfstep:Data.imfend;
    completed_tracks = cell(length(im_list),1);
    
    % Write the Experimental Summary
    write_expsummary(PTV_Data);
    
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
        Data.Track.Data  = Data;
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
        
    end
    
    %save processing parameters
    save(sprintf('%s%s%s',Data.Track.save_dir,'particle_TRACK_Val_PIV_parameters.mat'));
end
fprintf('---------------- Job Completed at %s ---------------------\n', datestr(now));
end


















































%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ID Functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [p_matrix,peaks,num_p]=particle_ID_MAIN_V1(im,particleIDprops,s_num)
%
% [p_matrix,peaks,num_p]=particle_ID_MAIN_V1(im,particleIDprops);
%
% This is a MAIN function for particle identification and segmentation from
% grayscale images.  It calls several differet sub-functions (contained
% below) depending the image parameters, specifically the particle image
% density.  For sparse-to-medium seeding densities, the 'blob' method is
% recommended.  For a densly seeded image, with many overlapped particle
% images, the 'dynamic' method is recommended.  The combined method is a
% beta-version and should be used with caution since it has not been 
% rigorly tested 
%
% INPUTS
%    im1,im2 - image to be processesed; should be uint8
%    particleIDprops - processing parameters related to particle ID
%        particleIDprops.method - selects the ID/segmentation method
%            'blob','dynamic','combined'
%        particleIDprops.v - base thresholding for the image
%        particleIDprops.contrast_ratio - ratio of max particle intensity to
%            possible pixel intensity, used to limit the intensity range a
%            particle can expand from its max value ('dynamic' only)
%            *typically set to 0, higher values can reduce processing time*
%        particleIDprops.save_dir - saveing filepath; make '0' if not saving
%        particleIDprops.slsh - '\' or '/' (PC or Linux)
%        particleIDprops.s_name - base saving name
%        particleIDprops.s_num - saving number
%
% OUTPUTS
%    p_matrix (2D array) - matrix of particle identification and extent
%    peaks (2D array) - matrix of identified image peaks (by image erosion)
%    num_p (num) - number of identified particles in p_matrix
%
%(v1) N. Cardwell - 10.29.2009
%(v1.1) N.Cardwell - 11.1.2009 (combined particleIDprops and saveprops into
%   a single input, particleIDprops 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% identify and segregate the particles %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%chose between different particle ID and segregation methods

switch lower(particleIDprops.method)
    case {'blob'}
        [p_matrix,peaks,num_p]=blob_segmentation(im,particleIDprops.v);
        
    case {'dynamic'}
        [p_matrix,peaks,num_p] = dynamic_threshold_segmentation_v3(im,...
            particleIDprops.v,particleIDprops.contrast_ratio);
    
    case {'combined'}
        [p_matrix,peaks,num_p]=combined_partID(im,particleIDprops.v);
    otherwise
        error('Unknown ID segmentation method\n')
end
%NOTE - currently the speed of each method is as follows:
%   'blob' (256x256) - 0.74sec per image pair - 1676 particles
%   'dynamic (256x256) - 9.07sec per image pair - 3592 particles
%   'combined (256x256) - 11.30sec per image pair - 4490 particles (hmm...)

if ~isempty(particleIDprops.save_dir)
    if exist(particleIDprops.save_dir,'dir')~=7
        fprintf('Making Save Directory %s \t',particleIDprops.save_dir)
        mkdir(particleIDprops.save_dir)
    end
    sname = sprintf('%%0%0.0fd',particleIDprops.Data.imzeros);
    save(fullfile(particleIDprops.save_dir,[particleIDprops.s_name sprintf(sname,s_num)]),...
        'particleIDprops','p_matrix','peaks','num_p');
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUB-FUNCTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [p_matrix,peaks,num_p]=blob_segmentation(im,v)
%
%[p_matrix,peaks,num_p]=mapparticles(im,v)
%
%Simple particle identification method which, after global thresholding,
%groups clusters of adjacent pixels together ('blobs') and labels them as a
%particle.  Very fast and works well for sparse-to-medium image seeding
%densities.  The method has no mechanism to separate overlapped particles
%other tham aggressive thresholding, which adversly affects the image.
%
%INPUTS
%   im - image to be segmented (matrix-uint8)
%   v  - intial threshold value (num)
%
%OUTPUTS
%   p_matrix (2D array) - matrix of particle identification and extent
%   peaks (2D array) - matrix of identified image peaks (by image erosion)
%   num_p (num) - number of identified particles in p_matrix
%
%N.Cardwell - 2.8.2008 (modified from the work by M.Brady)
%N.Cardwell - 10.27.2009 (added a base thresholding operation and a peak 
%   identification operation, removed section of code for the creation of
%   mapint and locxy - no longer needed with the TRACKING_V2 code)

%globally threshold the image
im_thresh=im;  im_thresh(im < v)=0;

%group adjacent pixels with nonzero intensity and assign each cluster a
%particle number and store in p_matrix, also returns total number of
%particles detected (num_p), searches by columns
[p_matrix,num_p]=bwlabel(im_thresh,8);

%determine the peak intensity pixel of each pixel blob
peaks=zeros(size(p_matrix));
for i=1:num_p
    [r_i,c_i]=find(p_matrix==i);  indx_i=sub2ind(size(peaks),r_i,c_i);
    I_i=im(indx_i);
    pixels_i=[r_i,c_i,I_i];
    pixels_i_sort=sortrows(pixels_i,3);
    peaks(pixels_i_sort(end,1),pixels_i_sort(end,2))=i;
end

end

function [p_matrix,peaks,num_p]=dynamic_threshold_segmentation_v3(im,v1,contrast_ratio)
%
%[p_matrix,peaks,num_p]=dynamic_threshold_segmentation_v3(im,v1,contrast_ratio)
%
%Uses an erosion/dilation process to identify peaks and then determine the
%extents of each particle.  Algorithm is very effective at separating
%overlapped particle images.  However, there is significant increase in
%processing time over the 'blob' method (apprximately 12X).
%
%INPUTS
%   im - image to be segmented (matrix-uint8)
%   v1 - intial threshold value (num)
%   contrast_ratio - may remove this...for now just set to zero
%
%OUTPUTS
%   p_matrix (2D array) - matrix of particle identification and extent
%   peaks (2D array) - matrix of identified image peaks (by image erosion)
%   num_p (num) - number of identified particles in p_matrix
%
%N.Cardwell (v3)   - 8.12.09
%N.Cardwell (v3.1) - 10.12.09 (replaced 'sub2ind' function call with a
%   direct calculation to increase speed)

%intially threshold the image using 'threshold'
im_thresh=im;  im_thresh(im_thresh<=v1)=0;

%use built in function to identify regional maxima (i.e peaks)
% tic;  fprintf('Locating image peaks...');
BW_max=imregionalmax(im_thresh);
[p_matrix,num_p]=bwlabel(BW_max,8);
% fprintf('DONE!----')
% fprintf(strcat('elapsed time=',num2str(toc),'seconds\n'));

%set the maximum intensity for each particle
Imax=zeros(num_p,1);
for i=1:num_p
    Imax(i,1)=max(im_thresh(p_matrix==i));
end;  clear i

%perform dilation on each particle until the contrast criterion is met for
%the boundary pixels of each particle; also has a check to make sure that
%the expanding pixels cannot grap "brighter" pixels
% tic;  fprintf('Expanding peaks........');
peaks=p_matrix;
p_matrix_temp=p_matrix;  s=size(p_matrix);
flags=ones(num_p,1);
%figure
while nnz(flags) > 0
%     imagesc(p_matrix_temp);  set(gca,'DataAspectRatio',[1 1 1]);
%     pause(0.5)
    for i=1:num_p

        %check to see if the particle has been flaged (ie. completly expanded)
        if flags(i,1)==1

            %initialize particle conditions
            l_row=size(p_matrix_temp,1);
            [r,c]=find(p_matrix_temp==i);  part_pixels=zeros(length(r),length(c));
            part_pixels(:,1)=r; part_pixels(:,2)=c;
            part_index=(part_pixels(:,2)-1).*l_row+part_pixels(:,1);
%             [part_pixels(:,1),part_pixels(:,2)]=find(p_matrix_temp==i);

%             part_index=find(p_matrix_temp==i);  part_pixels=zeros(length(part_index),2);
            

            %expand all particle pixels by one in each direction
            possible_pixels=[];            
            for j=1:size(part_pixels,1)
                p_pix_j=part_pixels(j,:);
                I_pix_j=im_thresh(p_pix_j(1),p_pix_j(2));
                if p_pix_j(2) > 1
%                     poss_pix=sub2ind(s,p_pix_j(1),p_pix_j(2)-1);
                    poss_pix=((p_pix_j(2)-1)-1)*l_row+p_pix_j(1);
                    if im_thresh(poss_pix) < I_pix_j
                        possible_pixels=[possible_pixels,poss_pix];
                    end
                end
                if p_pix_j(2) < s(2)
%                     poss_pix=sub2ind(s,p_pix_j(1),p_pix_j(2)+1);
                    poss_pix=((p_pix_j(2)-1)+1)*l_row+p_pix_j(1);
                    if im_thresh(poss_pix) < I_pix_j
                        possible_pixels=[possible_pixels,poss_pix];
                    end
                end
                if p_pix_j(1) > 1
%                     poss_pix=sub2ind(s,p_pix_j(1)-1,p_pix_j(2));
                    poss_pix=((p_pix_j(2)-1))*l_row+(p_pix_j(1)-1);
                    if im_thresh(poss_pix) < I_pix_j
                        possible_pixels=[possible_pixels,poss_pix];
                    end
                end
                if p_pix_j(1) < s(1)
%                     poss_pix=sub2ind(s,p_pix_j(1)+1,p_pix_j(2));
                    poss_pix=((p_pix_j(2)-1))*l_row+(p_pix_j(1)+1);
                    if im_thresh(poss_pix) < I_pix_j
                        possible_pixels=[possible_pixels,poss_pix];
                    end
                end
            end

            %remove all non-unique identifications in 'possible_pixels'
            possible_pixels=unique(possible_pixels);

            %check to see if the possible pixels are already part of the
            %particle or another particle (remove if true)
            check=p_matrix_temp(possible_pixels)~=0;
            possible_pixels=possible_pixels(check==0);

            %see if any of the border pixels satisy the contrast criterion
            %(if so then attach to the particle)
            check2=single(im_thresh(possible_pixels))./single(Imax(i,1)) > contrast_ratio;
            if nnz(check2)~=0
                %assign the border pixels to p_matrix final and border_pixels
                p_matrix_temp(possible_pixels(check2))=i;
            else
                flags(i,1)=0;
            end
        end
    end
end
% fprintf('DONE!----')
% fprintf(strcat('elapsed time=',num2str(toc),'seconds\n'));

p_matrix=p_matrix_temp;
end

function [p_matrix,peaks,num_p]=combined_partID(im,v)
%
%This function uses a combination of the 'blob' and 'dynamic threshold
%segmentation' algorithms in an attempt to reduce the computational cost of
%the dynamic while still retaining its ability to accuartely segment 
%overlaped particles.  Special care must be taken in the intital
%thresholding of the image to partially segment the blobs without removing
%too many of the low intensity particle images
%
%(beta) N.Cardwell - 10.28.2009

%call the mapparticles subfunction ('blob' particle ID analysis)
[p_matrix,peaks,num_p]=mapparticles(im,v);

%preallocate arrays
p_matrix_new=zeros(size(p_matrix));  num_p_new=0;  peaks_new=zeros(size(p_matrix));

%main program loop - to segment or not to segment, that is the question!
for i=1:num_p
    blob_i=double(p_matrix==i);
    
    %use 'regionprops to det. the extent of blob_i and segment the blob out
    STATS = regionprops(blob_i,'BoundingBox');
    B_Box=ceil(STATS.BoundingBox);
    if B_Box(3)==1;  width=0;  else width=B_Box(3)-1;  end
    if B_Box(4)==1;  height=0;  else height=B_Box(4)-1;  end
    part_i=im( B_Box(2):B_Box(2)+height , B_Box(1):B_Box(1)+width );
    part_i(part_i<v)=0;
    
    %determine if blob_i has multiple peaks, if so then pass it to the
    %dynamic function for additional segmentation, otherwise keep the
    %particle and move on
    BW_max=imregionalmax(part_i);
    if nnz(BW_max)~=1
        %crop out the blob and segment it
        im_crop=im(B_Box(2):B_Box(2)+height , B_Box(1):B_Box(1)+width);
        [p_matrix_crop,peaks_crop,num_p_crop]=dynamic_threshold_segmentation_v3(im_crop,0,0);
        
        %reset the 'particle number' of the segmented blod so it can be
        %directly inserted into the 'new' parameters
        for j=1:num_p_crop
            num_p_new=num_p_new+1;
            peaks_crop(peaks_crop==j) = num_p_new;
            p_matrix_crop(p_matrix_crop==j) = num_p_new;
        end
        
        %insert the segmented/cropped particle into the 'new' parameters
        peaks_new( B_Box(2):B_Box(2)+height ,...
            B_Box(1):B_Box(1)+width ) = peaks_crop;
        p_matrix_new( B_Box(2):B_Box(2)+height ,...
            B_Box(1):B_Box(1)+width ) = p_matrix_crop;
    else
        %assign the particle to the 'new' parameters
        num_p_new=num_p_new+1;
        [r_peak,c_peak]=find(BW_max==1);
        peaks_new(B_Box(2)+(r_peak-1) , B_Box(1)+(c_peak-1)) = num_p_new;
        p_matrix_new(B_Box(2):B_Box(2)+height , B_Box(1):B_Box(1)+width) = num_p_new;
    end
end

%reassign output variables
num_p=num_p_new;
peaks=peaks_new;
p_matrix=p_matrix_new;

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Sizing Functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [XYDiameter,mapsizeinfo,locxy] = particle_size_MAIN_V1(im,p_matrix,num_p,sizeprops)
%
% [XYDiameter,mapsizeinfo,locxy] = particle_size_MAIN_V1(im,p_matrix,...
%     num_p,sizeprops);
%
% PROGRAM DESCRIPTION
% This is a MAIN function for the sizing of previously identified particles
% For a given image, this function reads in detected particles and then
% calculates the diameter and centroid location of each particle using 
% various fitting algorithms
%
% INPUTS
%   im - input image
%   p_matrix (2D array) - matrix of particle identification and extent
%   num_p (num) - number of identified particles in p_matrix
%   sizeprops - processing parameters related to particle sizing
%       sizeprops.method - choose sizing method
%           1-Intensity Weighted Centroid
%           2-Three Point Gaussian Method
%           3-Standard Four Point Gaussian
%           4-Continuous Four Point Gaussian
%           5-Least Squares Gaussian Fit
%           6-Continuous Least Squares Gaussian Fit
%       sizeprops.p_area - minimum area (in pixels) to identify a particle
%       sizeprops.sigma - number of standard deviations used to compute the diameter
%           if a gaussian sizing method (ie - anything but the intensity 
%           weighted centroid method) is selected
%       sizeprops.errors - switch: 
%           1-return NaN's if the selected method fails
%           0-return Intensity-Weighted Centroid results if method fails
%       sizeprops.save_dir - saveing filepath; make '0' if not saving
%       sizeprops.slsh - '\' or '/' (PC or Linux)
%       sizeprops.s_name - base saving name
%       sizeprops.s_num - saving number
%
% OUTPUTS
%   XYDiameter = 6-column matrix; 1st column is x-centorid location, 2nd
%       column is y-centroid location, 3rd column is particle diameter, 4th
%       column is true max. intensity, 5th column is particle id#, 6th
%       column is sizing method
%   mapsizeinfo - (2 column array) defines the size (row col) of the each
%       sized particle
%   locxy - (2 column array) ROW/COL? location associated with the upper 
%       left pixel of the particle's rectangular projection
%
%S. Raben - 9.15.08
%(v1) N. Cardwell - 10.29.2009 (now accepts generalized inputs from particleID 
%       methods; also revised header info and added save command)
%(v1.1) N. Cardwell - 11.1.2009 (combined particleIDprops and saveprops into
%       a single input, particleIDprops)
%(v1.2) N.Cardwell - 1.7.2010 (fixed the portions of the code related to
%       saving the IWC estimate whenever the user selected method fails)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%identifies/locates all particles with (mapparticles)
%(adjacent pixels of greater than 0 intensity)
%   mapint - a cellular 2-D array of particle intensity profiles
%       locxy - a 2 column array of location associated with the upper
%       left pixel index of the particle intensity profile's
%       rectangular projection, the format of locxy is column then row 
%       (or x then y)
[mapint,locxy,num_p,mapsizeinfo]=mapparticles_v3(im,p_matrix,num_p,sizeprops.p_area);

%initialize particle properties array
particleprops=zeros(num_p,6);

%calculate the Intensity Weigthed Centorid (IWC) using (centroidfit)
%these results will be outputted if the other methods return NaN's and
%the user does not want errors in the output
%returns    *x_centroid - x (column) index location of the particles IWC
%           *y_centroid - y (row) index location of the particles IWC
%           *diameter   - calculated diameter of the particle
%           *I0         - intensity value of the brightest pixel
c = 1;
if sizeprops.method<7
    for i=1:num_p
        [particleprops(c,1),particleprops(c,2),particleprops(c,3),particleprops(c,4)]=...
            centroidfit(mapint{i},locxy(i,:));
        particleprops(c,1)=particleprops(c,1) - 1;
        particleprops(c,2)=particleprops(c,2) - 1;
        particleprops(c,5)=i;
        particleprops(c,6)=1;
        c=c+1;
    end
end

%calculate the centroid location and particle diameter using the 3-point
%approximation method (threeptgausfit)
%returns    *x_centroid - x (column) index location of the particles IWC
%           *y_centroid - y (row) index location of the particles IWC
%           *diameter   - calculated diameter of the particle
%           *I0         - calculated maximum particle intensity
if sizeprops.method==2
    for i=1:num_p
        [x_cg,y_cg,D,I0,E,Meth]=Gaussfit(mapint{i},sizeprops.method-1,sizeprops.sigma);
        x_c = locxy(i,1)+x_cg-1;
        y_c = locxy(i,2)+y_cg-1;

        if nnz(isnan([x_c,y_c,D,I0])>1)
            if sizeprops.errors==0
                %retain the IWC estimate
            else
                particleprops(i,:)=[x_c,y_c,max(D),max(I0),i,Meth+1];
            end
        else 
            particleprops(i,:)=[x_c,y_c,max(D),max(I0),i,Meth+1];
        end
    end
    
%calculate the centroid location and particle diameter using the 4-point
%or continuous 4-point approximation method (fourptgausfit)
%returns    *x_centroid - x (column) index location of the particles IWC
%           *y_centroid - y (row) index location of the particles IWC
%           *diameter   - calculated diameter of the particle
%           *I0         - calculated maximum particle intensity
elseif sizeprops.method==3 || sizeprops.method==4
    for i=1:num_p
        [x_cg,y_cg,D,I0,E,Meth]=Gaussfit(mapint{i},sizeprops.method,sizeprops.sigma);
        x_c = locxy(i,1)+x_cg-1;
        y_c = locxy(i,2)+y_cg-1;

        if nnz(isnan([x_c,y_c,D,I0])>1)
            if sizeprops.errors==0
                %retain the IWC estimate
            else
                particleprops(i,:)=[x_c,y_c,max(D),max(I0),i,Meth];
            end
        else 
            particleprops(i,:)=[x_c,y_c,max(D),max(I0),i,Meth];
        end
    end

%calculate the centroid location and particle diameter using the local
%least squares or continuous local least squares method (localleastsquares)
%returns    *x_centroid - x (column) index location of the particles IWC
%           *y_centroid - y (row) index location of the particles IWC
%           *diameter   - calculated diameter of the particle
%           *I0         - calculated maximum particle intensity
elseif sizeprops.method==5 || sizeprops.method==6
    for i=1:num_p
        [x_cl,y_cl,D_l,I0_l,E_l,Meth] = Leastsqrfit(mapint{i},sizeprops.method-2,sizeprops.sigma);
        x_c = locxy(i,1)+x_cl-1;
        y_c = locxy(i,2)+y_cl-1;
        
        if nnz(isnan([x_c,y_c,D_l,I0_l]))>1
            if sizeprops.errors==0
                %retain the IWC estimate
            else
                particleprops(i,:)=[x_c,y_c,max(D_l),max(I0_l),i,Meth+2];
            end
        else
            particleprops(i,:)=[x_c,y_c,max(D_l),max(I0_l),i,Meth+2];
        end
    end
    fprintf('\n\t%3.2f%% of particles sized with method %0.0f ',(sum(particleprops(:,6)==sizeprops.method)/num_p)*100,sizeprops.method)
end

XYDiameter=particleprops;

%section for saving sized particle information
if ~isempty(sizeprops.save_dir)
    if exist(sizeprops.save_dir,'dir')~=7
        fprintf('Making Save Directory %s \n',sizeprops.save_dir)
        mkdir(sizeprops.save_dir)
    end
    sname = sprintf('%%s%%0%0.0fd',sizeprops.Data.imzeros);
    save(fullfile(sizeprops.save_dir,sprintf(sname,sizeprops.s_name,sizeprops.s_num)),...
        'sizeprops','XYDiameter','mapsizeinfo','locxy');
end

% if ~isempty(particleIDprops.save_dir)
%     if exist(particleIDprops.save_dir,'dir')~=7
%         fprintf('Making Save Directory %s \t',particleIDprops.save_dir)
%         mkdir(particleIDprops.save_dir)
%     end
%     sname = sprintf('%%0%0.0fd',particleIDprops.Data.imzeros);
%     save(fullfile(particleIDprops.save_dir,[particleIDprops.s_name sprintf(sname,particleIDprops.s_num)]),...
%         'particleIDprops','p_matrix','peaks','num_p');
% end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUB-FUNCTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [mapint,locxy,num_p,mapsize]=mapparticles_v3(im,p_matrix,num_p,p_area)
%
%[mapint,locxy,num_p,mapsize]=mapparticles_v3(im,p_matrix,num_p,p_area)
%
% This function segments individual sections of an image using the original
% image as well as a pre-constructed label matrix.  The function also
% employs a spatial filter (p_area) which allows the user to only return
% sections above a user-defined "size"
%
%INPUTS
%   im (2D array,uint8) - original image in intensity units
%   p_matrix (2D array) - matrix of particle identification and extent
%   num_p (num) - number of identified particles in p_matrix
%   p_area (num) - minimum area in pixels to identify a particle 
%
%OUTPUTS
%   mapint - (2D cell array) individual particle intensity profiles
%   locxy - (2 column array) ROW/COL? location associated with the upper 
%       left pixel of the particle's rectangular projection
%   mapsize - (2 column array) defines the size (row col) of the each
%       array in 'mapint'
%   num_p - number of identified particles
%
%N.Cardwell - 8.13.2009 (v3) - totally revamped the entire particle mapping
%   code; updated to utilize more built-in functions

%convert all inputs to single precision
im=single(im);  p_matrix=single(p_matrix);  
num_p=single(num_p);  p_area=single(p_area);

%check each particle to assure that its area is >= 'p_area'; remove if this
%condition is violated and reset the particle index number
% if p_area < 0
%     fprintf('\nWARNING!  Minimum particle identification area is < 0\n');
% elseif p_area > 0
%     STATS=regionprops(p_matrix,'Area');
%     idx=find([STATS.Area] >= p_area);
%     dummy=zeros(size(p_matrix));
%     for i=1:length(idx)
%         dummy=dummy+(p_matrix==idx(i)).*i;
%         %imagesc(dummy,[0 44]);  set(gca,'DataAspectRatio',[1 1 1]);
%         %pause(0.05);
%     end
%     num_p=length(idx);  p_matrix=dummy;
%     clear dummy idx STATS
% end

%determine the bounding box and location for each identified particle
STATS=regionprops(p_matrix,im,'BoundingBox','PixelIdxList','PixelList','Area');

%populate locxy (upper-left pixel of the bounding box)
locxy=zeros(length(STATS),4);
for i=1:length(STATS);  locxy(i,:)=STATS(i,1).BoundingBox(1:4);  end
locxy=ceil(locxy);

%populate mapint (particle confined by the bounding box)
% mapint=cell(1,length(STATS)); 
mapint=cell(1,1); keep=zeros(1,1); c=1;
for i=1:length(STATS)
    %check each particle to assure that its area is >= 'p_area'; remove if this
    %condition is violated and reset the particle index number; Also remove
    %particles that are the size of only 1 pixel.
    if length(STATS(i,1).PixelIdxList)>1 && STATS(i,1).Area >=p_area
        mapint{1,c}=zeros(locxy(i,4),locxy(i,3));
        for j=1:length(STATS(i,1).PixelIdxList)
            pix_loc_row=STATS(i,1).PixelList(j,2)-locxy(i,2)+1;
            pix_loc_col=STATS(i,1).PixelList(j,1)-locxy(i,1)+1;
            mapint{1,c}(pix_loc_row,pix_loc_col)=im(STATS(i,1).PixelIdxList(j,1));
        end
        keep(c) = i;
        c=c+1;
    end
end

if keep(1) == 0
    mapsize = []; locxy = [];
    num_p = 0;
else
    %populate mapsize (size of each array in mapint) and also trim locxy
    mapsize=[locxy(keep,4),locxy(keep,3)];  locxy=locxy(keep,1:2);
    num_p = length(keep);
end
%code for plotting the results of this function 
%COMMENT OUT FOR NORMAL OPERATION
% figure
% for i=1:length(mapint)
%     subplot(1,2,1)
%     imagesc(im,[0 max(max(im))]);  set(gca,'DataAspectRatio',[1 1 1]); hold on
%     rectangle('Position',[locxy(i,:)-0.5,mapsize(i,2),mapsize(i,1)],'EdgeColor','w')
%     subplot(1,2,2)
%     imagesc(mapint{1,i},[0 max(max(im))]);  set(gca,'DataAspectRatio',[1 1 1]);
%     pause(0.5)
% end

end

function [x_centroid,y_centroid,diameter,I0] = centroidfit(mapint_i,locxy_i)

%
%[x_centroid,y_centroid,diameter]=centroidfit(mapint_i,locxy_i)
%
%given an input particle intensity profile, the function computes the x and
%y Intensity Weighted Centroid (IWC) location in index units relative to
%the upper left hand pixel of the entire input image
%
%function also estimates a particle's diameter given an input particle 
%intensity profile(mapint_i), relative locator to the entire image(locxy_i), 
%and the previously derived particle centroid location(x_centroid,y_centroid).
%
%The function extimates the diameter by meauring the total intensity of 
%the particle then calculating the portion of that associated with the 
%actual particle (86.4665%) from light scattering theory.  Starting with 
%the pixel closest to the particle center and working outward, the function 
%calculates a running intensity summation until 86.4665% of the total 
%intensity is reached.
%
%mapint_i - input particle intensity profile
%locxy_i  - index location of the upper left pixel in the particles square
%           projection, used to orient the IWC to the entire image
%
%N. Cardwell - 2.8.2008
%B. Drew - 7.31.2008

%The first section of the code calcualtes the centroid locations, the
%second calculates the particle diameter

%%%%%%%%%%%%%%%%%Section 1%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%determine extents of mapint_i (particle's square projection)
xintsize=size(mapint_i,2);
yintsize=size(mapint_i,1);

x_c=0;
y_c=0;

%calcualte the IWC x-index location
totalintensity=sum(sum(mapint_i));
for k=1:xintsize
    x_index=(k-0.5)+locxy_i(1,1);%+0.5
    x_c=x_c+sum(mapint_i(:,k))*x_index;
end
x_centroid=x_c/totalintensity+0.5;

%calcualte the IWC y-index location
for k=1:yintsize
    y_index=(k-0.5)+locxy_i(1,2);%+0.5
    y_c=y_c+sum(mapint_i(k,:))*y_index;
end
y_centroid=y_c/totalintensity+0.5;

%%%%%%%%%%%%%%%%%Section 2%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%preallocate distance array of pixel locations and distances from the centroid 
dist_cent=zeros(nnz(mapint_i),4);

%use the "regionprops" function to survey the particle intensity profile
bw_mapint_i=bwlabel(mapint_i);
[mapint_i_info]=regionprops(bw_mapint_i,'PixelIdxList','PixelList');

%populate the distance array
try
    % Check to see if mapint_i is all zeros
    if nnz(mapint_i) == 0
        diameter=NaN;
        I0=NaN;
        %keyboard
        return
    else
        %This loop will not run if size(mapint_i_info,1)>1 (ie - more than one
        %particle in mapint_i)
        for k=1:nnz(mapint_i);
            dist_cent(k,2)=mapint_i_info(1,1).PixelList(k,1)+(locxy_i(1,1)+0.5);%+1.5
            dist_cent(k,3)=mapint_i_info(1,1).PixelList(k,2)+(locxy_i(1,2)+0.5);%+1.5
            dist_cent(k,1)=(dist_cent(k,2)-x_centroid)^2+(dist_cent(k,3)-y_centroid)^2;
            dist_cent(k,4)=mapint_i_info(1,1).PixelIdxList(k);
        end
    end
catch
    %Return NaN's for the diameter and intensity and return to 'detectandsize'
    %    x_centroid=locxy_i(1);
    %    y_centroid=locxy_i(2);
    diameter=NaN;
    I0=NaN;
    %keyboard
    return
end

%sort the distance array from lowest to highest
sort_dist_cent=sortrows(dist_cent);

%calcuate required intensity amount, 86.4665% of total particle intensity
total_int=sum(sum(mapint_i));
req_int=0.864665*total_int;

%preallocate
diamc=0;
cum_int=0;
tag=1;
pp=0;

%sum closest pixel intensities until req_int is reached, also taken into
%account the "fraction pixel amount" needed to reach req_int(line 55)
while tag
    pp=pp+1;
    pixel_int=mapint_i(sort_dist_cent(pp,4));
    cum_int=cum_int + pixel_int;
    if cum_int<req_int
        diamc=diamc+1;
    else
        diamc=diamc+(req_int-(cum_int-pixel_int))/pixel_int;
        tag=0;
    end
end

%Uses the area of a circle to approximate the calculate the particle
%diameter.  If you particles are non-spherical, change this equation to
%ellipse/rectangle/whatever is most appropriate.  "regionprops" can calulate 
%other parameters to help with non=spherical particle diameter calculation, 
%such as orientation and the ratio of major and minor axis
diameter=2*sqrt(diamc/pi);

%The maximum intensity is currently set to the intensity of the brightest
%pixel. This could be modified to calculate the actual intensity if
%desired.
I0=max(max(mapint_i));
end

function [x_c,y_c,D,P,E,Meth] = Gaussfit(intmap,method,sigma)
%
%[x_c,y_c,D,P,E] = Gaussfit(intmap,method,sigma)
%
%this function estimates a particle's diameter given an input particle
%intensity profile(intmap)
%
%The function calculated the particle centroid and diameter by assuming the
%light scattering profile (particle intensity profile) is Gaussian in
%shape.  Four points are chosen in both the x&y dimension which includes
%the maximum intensity pixel of the particle profile.
%
%intmap   - input particle intensity profile
%method   - switch: = 1 for standard 3-pt or = 2 for standard 4-pt or = 3 for continuous 4-pt
%sigma    - number of standard deviations in one diameter
%
%x_c      - X Centroid
%y_c      - Y Centroid
%D        - Diamter, for 3pt method it will return 2 Diameters (X and Y)
%P        - Peak Value (again for 3pt, 2 numbers are returned
%E        - Eccentricity of the 3pt fit that is always run first

%B.Drew  - 7.31.2008   (4-pt Gaussian method taken from the Master's
%                        Thesis of M. Brady)
%S.Raben - 9.20.2008

if method <1 || method > 4
    error('Unknown Method in Gaussfit Function')
end
Meth = method;
%Check to make sure intmap has at least four nonzero points
if method == 1
    if nnz(intmap)<3
        x_c = 1;
        y_c = 1;
        D   = NaN;
        P   = NaN;
        E   = NaN;
        fprintf('Not Enough Points for a Three Point Gauss Fit\n')
        return
    end
else 
    if nnz(intmap)<4
        x_c = 1;
        y_c = 1;
        D   = NaN;
        P   = NaN;
        E   = NaN;
        fprintf('Not Enough Points for a Four Point Gauss Fit\n')
        return
    end
end

%find all of the max intensity locations
[r c]=find(max(intmap(:)) == intmap);


[x_ct,y_ct,Dt,Pt,E] = threeptgaussfit(intmap,[r c],sigma);
Meth = method;
if sum(isnan(Dt)) > 0
    method = 0;
    Meth = method;
end

if method > 2
    [x_c,y_c,D,P,Meth] = fourptgaussfit(intmap,[r c],sigma,method);

    if sum(isnan([x_c,y_c,D,P]))>0
        x_c = x_ct;
        y_c = y_ct;
        D = Dt;
        P = Pt;
        Meth = 1;
    end
else
    x_c = x_ct;
    y_c = y_ct;
    D = Dt;
    P = Pt;
end

% Three Point Method
    function [x_centroid,y_centroid,diameter,I0,eccentricity] = threeptgaussfit(I,locxy,sigma)
        %
        %[x_centroid,y_centroid,diameter]=threeptgausfit(I,locxy_i, sigma)
        %
        %this function estimates a particle's diameter given an input particle
        %intensity profile(I), relative locator to the entire image(locxy_i),
        %and previously derived particle centroid location(xcent,ycent)
        %
        %The function calculated the particle centroid and diameter by assuming the
        %light scattering profile (particle intensity profile) is Gaussian in
        %shape.  Three points are chosen in both the x&y dimension which includes
        %the maximum intensity pixel of the particle profile.
        %
        %I        - input intensity profile
        %locxy    - location on max(iums) intensities
        %sigma    - number of standard deviations in one diameter
        %
        %N.Cardwell - 2.27.2008 (3-pt Gaussian method taken from the Master's Thesis of M. Brady)
        %B.Drew     - 7.31.2008
        %S.Raben    - 9.20.08

        % Find extents of I
        [ymax_I xmax_I] = size(I);

        % Find the maxium value of I and all of its locations
        [max_int,int_L] = max(I(:));%#ok
        %        max_row = locxy(:,1);
        %        max_col = locxy(:,2);
        max_locxy(1) = round(median(locxy(:,1)));
        max_locxy(2) = round(median(locxy(:,2)));
        diameter_x = 0;
        diameter_y = 0;
        
        Better = 0;% if =0 then NO pixel selection just left 1 center and right 1.
        %This method cannot fit a Gaussian shape to a particle intensity profile
        %where the highest intensity pixel is on the edge of I (assumes a
        %Gaussian light scattering profile)
        if Better ~= 1
            LE = [max_locxy(1)   max_locxy(2)-1];
            RE = [max_locxy(1)   max_locxy(2)+1];
            TE = [max_locxy(1)-1 max_locxy(2)];
            BE = [max_locxy(1)+1 max_locxy(2)];

            %Given that this method assumes a Gaussian shaped intensity profile, the
            %function is designed to return NaN should any of the nearest neighbor
            %pixels to max_int equal zero or if the particles have equal intensity
            %values
            if LE(2) > 1 && RE(2) < xmax_I %edge_ckx == 0;
                left_int  = I(LE(1),LE(2));
                right_int = I(RE(1),RE(2));
                LR_denom = (2*(log(left_int)+log(right_int)-2*log(max_int)));
                if left_int==0 || right_int==0 || LR_denom==0
                    x_centroid = max_locxy(2);
                    diameter_x=NaN;
                    I0x=NaN;
                else
                    %calculate the x&y centroid location using the 3-pt Gaussian Method
                    x_centroid=max_locxy(2)+log(left_int/right_int)/LR_denom;
                end
            else
                x_centroid = max_locxy(2);
                diameter_x=NaN;
                I0x=NaN;
            end

            if TE(1) > 1 && BE(1) < ymax_I %edge_cky == 0;
                top_int    = I(TE(1),TE(2));
                bottom_int = I(BE(1),BE(2));
                TB_denom = (2*(log(top_int)+log(bottom_int)-2*log(max_int)));
                if top_int==0 || bottom_int==0 || TB_denom==0
                    y_centroid = max_locxy(1);
                    diameter_y=NaN;
                    I0y=NaN;
                else
                    %calculate the x&y centroid location using the 3-pt Gaussian Method
                    y_centroid=max_locxy(1)+log(top_int/bottom_int)/TB_denom;
                end
            else
                y_centroid = max_locxy(1);
                diameter_y=NaN;
                I0y=NaN;
            end

            %The diameter of the particle is directly related to the variance (where
            %Beta=1/(variance)^2 and diameter=sigma*variance.  Logic statements deals with
            %the possibility that one pixel adjacent to max_int is equal to max_int
            %while the other adjacent pixel does not
            if ~isnan(diameter_x)
                betax=abs(log(left_int/max_int)/((LE(2)-x_centroid)^2-(max_locxy(2)-x_centroid)^2));
                I0x=left_int/exp(-betax*(LE(2)-x_centroid)^2);
                diameter_x=sigma./sqrt((2*betax));
            end

            if ~isnan(diameter_y)
                betay=abs(log(top_int/max_int)/((TE(1)-y_centroid)^2-(max_locxy(1)-y_centroid)^2));
                I0y=top_int/exp(-betay*(TE(1)-y_centroid)^2);
                diameter_y=sigma./sqrt((2*betay));
            end

            % General Solution to the 3pt Gausian Estimator
        else

            [TE,BE,LE,RE] = Edgesearch(I,max_locxy);

            if numel(locxy(:,2)) > 1
                if TE(1) > 1 && BE(1) < ymax_I && I(TE(1),TE(2)) ~= 0 && I(BE(1),BE(2)) ~= 0
                    if I(TE(1),TE(2)) > I(BE(1),BE(2))
                        [TBLoc] = find(I(1:TE(1),TE(2)) < I(TE(1),TE(2)));%#ok matlab is just being stupid
                        maxlocTB = [max(TBLoc) TE(2)];
                    else
                        [TBLoc] = find(I(BE(1):end,BE(2)) < I(BE(1),BE(2)));
                        maxlocTB = [min(BE(1)-1+TBLoc) BE(2)];
                    end
                else
                    maxlocTB = max_locxy;
                end

                if LE(2) > 1 && RE(2) < xmax_I && I(LE(1),LE(2)) ~= 0 && I(RE(1),RE(2)) ~= 0
                    if I(LE(1),LE(2)) > I(RE(1),RE(2))
                        [LRLoc] = find(I(LE(1),1:LE(2)) < I(LE(1),LE(2)));%#ok matlab is just being stupid
                        maxlocLR = [LE(1) max(LRLoc)];
                    else
                        [LRLoc] = find(I(RE(1),RE(2):end) < I(RE(1),RE(2)));
                        maxlocLR = [RE(1) min(RE(2)-1+LRLoc)];
                    end
                else
                    maxlocLR = max_locxy;
                end

            else
                maxlocTB   = max_locxy;
                maxlocLR   = max_locxy;
            end
            try
                points_LR = [LE(1) LE(2) I(LE(1),LE(2)); RE(1) RE(2) I(RE(1),RE(2)); maxlocLR(1) maxlocLR(2) I(maxlocLR(1),maxlocLR(2))];
                points_TB = [TE(1) TE(2) I(TE(1),TE(2)); BE(1) BE(2) I(BE(1),BE(2)); maxlocTB(1) maxlocTB(2) I(maxlocTB(1),maxlocTB(2))];
            catch
                %keyboard
            end

            [LRD,LRI] = sort(points_LR(:,3),'ascend');
            [TBD,TBI] = sort(points_TB(:,3),'ascend');
            points_LR = points_LR(LRI,:);
            points_TB = points_TB(TBI,:);
            %Given that this method assumes a Gaussian shaped intensity profile, the
            %function is designed to return NaN should any of the nearest neighbor
            %pixels to max_int equal zero or if the particles have equal intensity
            %values
            if LE(2) ~= 1 && RE(2) ~= xmax_I %edge_ckx == 0;
                LR_denom = 2*((-points_LR(2,2) + points_LR(3,2))*log(points_LR(1,3))...
                    + (points_LR(1,2) - points_LR(3,2))*log(points_LR(2,3)) + (-points_LR(1,2) + points_LR(2,2))*log(points_LR(3,3)));
                if points_LR(1,3)==0 || points_LR(3,3)==0 || LR_denom==0
                    x_centroid = max_locxy(2);
                    diameter_x=NaN;
                    I0x=NaN;
                else
                    %calculate the x&y centroid location using the 3-pt Gaussian Method
                    %x_centroid=max_locxy(2)+log(left_int/right_int)/LR_denom;
                    x_centroid = ((points_LR(1,2)^2 - points_LR(3,2)^2)*log(points_LR(2,3)/points_LR(1,3))...
                        + (-points_LR(1,2)^2 + points_LR(2,2)^2)*log(points_LR(3,3)/points_LR(1,3)))/(LR_denom);
                end
            else
                x_centroid = max_locxy(2);
                diameter_x=NaN;
                I0x=NaN;
            end

            if TE(1) ~= 1 && BE(1) ~= ymax_I %edge_cky == 0;
                TB_denom = 2*((-points_TB(2,1) + points_TB(3,1))*log(points_TB(1,3))...
                    + (points_TB(1,1) - points_TB(3,1))*log(points_TB(2,3)) + (-points_TB(1,1) + points_TB(2,1))*log(points_TB(3,3)));
                if points_TB(1,3)==0 || points_TB(3,3)==0 || TB_denom==0
                    y_centroid = max_locxy(1);
                    %keyboard
                    diameter_y=NaN;
                    I0y=NaN;
                else
                    %calculate the x&y centroid location using the 3-pt Gaussian Method
                    %y_centroid=max_locxy(1)+log(top_int/bottom_int)/TB_denom;
                    y_centroid = ((points_TB(1,1)^2 - points_TB(3,1)^2)*log(points_TB(2,3)/points_TB(1,3))...
                        + (-points_TB(1,1)^2 + points_TB(2,1)^2)*log(points_TB(3,3)/points_TB(1,3)))/(TB_denom);
                end
            else
                y_centroid = max_locxy(1);
                diameter_y=NaN;
                I0y=NaN;
            end

            %The diameter of the particle is directly related to the variance (where
            %Beta=1/(variance)^2 and diameter=sigma*variance.  Logic statements deals with
            %the possibility that one pixel adjacent to max_int is equal to max_int
            %while the other adjacent pixel does not
            if ~isnan(diameter_x)
                %                     betax=abs(log(left_int/max_int)/((max_locxy(2)-x_centroid)^2-(LE(2)-x_centroid)^2));
                %                     I0x=left_int/exp(-betax*(LE(2)-x_centroid)^2);
                sigma_x = sqrt(abs(((-points_LR(1,2)+points_LR(2,2))*(points_LR(1,2)-points_LR(3,2))...
                    *(points_LR(2,2)-points_LR(3,2))/LR_denom)));
                I0x = points_LR(3,3)/exp((-(points_LR(3,2)-x_centroid)^2)/(2*sigma_x^2));
                diameter_x=sigma.*sigma_x;
            end

            if ~isnan(diameter_y)
                %                     betay=abs(log(top_int/max_int)/((max_locxy(1)-y_centroid)^2-(TE(1)-y_centroid)^2));
                %                     I0y=top_int/exp(-betay*(TE(1)-y_centroid)^2);
                sigma_y = sqrt(abs(((-points_TB(1,1)+points_TB(2,1))*(points_TB(1,1)-points_TB(3,1))...
                    *(points_TB(2,1)-points_TB(3,1))/TB_denom)));
                I0y = points_TB(3,3)/exp((-(points_TB(3,1)-y_centroid)^2)/(2*sigma_y^2));
                diameter_y=sigma.*sigma_y;
            end
        end

        %Check to make sure the diameter is not more than 1.5 times the size of the
        %interrogation window
        if diameter_x > 1.5*xmax_I
            diameter_x=NaN;
            I0x=NaN;
        end
        if diameter_y > 1.5*ymax_I
            diameter_y=NaN;
            I0y=NaN;
        end

        %Currently the diameter is approximated by a sphere of diameter equal to
        %the max size in the x/y direction.  For non-spherical particles, this
        %equation should be altered to an ellipse using diameterx and diametery as
        %the major and minor axis
        diameter=[diameter_x,diameter_y];
        I0=[I0x,I0y];
        eccentricity(1) = sqrt(max(diameter).^2 - min(diameter).^2)./max(diameter);
        eccentricity(2) = sqrt(max(diameter).^2 - min(diameter).^2)./min(diameter);

        %Validation to make sure the results are not negative
        if min([x_centroid,y_centroid,diameter,I0])<=0
            x_centroid=max_locxy(2);
            y_centroid=max_locxy(1);
            diameter=NaN;
            I0=NaN;
            return
        end
    end

% Four Point Method
    function [x_centroid,y_centroid,diameter,I0,Meth] = fourptgaussfit(I,locxy,sigma,method)
        %
        %[x_centroid,y_centroid,diameter]=fourptgausfit(mapint,locxy_i,method,sigma)
        %
        %this function estimates a particle's diameter given an input particle
        %intensity profile(mapint), relative locator to the entire image(locxy_i),
        %
        %The function calculated the particle centroid and diameter by assuming the
        %light scattering profile (particle intensity profile) is Gaussian in
        %shape.  Four points are chosen in both the x&y dimension which includes
        %the maximum intensity pixel of the particle profile.
        %
        %I        - input particle intensity profile
        %locxy    - index location of the upper left pixel in the particles square
        %           projection, used to orient the IWC to the entire image
        %method   - switch: =2 for standard 4-pt or =3 for continuous 4-pt
        %sigma    - number of standard deviations in one diameter
        %
        %B.Drew  - 7.31.2008   (4-pt Gaussian method taken from the Master's
        %                        Thesis of M. Brady)
        %S.Raben - 9.20.2008


        % Find extents of I
        [ymax_I xmax_I] = size(I);
        Meth = method;
        % Find the maxium value of I and all of its locations
        [max_int,int_L] = max(I(:));%#ok
        max_int_row = locxy(:,1);
        max_int_col = locxy(:,2);
        max_int_locxy(1) = round(median(locxy(:,1)));
        max_int_locxy(2) = round(median(locxy(:,2)));

        %Pick the four pixels to be used, starting with max_int
        points=zeros(4,3);
        points(1,:)=[max_int_locxy,max_int];

%         if numel(locxy(:,1)) > 1
%             max_try = max(I(I~=max_int));
%         else
%             max_try = max_int;
%         end
        sat_int = max_int;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        I2 = zeros(size(I)+2);
        I2(2:end-1,2:end-1) = I;
        I = I2;
        locxy = locxy + 1;
        xmax_I = xmax_I+2;
        ymax_I = ymax_I+2;
        points(1,1:2) = points(1,1:2) + 1;
        max_int_locxy(1) = max_int_locxy(1) + 1;
        max_int_locxy(2) = max_int_locxy(2) + 1;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        if min(locxy(:,1)) > 1 && max(locxy(:,1)) < ymax_I-2 && ...
                min(locxy(:,2)) > 1 && max(locxy(:,2)) < xmax_I-2
            %Try different combinations of points to find a suitable set
            %A suitable set of points is one where a circle cannot be drawn
            %which intersects all four, and no pixels are empty or
            %saturated

            [TE,BE,LE,RE] = Edgesearch(I,max_int_locxy);
            if numel(locxy(:,1))>1
                if I(TE(1)-1,TE(2))~=0
                    [pts] = circlecheck([TE;[TE(1)-1 TE(2)];LE;RE],I);
                elseif I(BE(1)+1,BE(2))~=0
                    [pts] = circlecheck([BE;[BE(1)+1 BE(2)];LE;RE],I);
                elseif I(LE(1),LE(2)-1)~=0
                    [pts] = circlecheck([TE;[LE(1) LE(2)-1];LE;RE],I);
                elseif I(RE(1),RE(2)+1)~=0
                    [pts] = circlecheck([TE;[RE(1) RE(2)+1];LE;RE],I);
                else
                    %keyboard
                end
            else
                pts = [TE;locxy;BE;LE];
            end
            points(:,1:2) = pts(:,2:-1:1);

            for pp = 1:4
                points(pp,3)   = I(pts(pp,1),pts(pp,2));
            end

            if sum(diff(points(:,1))) == 0 || sum(diff(points(:,2))) == 0
                error('Points are in a Staight line')
            end

        else
            %keyboard
            %If max_int is on the edge of the window, then try these point combinations
            if max_int_row==1 && max_int_col==1
                points(2,:)=[points(1,1)+1 , points(1,2)   , I(max_int_row,max_int_col+1)];
                points(3,:)=[points(1,1)   , points(1,2)+1 , I(max_int_row+1,max_int_col)];
                points(4,:)=[points(1,1)   , points(1,2)+2 , I(max_int_row+2,max_int_col)];
            elseif max_int_row==1 && max_int_col==xmax_I
                points(2,:)=[points(1,1)-1 , points(1,2)   , I(max_int_row,max_int_col-1)];
                points(3,:)=[points(1,1)   , points(1,2)+1 , I(max_int_row+1,max_int_col)];
                points(4,:)=[points(1,1)   , points(1,2)+2 , I(max_int_row+2,max_int_col)];
            elseif max_int_row==ymax_I && max_int_col==1
                points(2,:)=[points(1,1)+1 , points(1,2)   , I(max_int_row,max_int_col+1)];
                points(3,:)=[points(1,1)   , points(1,2)-1 , I(max_int_row-1,max_int_col)];
                points(4,:)=[points(1,1)   , points(1,2)-2 , I(max_int_row-2,max_int_col)];
            elseif max_int_row==ymax_I && max_int_col==xmax_I
                points(2,:)=[points(1,1)-1 , points(1,2)   , I(max_int_row,max_int_col-1)];
                points(3,:)=[points(1,1)   , points(1,2)-1 , I(max_int_row-1,max_int_col)];
                points(4,:)=[points(1,1)   , points(1,2)-2 , I(max_int_row-2,max_int_col)];
            elseif max_int_row==1
                points(2,:)=[points(1,1)+1 , points(1,2)   , I(max_int_row,max_int_col+1)];
                points(3,:)=[points(1,1)-1 , points(1,2)   , I(max_int_row,max_int_col-1)];
                points(4,:)=[points(1,1)   , points(1,2)+1 , I(max_int_row+1,max_int_col)];
            elseif max_int_row==ymax_I
                points(2,:)=[points(1,1)+1 , points(1,2)   , I(max_int_row,max_int_col+1)];
                points(3,:)=[points(1,1)-1 , points(1,2)   , I(max_int_row,max_int_col-1)];
                points(4,:)=[points(1,1)   , points(1,2)-1 , I(max_int_row-1,max_int_col)];
            elseif max_int_col==xmax_I
                points(2,:)=[points(1,1)-1 , points(1,2)   , I(max_int_row,max_int_col-1)];
                points(3,:)=[points(1,1)   , points(1,2)+1 , I(max_int_row+1,max_int_col)];
                points(4,:)=[points(1,1)   , points(1,2)-1 , I(max_int_row-1,max_int_col)];
            end
        end%end

        %If a suitable set of points could not be found:
        if numel(max(points(:,3))==sat_int)>1 || min(points(:,3))==0
            x_centroid=locxy(1);
            y_centroid=locxy(2);
            %keyboard
            diameter=NaN;
            I0=NaN;
            return
        end

        %This code is specific to the standard four point gaussian estimator, but
        %will be used to come up with guess values for the continuous 4-point
        %method

        [Isort,IsortI] = sort(points(:,3),'descend');
        points = points(IsortI,:);

        x1=points(1,1);
        x2=points(2,1);
        x3=points(3,1);
        x4=points(4,1);
        y1=points(1,2);
        y2=points(2,2);
        y3=points(3,2);
        y4=points(4,2);
        a1=points(1,3);
        a2=points(2,3);
        a3=points(3,3);
        a4=points(4,3);

        alpha(1) = (x4^2)*(y2 - y3) + (x3^2)*(y4 - y2) + ((x2^2) + (y2 - y3)*(y2 - y4))*(y3 - y4);
        alpha(2) = (x4^2)*(y3 - y1) + (x3^2)*(y1 - y4) - ((x1^2) + (y1 - y3)*(y1 - y4))*(y3 - y4);
        alpha(3) = (x4^2)*(y1 - y2) + (x2^2)*(y4 - y1) + ((x1^2) + (y1 - y2)*(y1 - y4))*(y2 - y4);
        alpha(4) = (x3^2)*(y2 - y1) + (x2^2)*(y1 - y3) - ((x1^2) + (y1 - y2)*(y1 - y3))*(y2 - y3);

        gamma(1) = (-x3^2)*x4 + (x2^2)*(x4 - x3) + x4*((y2^2) - (y3^2)) + x3*((x4^2) - (y2^2) + (y4^2)) + x2*(( x3^2) - (x4^2) + (y3^2) - (y4^2));
        gamma(2) = ( x3^2)*x4 + (x1^2)*(x3 - x4) + x4*((y3^2) - (y1^2)) - x3*((x4^2) - (y1^2) + (y4^2)) + x1*((-x3^2) + (x4^2) - (y3^2) + (y4^2));
        gamma(3) = (-x2^2)*x4 + (x1^2)*(x4 - x2) + x4*((y1^2) - (y2^2)) + x2*((x4^2) - (y1^2) + (y4^2)) + x1*(( x2^2) - (x4^2) + (y2^2) - (y4^2));
        gamma(4) = ( x2^2)*x3 + (x1^2)*(x2 - x3) + x3*((y2^2) - (y1^2)) - x2*((x3^2) - (y1^2) + (y3^2)) + x1*((-x2^2) + (x3^2) - (y2^2) + (y3^2));

        delta(1) = x4*(y2 - y3) + x2*(y3 - y4) + x3*(y4 - y2);
        delta(2) = x4*(y3 - y1) + x3*(y1 - y4) + x1*(y4 - y3);
        delta(3) = x4*(y1 - y2) + x1*(y2 - y4) + x2*(y4 - y1);
        delta(4) = x3*(y2 - y1) + x2*(y1 - y3) + x1*(y3 - y2);

        deno = 2*(log(a1)*delta(1) + log(a2)*delta(2) + log(a3)*delta(3) + log(a4)*delta(4));

        x_centroid = (log(a1)*alpha(1) + log(a2)*alpha(2) + log(a3)*alpha(3) + log(a4)*alpha(4))/deno;

        y_centroid = (log(a1)*gamma(1) + log(a2)*gamma(2) + log(a3)*gamma(3) + log(a4)*gamma(4))/deno;

        if a2 ~= a1
            betas = abs((log(a2)-log(a1))/((x2-x_centroid)^2+(y2-y_centroid)^2-(x1-x_centroid)^2-(y1-y_centroid)^2));
        elseif a3 ~= a1
            betas = abs((log(a3)-log(a1))/((x3-x_centroid)^2+(y3-y_centroid)^2-(x1-x_centroid)^2-(y1-y_centroid)^2));
        elseif a4 ~= a1
            betas = abs((log(a4)-log(a1))/((x4-x_centroid)^2+(y4-y_centroid)^2-(x1-x_centroid)^2-(y1-y_centroid)^2));
        else
            %keyboard
            I0 = NaN;
            diameter = NaN;
            return
        end

        I0=a1/exp(-betas*((x1-x_centroid)^2+(y1-y_centroid)^2));
        
        if sum(isnan([x_centroid,y_centroid,betas])) > 0
            %keyboard
        end
        %Check solutions for errors - set guess values if the 4-point
        %continuous method is being used
        if (min([x_centroid,y_centroid,betas])<=0 || max([x_centroid,y_centroid,betas])>=10^4 ...
                || max(isnan([x_centroid,y_centroid,betas]))~=0) && method==4
            %keyboard
            x_centroid=points(1,1);
            y_centroid=points(1,2);
            betas=.25;
            I0=points(1,3);
        elseif min([x_centroid,y_centroid,betas])<=0
            %keyboard
            x_centroid=locxy(1);
            y_centroid=locxy(2);
            diameter=NaN;
            I0=NaN;
            return
        elseif x_centroid > xmax_I || y_centroid > ymax_I || betas >= 10^3%max([x_centroid,y_centroid,betas])>=10^3
            %keyboard
            x_centroid=locxy(1);
            y_centroid=locxy(2);
            diameter=NaN;
            I0=NaN;
            return
        end

        %This code is specific to the continuous four point gaussian estimator:
        if method==4

            %Guess value for fsolve
            x0=[I0,betas,x_centroid,y_centroid];

            %Convert row/column measurements to center-of-pixel measurements
            points(:,1)=points(:,1)+0.5;
            points(:,2)=points(:,2)+0.5;
            options=optimset('MaxIter',400,'MaxFunEvals',2000,'LargeScale','off','Display','off');
            try
                %Call the fsolve function to find solutions to the equation in the
                %Powell_optimized_eq function
                [xvar,fval,exitflag]=fsolve(@Powell_optimized_eq,x0,options,points);%#ok
                I0=xvar(1);
                betas=xvar(2);
                x_centroid=xvar(3)-1;
                y_centroid=xvar(4)-1;
            catch %#ok
                %keyboard
                x_centroid=locxy(1);
                y_centroid=locxy(2);
                diameter=NaN;%#ok
                I0=NaN;
            end
        end

        %Calculate diameter
        diameter=sqrt(sigma^2/(2*betas));
        x_centroid = x_centroid - 1;
        y_centroid = y_centroid - 1;

        % Powel Optimized
        function F = Powell_optimized_eq(x,points)
            %This function is called by fsolve when the continuous four-pt method is
            %selected. X is a matrix containing the initial guesses [I0, betas, x_c,
            %y_c]. Points is a matrix containing the four selected points [x1 y1 I1;
            %...;x4 y4 I4]
            %
            %Adapted from M. Brady's 'fourptintgaussfit'
            %B.Drew - 7.18.2008

            x1=points(1,1);
            x2=points(2,1);
            x3=points(3,1);
            x4=points(4,1);
            y1=points(1,2);
            y2=points(2,2);
            y3=points(3,2);
            y4=points(4,2);
            a1=points(1,3);
            a2=points(2,3);
            a3=points(3,3);
            a4=points(4,3);

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %For unknown reasons, the fsolve function tries negative values of x(2),
            %causing these equations to return an error. Putting the
            %abs() function in front of all the x(2)'s fixes this error, and fsolve's
            %final x(2) value ends up positive.
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            F = [pi*x(1)/4/x(2)*((erf(sqrt(abs(x(2)))*(x1-x(3)))-erf(sqrt(abs(x(2)))*(x1+1-x(3))))*(erf(sqrt(abs(x(2)))*(y1-x(4)))-erf(sqrt(abs(x(2)))*(y1+1-x(4)))))-a1;
                pi*x(1)/4/x(2)*((erf(sqrt(abs(x(2)))*(x2-x(3)))-erf(sqrt(abs(x(2)))*(x2+1-x(3))))*(erf(sqrt(abs(x(2)))*(y2-x(4)))-erf(sqrt(abs(x(2)))*(y2+1-x(4)))))-a2;
                pi*x(1)/4/x(2)*((erf(sqrt(abs(x(2)))*(x3-x(3)))-erf(sqrt(abs(x(2)))*(x3+1-x(3))))*(erf(sqrt(abs(x(2)))*(y3-x(4)))-erf(sqrt(abs(x(2)))*(y3+1-x(4)))))-a3;
                pi*x(1)/4/x(2)*((erf(sqrt(abs(x(2)))*(x4-x(3)))-erf(sqrt(abs(x(2)))*(x4+1-x(3))))*(erf(sqrt(abs(x(2)))*(y4-x(4)))-erf(sqrt(abs(x(2)))*(y4+1-x(4)))))-a4];

        end
    end

% Find Unstaturated Edges
    function [Top,Bottom,Left,Right] = Edgesearch(Im,locmax)


        [max_r max_c] = find(max(Im(:))==Im);
        maxI          = max(Im(:));
        [ymax xmax]   = size(Im);

        if min(max_r)==1
            stop_ty = 1;
            Top = [round(mean(max_r)) round(mean(max_c))];
        else
            stop_ty = 0;
        end

        if max(max_r)==ymax
            stop_by = 1;
            Bottom = [round(mean(max_r)) round(mean(max_c))];
        else
            stop_by = 0;
        end

        if min(max_c)==1
            stop_lx = 1;
            Left = [round(mean(max_r)) round(mean(max_c))];
        else
            stop_lx = 0;
        end

        if max(max_c)==xmax
            stop_rx = 1;
            Right = [round(mean(max_r)) round(mean(max_c))];
        else
            stop_rx = 0;
        end

        %2 Logical loops to deal with saturated or equvalent pixels values surrounding
        %max_int.  Incrementally moves outward until suitable values are found or
        %the extents of I are reached
        %stop_x=0;
        %stop_y=0;
        ii=1;
        extent_left   = min(max_c);
        extent_right  = max(max_c);
        extent_top    = min(max_r);
        extent_bottom = max(max_r);

        %search in the x-dimension
        while stop_lx==0 && extent_left>=1
            %locate x-index and intensity values for the nearest neighbors of max_int
            left_locxy  = [locmax(1) locmax(2)-ii];
            left_int    = Im(left_locxy(1),left_locxy(2));
            %logical check to assure that the three pixels chosen are not equal
            if left_int==maxI
                ii=ii+1;
            else
                stop_lx=1;%#ok
            end

            %make sure there are data points to the left & right
            Left  = left_locxy;%#ok
        end
        ii = 1;
        while stop_rx==0 && extent_right<=xmax
            %locate x-index and intensity values for the nearest neighbors of max_int
            right_locxy = [locmax(1) locmax(2)+ii];
            right_int   = Im(right_locxy(1),right_locxy(2));
            %logical check to assure that the three pixels chosen are not equal
            if right_int==maxI
                ii=ii+1;
            else
                stop_rx=1;%#ok
            end

            %make sure there are data points to the left & right
            Right = right_locxy;%max(max_c)+i;%#ok
        end
        ii = 1;
        %search in the y-dimension
        while stop_ty==0 && extent_top>=1
            %locate y-index and intensity values for the nearest neighbors of max_int
            top_locxy    = [locmax(1)-ii locmax(2)];
            top_int      = Im(top_locxy(1),top_locxy(2));
            %logical check to assure that the three pixels chosen are not equal
            if top_int==maxI
                ii=ii+1;%#ok
            else
                stop_ty=1;%#ok
            end

            %make sure there are data points to the top & bottom
            Top    = top_locxy;%min(max_r)-j;
        end
        ii = 1;
        while stop_by==0 && extent_bottom<=ymax
            %locate y-index and intensity values for the nearest neighbors of max_int
            bottom_locxy = [locmax(1)+ii locmax(2)];
            bottom_int   = Im(bottom_locxy(1),bottom_locxy(2));
            %logical check to assure that the three pixels chosen are not equal
            if bottom_int==maxI
                ii=ii+1;%#ok
            else
                stop_by=1;%#ok
            end

            %make sure there are data points to the top & bottom
            Bottom = bottom_locxy;%max(max_r)+j;
        end
    end
    
% Check to make sure points arn't in a cicle or straight line.
    function [locations] = circlecheck(locations,Icir)
        [ymaxy,xmaxx] = size(Icir);
        C_x = mean(locations(:,2));
        C_y = mean(locations(:,1));
        Cir = sqrt(((locations(:,1)-C_x).^2)+((locations(:,2)-C_y).^2));

        while sum(abs(diff(Cir))) == 0
            if locations(1,1)-1 >= 1 && Icir(locations(1,1)-1,locations(1,2)) ~= 0
                locations(1,1) = locations(1,1)-1;
            elseif locations(2,1)+1 <= ymaxy  && Icir(locations(2,1)+1,locations(2,2)) ~= 0
                locations(2,1) = locations(2,1)+1;
            elseif locations(3,2)-1 >= 1 && Icir(locations(3,1),locations(3,2)-1) ~= 0
                locations(3,2) = locations(3,2)-1;
            elseif locations(4,2)+1 <= xmaxx && Icir(locations(4,1),locations(4,2)+1) ~= 0
                locations(4,2) = locations(4,2)+1;
            else
                %keyboard
            end                
                
            C_x = mean(locations(:,2));
            C_y = mean(locations(:,1));
            Cir = sqrt(((locations(:,2)-C_x).^2)+((locations(:,1)-C_y).^2));
        end
    end

end

function [x_c,y_c,D,P,E,Meth] = Leastsqrfit(I_in,method_in,sigma_in)

%Desicription
%
%
%
%







%Options for the lsqnonlin solver
options=optimset('MaxIter',1200,'MaxFunEvals',5000,'TolX',5e-6,'TolFun',5e-6,...
    'Display','off','DiffMinChange',1e-7,'DiffMaxChange',1);%,'LevenbergMarquardt','off');%,'LargeScale','off');

% Find the center Max and all of the saturated points
[locxy_in(:,1) locxy_in(:,2)] = find(I_in == max(I_in(:)));
max_locxy_in(1) = round(median(locxy_in(:,1)));
max_locxy_in(2) = round(median(locxy_in(:,2)));

% If there are not enough points don't use the method
if nnz(I_in) - numel(find(I_in == max(I_in(:)))) + 1 < 5
    x_c = 1;
    y_c = 1;
    D   = NaN;
    P   = NaN;
    E   = NaN;
    Meth = 0;
%     fprintf('Not Enough Points for Least Squres Fit\n')
    return
end

%Removes negitive values and puts in zeros
Ils = I_in;
Ils(I_in<0) = 0;

[x_c,y_c,D,P,E,Meth] = Leastsqrmethods(Ils,method_in,sigma_in,options,locxy_in,max_locxy_in);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Sub Functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [x_c,y_c,D,P,E,Meth] = Leastsqrmethods(I,method,sigma,options,locxy,max_locxy)

        % Find the maxium value of I and all of its locations
        [max_int,int_L] = max(I(:));%#ok
        max_row = locxy(:,1);
        max_col = locxy(:,2);
        E = 0;
        Meth = method;
        %Try making a 5x5 window around the maximum intensity
        %pixel to reduce the chance of grouping particles together.
        %If a 5x5 window isn't possible due to the particle's location on
        %the image, then use the original window.

        % Find extents of I
        [ymax_I xmax_I] = size(I);

        %This method cannot fit a Gaussian shape to a particle intensity profile
        %where the highest intensity pixel is on the edge of I (assumes a
        %Gaussian light scattering profile)
        if min(max_row)==1
            stop_y = 1;
            edge_cky = 1;            
        elseif max(max_row)==ymax_I
            stop_y = 1;
            edge_cky = 1;            
        else
            stop_y = 0;
            edge_cky = 0;
        end

        if min(max_col)==1
            stop_x = 1;
            edge_ckx = 1;            
        elseif max(max_col)==xmax_I
            stop_x = 1;
            edge_ckx = 1;
        else
            stop_x = 0;
            edge_ckx = 0;            
        end

        %2 Logical loops to deal with saturated or equvalent pixels values surrounding
        %max_int.  Incrementally moves outward until suitable values are found or
        %the extents of I are reached
        i=1;
        j=1;
        extent_left=min(max_col);
        extent_right=max(max_col);
        extent_top=min(max_row);
        extent_bottom=max(max_row);
        max_edge_left = [0 0];
        max_edge_top = [0 0];
        diameter_x = 0;
        diameter_y = 0;
        left_int = max_int;
        left_locxy = max_locxy;
        right_int = max_int;
        right_locxy = max_locxy;
        top_int = max_int;
        top_locxy = max_locxy;
        bottom_int = max_int;
        bottom_locxy = max_locxy;
        
        %search in the x-dimension
        while stop_x==0 && extent_left>=1 && extent_right<=xmax_I
            %locate x-index and intensity values for the nearest neighbors of max_int
            left_locxy  = [max_locxy(1) max_locxy(2)-i];
            left_int    = I(left_locxy(1),left_locxy(2));
            right_locxy = [max_locxy(1) max_locxy(2)+i];
            right_int   = I(right_locxy(1),right_locxy(2));
            %logical check to assure that the three pixels chosen are not equal
            if left_int==max_int || right_int==max_int
                i=i+1;
                max_edge_left = left_locxy;
            else
                stop_x=1;
            end

            %make sure there are data points to the left & right
            extent_left=min(max_col)-i;
            extent_right=max(max_col)+i;
        end

        %search in the y-dimension
        while stop_y==0 && extent_top>=1 && extent_bottom<=ymax_I
            %locate y-index and intensity values for the nearest neighbors of max_int
            top_locxy    = [max_locxy(1)-j max_locxy(2)];
            top_int      = I(top_locxy(1),top_locxy(2));
            bottom_locxy = [max_locxy(1)+j max_locxy(2)];
            bottom_int   = I(bottom_locxy(1),bottom_locxy(2));
            %logical check to assure that the three pixels chosen are not equal
            if top_int==max_int || bottom_int==max_int
                j=j+1;
                max_edge_top = top_locxy;
            else
                stop_y=1;
            end

            %make sure there are data points to the top & bottom
            extent_top=min(max_row)-j;
            extent_bottom=max(max_row)+j;
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
        % 3pt Gauss Fit - Start
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Given that this method assumes a Gaussian shaped intensity profile, the
        %function is designed to return NaN should any of the nearest neighbor
        %pixels to max_int equal zero or if the particles have equal intensity
        %values
        if edge_ckx == 0;
            LR_denom = (2*(log(left_int)+log(right_int)-2*log(max_int)));
            if left_int==0 || right_int==0 || LR_denom==0
                x_centroid = max_locxy(2);
                %keyboard
                diameter_x=NaN;
                I0x=NaN;
            else
                %calculate the x&y centroid location using the 3-pt Gaussian Method
                x_centroid=max_locxy(2)+log(left_int/right_int)/LR_denom;
            end
        else
            x_centroid = max_locxy(2);
            diameter_x=NaN;
            I0x=NaN;
        end

        if edge_cky == 0
            TB_denom = (2*(log(top_int)+log(bottom_int)-2*log(max_int)));
            if top_int==0 || bottom_int==0 || TB_denom==0
                y_centroid = max_locxy(1);
                %keyboard
                diameter_y=NaN;
                I0y=NaN;
            else
                %calculate the x&y centroid location using the 3-pt Gaussian Method
                y_centroid=max_locxy(1)+log(top_int/bottom_int)/TB_denom;
            end
        else
            y_centroid = max_locxy(1);
            diameter_y=NaN;
            I0y=NaN;
        end

        %The diameter of the particle is directly related to the variance (where
        %Beta=1/(variance)^2 and diameter=sigma*variance.  Logic statements deals with
        %the possibility that one pixel adjacent to max_int is equal to max_int
        %while the other adjacent pixel does not
        if ~isnan(diameter_x)
            if nnz(max_edge_left) < 1
                betax=log(left_int/max_int)/((max_locxy(2)-x_centroid)^2-(left_locxy(2)-x_centroid)^2);
                I0x=left_int/exp(-betax*(left_locxy(2)-x_centroid)^2);
            else
                betax=log(left_int/max_int)/((max_edge_left(2)-x_centroid)^2-(left_locxy(2)-x_centroid)^2);
                I0x=left_int/exp(-betax*(left_locxy(2)-x_centroid)^2);
            end
            diameter_x=sigma./sqrt((2*betax));
        end

        if ~isnan(diameter_y)
            if nnz(max_edge_top) < 1
                betay=log(top_int/max_int)/((max_locxy(1)-y_centroid)^2-(top_locxy(1)-y_centroid)^2);
                I0y=top_int/exp(-betay*(top_locxy(1)-y_centroid)^2);
            else
                betay=log(top_int/max_int)/((max_edge_top(1)-y_centroid)^2-(top_locxy(1)-y_centroid)^2);
                I0y=top_int/exp(-betay*(top_locxy(1)-y_centroid)^2);
            end
            diameter_y=sigma./sqrt((2*betay));
        end

        %Check to make sure the diameter is not more than 1.5 times the size of the
        %interrogation window
        if diameter_x > 1.5*xmax_I
            diameter_x=NaN;
            I0x=NaN;
        end
        if diameter_y > 1.5*ymax_I
            diameter_y=NaN;
            I0y=NaN;
        end

        if isnan(diameter_x) && isnan(diameter_y)
            x_c = x_centroid;
            y_c = y_centroid;
            D   = NaN;
            P   = NaN;
            E   = NaN;
            Meth = 0;
            return
        end

        %Currently the diameter is approximated by a sphere of diameter equal to
        %the max size in the x/y direction.  For non-spherical particles, this
        %equation should be altered to an ellipse using diameterx and diametery as
        %the major and minor axis
        diameter=[diameter_x,diameter_y];
        I0=[I0x,I0y];
        E(1) = sqrt(max(diameter).^2 - min(diameter).^2)./max(diameter);
        E(2) = sqrt(max(diameter).^2 - min(diameter).^2)./min(diameter);

        %Validation to make sure the results are not negative
        if min([x_centroid,y_centroid,diameter,I0])<=0
            x_c      = x_centroid; %x and y should not be negitive
            y_c      = y_centroid; 
            D        = NaN;
            P        = NaN;
            Meth     = 0; 
            return
        end
        if sum(isnan(diameter))>1 || sum(isnan(I0))>1
            x_c      = x_centroid; %x and y should not be negitive
            y_c      = y_centroid; 
            D        = NaN;
            P        = NaN;
            Meth     = 0; 
            return
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % 3pt Gauss Fit - End
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        %Initial guesses for lsqnonlin - [I0 beta x_c y_c]
        if method < 3
            x0=[I0;(8./diameter.^2);repmat(x_centroid,[1 2]);repmat(y_centroid,[1 2])].';
        else
            x0=[max(I0(diameter == max(diameter))),8/max(diameter)^2, x_centroid,y_centroid];
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        try
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % 1D Methods
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if method < 3
                %point Selection
                if left_locxy(2) > 2 && right_locxy(2) + 2 <= xmax_I
                    if numel(max_row) > 1
                        points_x = (I(max_locxy(1,1),[left_locxy(2)-2:left_locxy(2),right_locxy(2):right_locxy(2)+2])).';
                        points_locx = ([left_locxy(2)-2:left_locxy(2),right_locxy(2):right_locxy(2)+2]).';
                    else
                        points_x = (I(max_locxy(1,1),max_locxy(1,2)-2:max_locxy(1,2)+2)).';
                        points_locx = (max_locxy(1,2)-2:max_locxy(1,2)+2).';
                        %If there are brighter particles in the window than the chosen
                        %max_int, then use a 3x3 window
                        if numel(find(points_x>max_int))>0
                            points_x = (I(max_locxy(1,1),max_locxy(1,2)-1:max_locxy(1,2)+1)).';
                            points_locx = (max_locxy(1,2)-1:max_locxy(1,2)-1).';
                            %keyboard
                        end
                    end
                else
                    %keyboard
                    points = I;
                    points_locxy = max_locxy;
                end
                if top_locxy(1) > 2 && bottom_locxy(1) + 2 <= ymax_I
                    if numel(max_row) > 1
                        points_y = I([top_locxy(1,1)-2:top_locxy(1),bottom_locxy(1):bottom_locxy(1)+2],max_locxy(1,2));
                        points_locy = ([top_locxy(1,1)-2:top_locxy(1),bottom_locxy(1):bottom_locxy(1)+2]).';
                    else
                        points_y = (I(max_locxy(1,1)-2:max_locxy(1,1)+2,max_locxy(1,2))).';
                        points_locy = (max_locxy(1,1)-2:max_locxy(1,1)+2).';                        
                        %If there are brighter particles in the window than the chosen
                        %max_int, then use a 3x3 window
                        if numel(find(points_y>max_int))>0
                            points_y = (I(max_locxy(1,1)-1:max_locxy(1,1)+1,max_locxy(1,2))).';
                            points_locy = (max_locxy(1,1)-1:max_locxy(1,1)+1).';                        
                            %keyboard
                        end
                    end
                else
                    %keyboard
                    points = I;
                    points_locxy = max_locxy;
                end

                %Run a nonlinear least squares regression on the window -
                %method=1 will run the standard least squares equations, and
                %method=2 will use the continuous least squares equations
                xvars=lsqnonlin(@leastsquares1D,  x0(1,1:3)  ,[],[],options,points_x,points_locx,method);
                yvars=lsqnonlin(@leastsquares1D,x0(2,[1 2 4]),[],[],options,points_y,points_locy,method);
                x_c=xvars(3);
                y_c=yvars(3);

                %If a good solution could not be reached, try making the points
                %matrix smaller to buffer out any pixels that might be affected
                %by neighboring particles.
                if min([x_c,xvars(1),xvars(2)])<=0 || max(isnan([x_c,xvars(1),xvars(2)]))>0
                    %keyboard
                    points_x=points(2:size(points,1)-1,2:size(points,2)-1);
                    points_locx=points_locxy+1;
                    xvars=lsqnonlin(@leastsquares1D,x0,[],[],options,points_x,points_locx,method);
                    x_c=xvars(3);
                end
                if min([y_c,yvars(1),yvars(2)])<=0 || max(isnan([y_c,yvars(1),yvars(2)]))>0
                    %keyboard
                    points_y=points(2:size(points_y,1)-1,2:size(points_y,2)-1);
                    points_locy=points_locy+1;
                    yvars=lsqnonlin(@leastsquares1D,x0,[],[],options,points_y,points_locy,method);
                    y_c=yvars(3);
                end

                D(1)=sqrt(sigma^2/(2*xvars(2)));
                D(2)=sqrt(sigma^2/(2*yvars(2)));
                P(1)=xvars(1);
                P(2)=yvars(1);
                E(1) = sqrt(max(D).^2 - min(D).^2)./max(D);
                E(2) = sqrt(max(D).^2 - min(D).^2)./min(D);

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % 2D Methods
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            else
                % Point Selection
                if left_locxy(2) > 2 && right_locxy(2) + 2 <= xmax_I && top_locxy(1) > 2 && bottom_locxy(1) + 2 <= ymax_I
                    if numel(max_row) > 1
                        [xloc yloc] = meshgrid(left_locxy(2)-2:right_locxy(2)+2,top_locxy(1)-2:bottom_locxy(1)+2);
                        I2 = I(top_locxy(1)-2:bottom_locxy(1)+2,left_locxy(2)-2:right_locxy(2)+2);
                        points = I2(I2 ~=max_int);
                        points_locxy(:,1) = xloc(I2 ~=max_int);
                        points_locxy(:,2) = yloc(I2 ~=max_int);
                    else
                        [xloc yloc] = meshgrid(1:xmax_I,1:ymax_I); 
                        points = reshape(I(max_locxy(1)-2:max_locxy(1)+2,max_locxy(2)-2:max_locxy(2)+2),[25 1]);
                        points_locxy = [reshape(xloc(max_locxy(1)-2:max_locxy(1)+2,max_locxy(2)-2:max_locxy(2)+2),[25 1]),...
                            reshape(yloc(max_locxy(1)-2:max_locxy(1)+2,max_locxy(2)-2:max_locxy(2)+2),[25 1])];

                        %If there are brighter particles in the window than the chosen
                        %max_int, then use a 3x3 window
                        if numel(find(points>max_int))>0
                            %keyboard
                            points = reshape(I(max_locxy(1)-1:max_locxy(1)+1,max_locxy(2)-1:max_locxy(2)+1),[25 1]);
                            points_locxy = [reshape(xloc(max_locxy(1)-1:max_locxy(1)+1,max_locxy(2)-1:max_locxy(2)+1),[25 1]),...
                                reshape(yloc(max_locxy(1)-1:max_locxy(1)+1,max_locxy(2)-1:max_locxy(2)+1),[25 1])];
                        end
                    end
                else
%                     [xloc yloc] = meshgrid(1:xmax_I,1:ymax_I);
%                     points = I(top_locxy(1):bottom_locxy(1),left_locxy(2):right_locxy(2));
%                     points = reshape(points,[numel(points) 1]);
%                     xloc = xloc(top_locxy(1):bottom_locxy(1),left_locxy(2):right_locxy(2));
%                     yloc = yloc(top_locxy(1):bottom_locxy(1),left_locxy(2):right_locxy(2));
%                     points_locxy = [reshape(xloc,[numel(xloc) 1]), reshape(yloc,[numel(yloc) 1])];
                    [xloc yloc] = meshgrid(1:xmax_I,1:ymax_I);
                    points = I(I~=0);
                    xloc   = xloc(I~=0);
                    yloc   = yloc(I~=0);
                    points_locxy = [reshape(xloc,[numel(xloc) 1]), reshape(yloc,[numel(yloc) 1])];
                end
                % Need to have atleast 4 points to solve the nonlinear
                % leasquares problem.  Use the results from the 3pt Gauss
                % when this fails.
                if numel(points)<4
                    x_c  = x_centroid;
                    y_c  = y_centroid;
                    D    = diameter;
                    P    = I0;
                    Meth = 0;
                    return
                end                
                
                %Run a nonlinear least squares regression on the window -
                %method=3 will run the standard least squares equations, and
                %method=4 will use the continuous least squares equations
                xvars=lsqnonlin(@leastsquares2D,x0,[],[],options,points,points_locxy,method);
                x_c=xvars(3);
                y_c=xvars(4);

                %If a good solution could not be reached, try making the points
                %matrix smaller to buffer out any pixels that might be affected
                %by neighboring particles.
                if min([x_c,y_c,xvars(1)])<=0 || max(isnan([x_c,y_c,xvars(1),xvars(2)]))>0
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    %didn't check betas = xvars(2)
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                    
                    points=points(2:size(points,1)-1,2:size(points,2)-1);
                    points_locxy=points_locxy+1;
                    % Need to have atleast 4 points to solve the nonlinear
                    % leasquares problem.  Use the results from the 3pt Gauss
                    % when this fails.
                    if numel(points)>=4
                        xvars=lsqnonlin(@leastsquares2D,x0,[],[],options,points,points_locxy,method);
                        x_c=xvars(3);
                        y_c=xvars(4);
                    else
                        
                        x_c  = x_centroid;
                        y_c  = y_centroid;
                        D    = diameter;
                        P    = I0;
                        Meth = 0;
                        return
                    end
                end

                D=sqrt(sigma^2/(2*abs(xvars(2))));
                P=xvars(1);
            end

            %Return Nan's if there was an error
        catch errmess
            fprintf(errmess.message)
            keyboard
            x_c  = x_centroid;
            y_c  = y_centroid;
            D    = diameter;
            P    = I0;
            Meth = 0;
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Sub Function
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function F = leastsquares1D(x,mapint_i,locxy_i,method)
            %This function is called by lsqnonlin if the least squares or continuous
            %least squares method has been chosen. x contains initial guesses[I0, betas, x_c,
            %y_c]. mapint_i is a matrix containing pixel intensity values, and locxy_i
            %is a 1x2 vector containing the row/column coordinates of the top left
            %pixel in mapint_i
            %
            %F is the variable being minimized - the difference between the gaussian
            %curve and the actual intensity values.
            %
            %Adapted from M. Brady's 'leastsquaresgaussfit' and 'mapintensity'
            %B.Drew - 7.18.2008

            Io=x(1);
            betas=x(2);
            x_cent=x(3);

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %Just like in the continuous four-point method, lsqnonlin tries negative
            %values for x(2), which will return errors unless the abs() function is
            %used in front of all the x(2)'s.
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            %num1=pi/4;
            num1=0.5*Io*sqrt(pi);            
            num2=sqrt(abs(betas));
            S = size(mapint_i);

            if method==1
                gauss_int = zeros(S);
                mapint_mod = mapint_i(:);
                xp = zeros(S);
                for ii = 1:length(mapint_i)
                    xp(ii) = locxy_i(ii);
                    % map an intensity profile of a gaussian function:
                    gauss_int(ii)=Io*exp(-abs(betas)*((xp(ii))-x_cent)^2);
                end

            elseif method==2
                gauss_int = zeros(S);
                mapint_mod = zeros(S);
                xp1 = zeros(S);
                erfx1 = zeros(S);
                erfx2 = zeros(S);
                for ii = 1:length(mapint_i)
                    xp1(ii) = locxy_i(ii)-0.5;
                    erfx1(ii) = erf(num2*(xp1(ii)-x_cent));
                    erfx2(ii) = erf(num2*(xp1(ii)+1-x_cent));
                    % map an intensity profile of a gaussian function:
                    gauss_int(ii)=(num1/num2)*(erfx2(ii)-erfx1(ii));
                    mapint_mod(ii)=mapint_i(ii);
                end
            end

            % compare the Gaussian curve to the actual pixel intensities
            F=mapint_mod-gauss_int;
        end

        function F = leastsquares2D(x,mapint_i,locxy_i,method)
            %This function is called by lsqnonlin if the least squares or continuous
            %least squares method has been chosen. x contains initial guesses[I0, betas, x_c,
            %y_c]. mapint_i is a matrix containing pixel intensity values, and locxy_i
            %is a 1x2 vector containing the row/column coordinates of the top left
            %pixel in mapint_i
            %
            %F is the variable being minimized - the difference between the gaussian
            %curve and the actual intensity values.
            %
            %Adapted from M. Brady's 'leastsquaresgaussfit' and 'mapintensity'
            %B.Drew - 7.18.2008

            Io=x(1);
            betas=x(2);
            x_cent=x(3);
            y_cent=x(4);

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %Just like in the continuous four-point method, lsqnonlin tries negative
            %values for x(2), which will return errors unless the abs() function is
            %used in front of all the x(2)'s.
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            num1=(Io*pi)/4;
            num2=sqrt(abs(betas));

            if method==3
                gauss_int = zeros(size(mapint_i));
                mapint_mod = zeros(size(mapint_i));
                xp = zeros(size(mapint_i));
                yp = zeros(size(mapint_i));
                for ii = 1:length(mapint_i)
                    xp(ii) = locxy_i(ii,1);
                    yp(ii) = locxy_i(ii,2);
                end

                % map an intensity profile of a gaussian function:
                for rr = 1:length(xp)
                    gauss_int(rr)=Io*exp(-abs(betas)*(((xp(rr))-x_cent)^2 + ...
                        ((yp(rr))-y_cent)^2));
                    mapint_mod(rr)=mapint_i(rr);
                end

            elseif method==4
                S = size(mapint_i);
                gauss_int = zeros(S(1),S(2));
                mapint_mod = zeros(S(1),S(2));
                xp = zeros(size(mapint_i));
                yp = zeros(size(mapint_i));
                for ii = 1:length(mapint_i)
                    xp(ii) = locxy_i(ii,1)-0.5;
                    yp(ii) = locxy_i(ii,2)-0.5;
                    erfx1 = erf(num2*(xp(ii)-x_cent));
                    erfy1 = erf(num2*(yp(ii)-y_cent));
                    erfx2 = erf(num2*(xp(ii)+1-x_cent));
                    erfy2 = erf(num2*(yp(ii)+1-y_cent));
                    % map an intensity profile of a gaussian function:
                    gauss_int(ii)=(num1/abs(betas))*(erfx1*(erfy1-erfy2)+erfx2*(-erfy1+erfy2));
                    mapint_mod(ii)=mapint_i(ii);
                end
            end
            mapint_i=mapint_mod;
            % compare the Gaussian curve to the actual pixel intensities
            F=mapint_i-gauss_int;
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Estimation Functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [X2_est,Y2_est,Z2_est]=particle_estloc_MAIN_V2(completed_tracks,SIZE1,estlocprops)
%
% [X2_est,Y2_est,Z2_est]=particle_estloc_MAIN_V1(completed_tracks,SIZE1,...
%   estlocprops)
%
% PROGRAM DESCRIPTION
% This is a MAIN function for providing a location prediction for the
% motion of identfied particles between two images.  Currently the function
% utilizses three predictive methods: PIV-driven, PTV-driven, and a hybrid 
% PIV/PTV-driven mode.  It is important to note that the PIV mode is 
% designed to use information from the current image pair, while the PTV 
% mode must utilize information from a previous image pair (i.e. previous
% tracks).  The PTV location prediction will work well assuming that the
% data is sufficiently resolved in time, otherwise the PIV suggested.
%
% INPUT
%   completed_tracks (cell array - [X1 X2 Y1 Y2 Z1 Z2 d1 d2 I1 I2 p#1 p#2
%       match_probability])  particle tracks from the previous image pair
%   SIZE1 (structured array) - contains XYDiameter which is has information
%       on the location/size/intensity of particles in the 1st image
%       (X1 Y1 D1 I1 p#1 size-method)
%   estlocprops (structured array) - control parameters
%       .method - 'none' 'piv' 'ptv' 'piv-ptv'
%       .PIV_PTV_weight - weighting value for the combined method; (1-0)
%           where 1-PIV and 0-PTV; linear weighting scheme
%       .PIVprops - parameters related to PIV location prediction
%           .load_dir - full path of the PIV flowfields
%           .load_name - basename of the files
%           .slsh - '\' or '/'
%           .precision - number of zeros in the filename
%           .extension - usually '.plt'
%           .framelist - vector of the frame numbers [0 1 2....]
%           .frame1 - current frame number, set to [] initially
%       .PTVprops - parameters related to PTV location prediction
%           .predict_mode - 'static' or 'dynamic'; static mode used a fixed
%               search radius to id neighboring vectors and perform
%               weighted averaging; dynamic mode varies the
%               search/averaging radius until a set # of vectors is located
%           .r_weight - radius of the search/averaging window
%           .edgeval - value of the Gaussian weighting function at
%               r_weight; should vary BETWEEEN 0 and 1
%           .numvecs - min # of vectors for serach\averaging (dynamic only)
%           .max_iterations - prevents infinite loops (dynamic only)
%        .save_dir - saving filepath; make '0' if not saving
%        .slsh - '\' or '/' (PC or Linux)
%        .s_name - base saving name
%        .s_num - saving number
%       
% OUTPUT
%   X2_est, Y2_est, Z2_est - estimated location of the particles
%
%(v1) N.Cardwell - 11.18.2009
%(v2) N.Cardwell - 04.01.2010

%extract particle information from the SIZE1 structured array
X1=SIZE1.XYDiameter(:,1);  Y1=SIZE1.XYDiameter(:,2);  Z1=zeros(size(X1));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%MAIN LOGIC BLOCK FOR DETERMINING THE LOCATION ESTIMATION FOR IM2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nnz( strcmpi(estlocprops.method,{'none' 'piv' 'ptv' 'piv-ptv'}) ) == 0
    fprintf('ERROR - An incorrect estimation location method was provided');
    return
end

switch lower(estlocprops.method)
    case {'none'}
        X2_est=X1;   Y2_est=Y1;   Z2_est=Z1;
    
    case {'piv'}
        %call the PIVestloc subfunction
        [X2_est,Y2_est,Z2_est]=PIVestloc_V1(X1,Y1,Z1,estlocprops);
    
    case {'ptv'}
        %call the PTVestloc subfunction
        [X2_est,Y2_est,Z2_est]=PTVestloc_V1(X1,Y1,Z1,estlocprops,completed_tracks);
    
    case {'piv-ptv'}
        %call the both the PIVestloc and PTVestloc subfunctions
        [X2_estPIV,Y2_estPIV,Z2_estPIV]=PIVestloc_V1(X1,Y1,Z1,estlocprops);
        [X2_estPTV,Y2_estPTV,Z2_estPTV]=PTVestloc_V1(X1,Y1,Z1,estlocprops,completed_tracks);
        
        %apply the user-defined weighting scheme
        X2_est = X2_estPIV.*estlocprops.PIV_PTV_weight + ...
            X2_estPTV.*(1-estlocprops.PIV_PTV_weight);
        Y2_est = Y2_estPIV.*estlocprops.PIV_PTV_weight + ...
            Y2_estPTV.*(1-estlocprops.PIV_PTV_weight);
        Z2_est = Z2_estPIV.*estlocprops.PIV_PTV_weight + ...
            Z2_estPTV.*(1-estlocprops.PIV_PTV_weight);
end

%section for saving sized particle information
if ~isempty(estlocprops.save_dir)
    if exist(estlocprops.save_dir,'dir')~=7
        fprintf('Making Save Directory %s \n',estlocprops.save_dir)
        mkdir(estlocprops.save_dir)
    end
    sname = sprintf('%%s%%0%0.0fd',estlocprops.Data.imzeros);
    save(fullfile(estlocprops.save_dir,sprintf(sname,estlocprops.s_name,estlocprops.s_num)),...
        'estlocprops','X2_est','Y2_est','Z2_est');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%plotting code - COMMENT OUT FOR NORMAL OPERATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% dummy=cell2mat(completed_tracks);
% figure;
% quiver(dummy(:,1),dummy(:,3),dummy(:,2)-dummy(:,1),dummy(:,4)-dummy(:,3),0,'Color','k');
% set(gca,'DataAspectRatio',[1 1 1]);
% hold on
% 
% %quiver(X1,Y1,X2_estPIV-X1,Y2_estPIV-Y1,0,'Color','r');
% %quiver(X1,Y1,X2_estPTV-X1,Y2_estPTV-Y1,0,'Color','b');
% quiver(X1,Y1,X2_est-X1,Y2_est-Y1,0,'Color','m');
% legend({'completed_tracks' 'PIV_est' 'PTV_est' 'Est_loc'});

end

function [X2_est,Y2_est,Z2_est]=PTVestloc_V1(X1,Y1,Z1,PTVprops,completed_tracks)
%
% [X2_est,Y2_est,Z2_est]=PTVestloc_V1(X1,Y1,Z1,PTVprops,completed_tracks)
%
% PROGRAM DESCRIPTION
% This function provides a diplacement estimation for particles in an image
% given previous tracks of the particles path.  A spatially weighted
% average (Gaussian) is used to provide the displacment prediction, based
% on the surrounding tracks.  Two main modes of operation are offered:
% static and dynamic.  The static used a fixed search raduis to search for
% neighboring particles.  The dynamic mode succesively increases the search
% radius until a user-defined number of vectors have been identified, or
% the max number of user-defined iterations has occured.  
% 
% The dynamic is suggested for sparse fields or when the spatial vector
% density varies thoroughout the image.  Otherwise the static is suggested.
%
% INPUTS
%   X1,Y1,Z1 - original particle locations
%   PTVprops (structured array)
%       .predict_mode - 'static' or 'dynamic'
%       .r_weight - radius (pixels) of the initial search window and the
%           gaussian weighting function
%       .edgeval - value of the GWF at r_weight (>0 & <1)
%       .numvecs - number of vectors to satisfy the dynamic mode
%       .max_iterations - prevents infinite runs in dynamic mode
%
% OUTPUTS
%   X2_est,Y2_est,Z2_est - estimated location of the particles
%
%(v1) N.Cardwell - 11.17.2009

if isempty(completed_tracks)==1
    %no completed tracks to use for PTV location prediction
    X2_est=X1;   Y2_est=Y1;   Z2_est=Z1;
else
    %completed tracks available for PTV location prediction
    switch lower(PTVprops.predict_mode)
        
        %static method uses a fixed radius to ID neighboring particles
        case {'static'}
            %intialize arrays
            X2_est=zeros(size(X1));  Y2_est=zeros(size(Y1));  Z2_est=zeros(size(Z1));
            dummy=cell2mat(completed_tracks);
            
            %perform static location prediction for each new particle
            for i=1:size(X1,1)
                loc=[X1(i,1),Y1(i,1),Z1(i,1)];

                %determine particles within the search radius 'r_weight'
                dX=dummy(:,1)-loc(1);  dY=dummy(:,3)-loc(2);  dZ=dummy(:,5)-loc(3);
                distance=sqrt(dX.^2+dY.^2+dZ.^2);
                check = (distance <= PTVprops.r_weight & distance ~= 0);

                %if no particles found within r_weight, initialize w/ the
                %origianl location
                if nnz(check)==0
                    X2_est(i,1)=loc(1);  Y2_est(i,1)=loc(2); Z2_est(i,1)=loc(3);

                %otherwise use the spatial weighting function to predict w/    
                else
                    data=zeros(nnz(check),4);
                    data(:,1)=dummy(check,1);  
                    data(:,2)=dummy(check,3);
                    data(:,3)=dummy(check,5);
                    data(:,4)=dummy(check,2)-dummy(check,1);
                    data(:,5)=dummy(check,4)-dummy(check,3);
                    data(:,6)=dummy(check,6)-dummy(check,5);
                    [U_est,V_est]=twoVar_spatial_weighting(...
                        PTVprops.r_weight,PTVprops.edgeval,loc,data);
                    X2_est(i,1)=X1(i,1)+U_est;  Y2_est(i,1)=Y1(i,1)+V_est;
                end
            end
%             figure; quiver(dummy(:,1),dummy(:,3),dummy(:,2)-dummy(:,1),dummy(:,4)-dummy(:,3),0,'Color','g');
%             hold on;  quiver(X1,Y1,X2_est-X1,Y2_est-Y1,0,'Color','r');
%             set(gca,'DataAspectRatio',[1 1 1]);

        %dynamic method varies to ID radius to get a specific # of part
        case {'dynamic'}
            %intialize arrays
            X2_est=zeros(size(X1));  Y2_est=zeros(size(Y1));  Z2_est=Z1;
            dummy=cell2mat(completed_tracks);

            %perform static location prediction for each new particle
            for i=1:size(X1,1)
                loc=[X1(i,1),Y1(i,1),Z1(i,1)];
                dX=dummy(:,1)-loc(1);  dY=dummy(:,3)-loc(2);  dZ=dummy(:,5)-loc(3);
                distance=sqrt(dX.^2+dY.^2+dZ.^2);
                
                %determine particles within the search radius 'r_weight',
                %increase the r_weight until numvecs is satified
                number_o_vecs=0;  count=1;  r_weight_temp=PTVprops.r_weight;
                while (number_o_vecs<PTVprops.numvecs && count<PTVprops.max_iterations)
                    check = (distance <= r_weight_temp & distance ~= 0);
                    number_o_vecs=nnz(check);
                    if number_o_vecs<PTVprops.numvecs
                        r_weight_temp=r_weight_temp+1;
                    end
                    count=count+1;
                end
                data=zeros(nnz(check),6);
                data(:,1)=dummy(check,1);  data(:,2)=dummy(check,3);  data(:,3)=dummy(check,5); 
                data(:,4)=dummy(check,2)-dummy(check,1);
                data(:,5)=dummy(check,4)-dummy(check,3);
                data(:,6)=dummy(check,6)-dummy(check,5);
                [U_est,V_est,W_est]=twoVar_spatial_weighting(...
                    r_weight_temp,PTVprops.edgeval,loc,data);
                X2_est(i,1)=X1(i,1)+U_est;  Y2_est(i,1)=Y1(i,1)+V_est; Z2_est(i,1)=Z1(i,1)+W_est;
            end
%             figure; quiver(dummy(:,1),dummy(:,3),dummy(:,2)-dummy(:,1),dummy(:,4)-dummy(:,3),0,'Color','g');
%             hold on;  quiver(X1,Y1,X2_est-X1,Y2_est-Y1,0,'Color','r');
%             set(gca,'DataAspectRatio',[1 1 1]); 
    end
end

end

function [X2_est,Y2_est,Z2_est]=PIVestloc_V1(X1,Y1,Z1,PIVprops)
%
% [X2_est,Y2_est,Z2_est]=PIVestloc_V1(X1,Y1,Z1,PIVprops);
%
% PROGRAM DESCRIPTION
% This function provides a diplacement estimation for particles in an image
% given a provided flowfield (typically from PIV).  Processing of PIV 
% flowfields is a separate step and must be completed BEFORE implementing 
% this function.  The user is cautioned to make sure that the supplied PIV 
% velocity field has the same orientation as the identified/sized/and 
% tracked particle fields.  This program assumes that the user will be
% supplying fields from the PIVadvance2 code (written by A.Eckstein) and 
% rotates/transforms the PIV field automatically.  This transformation can
% be edited in the 'PIV loading and transformation block' of the code if 
% desired (specifically line 41).
%
% INPUTS
%   X1,Y1,Z1 - original particle locations
%   PIVprops - loading parameters related to the PIV processed flowfield
%       PIVprops.load_dir - full path of the PIV-file directory
%       PIVprops.load_name - base name of the PIV-file
%       PIVprops.precision - base precision of the PIV-file
%       PIVprops.extension - base extension of the PIV-file (EX: '.plt')
%       PIVprops.frame1 - frame number of the PIV-file used to estimate the
%           particle displacment from SIZE1 to SIZE2
%       PIVprops.slsh - '\' or '/'
%
% OUTPUTS
%   X2_est,Y2_est,Z2_est - estimated location of the particles
%
%(v1) N.Cardwell - 11.1.2009

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%PIV loading and transformation block
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%load in the PIV processed data for the image pair 
%(rearrange coordiante system to match the output of the sizing program)
if any(strcmpi(PIVprops.PIVprops.extension,{'.plt' '.dat'}))
    [X,Y,U,V]=read_pltmod_NC(PIVprops.Data.PIV.Data.outbase,PIVprops.Data.Track.PIVprops.load_dir,...
        PIVprops.PIVprops.frame1,PIVprops.PIVprops.frame1,PIVprops.Data.imzeros);
    PIV1={X,Y,flipud(U),-1.*flipud(V)};  clear X Y U V
elseif strcmpi(PIVprops.PIVprops.extension,'.mat')
    lname = sprintf('%%s%%s%%0%0.0fd',PIVprops.Data.imzeros);
    PIV1t = load(sprintf(lname,PIVprops.Data.Track.PIVprops.load_dir,PIVprops.PIVprops.outbase,PIVprops.PIVprops.frame1,PIVprops.PIVprops.extension));
%     PIV1 = {PIV1t.flowvarunsteady(end:-1:1,:,1,1),PIV1t.flowvarunsteady(end:-1:1,:,1,2),...
%         PIV1t.flowvarunsteady(end:-1:1,:,1,3),PIV1t.flowvarunsteady(end:-1:1,:,1,4)};
    PIV1 = {PIV1t.X(end:-1:1,:),PIV1t.Y(end:-1:1,:),PIV1t.U(end:-1:1,:),PIV1t.V(end:-1:1,:)};

else
    error('Unknown PIV extension')
end

%interpolate the PIV flowfield at each particle location in im1
UI=interp2(PIV1{1},PIV1{2},PIV1{3},X1,Y1,'spline');
VI=interp2(PIV1{1},PIV1{2},PIV1{4},X1,Y1,'spline');

check=isnan(UI); %replace NaNs with the original particle location (occurs when interp2 fails)
X2_est=zeros(size(X1));  Y2_est=zeros(size(Y1));  Z2_est=zeros(size(Z1));
for j=1:length(UI)
    if check(j)==1
        X2_est(j)=X1(j);  Y2_est(j)=Y1(j);  Z2_est(j)=Z1(j)+0;
    else
        X2_est(j)=X1(j)+UI(j);  Y2_est(j)=Y1(j)+VI(j);  Z2_est(j)=Z1(j)+0;
    end
end

% %create overlaid quiver plot of the PIV and interpolated velocity field
% figure;  quiver(PIV1{1},PIV1{2},PIV1{3},PIV1{4},0);
% set(gca,'DataAspectRatio',[1 1 1]);
% set(gca,'YDir','reverse');
% hold on
% quiver(X1,Y1,UI,VI,0,'r');
% axis([1 64 1 64])
% 
% %create scatter plot of position and estimated position
% figure;  scatter(X1,Y1,'.r');
% set(gca,'DataAspectRatio',[1 1 1]);
% %set(gca,'YDir','reverse');
% hold on
% scatter(X2_est,Y2_est,'+b');
% scatter(X2,Y2,'xg');
% legend({'loc1','loc2'})
% axis([1 64 1 64])

end

function [var1_est,var2_est,var3_est]=twoVar_spatial_weighting(r_weight,edgeval,loc,data)
%
% [var1_est,var2_est]=twoVar_spatial_weighting(r_weight,edgeval,loc,data)
%
% PROGRAM DESCRIPTION
% This function provides, for a given location, a weigthed estimate on 
% two variables based on the surrounding data, which may be random or 
% structured.  The weighting is based on a Gaussian distribution, of which 
% the user can vary.
%
% INPUTS
%   r_weight - radius of the spatial weighting function
%   edgeval - value of the spatial weighting function at r_weight
%       (between 0-1, but not equal to either 0 or 1); lower values weight the
%       center higher, higher values weight everything more evenly)
%   loc - location predicted point [x,y]
%   data - surrounding points and their values (two-variable)
%       [x,y,var1,var2]
%
% OUTPUT
%   var1_est - predicted value of var1 at 'loc'
%   var2_est - predicted value of var2 at 'loc'
%
%(v1) N.Cardwell - 11.16.2009
%(v2) S. Raben   - 05.25.2010

%compute the required standard deviation of the Gaussian curve so that the
%value ar r_weight equals edgeval
std_dev=sqrt(-(r_weight^2)/(2*log(edgeval)));

%compute the Gaussian weighting values at each point in data
gaus_weight_x2 = exp(- ((data(:,1)-loc(1)).^2) / (2*std_dev^2));
gaus_weight_y2 = exp(- ((data(:,2)-loc(2)).^2) / (2*std_dev^2));
gaus_weight_z2 = exp(- ((data(:,3)-loc(3)).^2) / (2*std_dev^2));

%predict the value of both variables at 'loc' (sum of the products over
%the sum of the weights)
var1_est=sum(gaus_weight_x2.*data(:,4))/sum(gaus_weight_x2);
var2_est=sum(gaus_weight_y2.*data(:,5))/sum(gaus_weight_y2);
var3_est=sum(gaus_weight_z2.*data(:,6))/sum(gaus_weight_z2);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%plotting code - COMMENT OUT FOR NORMAL OPERATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% figure;
% quiver(data(:,1),data(:,2),data(:,3),data(:,4),0);  hold on
% set(gca,'DataAspectRatio',[1 1 1]);
% axis([loc(1)-2*r_weight loc(1)+2*r_weight loc(2)-2*r_weight loc(2)+2*r_weight])
% scatter(loc(1),loc(2),'r');
% rectangle('Position',[loc(1)-r_weight,loc(2)-r_weight,r_weight*2,r_weight*2],...
%     'Curvature',[1,1],'LineStyle','--')
% scatter(loc(1)+var1_est,loc(2)+var2_est,'m','+');
% line([loc(1) loc(1)+var1_est],[loc(2) loc(2)+var2_est],'Color','m');
% 
% %compute the unweighted average and overlay on the plot
% avg_U=mean(data(:,3));  avg_V=mean(data(:,4));
% scatter(loc(1)+avg_U,loc(2)+avg_V,'g');
% line([loc(1) loc(1)+avg_U],[loc(2) loc(2)+avg_V],'Color','g');

% %set up the Gaussian weighting scheme
% mean_var=loc(1);
% x = (mean_var-r_weight : 2*r_weight/100 : mean_var+r_weight);
% std_dev=sqrt(-(r_weight^2)/(2*log(edgeval)));
% gaussian_weight_x=exp(- ((x-mean_var).^2) / (2*std_dev^2));
% plot(x,gaussian_weight_x.*r_weight,'b');

end

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Tracking Functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [tracks]=particle_track_MAIN_V1(X2_est,Y2_est,Z2_est,SIZE1,SIZE2,trackprops,valprops)
%
% [tracks]=particle_track_MAIN_V1(X2_est,Y2_est,Z2_est,SIZE1,SIZE2,...
%   trackprops,valprops)
%
% PROGRAM DESCRIPTION
% This is a MAIN function for tracking particles between two consecutive
% images, where the particles have already been identified and sized using
% the particle_ID_MAIN and particle_size_MAIN functions.  The tracking
% fucntion also assumes that a displacement estimate has already been
% computed, typically using the particle_estloc_MAIN function.
%
% INPUT
%   X2_est, Y2_est, Z2_est - estimated location of the particles
%   SIZE1/SIZE2 - structured array containing the sizing results
%       .XYDiameter = 6-column matrix; 1st column is x-centorid location, 
%           2nd column is y-centroid location, 
%           3rd column is particle diameter, 
%           4th column is true max. intensity, 
%           5th column is particle id#, 
%           6th column is sizing method
%       .mapsizeinfo - (2 column array) defines the size (row col) of the 
%           each sized particle
%       .locxy - (2 column array) ROW/COL? location associated with the 
%           upper left pixel of the particle's rectangular projection
%   trackprops (structured array) tracking control parameters
%       .s_radius - radius to search for neighboring particles
%       .weights - [0-1 0-1 0-1] sets the relative emphasis for the pair
%               matching algorithm: inter-particle distance, size diff, &
%               max-intensity diff
%       .save_dir - full path of the save directory, make '0' for no saving
%       .slsh - '\' or '/'
%       .s_name - base save name
%       .s_num - intialize to [], controlled by the program
%       .plotfig - set to '1' to plot the results, otherwise set to '0'
%   valprops (structured array) validation control parameters
%       .num_pass - number of validation passes
%       .method - (vector of strings) 'none' 'coeff' 'mean' 'median'
%           'relaxation' sets the method to ID bad vectors and subsequently
%           revise the displacement estimate
%       .C_cutoff - max matching coeff to accept a track (between 0-1)
%       .s_radius - radius to search for neighboring tracks
%       .MAD_U/MAD_V - maximum allowable mean/median absolute deviation
%           before a track is considered 'bad'
%
%                               SAMPLE VALPROPS
%   valprops.numpass=4;
%   valprops.method={'median' 'median' 'median' 'coeff'};
%   valprops.C_cutoff=[1 1 1 0.2];
%   valprops.s_radius=[15 15 15 15];
%   valprops.MAD_U=[2 1.5 1 1];
%   valprops.MAD_V=[2 1.5 1 1];
%
% OUTPUT
%   tracks - main output array of the matched particle pairs:
%       [X1 X2 Y1 Y2 Z1 Z2 d1 d2 I1 I2 p#1 p#2 match_probability]
%       *where a lower match_coefficient is better, should vary
%        between 0 and 1*
%
%(v1) N.Cardwell - 11.18.2009

%set the locations for the particles in image 1 and 2
X1=SIZE1.XYDiameter(:,1);  X2=SIZE2.XYDiameter(:,1);
Y1=SIZE1.XYDiameter(:,2);  Y2=SIZE2.XYDiameter(:,2);
d1=SIZE1.XYDiameter(:,3);  d2=SIZE2.XYDiameter(:,3);
I1=SIZE1.XYDiameter(:,4);  I2=SIZE2.XYDiameter(:,4);
Z1=zeros(size(X1));        Z2=zeros(size(X2));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%MAIN TRACKING BLOCK
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%track the particles in the image pair using the 3D weighted
%nearest neighbor tracking method
[tracks]=weighted_nearest_neighbor3D(X1,X2,X2_est,Y1,Y2,Y2_est,...
    Z1,Z2,Z2_est,d1,d2,I1,I2,trackprops.weights,trackprops.s_radius);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Validation of the determined tracks
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%decision block - validate the measured tracks?
if valprops.run==1
for i=1:valprops.numpass
    %intialize temporary valprops array
    valprops_i.method=valprops.method{i};
    valprops_i.C_cutoff=valprops.C_cutoff(i);
    valprops_i.s_radius=valprops.s_radius(i);
    valprops_i.MAD_U=valprops.MAD_U(i);
    valprops_i.MAD_V=valprops.MAD_V(i);
    
    switch lower(valprops.method{i})
        case {'none'}
            %accepts the measured tracks with no post processing

        case {'coeff'}
            %retain tracks below the coefficient threshold - discard those above
            tracks=tracks((tracks(:,13) <= valprops_i.C_cutoff),:);

        case {'mean','median'}
            %call validation mean/median subfunction
            [MAD_ratio,MAD_ratio_hdr]=PTVval_meanandmedian(tracks,valprops_i);

            %refresh the location estimation with the validation results
            X2_est( MAD_ratio(:,5) ,1) = MAD_ratio(:,6);
            Y2_est( MAD_ratio(:,5) ,1) = MAD_ratio(:,7);

            %rerun tracking with new estimated locations
            [tracks]=weighted_nearest_neighbor3D(X1,X2,X2_est,Y1,Y2,Y2_est,...
                Z1,Z2,Z2_est,d1,d2,I1,I2,trackprops.weights,trackprops.s_radius);

        case {'relaxation'}
            %TO BE COMPLETED - see paper by Ohmi and Lee (2000) on the new
            %relaxation method for more information

    end
end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Saving of determined tracks
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%save tracking results and relevant info if called for by the user
if ~isempty(trackprops.save_dir)
    if exist(trackprops.save_dir,'dir')~=7
        fprintf('Making Save Directory %s \n',trackprops.save_dir)
        mkdir(trackprops.save_dir)
    end
    sname = sprintf('%%s%%0%0.0fd.mat',trackprops.Data.imzeros);
    save(fullfile(trackprops.save_dir,sprintf(sname,trackprops.s_name,trackprops.s_num)),...
        'tracks','trackprops','valprops');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Plotting of the determined tracks
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if trackprops.plotfig==1
    %plot the tracking results as a connected scatter plot of particle
    %positions in IM1 and IM2
    im_bounds=[1 800 1 1400 -1 1];
    %figure;  
    scatter(tracks(:,1),tracks(:,3),'.b');
    set(gca,'DataAspectRatio',[1 1 1]);
    set(gca,'YDir','reverse');
    hold on
    scatter(tracks(:,2),tracks(:,4),'+r');  hold on
    for j=1:size(tracks,1)
        line([tracks(j,1);tracks(j,2)],[tracks(j,3);tracks(j,4)]);
        hold on
    end
    axis(im_bounds(1:4))
    legend({'image1','image2'},'Location','NorthOutside','Orientation','horizontal')
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUB-FUNCTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [tracks]=weighted_nearest_neighbor3D(X1,X2,X2_est,Y1,Y2,Y2_est,Z1,Z2,Z2_est,d1,d2,I1,I2,weights,s_radius)
%
% This function tracks discrete locations (particles) between two successive
% realizations.  Pairing is based on a improved nearest-neighbor algorithm 
% which makes use of additional information about the particles to improve
% the pair-matching effecitiveness, specifically the diameter and peak
% intensity value.  Additionally, the program also makes use of an
% estimation of particle displacement if available
%
% INPUTS
%   (1-D number vectors of equal length)
%   X1,Y1,Z1 - spatial location of the particles in image 1
%   d1,I1 - diameter and maximum intensity of the particles in image 1
%   X2,Y2,Z2 - spatial location of the particles in image 2
%   d2,I2 - diameter and maximum intensity of the particles in image 2
%   X2_est,Y2_est,Z2_est - estimated particle location in image 2 
%
%   weights - [0-1 0-1 0-1] sets the relative emphasis for the pair
%       matching algorithm: inter-particle distance, size diff, & max-intensity diff
%   s_radius - (num, pixles) 3D search radius to identify particle pairs
%
% OUTPUTS
%   tracks - main output array of the matched particle pairs:
%       [X1 X2 Y1 Y2 Z1 Z2 d1 d2 I1 I2 p#1 p#2 match_probability]
%       *where a lower match_coefficient is better, should vary
%        between 0 and 1*
%
%(v6) N. Cardwell - 7.21.2009
%(v7) S. Raben    - 4.2011

%compute relevant parameters for each image
num_p1=length(X1);  num_p2=length(X2);
d_diff_max=max([d1;d2])-min([d1;d2]);
I_diff_max=max([I1;I2])-min([I1;I2]);
% If there is no difference between values (i.e. all the values are the
% same) then set variable to inf so that they are removed from the
% probablity calcuation.
if d_diff_max == 0;
    d_diff_max = inf;
end
if I_diff_max == 0
    I_diff_max = inf;
end

%adjust the particle locations in image 1 to match the estimated locations
%provided in X2_est, Y2_est, and Z2_est (save org. variables for later use)
X1_org=X1;  Y1_org=Y1;  Z1_org=Z1;
X1=X2_est;  Y1=Y2_est;  Z1=Z2_est;

%for each particle in image 1, create a "comparison" array for all particles
%in image 2 contained within the search radius and centered on the particle
%location in image 1
compare=cell(num_p1,1);
for i=1:num_p1
    dX=X2-X1(i,1);  dY=Y2-Y1(i,1);  dZ=Z2-Z1(i,1);
    distance=sqrt(dX.^2+dY.^2+dZ.^2);
    
    %build comparison array using particle whose linear distance is within
    %the user defined search radius (s_radius)
    p_index=find(distance<=s_radius);  compare_i=zeros(length(p_index),3);
    compare_i(:,1)=ones(size(p_index)).*i;  compare_i(:,2)=p_index;
    
    %compute the match probability for each possible pairing
    Prob_D=(distance(p_index)./s_radius).*weights(1);
    Prob_d=(abs(d2(p_index)-d1(i,1))./d_diff_max).*weights(2);
    Prob_I=(abs(I2(p_index)-I1(i,1))./I_diff_max).*weights(3);
    compare_i(:,3)=(Prob_D+Prob_d+Prob_I)./sum(weights);
    
    %populate main 'compare' array 
    compare{i}=compare_i;
    

    clear dX dY dZ distance p_index compare_i Prob_D Prob_d Prob_I
end
clear i
compare = cell2mat(compare);
%sort the comparison array in ascending order based on the pairing
%probability
compare_sort=sortrows(compare,3);

%determine the best particle matches and successively match particles
%until the max number of pairs has been reached (or no more exist)
p_pairs=[];
max_matches=min(num_p1,num_p2);
found_matches_im1=[];  found_matches_im2=[];
c=0;  i=1;
while (c<=max_matches) && (i<=size(compare_sort,1))
    test1=found_matches_im1==compare_sort(i,1);
    test2=found_matches_im2==compare_sort(i,2);
    if (any(test1) || any(test2))
        i=i+1;
    else
        found_matches_im1=vertcat(found_matches_im1,compare_sort(i,1));
        found_matches_im2=vertcat(found_matches_im2,compare_sort(i,2));
        p_pairs=vertcat(p_pairs,compare_sort(i,:));
        i=i+1;
        c=c+1;
    end
end

%populate the 'tracks' arrays with the following structure: 
%   tracks=[X1, X2, Y1, Y2, Z1, Z2, d1, d2, I1, I2, p#1 p#2 match_probability]
try;tracks=[X1_org(p_pairs(:,1)), X2(p_pairs(:,2)), Y1_org(p_pairs(:,1)), Y2(p_pairs(:,2)), ...
    Z1_org(p_pairs(:,1)), Z2(p_pairs(:,2)), d1(p_pairs(:,1)), d2(p_pairs(:,2)), ...
    I1(p_pairs(:,1)), I2(p_pairs(:,2)), p_pairs(:,:)];
catch;keyboard;end
end

function [MAD_ratio,MAD_ratio_hdr]=PTVval_meanandmedian(tracks,valprops)
%
% [MAD_ratio,MAD_ratio_hdr]=PTVval_meanandmedian(tracks,valprops);
%
% PROGRAM DESCRIPTION
% This function validates each particle track by statistically comparing it
% to neighboring tracks within a user defined search radius.  Once deemed
% an incorrect track, an position estimate is provided based on the motion
% of the surrounding particles
% 
% INPUTS
%   tracks - main  array of the matched particle pairs:
%       [X1 X2 Y1 Y2 Z1 Z2 d1 d2 I1 I2 p#1 p#2 match_probability]
%   valprops - control parameters for the validation of tracks
%       valprops.method - 'mean' or 'median' ; controls if the validation
%           is performed with the mean or median statistical properties
%       valprops.s_radius - pixel value to search for neighboring tracks
%       valprops.MAD_U - allowable ratio of mean/median absolute deviation
%           of each track vs. the avg. MAD of the neighboring tracks
%       valprops.MAD_V - same as above, but for the V-component of vel.
%
% OUTPUTS
%   MAD_ratio - estimated location of the ID'ed bad vectors for next
%       tracking step (2-D array)
%           [MAD_ratio_U; MAD_ratio_V; bad_U(0|1); bad_V(0|1); p#1; 
%           est_locx; est_locy]
%
%(v1) N. Cardwell - 10.31.2009

%get info and preallocate arrays
num_tracks=size(tracks,1);
MAD_ratio=zeros(num_tracks,7);

%main program loop for all tracks
for i=1:num_tracks

    %remove the track currently being evaluated (replace w/ NaN)
    trk_i=tracks(i,:);   tracks_i=tracks;   tracks_i(i,:)=NaN;
    
    %locate neighboring particles within user specified s_radius
    dist_i=sqrt( (tracks_i(:,1)-trk_i(1)).^2 + (tracks_i(:,3)-trk_i(3)).^2 + (tracks_i(:,5)-trk_i(5)).^2 );
    [r_i,c_i]=find(dist_i <= valprops.s_radius);

    %decision block - which method to use to ID and recalc "bad" tracks?
    switch lower(valprops.method)
        %use MEAN absoulte deviation
        case {'mean'}
            MAD_ratio_hdr={'MeanADratio Ui' 'MeanADratio Vi' 'bad vec U' 'bad vec V'...
                'part index' 'new est_locx' 'new est_locy'};

            %compute the mean of the neighboring vectors
            Ui=tracks_i(r_i,2) - tracks_i(r_i,1);
            Vi=tracks_i(r_i,4) - tracks_i(r_i,3);
            Mean_Ui=mean( Ui );   Mean_Vi=mean( Vi );

            %compute the mean absolute deviation of all neighboring vectors
            MeanAD_Ui=mad(Ui,0);   MeanAD_Vi=mad(Vi,0);

            %compute the mean deviation for trk_i
            trk_i_U = trk_i(2)-trk_i(1);   trk_i_V = trk_i(4)-trk_i(3);
            MeanAD_trk_iU = abs(trk_i_U - Mean_Ui);
            MeanAD_trk_iV = abs(trk_i_V - Mean_Vi);

            %compute the ratio of the MeanAD of the neighboring tracks
            %with the MeanAD of trk_i
            MAD_ratio(i,1) = MeanAD_trk_iU/MeanAD_Ui;
            MAD_ratio(i,2) = MeanAD_trk_iV/MeanAD_Vi;

            %determine if the MeanAD_ratio (either U or V) exceeds the user
            %specified limits (MAD_U and MAD_V)
            MAD_ratio(i,3) = MAD_ratio(i,1) > valprops.MAD_U;
            MAD_ratio(i,4) = MAD_ratio(i,2) > valprops.MAD_V;

            if ( MAD_ratio(i,3)==1 || MAD_ratio(i,4)==1 )==1
                %grab the particle number (related to SIZE1)
                MAD_ratio(i,5) = trk_i(11);
                
                %use the mean U and V velocity to get the new est_loc X & Y
                MAD_ratio(i,6) = trk_i(1) + Mean_Ui;
                MAD_ratio(i,7) = trk_i(3) + Mean_Vi;
            end

            %use MEDIAN absolute deviation
        case {'median'}
            MAD_ratio_hdr={'MedianADratio Ui' 'MedianADratio Vi' 'bad vec U' 'bad vec V'...
                'part index' 'new est_locx' 'new est_locy'};
            
            %compute the median of the neighboring vectors
            Ui=tracks_i(r_i,2) - tracks_i(r_i,1);
            Vi=tracks_i(r_i,4) - tracks_i(r_i,3);
            Median_Ui=median( Ui );   Median_Vi=median( Vi );

            %compute the median absolute deviation of all neighboring vectors
            MedianAD_Ui=mad(Ui,1);   MedianAD_Vi=mad(Vi,1);

            %compute the median deviation for trk_i
            trk_i_U = trk_i(2)-trk_i(1);   trk_i_V = trk_i(4)-trk_i(3);
            MedianAD_trk_iU = abs(trk_i_U - Median_Ui);
            MedianAD_trk_iV = abs(trk_i_V - Median_Vi);

            %compute the ratio of the MeanAD of the neighboring tracks
            %with the MeanAD of trk_i
            MAD_ratio(i,1) = MedianAD_trk_iU/MedianAD_Ui;
            MAD_ratio(i,2) = MedianAD_trk_iV/MedianAD_Vi;

            %determine if the MeanAD_ratio (either U or V) exceeds the user
            %specified limits (MAD_U and MAD_V)
            MAD_ratio(i,3) = MAD_ratio(i,1) > valprops.MAD_U;
            MAD_ratio(i,4) = MAD_ratio(i,2) > valprops.MAD_V;

            if ( MAD_ratio(i,3)==1 || MAD_ratio(i,4)==1 )==1
                %grab the particle number (related to SIZE1)
                MAD_ratio(i,5) = trk_i(11);
                
                %use the mean U and V velocity to get the new est_loc X & Y
                MAD_ratio(i,6) = trk_i(1) + Median_Ui;
                MAD_ratio(i,7) = trk_i(3) + Median_Vi;
            end

    end
end

% %Plotting code for all tracks and validated tracks
% im_bounds=[1 64 1 64 -1 1];
% figure;  scatter(tracks(:,1),tracks(:,3),'.b');  hold on
% scatter(tracks(:,2),tracks(:,4),'+r');
% for j=1:size(tracks,1)
%     line([tracks(j,1);tracks(j,2)],[tracks(j,3);tracks(j,4)]);
% end
% scatter(MAD_ratio(MAD_ratio(:,5)~=0,6),MAD_ratio(MAD_ratio(:,5)~=0,7),'xg');
% new_vecs=MAD_ratio(MAD_ratio(:,5)~=0,5:7);
% for j=1:nnz(MAD_ratio(:,5))
%     line([tracks(tracks(:,11)==new_vecs(j,1),1) ; new_vecs(j,2) ],...
%         [tracks(tracks(:,11)==new_vecs(j,1),3) ; new_vecs(j,3)],'Color','g');
% end
% legend({'image1','image2','im1-im2','image2:val','im1-im2:val'},...
%     'Location','NorthOutside','Orientation','horizontal')
% set(gca,'DataAspectRatio',[1 1 1]);
% set(gca,'YDir','reverse');
% axis(im_bounds(1:4))

%compress output array
MAD_ratio=MAD_ratio(MAD_ratio(:,5)~=0,:);

end

