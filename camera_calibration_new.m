function [calibration_data,calibration_plane_data]=camera_calibration_new(InputJobFile)
% This function is used to generate a set of calibration coordinates from
% the stereo calibration plate of the form
%
%  x = f(X,Y,Z)
%  y = g(X,Y,Z)
%
% where (X,Y,Z) are the lab reference frame world coordinates (defined by
% the grid spacing in the calibration plate and its z traverse position)
% and (x,y) are the image coordinates of the grid points.  This function
% must first identify the image coordinates of the grid points and then
% correlate these points with the appropriate world coordinates.

% This is a list (vector) of integers specifying camera numbers.
camera_numbers_list =InputJobFile.Camera_Numbers_List;
% This is the number of calibration planes
calibration_planes_list = InputJobFile.Plane_Numbers_List;
% This is the spacing betweent the calibration planes in mm
plane_spacing = InputJobFile.Plane_Spacing;
% This is the starting z depth in mm
z_grid_start = InputJobFile.Z_Grid_Start;
% This is the x spacing between grid points in mm
x_grid_spacing=InputJobFile.X_Grid_Spacing;
% This is the y spacing between grid points in mm
y_grid_spacing=InputJobFile.Y_Grid_Spacing;
% This is the diameter of the grid points in mm
grid_point_diameter=InputJobFile.Grid_Point_Diameter ;
% This states whether the calibration plate has multiple levels.  If the
% calibration plate does then set this to true.  Otherwise for only one
% level plates, set this to false.
multiple_level_plate=InputJobFile.Multiple_Level_Plate ;
% This is the number of levels that the calbration plate has.  This is the
% number of levels that is visible to one camera at a time, thus the number
% does not double if there are both front and back cameras.
plate_level_number=InputJobFile.Plate_Level_Number;
% This is the spacing between the grid levels.  The first level is always
% equal to 0; the subsequent levels follow the (front)  coordinate system
% convention, ie a value of -10 corresponds to the second level being -10
% distance units from the first plane.
plate_level_spacing=InputJobFile.Plate_Level_Spacing;
% This is the shift in the x location of the grid points from the first level
% of the calibration grid.  This distance should always be taken to be
% positive.
x_plate_level_shift=InputJobFile.X_Plate_Level_Shift;
% This is the shift in the y location of the grid points from the first level
% of the calibration grid.  This distance should always be taken to be
% positive.
y_plate_level_shift=InputJobFile.Y_Plate_Level_Shift;
% This is the thickness of the plate from the front top plane to the back
% top plane
front_to_back_plate_thickness=InputJobFile.Front_To_Back_Plate_Thickness;
% This stores the 'Front' or 'Back orientation for each camera
orientation={'front','back'};
cam1_orient=orientation{InputJobFile.targetsidecam1};
cam2_orient=orientation{InputJobFile.targetsidecam2};
% This is a binary variable stating whether the grid points are white
% points on a black background (if true).  If the variable is false, then
% this specifies black points on a white background.
grid_points_bright=InputJobFile.Grid_Points_Bright;
% This is the unit system ( . . . in mm???)
length_unit=InputJobFile.Length_Unit ;
% This is the camera resolution in pixels
x_pixel_number = InputJobFile.X_Pixel_Number;
y_pixel_number = InputJobFile.Y_Pixel_Number ;

% This is the file extension of the images.
image_extension = InputJobFile.ImageExtension;

% This is the directory that contains all of the camera-specific subdirectories.
% calibration_image_repository = InputJobFile.Calibration_Image_Repository;

% This is a string in stdio printf format specifying the image name. 
% This string is formatted as "base_name_%01d_z_%01d". The base name and 'z', can be changed
% as can the integer format specifiers.
% The integers specify camera number and z-plane number.
% The string should be used together with printf (or sprintf, etc.)
% like "sprintf(image_name_format, 1, 1)" for camera 1, plane 1
% image_name_format = InputJobFile.Image_Name_Format;

% This is a string in stdio printf format specifying the image directory name. 
% It is formatted as "base_name_%02d". The base name and the integer format 
% specifier can be changed .
% The integer specifies camera number.
% The string should be used together with printf (or sprintf, etc.)
% like "sprintf(image_directory_name_format, 1)" for camera 1, etc.
%image_directory_name_format = InputJobFile.Image_Directory_Name_Format;

% This is the name of the calibration output file.
calibration_output_file_name = 'calibration_test_data';%InputJobFile.Calibration_Output_File_Name;

% This is the directory in which the calibration output file will be saved.
calibration_output_directory = 'D:\Stereo Gui design\quickfix1\Prana_Stereo';%InputJobFile.Calibration_Output_Directory;

%%% This is the end of the inputs. %%%

% This is the number of cameras used.
number_of_cameras = length(camera_numbers_list);

% This is the number of calibration planes used.
number_of_planes = length(calibration_planes_list);

% This initializes the image filename cell array
image_filenames = cell(number_of_cameras, number_of_planes);
plane_spacing_temp=zeros(number_of_cameras, number_of_planes);
% This initializes the calibration plate orientation cell array
plate_orientation = cell(number_of_cameras, number_of_planes);
%keyboard;
% This is generates a cell array containing the locations of the image
% files for each camera ii and for each depth jj
for ii=1 : number_of_cameras;
	
	% Determine the number of the current camera.
	%current_camera_number = camera_numbers_list(ii);
	
    % Loop over all the planes
    for jj=1 : number_of_planes;
        
% %         % Determine the number of the current calibration plane.
% %         current_plane_number = calibration_planes_list(jj);
% %         		
% % 		% This is the directory containing the images of the ii'th camera.
% % 		%camera_directory_name = ['Camera_' sprintf('%02d',ii_temp)];
% % 		camera_directory_name = sprintf(image_directory_name_format, current_camera_number);
% % 		
% % 		% This is the name of the the image from the ii'th camera for the jj'th calibration plane z-position.
% % % 		image_file_name = ['Z_camera_' sprintf('%01d',ii_temp) '_z_' sprintf('%01d',jj) image_extension];
% % 		%keyboard;
% %         %image_file_name = sprintf([image_name_format image_extension], current_camera_number, current_plane_number);
% % 		image_file_name = sprintf([image_name_format image_extension],current_plane_number);
% % 
% % 		% This is the full path to the image file.
% %         filename_string = fullfile(calibration_image_repository, camera_directory_name, image_file_name);
% % 	
        % Image file name.
        image_filenames{ii,jj}=InputJobFile.CalImageList{ii,jj};%filename_string;
        %This checks for the plane_spacing, if plane_spacing is uniform in
        %which case it will have one element it assigns uniform plane
        %spacing otherwise stores the different plana_spacings between
        %calibration planes.
        if numel(plane_spacing)==1 && jj>1
            plane_spacing_temp(ii,jj)=(jj-1)*plane_spacing;
        elseif numel(plane_spacing)>1 && jj>1
            plane_spacing_temp(ii,jj)=plane_spacing(jj-1);
        end
                
        % This saves the calibration plate orientation data into the cell array
        if ii==1;
            plate_orientation{ii,jj}=cam1_orient; % takes cam1 orientation
        elseif ii==2;
            plate_orientation{ii,jj}=cam2_orient; % takes cam2 orientation
        elseif ii==3;
            plate_orientation{ii,jj}='front';
        elseif ii==4;
            plate_orientation{ii,jj}='back';
        end;
    end;
end;
% clearing variables and reassigning plane_spacing with an array containing
% plane _spacing for each camera each plane.
clear plane_spacing;
plane_spacing=plane_spacing_temp;
clear plane_spacing_temp;
% keyboard;
% This is the filename to save the calibration data structure to
calibration_output_file_path = fullfile(calibration_output_directory, calibration_output_file_name);

% This creates a structure to save these parameters to
calibration_data=struct;

% This saves the list of camera numbers.
calibration_data.camera_numbers_list = camera_numbers_list;

% This saves the list of plane numbers.
calibration_data.calibration_planes_list = calibration_planes_list;

% This saves the camera number
calibration_data.camera_number = number_of_cameras;

% This saves the number of calibration planes
calibration_data.plane_number = number_of_planes;

% This saves the spacing betweent the calibration planes in mm
calibration_data.plane_spacing = plane_spacing;

% This saves the starting z-depth in mm
calibration_data.z_grid_start = z_grid_start;

% This saves the x spacing between grid points in mm
calibration_data.x_grid_spacing = x_grid_spacing;

% This saves the y spacing between grid points in mm
calibration_data.y_grid_spacing = y_grid_spacing;

% This saves the diameter of the grid points in mm
calibration_data.grid_point_diameter = grid_point_diameter;

% This saves whether the calibration plate has multiple levels
calibration_data.multiple_level_plate = multiple_level_plate;

% This saves the number of calibration plate levels
calibration_data.plate_level_number = plate_level_number;

% This saves the spacing between the different plate levels
calibration_data.plate_level_spacing = plate_level_spacing;

% This saves the shift in the x location of the grid points between the
% different calibration plate levels
calibration_data.x_plate_level_shift = x_plate_level_shift;

% This saves the shift in the y location of the grid points between the
% different calibration plate levels
calibration_data.y_plate_level_shift = y_plate_level_shift;

% This saves the thickness of the plate from the front top plane to the back
% top plane
calibration_data.front_to_back_plate_thickness = front_to_back_plate_thickness;

% This saves a binary variable stating whether the grid points are white
% points on a black background (if true).  If the variable is false, then
% this specifies black points on a white background.
calibration_data.grid_points_bright = grid_points_bright;

% This saves the unit system
calibration_data.length_unit = length_unit;

% This is the camera resolution in pixels
calibration_data.x_pixel_number = x_pixel_number;
calibration_data.y_pixel_number = y_pixel_number;

% This saves the calibration image locations
calibration_data.image_filenames = image_filenames;

% This saves the calibration plate orientation to the camera
calibration_data.plate_orientation = plate_orientation;

% This saves the filename to save the data structure to the data structure
calibration_data.filename_save = calibration_output_file_path;

% This initializes the coordinate data cell array in the structure
calibration_data.coordinate_data = cell(number_of_cameras, number_of_planes);

% This creates the calibration function
[calibration_data,calibration_plane_data]=generate_calibration(calibration_data);


function [calibration_data,calibration_plane_data]=generate_calibration(calibration_data)
% This function generates the calibration data using magic, et cetera . . .


% This extracts the list of camera numbers from the structure.
camera_numbers_list = calibration_data.camera_numbers_list;

% This counts the number of cameras.
number_of_cameras = length(camera_numbers_list);

% This extracts the calibration plane number list from the structure
calibration_planes_list = calibration_data.calibration_planes_list;

% Determine the number of calibration planes.
number_of_planes = length(calibration_planes_list);

% This is the filename to save the data structure to the data structure
filename_save = calibration_data.filename_save;

% This initializes the calibration_plane_data cell array
calibration_plane_data=cell(number_of_cameras,number_of_planes);
% This iterates through the cameras
for ii=1:number_of_cameras;
        
    % This iterates through the calibration planes
    for jj=1:number_of_planes;        
        
        % This lets the user input the coordinate system for the current
        % image
        calibration_plane_data_temp = select_grid_coordinates(calibration_data, ii, jj);
        % This saves the grid data to a cell array
        calibration_plane_data{ii,jj}=calibration_plane_data_temp;
    end;
end;


% This saves the calibration data structure to the harddrive
%save(filename_save, 'calibration_data', 'calibration_plane_data');

% This loads the file back in.
%load(filename_save);

% This iterates through the cameras
for ii = 1:number_of_cameras;
    
    % This iterates through the calibration planes
    for jj = 1:number_of_planes;
                
        % This calculates the image and world coordinates of the grid points
        calibration_data = calculate_grid_coordinates(calibration_data,calibration_plane_data, ii, jj);
    end;
%     % This saves the calibration data structure to the harddrive
%     save(filename_save,'calibration_data','calibration_plane_data');
end;

% This initializes the full calibration coordinate data cell arrays in the
% calibration data structure
calibration_data.x_image_full=cell(number_of_cameras,1);
calibration_data.y_image_full=cell(number_of_cameras,1);
calibration_data.x_world_full=cell(number_of_cameras,1);
calibration_data.y_world_full=cell(number_of_cameras,1);
calibration_data.z_world_full=cell(number_of_cameras,1);
% This iterates through the cameras collecting all the plane data into
% single vectors
for ii=1:number_of_cameras;
    % This initializes the full calibration data vectors
    x_image_full=[];
    y_image_full=[];
    x_world_full=[];
    y_world_full=[];
    z_world_full=[];
    % This iterates through the calibration planes extracting the
    % calibration data from each plane and adding it to the full data
    % vectors
    for jj=1 : number_of_planes;
        % This is the image coodinate data
        x_image_temp=calibration_data.x_image{ii,jj};
        y_image_temp=calibration_data.y_image{ii,jj};
        % This is the world coordinate data
        x_world_temp=calibration_data.x_world{ii,jj};
        y_world_temp=calibration_data.y_world{ii,jj};
        z_world_temp=calibration_data.z_world{ii,jj};
        % This adds the image coordinate data to the full data vector
        x_image_full=[x_image_full;x_image_temp];
        y_image_full=[y_image_full;y_image_temp];
        % This adds the world coordinate data to the full data vector
        x_world_full=[x_world_full;x_world_temp];
        y_world_full=[y_world_full;y_world_temp];
        z_world_full=[z_world_full;z_world_temp];
    end;
    % This adds the full image coordinate data vector to the calibration
    % data structure
    calibration_data.x_image_full{ii}=x_image_full;
    calibration_data.y_image_full{ii}=y_image_full;
    % This adds the full world coordinate data vector to the calibration
    % data structure
    calibration_data.x_world_full{ii}=x_world_full;
    calibration_data.y_world_full{ii}=y_world_full;
    calibration_data.z_world_full{ii}=z_world_full;
end;
% This saves the calibration data structure to the harddrive
%save(filename_save,'calibration_data','calibration_plane_data');

% % This creates camera calibration files for each camera
% save_calibration_data(calibration_data);



function save_calibration_data(calibration_data)
% This function saves the calibration data into a structure containing only
% the image and world coordinates in the form of a 'mat' file with the
% variables x, y, X, Y, Z.  This can then be read into 
% reprojection_tomography_06 for creating a calibration function.

% This extracts the camera number from the structure
camera_number=calibration_data.camera_number;

% This iterates through the cameras creating the calibration data files
for ii=1:camera_number;
    % This extracts the image and world coordinates from the data structure
    x=calibration_data.x_image_full{ii,1};
    y=calibration_data.y_image_full{ii,1};
    X=calibration_data.x_world_full{ii,1};
    Y=calibration_data.y_world_full{ii,1};
    Z=calibration_data.z_world_full{ii,1};
    % This is the filename to save the calibration data as
    %filename_save=['/mnt/current_storage/Projects/Tomo_PIV/101107_Vortex_Ring_Processing/Camera_Position_01/Calibration_01/cam_',num2str(ii),'_calibration.mat'];
    % This saves the calibration data
    save(filename_save,'x','y','X','Y','Z');
end;



function calibration_data=calculate_grid_coordinates(calibration_data,calibration_plane_data,current_camera_number,current_plane_number);
% This function loads the current image and lets the user select points to
% help generate the calibration coordinates.

% This is the number of calibration plate grid levels
plate_level_number=calibration_data.plate_level_number;

% This is the current calibration_plane_data structure
calibration_plane_data_current=calibration_plane_data{current_camera_number,current_plane_number};

% This is the list of camera numbers
camera_numbers_list = calibration_data.camera_numbers_list;

% This is the list of calibration plane numbers
calibration_planes_list = calibration_data.calibration_planes_list;

% This iterates through the different calibration plate grid levels
% calculating the location of the grid points
for ii = 1 : plate_level_number;
    % This is the figure string to display
    figure_string=['Camera ',num2str(camera_numbers_list(current_camera_number)),' Plane ', num2str(calibration_planes_list(current_plane_number)),': Locating the n = ',num2str(ii),' grid system.'];
    % This sets the figure title to tell the use to select the current grid plane
    set(1,'Name',figure_string);
    % Based upon the location of the origin and the user specification of the
    % coordinate system, this function identifies the other grid points in the
    % image that belong to the same coordinate system
    grid_data=locate_grid_points(calibration_plane_data_current.axis_data{ii},calibration_data,current_camera_number,current_plane_number);
    % This saves the current grid point location data to the structure
    calibration_plane_data_current.grid_data{ii}=grid_data;
end;

% This calculates the current grid point coordinates in world and image
% cooridnates
calibration_data=calculate_full_coordinates(calibration_data,calibration_plane_data_current,current_camera_number,current_plane_number);



function calibration_plane_data=select_grid_coordinates(calibration_data,current_camera_number,current_plane_number);
% This function loads the current image and lets the user select points to
% help generate the calibration coordinates.

% This is the list of camera numbers
camera_numbers_list = calibration_data.camera_numbers_list;

% This is the list of calibration plane numbers
calibration_planes_list = calibration_data.calibration_planes_list;

% This is the number of calibration plate grid levels
plate_level_number=calibration_data.plate_level_number;

% This initializes a structure to contain the multiple levels of grid
% planes
calibration_plane_data=struct;
calibration_plane_data.axis_data=cell(plate_level_number,1);
calibration_plane_data.grid_data=cell(plate_level_number,1);

% This iterates through the different calibration plate grid levels letting
% the user select the grid points of this system.
for ii=1:plate_level_number;
    % This is the figure string to display
    figure_string=['Camera ',num2str(camera_numbers_list(current_camera_number)),' Plane ',num2str(calibration_planes_list(current_plane_number)),': Select the n = ',num2str(ii),' grid system.'];
    % This sets the figure title to tell the use to select the current grid plane
%     set(1,'Name',figure_string);
    figure(1);
    set(gcf,'Name',figure_string);
    % This has the user select the grid point origin and an x and y axis.
    axis_data=specify_coordinate(calibration_data,current_camera_number,current_plane_number);
    % This saves the current grid plane data to the structure
    calibration_plane_data.axis_data{ii}=axis_data;
end;



function calibration_data=calculate_full_coordinates(calibration_data,calibration_plane_data,current_camera_number,current_plane_number)
% This function calculates the world and image coordinates for each
% calibration plate location.  This function takes into account multiple
% level calibration grids, camera orientation, and the current calibration
% plate location.

% This is the spacing betweent the calibration planes for particular camera
% number and particular plane
plane_spacing=calibration_data.plane_spacing(current_camera_number,current_plane_number);
% This is the starting z-depth
z_grid_start=calibration_data.z_grid_start;
% This is the orientation of the current calibration image
plate_orientation=calibration_data.plate_orientation{current_camera_number,current_plane_number};
% This is the number of calibration plate grid levels
plate_level_number=calibration_data.plate_level_number;
% This is the x spacing between grid points
x_grid_spacing=calibration_data.x_grid_spacing;
% This is the y spacing between grid points
y_grid_spacing=calibration_data.y_grid_spacing;
% This is the spacing between the different plate levels
plate_level_spacing=calibration_data.plate_level_spacing;
% This the shift in the x location of the grid points between the
% different calibration plate levels
x_plate_level_shift=calibration_data.x_plate_level_shift;
% This is the shift in the y location of the grid points between the
% different calibration plate levels
y_plate_level_shift=calibration_data.y_plate_level_shift;
% This is the thickness of the plate from the front top plane to the back
% top plane
front_to_back_plate_thickness=calibration_data.front_to_back_plate_thickness;

% This initializes the image and world coordinate vectors
x_image=[];
y_image=[];
x_world=[];
y_world=[];
z_world=[];
% This calculates the world coordinates if the camera is facing the front
% side of the plate; if not the world coodinates are calculated for the
% back side of the plate.
if strcmp(plate_orientation,'front');
    % This iterates through the calibration grid levels calculating the Z
    % world coodinates
    for ii=1:plate_level_number;
        % This extracts the current plate level grid point locations
        grid_data=calibration_plane_data.grid_data{ii};
        % These are the subpixel grid locations in image coordinates of the
        % current grid level
        x_image_temp=grid_data.x_subpixel_grid;
        y_image_temp=grid_data.y_subpixel_grid;
        % These are the grid locations in index form
        x_world_temp=x_grid_spacing*grid_data.x_index_grid+x_plate_level_shift(ii);
        y_world_temp=y_grid_spacing*grid_data.y_index_grid+y_plate_level_shift(ii);
        % This calculates the Z world coordinate of the current calbration
        % level/position
        z_world_temp=(z_grid_start+plane_spacing)*ones(size(x_image_temp))+plate_level_spacing(ii);

        % This adds the current calibration points to the complete set of
        % calibration points
        x_image=[x_image;x_image_temp];
        y_image=[y_image;y_image_temp];
        x_world=[x_world;x_world_temp];
        y_world=[y_world;y_world_temp];
        z_world=[z_world;z_world_temp];
    end;
    % This saves the current calibration data to the calibration data
    % structure
    calibration_data.x_image{current_camera_number,current_plane_number}=x_image;
    calibration_data.y_image{current_camera_number,current_plane_number}=y_image;
    calibration_data.x_world{current_camera_number,current_plane_number}=x_world;
    calibration_data.y_world{current_camera_number,current_plane_number}=y_world;
    calibration_data.z_world{current_camera_number,current_plane_number}=z_world;
elseif strcmp(plate_orientation,'back');
    % This iterates through the calibration grid levels calculating the Z
    % world coodinates
    for ii=1:plate_level_number;
        % This extracts the current plate level grid point locations
        grid_data=calibration_plane_data.grid_data{ii};
        % These are the subpixel grid locations in image coordinates of the
        % current grid level
        x_image_temp=grid_data.x_subpixel_grid;
        y_image_temp=grid_data.y_subpixel_grid;
        % These are the grid locations in index form
        x_world_temp=x_grid_spacing*grid_data.x_index_grid+x_plate_level_shift(ii);
        y_world_temp=y_grid_spacing*grid_data.y_index_grid+y_plate_level_shift(ii);
        % This calculates the Z world coordinate of the current calbration
        % level/position
        z_world_temp=(z_grid_start+plane_spacing-(front_to_back_plate_thickness-abs(diff(plate_level_spacing))))*ones(size(x_image_temp))+plate_level_spacing(ii);
        % the z_world_temp is calculated with taking into account the front
        % to back plane thickness for a zigzag type multilevel target:
        % meaning one which the same calibration dots on the inner plane on
        % one side and outer plane on the other side.
        % the case in which same dots are on inner side on both sides of
        % the target needs to be taken care of.
        % This adds the current calibration points to the complete set of
        % calibration points
        x_image=[x_image;x_image_temp];
        y_image=[y_image;y_image_temp];
        x_world=[x_world;x_world_temp];
        y_world=[y_world;y_world_temp];
        z_world=[z_world;z_world_temp];
    end;
    % This saves the current calibration data to the calibration data
    % structure
    calibration_data.x_image{current_camera_number,current_plane_number}=x_image;
    calibration_data.y_image{current_camera_number,current_plane_number}=y_image;
    calibration_data.x_world{current_camera_number,current_plane_number}=x_world;
    calibration_data.y_world{current_camera_number,current_plane_number}=y_world;
    calibration_data.z_world{current_camera_number,current_plane_number}=z_world;
end;
        
        


function grid_data=locate_grid_points(axis_data,calibration_data,current_camera_number,current_plane_number);
% This funtion attempts to locate the grid points that are in the same
% coordinate system as the user specified origin, and x and y axis.  The
% function uses the distance to the first grid points to estimate the
% location of nearby grid points.  In each succesive level the estimated
% distance to the next iteration is refined by previous data.  The
% subpixel grid locations are found by a minimization algorithm.

% This loads the origin data of the current image
origin_vector=axis_data.origin_vector;
% This loads the x-axis vector of the current image
x_axis_vector=axis_data.x_axis_vector;
% This loads the y-axis vector of the current image
y_axis_vector=axis_data.y_axis_vector;
% This extracts the x spacing between grid points
x_grid_spacing=calibration_data.x_grid_spacing;
% This extracts the y spacing between grid points
y_grid_spacing=calibration_data.y_grid_spacing;
% This extracts the diameter of the grid points
grid_point_diameter=calibration_data.grid_point_diameter;

% This estimates the magnification of the calibration image (in length_unit
% per pixel)
x_subpixel_magnification=x_grid_spacing/norm(x_axis_vector);
y_subpixel_magnification=y_grid_spacing/norm(y_axis_vector);
% This is the average of the two magnifications (because they are likely
% the same and thus the average is more accurate . . . and because I don't
% feel like dealing with different magnifications in each axis)
mean_subpixel_magnification=mean([x_subpixel_magnification,y_subpixel_magnification]);
% This is the estimate of the diameter of the grid points in pixels
subpixel_grid_point_diameter=grid_point_diameter/mean_subpixel_magnification;
% This is the correlation window size
L=subpixel_grid_point_diameter*sqrt(2*pi)/2;

% This is the centroid window jj size
Nx=0.40*norm(x_axis_vector);
% This is the centroid window ii size
Ny=0.40*norm(y_axis_vector);

% This loads the current image
I=imread(calibration_data.image_filenames{current_camera_number,current_plane_number});
% This diplays a high contrast version of the current image
figure(1);
imshow(imadjust(I));

% These are the vectors of known coordinate locations in pixels
X_Subpixel_Grid=[origin_vector(1);origin_vector(1)+x_axis_vector(1);origin_vector(1)+y_axis_vector(1)];
Y_Subpixel_Grid=[origin_vector(2);origin_vector(2)+x_axis_vector(2);origin_vector(2)+y_axis_vector(2)];
% These are the vectors of known coordinate locations in units from the
% origin
X_Index_Grid=[0;1;0];
Y_Index_Grid=[0;0;1];
% These are the x-axis and y-axis vectors from the known point in pixels
Grid_X_Axis=[x_axis_vector;x_axis_vector;x_axis_vector];
Grid_Y_Axis=[y_axis_vector;y_axis_vector;y_axis_vector];

% This plots the origin grid point
hold on;
plot(X_Subpixel_Grid(1),Y_Subpixel_Grid(1),'s','Color',[0,0.7,1],'Markersize',18,'Linewidth',2);
plot(X_Subpixel_Grid(1),Y_Subpixel_Grid(1),'+','Color',[0,0.7,1],'Markersize',18,'Linewidth',2);
text(X_Subpixel_Grid(1)+0.15*norm(x_axis_vector),Y_Subpixel_Grid(1)+0.15*norm(y_axis_vector),'(0,0)','Color','WHITE','FontSize',10);
hold off;
% This plots the user defined coordinate grid points
hold on;
plot(X_Subpixel_Grid(2),Y_Subpixel_Grid(2),'sYELLOW','Markersize',18,'Linewidth',2);
text(X_Subpixel_Grid(2)+0.15*norm(x_axis_vector),Y_Subpixel_Grid(2)+0.15*norm(y_axis_vector),'(1,0)','Color','WHITE','FontSize',10);
plot(X_Subpixel_Grid(3),Y_Subpixel_Grid(3),'sYELLOW','Markersize',18,'Linewidth',2);
text(X_Subpixel_Grid(3)+0.15*norm(x_axis_vector),Y_Subpixel_Grid(3)+0.15*norm(y_axis_vector),'(0,1)','Color','WHITE','FontSize',10);
hold off;

% This is the number of points from the origin to search
point_index=0;
% This iteratively expands the set of known coordinate points about the
% origin
while true;
    
    % This increments the point index
    point_index=point_index+1;

    % This creates a grid of searchable coordinates around the currently
    % known coordinates
    [X_Index_Search,Y_Index_Search]=meshgrid(-point_index:point_index,-point_index:point_index);
    
    % This linearizes the grid
    X_Index_Search=X_Index_Search(:);
    Y_Index_Search=Y_Index_Search(:);
    
    % This is a matrix of X inidices for checking for intersections (grid
    % points with known coordinates that are part of the current search)
    [X_Index_A,X_Index_B]=meshgrid(X_Index_Grid,X_Index_Search);
    % This is a matrix of Y inidices for checking for intersections (grid
    % points with known coordinates that are part of the current search)
    [Y_Index_A,Y_Index_B]=meshgrid(Y_Index_Grid,Y_Index_Search);
    
    % This is a set of logical indices into X_Index_Search and
    % Y_Index_Search that correspond to known grid points
    Index_Search_Intersect=logical(sum((X_Index_A==X_Index_B)&(Y_Index_A==Y_Index_B),2));
    
    % This removes the grid points with known locations from the index
    % search vectors
    X_Index_Search(Index_Search_Intersect)=[];
    Y_Index_Search(Index_Search_Intersect)=[];

    % These are the indices of the extrema values of the index vector
    Index_Search_Extrema_Index=find((X_Index_Search==-point_index)+(X_Index_Search==point_index)+(Y_Index_Search==-point_index)+(Y_Index_Search==point_index));
    
    % This is a logical vector that returns the number of the exterior
    % search points that are within the image.  If none of the search
    % points are within the window, the loop counting point_index will
    % break and the next plane will be searched.
    Grid_Search_Exterior_In_Image=true(length(Index_Search_Extrema_Index),1);
    
%     % This is a logical vector that returns the number of search points
%     % that were within the image.  If none of the search points are within
%     % the window, the loop counting point_index will break.
%     Grid_Search_In_Image=true(length(X_Index_Search),1);

    % This iterates through the search indices find grid point locations
    for ii=1:length(X_Index_Search);
        % Commented plotting the different iterations of searched dots in
        % the target, just plots the first estimated positiona and the
        % final estimated postion
        
        % This is a rough estimate of the location of the grid point based
        % on the user entered coordinate system
        X_First_Location_Search=X_Subpixel_Grid(1)+x_axis_vector(1)*X_Index_Search(ii)+y_axis_vector(1)*Y_Index_Search(ii);
        Y_First_Location_Search=Y_Subpixel_Grid(1)+x_axis_vector(2)*X_Index_Search(ii)+y_axis_vector(2)*Y_Index_Search(ii);
        
%         hold on;
%         plot(X_First_Location_Search,Y_First_Location_Search,'sBLUE','Markersize',18,'Linewidth',2);
%         hold off;
        
        % This is a vector of the distances to known grid points
        D=sqrt((X_First_Location_Search-X_Subpixel_Grid).^2+(Y_First_Location_Search-Y_Subpixel_Grid).^2);
        % This finds the minimum distance to the nearest known grid point
        [~,D_Min_Index]=min(D);
        
        % This is the average axis vector between the origin and the
        % current grid point - this can be used to make a refined guess at
        % the grid point location
        X_Axis_Mean=(x_axis_vector+Grid_X_Axis(D_Min_Index,:))/2;
        Y_Axis_Mean=(y_axis_vector+Grid_Y_Axis(D_Min_Index,:))/2;
        
        % This is a rough estimate of the location of the grid point based
        % on the user entered coordinate system and the previously known
        % grid points
        X_Second_Location_Search=X_Subpixel_Grid(1)+X_Axis_Mean(1)*X_Index_Search(ii)+Y_Axis_Mean(1)*Y_Index_Search(ii);
        Y_Second_Location_Search=Y_Subpixel_Grid(1)+X_Axis_Mean(2)*X_Index_Search(ii)+Y_Axis_Mean(2)*Y_Index_Search(ii);
            
%         hold on;
%         plot(X_Second_Location_Search,Y_Second_Location_Search,'sMAGENTA','Markersize',18,'Linewidth',2);
%         hold off;
        
        % This is a vector of the distances to known grid points
        D=sqrt((X_Second_Location_Search-X_Subpixel_Grid).^2+(Y_Second_Location_Search-Y_Subpixel_Grid).^2);
        % This finds the minimum distance to the nearest known grid point
        [~,D_Min_Index]=min(D);
        
        % This is the refined estimate of the location of the grid point
        % based on previously known grid point locations
        X_Third_Location_Search=X_Subpixel_Grid(D_Min_Index)+Grid_X_Axis(D_Min_Index,1)*(X_Index_Search(ii)-X_Index_Grid(D_Min_Index))+Grid_Y_Axis(D_Min_Index,1)*(Y_Index_Search(ii)-Y_Index_Grid(D_Min_Index));
        Y_Third_Location_Search=Y_Subpixel_Grid(D_Min_Index)+Grid_X_Axis(D_Min_Index,2)*(X_Index_Search(ii)-X_Index_Grid(D_Min_Index))+Grid_Y_Axis(D_Min_Index,2)*(Y_Index_Search(ii)-Y_Index_Grid(D_Min_Index));
        
%         hold on;
%         plot(X_Third_Location_Search,Y_Third_Location_Search,'sRED','Markersize',18,'Linewidth',2);
%         hold off;
        
%         % If the current index is an exterior point then this checks
%         % whether the point is on the image and if not continues to the
%         % next loop
%         if any(ii==Index_Search_Extrema_Index);
%             if (X_Third_Location_Search<1+L)||(X_Third_Location_Search>size(I,2)-L)||(Y_Third_Location_Search<1+L)||(Y_Third_Location_Search>size(I,1)-L);
%                 % This states that the current grid search point was not inside
%                 % the image
%                 Grid_Search_Exterior_In_Image(ii)=false;
%                 % This skips to the next loop
%                 continue;
%             end;
%         end;

        L=subpixel_grid_point_diameter*sqrt(2*pi)/4;
       
        % This checks whether the current point is exterior to the image
        % and if so also checks whether it is an external index search
        % point.  If it is exterior to the image, the function continues to
        % the next loop.
        if (X_Third_Location_Search<1+L)||(X_Third_Location_Search>size(I,2)-L)||(Y_Third_Location_Search<1+L)||(Y_Third_Location_Search>size(I,1)-L);
            % This is the index of the exterior grid point index of the
            % current search point
            Exterior_Grid_Point_Index=find(ii==Index_Search_Extrema_Index);
            % If this point exists as an exterior grid point, then it is
            % set as outside the image
            if not(isempty(Exterior_Grid_Point_Index));
                % This states that the current grid search point was not inside
                % the image
                Grid_Search_Exterior_In_Image(Exterior_Grid_Point_Index)=false;
            end;
            % This skips to the next loop
            continue;
        end;
            
%         % This checks whether the current search point is on the image and
%         % if not continues to the next loop
%         if (X_Third_Location_Search<1+L)||(X_Third_Location_Search>size(I,2)-L)||(Y_Third_Location_Search<1+L)||(Y_Third_Location_Search>size(I,1)-L);
%             % This states that the current grid search point was not inside
%             % the image
%             Grid_Search_In_Image(ii)=false;
%             % This skips to the next loop
%             continue;
%         end;
        
        % This calculates the centroid of the region about the estimated
        % location of the current grid point
        [X_Centroid_Search,Y_Centroid_Search]=calculate_centroid(I,X_Third_Location_Search,Y_Third_Location_Search,Nx,Ny);
        if isnan(X_Centroid_Search) || isnan(Y_Centroid_Search)
            continue;
        end
        hold on;
        plot(X_Centroid_Search,Y_Centroid_Search,'s','Color',[1,0.7,0],'Markersize',18,'Linewidth',2);
        hold off;
        
        % This is a vector of the distances to known grid points
        D=sqrt((X_Centroid_Search-X_Subpixel_Grid).^2+(Y_Centroid_Search-Y_Subpixel_Grid).^2);
        % This finds the minimum distance to the nearest known grid point
        [~,D_Min_Index]=min(D);
        
        % This estimates the magnification of the calibration image (in length_unit
        % per pixel)
        x_subpixel_magnification=x_grid_spacing/norm(Grid_X_Axis(D_Min_Index));
        y_subpixel_magnification=y_grid_spacing/norm(Grid_Y_Axis(D_Min_Index));
        % This is the min of the two magnifications in case the axis vector
        % is nearly zero in one direction to the nearest grid point
        mean_subpixel_magnification=min([x_subpixel_magnification,y_subpixel_magnification]);
        % This is the estimate of the diameter of the grid points in pixels
        subpixel_grid_point_diameter=grid_point_diameter/mean_subpixel_magnification;
        % This is the correlation window size
        L=subpixel_grid_point_diameter*sqrt(2*pi)/2;

        % This calculates the subpixel location of the current grid point
        [X_Subpixel_Search,Y_Subpixel_Search,~]=subpix_region(I,X_Centroid_Search,Y_Centroid_Search,subpixel_grid_point_diameter/2);
        
        % This checks whether the subpixel approximation is within a
        % specified distance of the centroid estimate
        if norm([X_Subpixel_Search-X_Centroid_Search,Y_Subpixel_Search-Y_Centroid_Search])>0.5*subpixel_grid_point_diameter;
            
            hold on;
            plot(X_Subpixel_Search,Y_Subpixel_Search,'s','Color',[0.1,0.7,0],'Markersize',18,'Linewidth',2);
            hold off;
            
            % This skips to the next loop
            continue;
        end;

        % This adds the current grid point to the vectors of known grid
        % point locations
        X_Subpixel_Grid(end+1)=X_Subpixel_Search;
        Y_Subpixel_Grid(end+1)=Y_Subpixel_Search;
        % This adds the vectors of current known coordinate locations in 
        % units from the origin
        X_Index_Grid(end+1)=X_Index_Search(ii);
        Y_Index_Grid(end+1)=Y_Index_Search(ii);
        
        % This is a vector of the distances to known grid points
        D=(X_Subpixel_Search-X_Subpixel_Grid).^2+(Y_Subpixel_Search-Y_Subpixel_Grid).^2;
        % This sorts the distances to the nearest known grid points
        [~,D_Index_Sort]=sort(D,'ascend');
        
        % This initializes the x and y axis vectors as nulls
        X_Near_X_Axis=[];
        Y_Near_X_Axis=[];
        X_Near_Y_Axis=[];
        Y_Near_Y_Axis=[];
        % This iterates through the points nearest to the current grid
        % point trying to find a pair from which x and y axis vectors may
        % be established
        for jj=1:length(D);
            % These are the current indices of the grid point
            X_Index_Current=X_Index_Grid(D_Index_Sort(jj));
            Y_Index_Current=Y_Index_Grid(D_Index_Sort(jj));
            % This is the index of the grid point correponding to a point
            % with indices (ii,jj+1) to establish an x-axis vector (if it
            % exists)
            Index_X_Axis_Positive=find((X_Index_Grid==X_Index_Current+1)&(Y_Index_Grid==Y_Index_Current));
            % This is the index of the grid point correponding to a point
            % with indices (ii,jj-1) to establish an x-axis vector (if it
            % exists)
            Index_X_Axis_Negative=find((X_Index_Grid==X_Index_Current-1)&(Y_Index_Grid==Y_Index_Current));
            % This is the index of the grid point correponding to a point
            % with indices (ii+1,jj) to establish an y-axis vector (if it
            % exists)
            Index_Y_Axis_Positive=find((X_Index_Grid==X_Index_Current)&(Y_Index_Grid==Y_Index_Current+1));
            % This is the index of the grid point correponding to a point
            % with indices (ii-1,jj) to establish an y-axis vector (if it
            % exists)
            Index_Y_Axis_Negative=find((X_Index_Grid==X_Index_Current)&(Y_Index_Grid==Y_Index_Current-1));
            % If a positive or negative index grid point exists, then the 
            % x-axis vector is established
            if not(isempty(Index_X_Axis_Positive))&&(isempty(X_Near_X_Axis));
                % This is the x-axis vector
                X_Near_X_Axis=X_Subpixel_Grid(Index_X_Axis_Positive)-X_Subpixel_Grid(D_Index_Sort(jj));
                Y_Near_X_Axis=Y_Subpixel_Grid(Index_X_Axis_Positive)-Y_Subpixel_Grid(D_Index_Sort(jj));
            elseif not(isempty(Index_X_Axis_Negative))&&(isempty(X_Near_X_Axis));
                % This is the x-axis vector
                X_Near_X_Axis=-X_Subpixel_Grid(Index_X_Axis_Negative)+X_Subpixel_Grid(D_Index_Sort(jj));
                Y_Near_X_Axis=-Y_Subpixel_Grid(Index_X_Axis_Negative)+Y_Subpixel_Grid(D_Index_Sort(jj));
            end;
            % If a positive or negative index grid point exists, then the 
            % y-axis vector is established
            if not(isempty(Index_Y_Axis_Positive))&&(isempty(X_Near_Y_Axis));
                % This is the y-axis vector
                X_Near_Y_Axis=X_Subpixel_Grid(Index_Y_Axis_Positive)-X_Subpixel_Grid(D_Index_Sort(jj));
                Y_Near_Y_Axis=Y_Subpixel_Grid(Index_Y_Axis_Positive)-Y_Subpixel_Grid(D_Index_Sort(jj));
            elseif not(isempty(Index_Y_Axis_Negative))&&(isempty(X_Near_Y_Axis));
                % This is the y-axis vector
                X_Near_Y_Axis=-X_Subpixel_Grid(Index_Y_Axis_Negative)+X_Subpixel_Grid(D_Index_Sort(jj));
                Y_Near_Y_Axis=-Y_Subpixel_Grid(Index_Y_Axis_Negative)+Y_Subpixel_Grid(D_Index_Sort(jj));
            end;
            % This breaks the loop if both the x and y axis vectors are
            % defined
            if not(isempty(X_Near_X_Axis))&&not(isempty(X_Near_Y_Axis));
                % Thie breaks the loop
                break;
            end;
        end;
        % This adds the x-axis and y-axis vectors to the current grid point
        % vector
        Grid_X_Axis(end+1,:)=[X_Near_X_Axis,Y_Near_X_Axis];
        Grid_Y_Axis(end+1,:)=[X_Near_Y_Axis,Y_Near_Y_Axis];
        
        % This is the text string for plotting the current point
        % coordinates
        text_string=['(',num2str(X_Index_Search(ii)),',',num2str(Y_Index_Search(ii)),')'];
        % This plots the current grid point
        hold on;
        plot(X_Subpixel_Search,Y_Subpixel_Search,'sYELLOW','Markersize',18,'Linewidth',2);
        text(X_Subpixel_Search+0.15*norm(x_axis_vector),Y_Subpixel_Search+0.15*norm(y_axis_vector),text_string,'Color','WHITE','FontSize',10);
        hold off;
        drawnow;
        pause(0.1);
        
    end;
    
    % If all the exterior search points were outside the image frame, then 
    % no additional grid points may be found and the loop quits
    if all(not(Grid_Search_Exterior_In_Image));
        % This breaks the infinite loop
        break;
    end;
    
%     % If all the search points were outside the image frame, then no
%     % additional grid points may be found and the loop quits
%     if all(not(Grid_Search_In_Image));
%         % This breaks the infinite loop
%         break;
%     end;

end;

% This creates a structure to output the data to
grid_data=struct;
% This saves the subpixel grid location data to the grid_data structure
grid_data.x_subpixel_grid=X_Subpixel_Grid;
grid_data.y_subpixel_grid=Y_Subpixel_Grid;
% This saves the grid location indexing data to the grid_data structure
grid_data.x_index_grid=X_Index_Grid;
grid_data.y_index_grid=Y_Index_Grid;
% This saves the grid axis vectors data to the grid_data  structure
grid_data.grid_x_axis=Grid_X_Axis;
grid_data.grid_y_axis=Grid_Y_Axis;



function axis_data=specify_coordinate(calibration_data,current_camera_number,current_plane_number);
% This lets the user select an origin grid point for the coordinate system
% on the calibration grid as well as an x-axis coordinate point and a
% y-axis coordinate point.

% This extracts the x spacing between grid points
x_grid_spacing=calibration_data.x_grid_spacing;
% This extracts the y spacing between grid points
y_grid_spacing=calibration_data.y_grid_spacing;
% This extracts the diameter of the grid points
grid_point_diameter=calibration_data.grid_point_diameter;

% This loads the current image
I=imread(calibration_data.image_filenames{current_camera_number, current_plane_number});
% This diplays a high contrast version of the current image
figure(1);
imshow(imadjust(I));

% This displays a title telling the user to select the origin grid point
title('Select origin grid point.');
% This lets the user generate a grid point with the mouse
[X_Mouse_Origin,Y_Mouse_Origin,~]=ginput(1);
% This plots the user selected point
hold on;
h1_origin=plot(X_Mouse_Origin,Y_Mouse_Origin,'o','Color',[1,0,0],'Markersize',18,'Linewidth',2);
h2_origin=plot(X_Mouse_Origin,Y_Mouse_Origin,'+','Color',[1,0,0],'Markersize',18,'Linewidth',2);
hold off;
drawnow;

% This displays a title telling the user to select the next grid point to
% the origin in the x-direction
title('Select first positive x grid point.');
% This lets the user generate a grid point with the mouse
[X_Mouse_X_Axis,Y_Mouse_X_Axis,~]=ginput(1);
% This plots the user selected point
hold on;
h1_x_axis=plot(X_Mouse_X_Axis,Y_Mouse_X_Axis,'o','Color',[1,0.7,0],'Markersize',18,'Linewidth',2);
h2_x_axis=plot(X_Mouse_X_Axis,Y_Mouse_X_Axis,'+','Color',[1,0.7,0],'Markersize',18,'Linewidth',2);
hold off;
drawnow;

% This displays a title telling the user to select the next grid point to
% the origin in the y-direction
title('Select first positive y grid point.');
% This lets the user generate a grid point with the mouse
[X_Mouse_Y_Axis,Y_Mouse_Y_Axis,~]=ginput(1);
% This plots the user selected point
hold on;
h1_y_axis=plot(X_Mouse_Y_Axis,Y_Mouse_Y_Axis,'o','Color',[1,0.7,0],'Markersize',18,'Linewidth',2);
h2_y_axis=plot(X_Mouse_Y_Axis,Y_Mouse_Y_Axis,'+','Color',[1,0.7,0],'Markersize',18,'Linewidth',2);
hold off;
drawnow;

% This estimates the magnification of the calibration image (in length_unit
% per pixel)
x_mouse_magnification=x_grid_spacing/norm([X_Mouse_X_Axis-X_Mouse_Origin,Y_Mouse_X_Axis-Y_Mouse_Origin]);
y_mouse_magnification=y_grid_spacing/norm([X_Mouse_Y_Axis-X_Mouse_Origin,Y_Mouse_Y_Axis-Y_Mouse_Origin]);
% This is the average of the two magnifications (because they are likely
% the same and thus the average is more accurate . . . and because I don't
% feel like dealing with different magnifications in each axis)
mean_mouse_magnification=mean([x_mouse_magnification,y_mouse_magnification]);
% This is the estimate of the diameter of the grid points in pixels
mouse_grid_point_diameter=grid_point_diameter/mean_mouse_magnification;

% This is the x-distance between points (so far)
x_mouse_spacing=abs(X_Mouse_X_Axis-X_Mouse_Origin);
% This is the y-distance between points (so far)
y_mouse_spacing=abs(Y_Mouse_Y_Axis-Y_Mouse_Origin);

% This is the centroid window jj size
Nx=0.40*x_mouse_spacing;
% This is the centroid window ii size
Ny=0.40*y_mouse_spacing;

% This produces a warning if the centroid window might not fully cover the
% grid point
if (mouse_grid_point_diameter>Nx)||(mouse_grid_point_diameter>Ny);
    % This displays the warning
    warning('Calibration:centroid_window_size','The grid point centroid may not be accurately calculated due to a small window size.');
end;

% This calculates the centroid of the origin grid point
[X_Centroid_Origin,Y_Centroid_Origin]=calculate_centroid(I,X_Mouse_Origin,Y_Mouse_Origin,Nx,Ny);
% This calculates the centroid of the x grid point
[X_Centroid_X_Axis,Y_Centroid_X_Axis]=calculate_centroid(I,X_Mouse_X_Axis,Y_Mouse_X_Axis,Nx,Ny);
% This calculates the centroid of the y grid point
[X_Centroid_Y_Axis,Y_Centroid_Y_Axis]=calculate_centroid(I,X_Mouse_Y_Axis,Y_Mouse_Y_Axis,Nx,Ny);

% This calculates the subpixel location of the origin grid point
[X_Subpixel_Origin,Y_Subpixel_Origin,~]=subpix_region(I,X_Centroid_Origin,Y_Centroid_Origin,mouse_grid_point_diameter/2);
% This calculates the subpixel location of the x grid point
[X_Subpixel_X_Axis,Y_Subpixel_X_Axis,~]=subpix_region(I,X_Centroid_X_Axis,Y_Centroid_X_Axis,mouse_grid_point_diameter/2);
% This calculates the subpixel location of the y grid point
[X_Subpixel_Y_Axis,Y_Subpixel_Y_Axis,~]=subpix_region(I,X_Centroid_Y_Axis,Y_Centroid_Y_Axis,mouse_grid_point_diameter/2);

% This creates a structure to store the data in
axis_data=struct;
% This saves the origin data
axis_data.origin_vector=[X_Subpixel_Origin,Y_Subpixel_Origin];
% This saves the x-axis vector
axis_data.x_axis_vector=[X_Subpixel_X_Axis-X_Subpixel_Origin,Y_Subpixel_X_Axis-Y_Subpixel_Origin];
% This saves the y-axis vector
axis_data.y_axis_vector=[X_Subpixel_Y_Axis-X_Subpixel_Origin,Y_Subpixel_Y_Axis-Y_Subpixel_Origin];

% This makes the first plots of the mouse grid points invisible
set(h1_origin,'Visible','off');
set(h2_origin,'Visible','off');
set(h1_x_axis,'Visible','off');
set(h2_x_axis,'Visible','off');
set(h1_y_axis,'Visible','off');
set(h2_y_axis,'Visible','off');
% This plots the centroid estimates of the grid point locations
figure(1);
hold on;
plot(X_Subpixel_Origin,Y_Subpixel_Origin,'o','Color',[1,0,0],'Markersize',18,'Linewidth',2);
plot(X_Subpixel_Origin,Y_Subpixel_Origin,'+','Color',[1,0,0],'Markersize',18,'Linewidth',2);
plot(X_Subpixel_X_Axis,Y_Subpixel_X_Axis,'o','Color',[1,0.7,0],'Markersize',18,'Linewidth',2);
plot(X_Subpixel_X_Axis,Y_Subpixel_X_Axis,'+','Color',[1,0.7,0],'Markersize',18,'Linewidth',2);
plot(X_Subpixel_Y_Axis,Y_Subpixel_Y_Axis,'o','Color',[1,0.7,0],'Markersize',18,'Linewidth',2);
plot(X_Subpixel_Y_Axis,Y_Subpixel_Y_Axis,'+','Color',[1,0.7,0],'Markersize',18,'Linewidth',2);
hold off;
drawnow;
pause(1);



function [X_Centroid,Y_Centroid]=calculate_centroid(I,X,Y,Nx,Ny);
% This function returns the centroid of the image I centered about the
% point (X,Y) with a window of size Ny x Nx.

% These are the index ranges to extract the window from I
ii_min=round(Y-Ny/2);
ii_max=round(Y+Ny/2);
jj_min=round(X-Nx/2);
jj_max=round(X+Nx/2);
% This checks whether the points are beyond the edge of the image and
% changes them appropriately
if ii_min<1;
    ii_min=1;
end;
if ii_max>size(I,1);
    ii_max=size(I,1);
end;
if jj_min<1;
    jj_min=1;
end;
if jj_max>size(I,2);
    jj_max=size(I,2);
end;
% This extracts the region from I
I_Win=double(I(ii_min:ii_max,jj_min:jj_max));
% This is the sum of values in the window
I_Win_Sum=double(sum(I_Win(:)));
% These are the pixel coordinates over the window
[jj_matrix,ii_matrix]=meshgrid(jj_min:jj_max,ii_min:ii_max);
% This is the x value of the centroid
X_Centroid=sum(jj_matrix(:).*I_Win(:))/I_Win_Sum;
% This is the y value of the centroid
Y_Centroid=sum(ii_matrix(:).*I_Win(:))/I_Win_Sum;



function [XC_Sub,YC_Sub,res]=subpix_region(I,XC,YC,R);
% This function extracts the subpixel coordinates of each of the grid
% points.

% This is the initial horizontal scaling factor
sigma_x=2;
% This is the initial vertical scaling factor
sigma_y=2;
% This is the initial region angle
phi=(0*pi/180);
% This is the initial guassian exponent
gamma=8;
% This initializes the fitting parameter vector
P0=zeros(6,1);
% This fills in the currently known initial parameter values
P0(3)=phi;
P0(4)=sigma_x;
P0(5)=sigma_y;
P0(6)=gamma;
% This is the length of a the side of a square that has twice the area as a
% circle of radius R
L=R*sqrt(2*pi);
% This initializes the subpixel coordinate vectors
XC_Sub=zeros(size(XC));
YC_Sub=zeros(size(YC));
% This initializes the residual vector
res=zeros(size(XC));
% This sets the options for the minimization function
options=optimset('Display','off');
% This iterates through the grid points extracting an ROI about each grid
% point
for ii=1:length(XC);
    % These are the initial estimates of the region center
    P0(1)=XC(ii);
    P0(2)=YC(ii);
    % This calculates the current region of interest, its extrema, and the
    % median light and dark regions about the current grid point
    [I_ROI,xmin,xmax,ymin,ymax,v_light,v_dark]=grid_region_of_interest(I,XC(ii),YC(ii),median(L)); 
    % This calculates the best fit of the gaussian function to the current
    % grid point region
    [P,res_temp]=fminsearch(@(P)grid_residual(P,I_ROI,xmin,xmax,ymin,ymax,v_light,v_dark,R(ii)),P0,options);
    % This extracts the subpixel coordinates of the grid point center
    XC_Sub(ii)=P(1);
    YC_Sub(ii)=P(2);
    % This saves the residual of the fit to the residual vector
    res(ii)=res_temp;
    % This is a pause for breaking purposes
    pause(0.001);
end;



function [I_ROI,xmin,xmax,ymin,ymax,v_light,v_dark]=grid_region_of_interest(I,XC,YC,L);
% This function extracts the region of interest about the ii-th grid point
% and returns various parameters about the grid point.

% These are the pixel extrema about the current grid point to exract
xmin=round(XC-L/2);
xmax=round(XC+L/2);
ymin=round(YC-L/2);
ymax=round(YC+L/2);
%Checking limits such that it is not outside the image
if xmin==0
    xmin=1;
end
if ymin==0
    ymin=1;
end
if xmax>size(I,2);
    xmax=size(I,2);
end
if ymax>size(I,1);
    ymax=size(I,1);
end

% This extracts the current ROI
try
    I_ROI=I(ymin:ymax,xmin:xmax);
catch
    keyboard; %if still fails see what is the error
end
% This is the sorted range of intensity values
v_sort=double(sort(I_ROI(:)));
v_index=linspace(-0.5,0.5,length(v_sort))';
% The distribution of intensity values should be roughly bimodal, thus
% the list of sorted intensity values should have a short quickly
% increasing section, a slowly increasing section, a quickly increasing
% section, a slowly increasing section, and finally a short quickly
% increasing section (plot it if you don't believe me . . . ) so to
% find the peaks of the two modes fifth order polynomial is fit to the
% sorted pixel values and the first and last inflection points are
% calculated
p=polyfit(v_index,v_sort,5);
% The inflection points occur where the second derivative equals zero
% or alternatively at
r=roots([20*p(1),12*p(2),6*p(3),2*p(4)]);
% The two inflection points of interest are the minimum and maximum
% values
r_min=min(r);
r_max=max(r);
% These are the indicies of the median dark and light regions
[~,dark_index]=min(abs(v_index-r_min));
[~,light_index]=min(abs(v_index-r_max));
% These are the median dark and light values
v_dark=v_sort(dark_index);
v_light=v_sort(light_index);



function r=grid_residual(P,I_ROI,xmin,xmax,ymin,ymax,v_light,v_dark,R);
% This function returns the sum of the square of the difference between the
% grid point image and the fit "gaussian" profile.

% This is the center point
XC=P(1);
YC=P(2);
% This is the rotation angle
phi=P(3);
% These are the x and y scale factors
sigma_x=P(4);
sigma_y=P(5);
% This is the gaussian exponent
gamma=P(6);
% This generates a grid of points about the grid point center
[X,Y]=meshgrid(xmin:xmax,ymin:ymax);
% This applies a rotation of the coordinate system (incase the grid points
% are being viewed from an angle)
X=(X-XC)*cos(phi)-(Y-YC)*sin(phi)+XC;
Y=(X-XC)*sin(phi)+(Y-YC)*cos(phi)+YC;
% This calculates the "radius" from the center of the circle
rho=sqrt((sigma_x*(X-XC)).^2+(sigma_y*(Y-YC)).^2);
% This calculates the "gaussian" function
G=(v_light-v_dark)*abs(exp(-(rho/(2*R)).^gamma))+v_dark;
% This is the error matrix
E=(double(I_ROI)-G).^2;
% This is the residual
r=sum(E(:));



