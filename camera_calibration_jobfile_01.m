function JOBFILE = camera_calibration_jobfile_01
% This function generates a job file containing the parameters needed to create a multi-camera calibration file.
% The output of this function should be used as an input to the function camera_calibration_01.m
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

% This is the directory that contains all of the camera-specific subdirectories.
calibration_image_repository = '';

% This is the file extension of the images.
image_extension = '.tif';

% This is the naming format of the calibration images.
image_name_format = 'z_%02d';

% This is the naming format of the calibration image directories.
image_directory_name_format = 'Camera_%02d';

% This is the name of the camera calibration output file.
calibration_output_file_name = '';

% This is the directory in which to save the camera calibration output file.
calibration_output_directory = '.';

% This is a vector containing the camera numbers, i.e., [1, 2, 3] or [1, 3]
camera_numbers_list = [4, 2];
% This is the number of calibration planes
plane_numbers_list = [1,2,3,4,5,6,7];
% This is the spacing betweent the calibration planes in mm
plane_spacing = 1;
% This is the starting z depth in mm
z_grid_start = -3;
% This is the x spacing between grid points in mm
x_grid_spacing = 15;
% This is the y spacing between grid points in mm
y_grid_spacing = 15;
% This is the diameter of the grid points in mm
grid_point_diameter = 3.2;
% This states whether the calibration plate has multiple levels.  If the
% calibration plate does then set this to true.  Otherwise for only one
% level plates, set this to false.
multiple_level_plate = 'true';
% This is the number of levels that the calbration plate has.  This is the
% number of levels that is visible to one camera at a time, thus the number
% does not double if there are both front and back cameras.
plate_level_number = 2;
% This is the spacing between the grid levels.  The first level is always
% equal to 0; the subsequent levels follow the (front)  coordinate system
% convention, ie a value of -10 corresponds to the second level being -10
% distance units from the first plane.
plate_level_spacing = [0,-3];
% This is the shift in the x location of the grid points from the first level
% of the calibration grid.  This distance should always be taken to be
% positive.
x_plate_level_shift = [0,7.5];
% This is the shift in the y location of the grid points from the first level
% of the calibration grid.  This distance should always be taken to be
% positive.
y_plate_level_shift = [0,7.5];
% This is the thickness of the plate from the front top plane to the back
% top plane
front_to_back_plate_thickness = 14.3;
% This is a binary variable stating whether the grid points are white
% points on a black background (if true).  If the variable is false, then
% this specifies black points on a white background.
grid_points_bright = true;

% This is the unit system ( . . . in mm???)
length_unit = 'mm';

% This is the camera resolution in pixels
x_pixel_number = 1024;
y_pixel_number = 1024;

% Populate the fields of the output structure.
JOBFILE.Camera_Numbers_List					= camera_numbers_list;
JOBFILE.Plane_Numbers_List 					= plane_numbers_list;
JOBFILE.Plane_Spacing 						= plane_spacing;
JOBFILE.Z_Grid_Start   						= z_grid_start;
JOBFILE.X_Grid_Spacing 						= x_grid_spacing;
JOBFILE.Y_Grid_Spacing 						= y_grid_spacing;
JOBFILE.Grid_Point_Diameter 				= grid_point_diameter;
JOBFILE.Multiple_Level_Plate 				= multiple_level_plate;
JOBFILE.Plate_Level_Number 					= plate_level_number;
JOBFILE.Plate_Level_Spacing 				= plate_level_spacing;
JOBFILE.X_Plate_Level_Shift 				= x_plate_level_shift;
JOBFILE.Y_Plate_Level_Shift 				= y_plate_level_shift;
JOBFILE.Front_To_Back_Plate_Thickness 		= front_to_back_plate_thickness;
JOBFILE.Grid_Points_Bright 					= grid_points_bright;
JOBFILE.Length_Unit 						= length_unit;
JOBFILE.X_Pixel_Number 						= x_pixel_number;
JOBFILE.Y_Pixel_Number 						= y_pixel_number;
JOBFILE.ImageExtension						= image_extension;
JOBFILE.Calibration_Image_Repository 		= calibration_image_repository;
JOBFILE.Image_Name_Format 					= image_name_format;
JOBFILE.Image_Directory_Name_Format 	  	= image_directory_name_format;
JOBFILE.Calibration_Output_File_Name      	= calibration_output_file_name;
JOBFILE.Calibration_Output_Directory        = calibration_output_directory;

end
