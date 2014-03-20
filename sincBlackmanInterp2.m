function ZI = sincBlackmanInterp2(Z, XI, YI, KERNELRADIUS, METHOD)
% ZI = sincBlackmanInterp2(Z, XI, YI, KERNELRADIUS)
% Sinc-interpolation with the option for blackman apodization. 
%
% INPUTS
%   X = Matrix of column positions of the original image
%   Y = Matrix of row positions of the original image
%   Z = Matrix of values of the original signal sampled on the [Y, X] grid
%   XI = Matrix of non-integer column positions at which to resample the original signal
%   YI = Matrix of non-integer row positions at which to resample the original signal
%   KERNELRADIUS = Radius of the Sinc interpolation kernel. 
%   METHOD = String specifying the interpolation type. Set METHOD =
%   'blackman' for blackman apodization and 'sinc' for no apodization. 
% 
% OUTPUTS
%   ZI = Matrix of values of the original signal resampled to the [YI, XI] grid.
%
% SEE ALSO
%   sinc
% 

% Default to blackman interpolation
if nargin < 5
    METHOD = 'blackman';
end

%find type of Z so we can make the output match
Zclass = class(Z);

% Image height (number of rows) and width (number of columns)
[height, width] = size(Z);

% Determine the class of the image.
imageClass = class(Z);

% Convert data types to singles
xi = single(XI);
yi = single(YI);

% Clear variables to save memory
clear XI YI

% Number of points to interpolate
nInterpPoints = numel(xi);

% Vector of integers corresponding to the kernel locations
kernelLocs = int16(-KERNELRADIUS : KERNELRADIUS);

% Length of the kernel
kernelLength = length(kernelLocs);

% Reshape coordinates into vectors
XIv = reshape(xi, nInterpPoints, 1);
YIv = reshape(yi, nInterpPoints, 1);

% Mirror the border of the image so that the kernel doesn't reach outside of the valid image
zPadded = padarray(uint16(Z), [KERNELRADIUS KERNELRADIUS], 0);

% Calculate size of padded image
[paddedHeight, paddedWidth] = size(zPadded);

% Save the padded image as a vector (vector operations are faster than
% matrix operations)
zPaddedVector = reshape(zPadded, numel(zPadded), 1);

% Clear a variable to save memory
clear zPadded

% Determine the integer locations of the anchor pixels
yAnchors = int16(round(YIv));
xAnchors = int16(round(XIv));

% Replicate the interpolated coordinate vectors and the anchor pixel
% vectors to save time on loops later. These operations are equivalent to
% "repmat" but faster.
repYI = YIv( : , ones(1, kernelLength) ); % Replicate the Interpolated row positions. This way is faster than REPMAT.
repXI = XIv( : , ones(1, kernelLength) ); % Replicate the Interpolated column positions

% Clear some variables to save memory
clear X Y YI XI YIv XIv Z yi xi;

% Calculate the (possibly non-integer) image-coordinates of the positions at which the
% interpolation function will be evaluated.
interpolantRows = single(kernelLocs( ones(1, nInterpPoints), :) + yAnchors( : , ones( 1, kernelLength ) ) );
interpolantCols = single(kernelLocs( ones(1, nInterpPoints), :) + xAnchors( : , ones( 1, kernelLength ) ) );

% Calculate the values of the intepolant at the interpolation coordinates
% (this line does both X- and Y- coordinates)
Interpolants = sinc([repYI - interpolantRows, repXI - interpolantCols]);

% Pull the X- and Y- interpolants out of the combined matrix
interpolantsY = Interpolants(:, 1 : kernelLength);
interpolantsX = Interpolants(:, kernelLength + 1 : end);

% Clear the Interpolants variable to save memory
clear Interpolants; 

% Calculate the interpolant apodization function
if strcmp(METHOD, 'blackman') % If blackman window was specified
    apodizationY =  0.42 + 0.5 * cos( pi *( interpolantRows - repYI ) / KERNELRADIUS ) + 0.08 * cos (2 * pi * ( interpolantRows- repYI ) / KERNELRADIUS ); % Blackman window ( vertical )
    apodizationX = 0.42 + 0.5 * cos( pi *( interpolantCols - repXI ) / KERNELRADIUS ) + 0.08 * cos (2 * pi * ( interpolantCols - repXI ) / KERNELRADIUS ); % Blackman window ( horizontal)
else % If no apodization function was specified then don't apodize the interpolant
    apodizationY =  ones(size(interpolantRows));
    apodizationX = ones(size(interpolantCols));
end

% Clear some variables to save memory
clear interpolantRows interpolantCols repXI repYI

% Replicate the interpolant and apodization matrices so that they can me
% multiplied by a single array operation
repInterpolantsX = interpolantsX(:, (1 : kernelLength)' * ones(1, kernelLength));
repApodizationX = apodizationX(:, (1 : kernelLength)' * ones(1, kernelLength));

% Clear some variables to save memory
clear interpolantsX apodizationX

% Initialize the vertical interpolant and apodization matrices
repInterpolantsY = zeros(size(repInterpolantsX), 'single');
repApodizationY = zeros(size(repApodizationX), 'single');

% Loop over the columns of the vertical interpolant and apodization
% matrices to populate them. This looping is necessary (or at least I
% haven't thought of another way to do it) because the elements of these
% matrices must appear in the order (a11, a11, a11, ... a12, a12, a12, ...)
% so that multiplying them by the corresponding horizontal interpolants and
% apodization functions results in a re-shaped version of real vector
% multiplication. 
for k = 1 : kernelLength;
    repInterpolantsY(:, kernelLength * (k - 1) + 1 : kernelLength * k ) = interpolantsY(:, k * ones(1, kernelLength ) ) ;
    repApodizationY(:, kernelLength * (k - 1) + 1 : kernelLength * k) = apodizationY(:, k * ones(1, kernelLength ) ) ;
end

% Clear some variables to save memory
clear interpolantsY apodizationY;

% Calculate the pixel intensity weighting function evaluated at each pixel
Weights = repInterpolantsY .* repInterpolantsX .* repApodizationY .* repApodizationX;

% Clear some variables to save memory
clear repInterpolantsY repInterpolantsX repApodizationY repApodizationX 

% Determine which x- and y- anchor positions are valid (i.e. which ones
% fall within the domain of the original image). This step is a memory-saving step: 
%the data class of these new variables is Logical, which is small in memory. 
validY = yAnchors >= 1 & yAnchors <= height;
validX = xAnchors >= 1 & xAnchors <= width;

% Source columns matrix 
kernelAdds = int16(0 : 2 * KERNELRADIUS); % Integers specifying the number of elements from the left-side of the interpolation kernel 
sourceColumns = kernelAdds( ones(1, nInterpPoints), :) + xAnchors( : , ones( 1, kernelLength ) ); % List of the source columns
sourceColumnsRep = uint16(sourceColumns( : , (1 : kernelLength)' * ones(1, kernelLength))); % Replicate the source columns matrix

% Clear a vairable to save memory
clear sourceColumns;

% Initialize the source rows matrix
sourceRowsRep = zeros(size(sourceColumnsRep), 'uint16');

% Determine row positions of source pixels
for k = 1 : kernelLength
    sourceRowsRep(:, kernelLength * (k - 1) + 1 : kernelLength * k) = yAnchors(:, ones(1, kernelLength) ) + (k -1);
end

% Clear some variables to save space. 
clear yAnchors xAnchors

% Set coordinates of points outside of the domain of the original signal to
% the values at the edges of the signal 
sourceRowsRep(sourceRowsRep  < 1) = 1;
sourceColumnsRep(sourceColumnsRep < 1) = 1;
sourceRowsRep(sourceRowsRep > paddedHeight ) = paddedHeight;
sourceColumnsRep(sourceColumnsRep > paddedWidth ) = paddedWidth;

% Values of the of source pixels. The list of incides here is an operation turns the 2-element [m, n]
% coordinates of an array into a single integer index.
% The operation here is equivalent to the matlab built-in command SUB2IND, 
% i.e. sourceIndices = (sourceColumns - 1) * size(zPadded, 1) + sourceRows;
% sourceIndices =  (double(sourceColumnsRep ) - 1) * size(zPadded, 1) + double(sourceRowsRep);
sourcePixels = zPaddedVector((double(sourceColumnsRep ) - 1) * paddedHeight + double(sourceRowsRep));

% Clear some variables to save memory
clear sourceRows sourceColumns

% Multiply the values of the source pixels by the weighting function
% evaluated at each pixel. Doing this as an array operation is faster than
% doing it in a loop because it takes a long time to loop over a million
% goddamn pixels. 
Zinterp = cast( sum( Weights .* single(sourcePixels) , 2 ) , Zclass );

% Clear a variable to save space
clear sourcePixels

% Set the values of the interpolated function whose source pixel values
% were outside of the domain of the original image to zero.
Zinterp( ~validY ) = 0;
Zinterp( ~validX ) = 0;

% Reshape the vector-version of the interpolated signal to a matrix of the
% same size as the original signal
ZI = reshape( Zinterp, [ height width ] );

end % End of function

% % %







