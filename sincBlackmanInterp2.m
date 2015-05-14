function ZI = sincBlackmanInterp2(Z, XI, YI, KERNELSIZE, METHOD)
% ZI = sincBlackmanInterp2(Z, XI, YI, KERNELSIZE, METHOD)
% Sinc-interpolation with the option for blackman apodization. 
%
% INPUTS
%   X = Matrix of column positions of the original image
%   Y = Matrix of row positions of the original image
%   Z = Matrix of values of the original signal sampled on the [Y, X] grid
%   XI = Matrix of non-integer column positions at which to resample the original signal
%   YI = Matrix of non-integer row positions at which to resample the original signal
%   KERNELSIZE = symmetric dimension of the kernel size. Should be even; 
%       odd inputs will be forced to the next largest even integer 
%   METHOD = String specifying the interpolation type. Set METHOD =
%   'blackman' for blackman apodization and 'sinc' for no apodization. 
% 
% OUTPUTS
%   ZI = Matrix of values of the original signal resampled to the [YI, XI] grid.
%
% SEE ALSO
%   sinc
% 

%     This file is part of prana, an open-source GUI-driven program for
%     calculating velocity fields using PIV or PTV.
%     Copyright (C) 2012  Virginia Polytechnic Institute and State
%     University
% 
%     prana is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     (at your option) any later version.
% 
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
% 
%     You should have received a copy of the GNU General Public License
%     along with this program.  If not, see <http://www.gnu.org/licenses/>.


% Default to blackman interpolation
if nargin < 5
    METHOD = 'blackman';
end

% Image height (number of rows) and width (number of columns)
[height, width] = size(Z);

% Determine the data class of the image.
imageClass = class(Z);

% Number of points to interpolate
nInterpPoints = numel(XI);

[nInterpRows, nInterpColumns] = size(XI);

% Reshape the grid of interpolation points into vectors.
% Convert to single precision to save memory.
XIv = single(reshape(XI, nInterpPoints, 1));
YIv = single(reshape(YI, nInterpPoints, 1));

% Determine the integer locations of the anchor pixels.
% Floor them so it's easy to specify the pixel locations locations of an even kernel.
xAnchors = int16(floor(XIv));
yAnchors = int16(floor(YIv));

% clear variables to save memory
clear XI YI

% Force the kernel size to an even integer
kernelSize = KERNELSIZE + 2 * mod(KERNELSIZE/2, 1);

% Mirror the border of the image so that the kernel doesn't reach outside of the valid image
zPadded = padarray(Z, [kernelSize/2, kernelSize/2], 0);

% Calculate size of padded image
[paddedHeight, paddedWidth] = size(zPadded);

% Save the padded image as a vector (vector operations are faster than
% matrix operations)
zPaddedVector = reshape(zPadded, numel(zPadded), 1);

% clear a variable to save memory
clear zPadded

% Vector of integers corresponding to the kernel locations
% The -1 is to specify that there are always
% (kernelSize/2 -1) points to the left of the floored anchor pixel
% and (kernelSize/2) points to the right of the floored anchor pixel.
kernelLocs = int16( -(kernelSize/2 - 1) : (kernelSize/2));

% Replicate the interpolated coordinate vectors and the anchor pixel
% vectors to save time on loops later. These operations are equivalent to
% "repmat" but faster.
repYI = YIv( : , ones(1, kernelSize) ); % Replicate the Interpolated row positions. This way is faster than REPMAT.
repXI = XIv( : , ones(1, kernelSize) ); % Replicate the Interpolated column positions

% clear some variables to save memory
clear YI XI YIv XIv Z;

% Calculate the (possibly non-integer) image-coordinates of the positions at which the
% interpolation function will be evaluated.
% interpolantRows = repmat(kernelLocs, nInterpPoints, 1) + repmat(yAnchors, 1, kernelSize);
interpolantRows = single(kernelLocs( ones(1, nInterpPoints), :) + yAnchors( : , ones( 1, kernelSize ) ) );

% interpolantCols = repmat(kernelLocs, nInterpPoints, 1) + repmat(xAnchors, 1, kernelSize);
interpolantCols = single(kernelLocs( ones(1, nInterpPoints), :) + xAnchors( : , ones( 1, kernelSize ) ) );

% Calculate the values of the intepolant at the interpolation coordinates
% (this line does both X- and Y- coordinates)
% Interpolants = sinc([repYI - interpolantRows, repXI - interpolantCols]);
Interpolants = sinc([interpolantRows - repYI, interpolantCols - repXI]);

% Pull the X- and Y- interpolants out of the combined matrix
interpolantsY = Interpolants(:, 1 : kernelSize);
interpolantsX = Interpolants(:, kernelSize + 1 : end);

% clear the Interpolants variable to save memory
clear Interpolants; 

% Calculate the interpolant apodization function
if strcmp(METHOD, 'blackman') % If blackman window was specified
    apodizationY =  0.42 + 0.5 * cos( pi *( interpolantRows - repYI ) / (kernelSize/2) ) + 0.08 * cos (2 * pi * ( interpolantRows - repYI ) / (kernelSize/2) ); % Blackman window ( vertical )
    apodizationX = 0.42 + 0.5 * cos( pi *( interpolantCols - repXI ) / (kernelSize/2) ) + 0.08 * cos (2 * pi * ( interpolantCols - repXI ) / (kernelSize/2) ); % Blackman window ( horizontal)
else % If no apodization function was specified then don't apodize the interpolant
    apodizationY =  ones(size(interpolantRows));
    apodizationX = ones(size(interpolantCols));
end

% clear some variables to save memory
clear interpolantRows interpolantCols repXI repYI

% Replicate the interpolant and apodization matrices so that they can be
% multiplied by a single array operation
% repInterpolantsX2 = repmat(interpolantsX, 1, kernelSize);
% repApodizationX2 = repmat(apodizationX, 1, kernelSize);
repInterpolantsX = interpolantsX(:, (1 : kernelSize)' * ones(1, kernelSize));
repApodizationX = apodizationX(:, (1 : kernelSize)' * ones(1, kernelSize));

% clear some variables to save memory
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
for k = 1 : kernelSize;
    repInterpolantsY(:, kernelSize * (k - 1) + 1 : kernelSize * k ) = interpolantsY(:, k * ones(1, kernelSize ) ) ;
    repApodizationY(:, kernelSize * (k - 1) + 1 : kernelSize * k) = apodizationY(:, k * ones(1, kernelSize ) ) ;
end

% clear some variables to save memory
clear interpolantsY apodizationY;

% Calculate the pixel intensity weighting function evaluated at each pixel
Weights = repInterpolantsY .* repInterpolantsX .* repApodizationY .* repApodizationX;

% clear some variables to save memory
clear repInterpolantsY repInterpolantsX repApodizationY repApodizationX 

% Determine which x- and y- anchor positions are valid (i.e. which ones
% fall within the domain of the original image). This step is a memory-saving step: 
%the data class of these new variables is Logical, which is small in memory. 
validY = yAnchors >= 1 & yAnchors <= height;
validX = xAnchors >= 1 & xAnchors <= width;

% Source columns matrix 
% Integers specifying the number of elements from the left-side of the interpolation kernel
% kernelAdds = int16(0 : (kernelSize - 1));
kernelAdds = int16(1 : kernelSize);

% List of the source columns
% sourceColumns = repmat(kernelAdds, nInterpPoints, 1) + repmat(xAnchors, 1, kernelSize)
sourceColumns = kernelAdds( ones(1, nInterpPoints), :) + xAnchors( : , ones( 1, kernelSize ) ); 

% Replicate the source columns matrix
% sourceColumnsRep = repmat(sourceColumns, 1, kernelSize)
sourceColumnsRep = uint16(sourceColumns( : , (1 : kernelSize)' * ones(1, kernelSize))); 

% clear a vairable to save memory
clear sourceColumns;

% Initialize the source rows matrix
sourceRowsRep = zeros(size(sourceColumnsRep), 'uint16');

% Determine row positions of source pixels
for k = 1 : kernelSize
%     sourceRowsRep(:, kernelSize * (k - 1) + 1 : kernelSize * k) = yAnchors(:, ones(1, kernelSize) ) + (k - 1);
    sourceRowsRep(:, kernelSize * (k - 1) + 1 : kernelSize * k) = yAnchors(:, ones(1, kernelSize) ) + kernelAdds(k);
end

% clear some variables to save space. 
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

% clear some variables to save memory
clear sourceRows sourceColumns

% Multiply the values of the source pixels by the weighting function
% evaluated at each pixel. Doing this as an array operation is faster than
% doing it in a loop because it takes a long time to loop over a million
% goddamn pixels. 
Zinterp = sum( Weights .* single(sourcePixels) , 2 );

% clear a variable to save space
clear sourcePixels

% Set the values of the interpolated function whose source pixel values
% were outside of the domain of the original image to zero.
Zinterp( ~validY ) = 0;
Zinterp( ~validX ) = 0;

% Reshape the vector-version of the interpolated signal to a matrix of the
% same size as the original signal. Save the output image to the output
% variable with the same class as the input image.
ZI = reshape( cast(Zinterp, imageClass), [ nInterpRows, nInterpColumns ] );

end % End of function

% % %







