function Q = imreadq(varargin)
% IMREADQ  Read image from graphics file and return quaternion image array.
%
% This function takes the same parameters as the Matlab function IMREAD,
% but it returns a quaternion array, with elements of type uint8 or uint16
% depending on the pixel data in the image file (8-bit or 16-bit). The way
% the image pixel data is converted to quaternions depends on the type of
% image stored in the file as follows:
%
% 1. RGB colour image in file. The data is returned as a pure quaternion
%    array with the RGB components in the XYZ components of the quaternion.
%
% No other cases are handled at present. For more precise control over the
% conversion from image to quaternion array, use the Matlab function IMREAD
% and then code the conversion explicitly.

% Copyright © 2009 Stephen J. Sangwine and Nicolas Le Bihan.
% See the file : Copyright.m for further details.

nargoutchk(0, 1) % The number of input parameters is checked by imread, so
                 % we don't check it here.
                 
% Before we read the image, check whether the file contains suitable data.
% By doing this, we are able to provide better diagnostic output than if we
% were to just read the file and try to cope with incorrect data.
                 
info = imfinfo(varargin{1}); % The first parameter should be the filename.

if ~strcmp(info.ColorType, 'truecolor')
   error(['Image file does not contain RGB data, found: ', info.ColorType]) 
end

A = imread(varargin{:}); % Read the specified image (call to standard
                         % Matlab function).
                      
% Now convert the data in A into quaternion format depending on what it is.

N = size(A, 3); % Find out how many image components are in the array A.

if N == 3 % The image is RGB.
    Q = quaternion(A(:,:,1), A(:,:,2), A(:,:,3));
else
    % TODO This error check should now be redundant, since we have checked
    % above for 'truecolor'.
    error('The number of components in the image is not 3. Not handled.');
end

end

% $Id: imreadq.m 1113 2021-02-01 18:41:09Z sangwine $
