function h = image(Q)
% IMAGE  Display pure quaternion array as image.
% (Quaternion overloading of standard Matlab function.)
%
% This function supports only one parameter, unlike the overloaded Matlab
% function of the same name, and the parameter supplied must be a pure
% quaternion array. It will be displayed as an image in a graphics window.
% For more sophisticated requirements, convert the quaternion array to a
% suitable double array with three planes and use the Matlab function.

% Copyright © 2009, 2010, 2011 Stephen J. Sangwine and Nicolas Le Bihan.
% See the file : Copyright.m for further details.

% TODO This function should implement automatic conversion to normalised
% pixel values, when the input is uint8 or uint16, by dividing by 255 or
% 65535 respectively.

narginchk(1, 1), nargoutchk(0, 1)

if ~isempty(Q.w)
    warning(['Quaternion array is not pure, ignoring scalar part.', ...
             ' Largest scalar part is ', num2str(max(max(abs(scalar(Q)))))]);
    Q = vector(Q);
end

C = class(Q.x);

if strcmp(C, 'double') %#ok<*STISA> Because we want to call class once only.
    
    Qname = inputname(1);
    if isempty(Qname), Qname = 'unnamed input parameter'; end
    
    % Clip the pixel values into the range 0.0 -> 1.0, otherwise any values
    % outside this range will raise an error. This is not needed in the
    % uint8 and uint16 cases because the numeric values cannot be outside
    % the range of uint8 or uint16.
    
    Q.x = clip(Q.x, 'x');
    Q.y = clip(Q.y, 'y');
    Q.z = clip(Q.z, 'z');
end

if strcmp(C, 'uint8') || strcmp(C, 'uint16') || strcmp(C, 'double')
    if nargout == 0
        image(cat(3, Q.x, Q.y, Q.z));
    else
        h = image(cat(3, Q.x, Q.y, Q.z));
    end
else
    error(['Quaternion array has unknown component class: ', C,...
           ' - expects uint8/16 or double'])
end

% TODO Maybe better to implement a separate function for gamut mapping, and
% call it from here using a second parameter, with an implicit choice of
% one of the algorithms. Simple clipping is one possible approach, as coded
% below. The code below was altered 24 September 2021 to make it a nested
% function so that the warning messages could output the variable name of
% the input parameter to the main function, for better diagnostics.

function y = clip(x, S)
% Clips the values in x that are less than 0.0 or greater than 1.0.

y = x;

L = x < 0.0;
if any(any(L))
    disp([num2str(sum(sum(L))) ' pixels of ', Qname, ' have ', S, ...
          ' components less than 0.0.']);
    disp(['Worst case was: ' num2str(min(min(x(L))))])
end
y(L) = 0.0;

U = x > 1.0;
if any(any(U))
    disp([num2str(sum(sum(U))) ' pixels of ', Qname, ' have ', S, ...
          ' components greater than 1.0.']);
    disp(['Worst case was: ' num2str(max(max(x(U))))])
end
y(U) = 1.0;
end

end

% $Id: image.m 1131 2021-12-21 21:52:48Z sangwine $
