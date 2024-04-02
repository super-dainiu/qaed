function n = ndims(o)
% NDIMS   Number of array dimensions.
% (Octonion overloading of standard Matlab function.)

% Copyright © 2011 Stephen J. Sangwine and Nicolas Le Bihan.
% See the file : Copyright.m for further details.

n = ndims(o.a); % Call the quaternion ndims on the first quaternion part.
                % (The second quaternion must have the same ndims result.)
end

% $Id: ndims.m 1113 2021-02-01 18:41:09Z sangwine $
