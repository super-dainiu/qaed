function n = length(o)
% LENGTH   Length of vector.
% (Octonion overloading of standard Matlab function.)

% Copyright © 2011 Stephen J. Sangwine and Nicolas Le Bihan.
% See the file : Copyright.m for further details.

narginchk(1, 1), nargoutchk(0, 1)

n = length(o.a); % This calls the quaternion length function.

end

% $Id: length.m 1113 2021-02-01 18:41:09Z sangwine $
