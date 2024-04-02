function p = imag(q)
% IMAG   Imaginary part of a quaternion.
% (Quaternion overloading of standard Matlab function.)
%
% This function returns the quaternion that is the imaginary
% part of q. If q is a real quaternion, it returns zero.

% Copyright © 2005 Stephen J. Sangwine and Nicolas Le Bihan.
% See the file : Copyright.m for further details.

narginchk(1, 1), nargoutchk(0, 1)

p = overload(mfilename, q);

end

% $Id: imag.m 1113 2021-02-01 18:41:09Z sangwine $
