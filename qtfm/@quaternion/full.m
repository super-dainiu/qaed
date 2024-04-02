function b = full(q)
% FULL  Convert sparse matrix to full matrix.
% (Quaternion overloading of standard Matlab function.)

% Copyright © 2013 Stephen J. Sangwine and Nicolas Le Bihan.
% See the file : Copyright.m for further details.

narginchk(1, 1), nargoutchk(0, 1)

b = overload(mfilename, q);

end

% $Id: full.m 1113 2021-02-01 18:41:09Z sangwine $
