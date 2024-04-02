function d = sym(q)
% SYM Construct symbolic quaternion.
% (Quaternion overloading of standard Matlab function.)

% Copyright © 2021 Stephen J. Sangwine and Nicolas Le Bihan.
% See the file : Copyright.m for further details.

narginchk(1, 1), nargoutchk(0, 1)

d = overload(mfilename, q);

end

% $Id: sym.m 1119 2021-03-08 14:52:23Z sangwine $
