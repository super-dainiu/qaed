function d = sym(q)
% SYM Construct symbolic octonion.
% (Octonion overloading of standard Matlab function.)

% Copyright © 2021 Stephen J. Sangwine and Nicolas Le Bihan.
% See the file : Copyright.m for further details.

narginchk(1, 1), nargoutchk(0, 1)

d = overload(mfilename, q);

end

% $Id: sym.m 1122 2021-08-13 19:09:07Z sangwine $

% Created automatically from the quaternion
% function of the same name on 13-Aug-2021.