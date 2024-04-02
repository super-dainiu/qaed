function b = repmat(a, m, n)
% REPMAT Replicate and tile an array.
% (Octonion overloading of standard Matlab function.)

% Copyright © 2013 Stephen J. Sangwine and Nicolas Le Bihan.
% See the file : Copyright.m for further details.

narginchk(2, 3), nargoutchk(0, 1) 

if nargin == 2
    b = overload(mfilename, a, m);
else
    b = overload(mfilename, a, m, n);
end

end

% $Id: repmat.m 1113 2021-02-01 18:41:09Z sangwine $
