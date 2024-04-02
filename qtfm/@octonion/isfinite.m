function tf = isfinite(A)
% ISFINITE  True for finite elements.  
% (Octonion overloading of standard Matlab function.)

% Copyright © 2012 Stephen J. Sangwine and Nicolas Le Bihan.
% See the file : Copyright.m for further details.

narginchk(1, 1), nargoutchk(0, 1)

tf = isfinite(A.a) & isfinite(A.b);

end

% $Id: isfinite.m 1113 2021-02-01 18:41:09Z sangwine $
