function tf = isnan(A)
% ISNAN  True for Not-a-Number.
% (Octonion overloading of standard Matlab function.)

% Copyright © 2015 Stephen J. Sangwine and Nicolas Le Bihan.
% See the file : Copyright.m for further details.

narginchk(1, 1), nargoutchk(0, 1)

tf = isnan(A.a) | isnan(A.b);

end

% $Id: isnan.m 1113 2021-02-01 18:41:09Z sangwine $
