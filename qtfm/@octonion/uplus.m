function u = uplus(a)
% +  Unary plus.
% (Octonion overloading of standard Matlab function.)

% Copyright © 2011 Stephen J. Sangwine and Nicolas Le Bihan.
% See the file : Copyright.m for further details.

narginchk(1, 1), nargoutchk(0, 1) 

u = a; % Since + does nothing, we can just return a.

end

% $Id: uplus.m 1113 2021-02-01 18:41:09Z sangwine $
