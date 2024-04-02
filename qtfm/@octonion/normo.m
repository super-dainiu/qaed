function a = normo(o)
% NORMO Norm of an octonion, the sum of the squares of the components.

% Copyright © 2011 Stephen J. Sangwine and Nicolas Le Bihan.
% See the file : Copyright.m for further details.

narginchk(1, 1), nargoutchk(0, 1) 

a = normq(o.a) + normq(o.b);

end

% $Id: normo.m 1113 2021-02-01 18:41:09Z sangwine $
