function p = vector(o)
% VECTOR   Octonion vector part.

% Copyright © 2011 Stephen J. Sangwine and Nicolas Le Bihan.
% See the file : Copyright.m for further details.

narginchk(1, 1), nargoutchk(0, 1)

p = o; p.a = v(p.a); % Copy, then take the vector part of the a quaternion.

end

% $Id: vector.m 1113 2021-02-01 18:41:09Z sangwine $
