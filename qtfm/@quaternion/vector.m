function p = vector(q)
% VECTOR   Quaternion vector part. Synonym of V.

% Copyright © 2005, 2010, 2020 Stephen J. Sangwine and Nicolas Le Bihan.
% See the file : Copyright.m for further details.

narginchk(1, 1), nargoutchk(0, 1)

p = q; p.w = cast([], class(p.x)); % Copy and set the scalar part to empty.

end

% $Id: vector.m 1113 2021-02-01 18:41:09Z sangwine $
