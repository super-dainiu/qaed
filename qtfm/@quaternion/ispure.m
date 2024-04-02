function tf = ispure(q)
% ISPURE   Tests whether a quaternion is pure.
% tf = ispure(q) returns true if q has an empty scalar part, 0 otherwise.
% Note that if q has a scalar part which is zero, ispure(q) returns false.
% Also, ispure(q) returns true if q is an empty quaternion, since the
% scalar part is empty (that is, pure and full empty quaternions cannot be
% distinguished from each other).

% Copyright © 2005, 2010 Stephen J. Sangwine and Nicolas Le Bihan.
% See the file : Copyright.m for further details.

narginchk(1, 1), nargoutchk(0, 1)

tf = isempty(q.w);

end

% $Id: ispure.m 1113 2021-02-01 18:41:09Z sangwine $
