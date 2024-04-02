function p = part(q, n)
% PART  Extracts the n-th component of a quaternion.
% This may be empty if the quaternion is pure.

% Copyright © 2013 Stephen J. Sangwine and Nicolas Le Bihan.
% See the file : Copyright.m for further details.

narginchk(2, 2), nargoutchk(0, 1)

if ~isnumeric(n)
    error('Second parameter must be numeric.')
end

switch n
    case 1
        p = q.w; % This could be empty, of course.
    case 2
        p = q.x;
    case 3
        p = q.y;
    case 4
        p = q.z;
    otherwise
        error('Second parameter must be an integer in 1:4.')
end

end

% $Id: part.m 1113 2021-02-01 18:41:09Z sangwine $
