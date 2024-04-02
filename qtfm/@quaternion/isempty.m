function tf = isempty(q)
% ISEMPTY True for empty matrix.
% (Quaternion overloading of standard Matlab function.)

% Copyright © 2005 Stephen J. Sangwine and Nicolas Le Bihan.
% See the file : Copyright.m for further details.

% narginchk(1, 1), nargoutchk(0, 1)

% It is sufficient to check the x component, because if x is empty so must
% be the y and z components. We must not check the scalar component,
% because this is empty for a non-empty pure quaternion.
     
tf = isempty(q.x);

end

% $Id: isempty.m 1113 2021-02-01 18:41:09Z sangwine $

