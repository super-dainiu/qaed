function tf = isreal(A)
% ISREAL True for real (octonion) array.
% (Octonion overloading of standard Matlab function.)

% Copyright © 2011 Stephen J. Sangwine and Nicolas Le Bihan.
% See the file : Copyright.m for further details.

% This function returns true if both quaternion components of A are real,
% that is, A is an octonion with real coefficients (a real octonion).

narginchk(1, 1), nargoutchk(0, 1)

tf = isreal(A.a) && isreal(A.b);

end

% $Id: isreal.m 1113 2021-02-01 18:41:09Z sangwine $
