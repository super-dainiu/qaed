function X = sqrtm(A)
% Matrix square root.
% (Quaternion overloading of standard Matlab function.)

% Copyright © 2008, 2010 Stephen J. Sangwine and Nicolas Le Bihan.
% See the file : Copyright.m for further details.

X = overloadm(mfilename, A);

end

% TODO Implement a more accurate dedicated algorithm for this function.

% $Id: sqrtm.m 1113 2021-02-01 18:41:09Z sangwine $
