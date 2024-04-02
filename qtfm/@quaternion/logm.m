function L = logm(A)
% Matrix logarithm.
% (Quaternion overloading of standard Matlab function.)

% Copyright © 2008, 2010 Stephen J. Sangwine and Nicolas Le Bihan.
% See the file : Copyright.m for further details.

L = overloadm(mfilename, A);

% TODO Tests have shown that this function sometimes returns a totally
% incorrect result, reason unknown. The problem seems to be more common
% with larger matrices (starting at around 8-by-8). As a result, the
% following will not be close to a zero matrix: logm(expm(A))-A

% TODO Implement a more accurate dedicated algorithm for this function.

end

% $Id: logm.m 1113 2021-02-01 18:41:09Z sangwine $
