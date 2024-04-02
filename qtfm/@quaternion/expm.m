function Y = expm(X)
% Matrix exponential.
% (Quaternion overloading of standard Matlab function.)

% Copyright © 2008, 2010 Stephen J. Sangwine and Nicolas Le Bihan.
% See the file : Copyright.m for further details.

Y = overloadm(mfilename, X);

% TODO Implement a more accurate dedicated algorithm for this function. A
% possible candidate is perhaps given in one of the two articles referenced
% by the Matlab documentation page on expm (Moler and van Loan, 1978/2003;
% or Higham, 2005).

end

% $Id: expm.m 1113 2021-02-01 18:41:09Z sangwine $
