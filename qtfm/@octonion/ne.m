function r = ne(a, b)
% ~=  Not equal.
% (Octonion overloading of standard Matlab function.)

% Copyright © 2011 Stephen J. Sangwine and Nicolas Le Bihan.
% See the file : Copyright.m for further details.

narginchk(2, 2), nargoutchk(0, 1) 

r = ~(a == b); % Use the octonion equality operator and invert the result.

end

% $Id: ne.m 1113 2021-02-01 18:41:09Z sangwine $
