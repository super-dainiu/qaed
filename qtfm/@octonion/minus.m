function o = minus(l, r)
% -   Minus.
% (Octonion overloading of standard Matlab function.)

% Copyright © 2011 Stephen J. Sangwine and Nicolas Le Bihan.
% See the file : Copyright.m for further details.

narginchk(2, 2), nargoutchk(0, 1) 

o = plus(l, -r); % We use uminus to negate the right argument.

end

% $Id: minus.m 1113 2021-02-01 18:41:09Z sangwine $
