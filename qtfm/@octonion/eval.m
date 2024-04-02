function a = eval(q)
% EVAL   Evaluate components of octonion.
% (Octonion overloading of standard Matlab function.)

% Copyright © 2020 Stephen J. Sangwine and Nicolas Le Bihan.
% See the file : Copyright.m for further details.

narginchk(1, 1), nargoutchk(0, 1)

if ~isa(x(q.a), 'sym')
    error('Argument must be octonion with sym components')
end

a = overload(mfilename, q);

end

% $Id: eval.m 1127 2021-08-30 19:35:41Z sangwine $

% Created automatically from the quaternion
% function of the same name on 13-Aug-2021.