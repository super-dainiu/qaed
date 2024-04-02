function a = eval(q)
% EVAL   Evaluate components of quaternion.
% (Quaternion overloading of standard Matlab function.)

% Copyright © 2020 Stephen J. Sangwine and Nicolas Le Bihan.
% See the file : Copyright.m for further details.

narginchk(1, 1), nargoutchk(0, 1)

if ~isa(q.x, 'sym')
    error('Argument must be quaternion with sym components')
end

a = overload(mfilename, q);

end

% $Id: eval.m 1113 2021-02-01 18:41:09Z sangwine $
