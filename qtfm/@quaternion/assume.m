function assume(expr, varargin)
% ASSUME  Assume condition on quaternion.
% (Quaternion overloading of standard Matlab function.)

% Copyright © 2020 Stephen J. Sangwine and Nicolas Le Bihan.
% See the file : Copyright.m for further details.

narginchk(2, Inf), nargoutchk(0, 0)

if ~isa(expr, 'quaternion') || ~isa(expr.x, 'sym')
    error('First argument must be a quaternion with symbolic components.')
end

overload(mfilename, expr, varargin{:});

end

% $Id: assume.m 1113 2021-02-01 18:41:09Z sangwine $
