function S = simplify(expr, varargin)
% SIMPLIFY  Algebraic simplification of symbolic quaternion.
% (Quaternion overloading of standard Matlab function.)

% Copyright © 2020, 2021 Stephen J. Sangwine and Nicolas Le Bihan.
% See the file : Copyright.m for further details.

narginchk(1, inf), nargoutchk(0, 1)

if ~isa(expr, 'quaternion') || ~isa(expr.x, 'sym')
    error('First argument must be a quaternion with symbolic components.')
end

S = overload(mfilename, expr, varargin{:});

% S = elide(S);

end

% $Id: simplify.m 1127 2021-08-30 19:35:41Z sangwine $
