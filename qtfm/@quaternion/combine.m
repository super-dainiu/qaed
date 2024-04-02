function Y = combine(P, varargin)
% COMBINE  Combine terms of identical algebraic structure.
% (Quaternion overloading of standard Matlab function.)

% Copyright © 2020 Stephen J. Sangwine and Nicolas Le Bihan.
% See the file : Copyright.m for further details.

narginchk(1, 4), nargoutchk(0, 1)

if ~isa(expr, 'quaternion') || ~isa(expr.x, 'sym')
    error('First argument must be a quaternion with symbolic components.')
end

Y = overload(mfilename, P, varargin{:});

end

% $Id: combine.m 1113 2021-02-01 18:41:09Z sangwine $
