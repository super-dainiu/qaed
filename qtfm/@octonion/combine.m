function Y = combine(P, varargin)
% COMBINE  Combine terms of identical algebraic structure.
% (Octonion overloading of standard Matlab function.)

% Copyright © 2020 Stephen J. Sangwine and Nicolas Le Bihan.
% See the file : Copyright.m for further details.

narginchk(1, 4), nargoutchk(0, 1)

if ~isa(expr, 'octonion') || ~isa(expr.x, 'sym')
    error('First argument must be an octonion with symbolic components.')
end

Y = overload(mfilename, P, varargin{:});

end

% $Id: combine.m 1122 2021-08-13 19:09:07Z sangwine $

% Created automatically from the quaternion
% function of the same name on 13-Aug-2021.