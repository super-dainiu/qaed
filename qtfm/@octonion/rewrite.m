function S = rewrite(expr, varargin)
% REWRITE  Rewrite component expressions in terms of another function.
% (Octonion overloading of standard Matlab function.)

% Copyright © 2020 Stephen J. Sangwine and Nicolas Le Bihan.
% See the file : Copyright.m for further details.

narginchk(1, inf), nargoutchk(0, 1)

if ~isa(expr, 'octonion') || ~isa(expr.x, 'sym')
    error('First argument must be an octonion with symbolic components.')
end

S = overload(mfilename, expr, varargin{:});

end

% $Id: rewrite.m 1124 2021-08-13 19:24:10Z sangwine $

% Created automatically from the quaternion
% function of the same name on 13-Aug-2021.