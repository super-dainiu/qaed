function C = collect(P, expr)
% COLLECT  Collect coefficients.
% (Octonion overloading of standard Matlab function.)

% Copyright © 2020, 2022 Stephen J. Sangwine and Nicolas Le Bihan.
% See the file : Copyright.m for further details.

narginchk(1, 2), nargoutchk(0, 1)

if ~isa(P, 'octonion') || ~isa(x(P.a), 'sym')
    error('First argument must be an octonion with symbolic components.')
end

if nargin == 1
    C = overload(mfilename, P);
else
    C = overload(mfilename, P, expr);
end

end

% $Id: collect.m 1136 2022-02-11 20:53:24Z sangwine $

% Created automatically from the quaternion
% function of the same name on 13-Aug-2021.