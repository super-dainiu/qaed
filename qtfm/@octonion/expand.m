function S = expand(expr, varargin)
% EXPAND  Expand expressions in octonion components.
% (Octonion overloading of standard Matlab function.)

% Copyright © 2020 Stephen J. Sangwine and Nicolas Le Bihan.
% See the file : Copyright.m for further details.

narginchk(1, inf), nargoutchk(0, 1)

S = overload(mfilename, expr, varargin{:});

end

% $Id: expand.m 1122 2021-08-13 19:09:07Z sangwine $

% Created automatically from the quaternion
% function of the same name on 13-Aug-2021.