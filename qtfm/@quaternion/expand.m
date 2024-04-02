function S = expand(expr, varargin)
% EXPAND  Expand expressions in quaternion components.
% (Quaternion overloading of standard Matlab function.)

% Copyright © 2020 Stephen J. Sangwine and Nicolas Le Bihan.
% See the file : Copyright.m for further details.

narginchk(1, inf), nargoutchk(0, 1)

S = overload(mfilename, expr, varargin{:});

end

% $Id: expand.m 1113 2021-02-01 18:41:09Z sangwine $
