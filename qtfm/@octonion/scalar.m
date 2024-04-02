function r = scalar(o)
% SCALAR   Octonion scalar part.
%
% This function returns zero in the case of pure octonions,
% whereas the function s gives empty if o is pure.

% Copyright © 2011 Stephen J. Sangwine and Nicolas Le Bihan.
% See the file : Copyright.m for further details.

narginchk(1, 1), nargoutchk(0, 1)

r = scalar(o.a); % The scalar part is in the 'a' quaternion.

end

% $Id: scalar.m 1113 2021-02-01 18:41:09Z sangwine $
