function unimplemented(name)
% Output error message about unimplemented function.

% Copyright © 2005 Stephen J. Sangwine and Nicolas Le Bihan.
% See the file : Copyright.m for further details.

disp(['Unimplemented quaternion function: ', name]);
error('Call to unimplemented function.');

end

% $Id: unimplemented.m 1113 2021-02-01 18:41:09Z sangwine $
