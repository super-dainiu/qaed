function ind = subsindex(~) %#ok<STOUT>
% Unimplemented subsindex function for quaternions.
%
% Implementation note: if the user defines a variable named
% e.g. s, x, y, z, and then enters s(a) (where a is a quaternion)
% expecting to get the scalar part of a, Matlab will try to use
% the quaternion a to index the user's variable s.

% Copyright © 2005 Stephen J. Sangwine and Nicolas Le Bihan.
% See the file : Copyright.m for further details.

error('Subsindex is not implemented for quaternions. Try help quaternion/subsindex.')

end

% $Id: subsindex.m 1113 2021-02-01 18:41:09Z sangwine $
