function r = ldivide(~, ~) %#ok<STOUT>
% Quaternion left elementwise division is not implemented. This is because
% Matlab defines A ./ B to be the matrix with elements A(i, j)/B(i, j) and
% A .\ B to be the matrix with elements B(i, j)/A(i, j). This is not
% consistent with the ideas of left and right division in a non-commutative
% algebra, therefore it is not supported for quaternion arrays. The ./
% operator implements elementwise division on the right, as expected.
% To obtain elementwise division on the left, multiply on the left by the
% elementwise inverse of B. This can be computed as conj(B) ./ normq(B).

% Copyright © 2005 Stephen J. Sangwine and Nicolas Le Bihan.
% See the file : Copyright.m for further details.

help quaternion/ldivide;
error('Elementwise left division is not implemented for quaternions.');

end

% $Id: ldivide.m 1154 2022-10-13 09:27:42Z sangwine $
