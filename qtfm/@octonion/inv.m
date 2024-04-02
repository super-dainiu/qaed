function R = inv(a)
% INV   Inverse of an octonion (matrix).
% (Octonion overloading of standard Matlab function.)

% TODO For a single octonion, this code works, but for matrices, we could
% add 'L' and 'R' parameters, and call linv or rinv as appropriate, with a
% warning if neither is specified. A third option would be to use [L, R] as
% the result, and compute both. Calling code can use [L, ~] or [~, R] to
% compute one only.

% Copyright © 2013, 2022 Stephen J. Sangwine and Nicolas Le Bihan.
% See the file : Copyright.m for further details.

[r, c] = size(a);

if r == 1 && c == 1
    R = conj(a) ./ normo(a);
else
    error(['Octonion matrices have left and right inverses.', ...
           ' See the linv and rinv functions for left and right inverses.'])
end

end

% $Id: inv.m 1141 2022-09-07 14:38:05Z sangwine $
