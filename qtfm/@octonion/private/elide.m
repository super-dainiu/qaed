function S = elide(E)
% ELIDE  Remove redundant parts of symbolic octonion.

% Copyright © 2021 Stephen J. Sangwine and Nicolas Le Bihan.
% See the file : Copyright.m for further details.

narginchk(1, 1), nargoutchk(0, 1)

if ~isa(E, 'octonion') || ~isa(x(E.a), 'sym')
    error('Argument must be an octonion with symbolic components.')
end

% Two cases are covered:
%
% 1. The vector part is all zero. Return the scalar part as a scalar,
%    whether zero or non-zero (we don't need to check).
% 2. The vector part is non-zero but the scalar part is zero. Return the
%    vector part as a pure octonion.

if isAlways(all(x(E.a(:)) == 0) && ... % TODO Why don't we just compare the
            all(y(E.a(:)) == 0) && ... % vector part to zero using vector(E)?
            all(z(E.a(:)) == 0) && ...
            all(E.b == 0))
    S = scalar(E); % The scalar function is used because it returns zero in
                   % the case of a empty scalar part.
elseif ~isempty(s(E.a)) && isAlways(all(s(E.a(:)) == 0))
    % The scalar part is present, but is entirely zero, so we can elide it.
    S = vector(E); % Discard the scalar part and return a pure result.
else
    S = E; % Default case, there is no elision to do.
end

end

% $Id: elide.m 1126 2021-08-17 21:16:34Z sangwine $

% Created automatically from the quaternion
% function of the same name on 13-Aug-2021.