function S = elide(E)
% ELIDE  Remove redundant parts of symbolic quaternion.

% Copyright © 2021 Stephen J. Sangwine and Nicolas Le Bihan.
% See the file : Copyright.m for further details.

narginchk(1, 1), nargoutchk(0, 1)

if ~isa(E, 'quaternion') || ~isa(E.x, 'sym')
    error('Argument must be a quaternion with symbolic components.')
end

% Two cases are covered:
%
% 1. The vector part is all zero. Return the scalar part as a scalar,
%    whether zero or non-zero (we don't need to check).
% 2. The vector part is non-zero but the scalar part is zero. Return the
%    vector part as a pure quaternion.

if isAlways(all(E.x(:) == 0)) && ... % TODO This is probably the most time-
   isAlways(all(E.y(:) == 0)) && ... % consuming part of the function, we
   isAlways(all(E.z(:) == 0))        % apply isAlways separately to each
                                     % part in order to skip the second and
                                     % third tests if the first fails.
    S = scalar(E); % The scalar function is used because it returns zero in
                   % the case of an empty scalar part.
elseif ~isempty(E.w) && isAlways(all(E.w(:) == 0))
    % The scalar part is present, but is entirely zero, so we can elide it.
    S = vector(E); % Discard the scalar part and return a pure result.
else
    S = E; % Default case, there is no elision to do.
end

end

% $Id: elide.m 1127 2021-08-30 19:35:41Z sangwine $
