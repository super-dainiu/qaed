function check(L, E)
% Test function to check a logical condition.

narginchk(2, 2), nargoutchk(0, 0)

if ~islogical(L)
    error('First parameter must be logical.');
end

if ~all(L)
    error(E);
end

% $Id: check.m 1113 2021-02-01 18:41:09Z sangwine $

