function S = subs(expr, old, new)
% SUBS  Symbolic substitution in quaternion expressions
% (Quaternion overloading of standard Matlab function.)

% Copyright © 2022 Stephen J. Sangwine and Nicolas Le Bihan.
% See the file : Copyright.m for further details.

narginchk(3, 3), nargoutchk(0, 1)

assert(isa(expr, 'sym'), ...
    ['First parameter must have class sym, given ', class(expr), '.'])

assert(isscalar(expr), 'Can handle only one equation at present.') % TODO

assert(isa(old, 'quaternion'), 'Second parameter must be a quaternion.')
assert(isa(new, 'quaternion'),  'Third parameter must be a quaternion.')

% The code below for splitting the expression was copied from SOLVE.
% See the comments there for discussion of its limitations.
% TODO Consider factoring this out into a private function.

% Check that the expression we are given is consistent with being a
% quaternion expression. We can't directly obtain the operators that join
% the children of the expression (why not, Mathworks?).

subexpr = children(expr);

if length(subexpr) ~= 4 || ~isequaln(subexpr{1} & subexpr{2} & subexpr{3} & subexpr{4}, expr)
    error('Input expression does not have correct form for a quaternion expression or equation')
end

S = subs(expr, {scalar(old), x(old), y(old), z(old)}, ...
               {scalar(new), x(new), y(new), z(new)});
end

% $Id: subs.m 1155 2022-10-13 09:55:58Z sangwine $
