function S = subs(expr, old, new)
% SUBS  Symbolic substitution in octonion expressions
% (Octonion overloading of standard Matlab function.)

% Copyright © 2022 Stephen J. Sangwine and Nicolas Le Bihan.
% See the file : Copyright.m for further details.

narginchk(3, 3), nargoutchk(0, 1)

assert(isa(expr, 'sym'), ...
    ['First parameter must have class sym, given ', class(expr), '.'])

assert(isscalar(expr), 'Can handle only one equation at present.') % TODO

assert(isa(old, 'octonion'), 'Second parameter must be an octonion.')
assert(isa(new, 'octonion'),  'Third parameter must be an octonion.')

% The code below for splitting the expression was copied from SOLVE.
% See the comments there for discussion of its limitations.
% TODO Consider factoring this out into a private function.

% Check that the expression we are given is consistent with being an
% octonion expression. We can't directly obtain the operators that join
% the children of the expression (why not, Mathworks?).

subexpr = children(expr);

if length(subexpr) ~= 8 || ...
        ~isequaln(subexpr{1} & subexpr{2} & subexpr{3} & subexpr{4} & ...
                  subexpr{5} & subexpr{6} & subexpr{7} & subexpr{8}, eqn)
    error('Input eqn does not have correct form for an octonion equation')
end

S = subs(expr, ...
    {part(old, 1), part(old, 2), part(old, 3), part(old, 4), ...
     part(old, 5), part(old, 6), part(old, 7), part(old, 8)}, ...
    {part(new, 1), part(new, 2), part(new, 3), part(new, 4), ...
     part(new, 5), part(new, 6), part(new, 7), part(new, 8)});
end

% $Id: subs.m 1155 2022-10-13 09:55:58Z sangwine $
