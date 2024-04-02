function result = solve(eqn, var)
% SOLVE Solve symbolic octonion equation(s).
% (Octonion overloading of standard Matlab function.)

% Copyright © 2022 Stephen J. Sangwine and Nicolas Le Bihan.
% See the file : Copyright.m for further details.

% TODO A major issue here is parentheses. An octonion expression must be
% parenthesised, but when evaluated into eight expressions involving
% complex numbers, the parentheses are no longer needed, and must be
% removed, otherwise the MATLAB SOLVE may not be able to make any progress.
% This can be done with EXPAND. Is it possible to detect that an expression
% contains parentheses, and warn the user that EXPAND may be needed?

% TODO Support the parameter profiles with Name/Value pairs, assuming it
% makes sense and that we can pass them to the individual solve calls on
% the components of the quaternion equation.

narginchk(1, 2), nargoutchk(0, 1)

assert(isa(eqn, 'sym'), ['First parameter must have class sym, given ', class(eqn), '.'])

assert(isscalar(eqn), 'Can handle only one equation at present.') % TODO

% The second parameter must be an octonion, since if it isn't, and the
% first would pass the test above, this function would not be called.

assert(isa(var, 'octonion'))

% We need to test here that eqn is a symbolic expression representing an
% octonion equation. This is tricky because it is actually an expression
% involving the components of an octonion, and contains no octonions. So
% we want to be sure it has 7 or 8 subexpressions combined by &.
% Unfortunately we don't seem to be able to access the operands in the
% expression (to check for &), so we resort to rebuilding the expression
% from its subexpressions and testing for equality. Clunky workaround!

% TODO Is there any case where the eqn could have 7 components? (I.e. a
% pure octonion expression.)

subexpr = children(eqn);

if length(subexpr) ~= 8 || ...
        ~isequaln(subexpr{1} & subexpr{2} & subexpr{3} & subexpr{4} & ...
                  subexpr{5} & subexpr{6} & subexpr{7} & subexpr{8}, eqn)
    error('Input eqn does not have correct form for an octonion equation')
end

r = solve([subexpr{:}], ...
    [part(var, 1), part(var, 2), part(var, 3), part(var, 4), ...
     part(var, 5), part(var, 6), part(var, 7), part(var, 8)]);

% r will be a struct, with fields named according to the variables which
% make up the octonion variable var. We need to convert this struct into
% an octonion. The trick for doing this is taken from Matlab Central
% answers 229604 by 'Guillaume' 13 Jul 2015.

rnames = fieldnames(r); % This will be a cell array of names.

f1 = rnames{1};
f2 = rnames{2};
f3 = rnames{3};
f4 = rnames{4};
f5 = rnames{5};
f6 = rnames{6};
f7 = rnames{7};
f8 = rnames{8};

result = octonion(r.(f1), r.(f2), r.(f3), r.(f4), ...
                  r.(f5), r.(f6), r.(f7), r.(f8));

end

% $Id: solve.m 1156 2022-10-13 16:12:30Z sangwine $
