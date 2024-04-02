function result = solve(eqn, var)
% SOLVE Solve symbolic quaternion equation(s).
% (Quaternion overloading of standard Matlab function.)

% Copyright © 2022 Stephen J. Sangwine and Nicolas Le Bihan.
% See the file : Copyright.m for further details.

% TODO Support the parameter profiles with Name/Value pairs, assuming it
% makes sense and that we can pass them to the individual solve calls on
% the components of the quaternion equation.

narginchk(1, 2), nargoutchk(0, 1)

assert(isa(eqn, 'sym'), ['First parameter must have class sym, given ', class(eqn), '.'])

assert(isscalar(eqn), 'Can handle only one equation at present.') % TODO

% The second parameter must be a quaternion, since if it isn't, and the
% first would pass the test above, this function would not be called.

assert(isa(var, 'quaternion'))

% We need to test here that eqn is a symbolic expression representing a
% quaternion equation. This is tricky because it is actually an expression
% involving the components of a quaternion, and contains no quaternions. So
% we want to be sure it has 3 or 4 subexpressions combined by &.
% Unfortunately we don't seem to be able to access the operands in the
% expression (to check for &), so we resort to rebuilding the expression
% from its subexpressions and testing for equality. Clunky workaround!

% TODO Is there any case where the eqn could have 3 components? (I.e. a
% pure quaternion expression.)

subexpr = children(eqn);

if length(subexpr) ~= 4 || ~isequaln(subexpr{1} & subexpr{2} & subexpr{3} & subexpr{4}, eqn)
    error('Input eqn does not have correct form for a quaternion equation')
end

r = solve([subexpr{:}], [scalar(var), x(var), y(var), z(var)]);

% r will be a struct, with fields named according to the variables which
% make up the quaternion variable var. We need to convert this struct into
% a quaternion. The trick for doing this is taken from Matlab Central
% answers 229604 by 'Guillaume' 13 Jul 2015.

rnames = fieldnames(r); % This will be a cell array of names.

f1 = rnames{1};
f2 = rnames{2};
f3 = rnames{3};
f4 = rnames{4};

result = quaternion(r.(f1), r.(f2), r.(f3), r.(f4));

end

% $Id: solve.m 1152 2022-09-29 20:10:33Z sangwine $
