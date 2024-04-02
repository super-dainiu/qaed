% Symbolic computation demonstration.
%
% This script demonstrates the symbolic capability of QTFM version 3.
%
% The problem studied is a simple one: to show that quaternions with
% complex components (biquaternions) do not commute with their complex
% conjugates, a surprising fact, since they do commute with their
% quaternion or Hamilton conjugates.

% Copyright  : Steve Sangwine, June 2020

if isempty(ver('symbolic'))
    disp('Symbolic Toolbox is not present, unable to complete script.')
    return
end

disp('Firstly, demonstrate numerically that a biquaternion does')
disp('not commute with its complex conjugate, by a counterexample.')
disp(' ')
 
q = complex(qi + qj + qk, qi - qj + qk);

disp(q .* conj(q, 'complex'))
disp(conj(q, 'complex') .* q)

disp(' ')
disp('... but it does commute with its quaternion conjugate:')
disp(' ')

disp(q .* conj(q))
disp(conj(q) .* q)

disp(' ')
disp('Now demonstrate the same symbolically (see code):')

syms w x y z % Declare four variables, assumed complex by default.

q = quaternion(w, x, y, z); % Make a symbolic biquaternion.

disp('Complex conjugate case - the warning indicates that this does not always hold true')
disp(' ')

isAlways(q .* conj(q, 'complex') == conj(q, 'complex') .* q)

disp(' ')
disp('Quaternion conjugate case')
disp(' ')

isAlways(q .* conj(q) == conj(q) .* q)

disp('')
disp('Now demonstrate the complex conjugate using assumptions.')
disp('')

assume(q, 'real') % Notice this is overloaded so we don't have to assume
                  % for each of w x y z.
                 
isAlways(q .* conj(q, 'complex') == conj(q, 'complex') .* q)

disp('Another method is to compute the two versions and display them.')
disp('Notice that there are sign differences between the two versions.')

assume(q, 'clear') % This removes the assumptions on q, so its components
                   % are now assumed complex again.

simplify(q .* conj(q, 'complex'))
simplify(conj(q, 'complex') .* q)

% $Id: symbolic_demonstration_1.m 1154 2022-10-13 09:27:42Z sangwine $