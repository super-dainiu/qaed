% Symbolic computation demonstration.
%
% This script demonstrates the symbolic capability of QTFM version 3. The
% problem covered is the solution of Sylvester's equation. Note that an
% explicit solution to this equation is known, due to Janovska and Opfer:
%
% Drahoslava Janovska and Gerhard Opfer, 'Linear equations in quaternionic
% variables', Mitteilungen der Mathematischen Gesellschaft in Hamburg,
% 27:223–234, 2008. Available at: https://www.math.uni-hamburg.de/mathges/
%
% so we can compare the solution produced by SOLVE with the known solution.
% QTFM includes a SYLVESTER function, which implements the Janovska and
% Opfer solution.

% Copyright  : Steve Sangwine, October 2022

if isempty(ver('symbolic'))
    disp('Symbolic Toolbox is not present, unable to complete script.')
    return
end

syms aw ax ay az; a = quaternion(aw, ax, ay, az);
syms bw bx by bz; b = quaternion(bw, bx, by, bz);
syms cw cx cy cz; c = quaternion(cw, cx, cy, cz);
syms xw xx xy xz; x = quaternion(xw, xx, xy, xz);
syms yw yx yy yz; y = quaternion(yw, yx, yy, yz);

eqn = a * x + x * b == c; % Sylvester's equation.

% Janovska and Opfer's pair of solutions, (2.2) and (2.3) in their paper:

ia = a^-1;
ib = b^-1;

y1 = (2 * scalar(b) + a + normq(b) * ia)^-1 * (c + ia * c * conj(b));
y2 = (c + conj(a) * c * ib) * (2 * scalar(a) + b + normq(a) * ib)^-1;

% Show that these two solutions are equivalent (the reason for two
% solutions is numerical, to handle cases where the norm of either of a or
% b is small). We can do this by subtraction and simplification to a
% quaternion zero.

disp(['Demonstrate that the two solutions given by Janovska and Opfer', ...
      ' are equivalent (subtract one from the other and simplify):'])

simplify(y1 - y2)

disp(' ')

%% Now check Janovska and Opfer's solutions by substitution:

z1 = subs(eqn, x, y1);
z2 = subs(eqn, x, y2);

disp('Result of substitution, Janovska/Opfer solutions:')
[simplify(z1), simplify(z2)] %#ok<NOPTS> 

%% Find the solution using QTFM's SOLVE (which uses MATLAB's SOLVE):

r  = solve(eqn, x);
z3 =  subs(eqn, x, r);

disp('Result of substitution, SOLVE solution:')

simplify(expand(z3)) % Without EXPAND, this does not simplify to SYMTRUE

% $Id: symbolic_demonstration_2.m 1154 2022-10-13 09:27:42Z sangwine $