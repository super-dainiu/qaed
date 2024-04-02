function test_symbolic
% Test code for the quaternion and octonion symbolic functions

% Copyright © 2020, 2021 Stephen J. Sangwine and Nicolas Le Bihan.
% See the file : Copyright.m for further details.

disp('Testing quaternion symbolic functionality ...')

% Start with the constructor to make sure that all the ways of building a
% quaternion with symbolic components will work correctly.

q = quaternion('a'); % Make a quaternion with symbolic zero vector part.
p = quaternion(a);   % This should be identical to q, since it has the same
                     % scalar part.
if p ~= q
    error('Symbolic test 1 fails.')
end

if vector(p) ~= 0
    error('Symbolic test 2 fails.')
end

syms a b c d

q = quaternion(a, b, c, d);
p = quaternion(a, vector(q));

if p ~= q
    error('Symbolic test 3 fails.')
end

q = quaternion('t', vector(q));
p = quaternion(t, b, c, d);

if p ~= q
    error('Symbolic test 4 fails.')
end

clear t

p = quaternion('x', 'y', 'z');
q = quaternion(x, y, z);

if p - q ~= 0 % Do a subtraction here rather than comparing p with q so we
              % test a bit more code.
    error('Symbolic test 5 fails.')
end

clear x y z

p = quaternion('w', 'x', 'y', 'z');
q = quaternion(w, x, y, z);

if p - q ~= 0 % Do a subtraction here rather than comparing p with q so we
              % test a bit more code.
    error('Symbolic test 6 fails.')
end

% A couple of further tests on the constructor.

if ~eval(cast(qi, 'sym') == quaternion(sym(1), sym(0), sym(0)))
   error('Symbolic test 7 fails.') 
end

if ~eval(cast(unit(quaternion(1,1,1)), 'sym') == ...
         unit(quaternion(sym(1), sym(1), sym(1))))
     error('Symbolic test 8 fails.')
end

% Test the inverse, and also simplification down to unity.

p = inv(q);
r =  p .* q;

if eval(simplify(r) ~= 1)
    error('Symbolic test of inverse fails on single quaternion.')
end

% Check that we can make sizes greater than [1,1], and invert a small
% matrix.

clear w x y z

Q = quaternion(sym('w', [2,2]), sym('x', [2,2]), ...
               sym('y', [2,2]), sym('z', [2,2]));
           
check(isAlways(simplify(Q * inv(Q)) == sym(eye(2))), ...
    'Symbolic 2-by-2 inverse fails'); %#ok<MINV>

% Now check out theorem proving. This is a much more demanding test of the
% toolbox code, but tests are a bit harder to devise.

% Perpendicular vectors negate on re-ordering. That is uv = -vu provided u
% and v are orthogonal. Note that we use the elementwise multiply on the
% LHS and the matrix multiply on the RHS, in order to verify that these two
% products both give the same result.

syms w x y z
v = quaternion(b, c, d);
u = quaternion(x, y, z);
assume(scalar_product(u, v) == 0); % Declares u perpendicular to v
if ~isAlways(u .* v == -v * u)
    error('Symbolic test fails: perpendicular vectors do not negate on re-ordering.')
end

% Test that a unit pure quaternion squares to minus 1. We can do this test
% in more than one way.

if simplify(unit(u).^2) ~= -1
   error('Symbolic test fails: unit pure quaternion squared fails to simplify to -1.') 
end

assume(normq(v) == 1);
if ~isAlways(v^2 == -1) % NB We use the matrix square here, but elementwise
                        % square in the previous test, deliberately.
   error('Symbolic test fails: unit pure quaternion fails to square to -1.') 
end

% TODO Insert some more demanding quaternion tests here.
% ******************************************************

disp('Passed')

clear variables

disp('Testing octonion symbolic functionality ...')

% Start with the constructor to make sure that all the ways of building an
% octonion with symbolic components will work correctly.

q = octonion('a'); % Make an octonion with symbolic zero vector part.
p = octonion(a);   % This should be identical to q, since it has the same
                   % scalar part.
if p ~= q
    error('Symbolic test 1 fails.')
end

if vector(p) ~= 0
    error('Symbolic test 2 fails.')
end

syms a b c d w x y z

q = octonion(w, x, y, z, a, b, c, d);
p = octonion(w, vector(q));

if p ~= q
    error('Symbolic test 3 fails.')
end

q = octonion('t', vector(q));
p = octonion(t, x, y, z, a, b, c, d);

if p ~= q
    error('Symbolic test 4 fails.')
end

clear t x y z a b c d

p = octonion('x', 'y', 'z', 'a', 'b', 'c', 'd');
q = octonion(x, y, z, a, b, c, d);

if p - q ~= 0 % Do a subtraction here rather than comparing p with q so we
              % test a bit more code.
    error('Symbolic test 5 fails.')
end

clear w x y z a b c d

p = octonion('w', 'x', 'y', 'z', 'a', 'b', 'c', 'd');
q = octonion(w, x, y, z, a, b, c, d);

if p - q ~= 0 % Do a subtraction here rather than comparing p with q so we
              % test a bit more code.
    error('Symbolic test 6 fails.')
end

% A couple of further tests on the constructor.

if ~eval(cast(oi, 'sym') == ...
         octonion(sym(1), sym(0), sym(0), sym(0), sym(0), sym(0), sym(0)))
   error('Symbolic test 7 fails.') 
end

if ~eval(cast(unit(octonion(1,1,1,1,1,1,1)), 'sym') == ...
         unit(octonion(sym(1), sym(1), sym(1), sym(1), sym(1), sym(1), sym(1))))
     error('Symbolic test 8 fails.')
end

% Check that we can make an octonion from two quaternions.

q = quaternion(w, x, y, z);
p = quaternion(a, b, c, d);
check(isAlways(octonion(q, p) == octonion(w, x, y, z, a, b, c, d)), ...
    'Symbolic test 9 fails');

% Test the inverse, and also simplification down to unity.

p = inv(q);
r =  p .* q;

if eval(simplify(r) ~= 1)
    error('Symbolic test of inverse fails on single octonion.')
end

% Check that we can make sizes greater than [1,1], and invert a small
% matrix (we have to use linv or rinv because octonion matrices have two
% inverses).

clear w x y z a b c d

Q = octonion(sym('w', [2,2]), sym('x', [2,2]), ...
             sym('y', [2,2]), sym('z', [2,2]), ...
             sym('a', [2,2]), sym('b', [2,2]), ...
             sym('c', [2,2]), sym('d', [2,2]));
% We put in the test below in order to do something with Q, since the rinv
% won't currently work.
P = Q.';
check(isAlways(Q == P.'), 'Symbolic test of transpose fails.');
         
% TODO The test below with rinv doesn't currently work because the poly
% function has been removed from Matlab for symbolics. Substituting the
% recommended charpoly function results in the numerical code running very
% much slower. So this needs more detailed study, and can wait for now (Aug
% 2021). Clearly from Tian's method in rinv/linv, even a symbolic 2-by-2
% inverse is going to be a big challenge and it may not be worth trying it.

% *********************************************************
% check(isAlways(simplify(Q * rinv(Q)) == sym(eye(2))), ...
%     'Symbolic 2-by-2 inverse fails');
% *********************************************************

% Now check out theorem proving. This is a much more demanding test of the
% toolbox code, but tests are a bit harder to devise.

% Perpendicular vectors negate on re-ordering. That is uv = -vu provided u
% and v are orthogonal. Note that we use the elementwise multiply on the
% LHS and the matrix multiply on the RHS, in order to verify that these two
% products both give the same result.

syms x1 y1 z1 a1 b1 c1 d1
syms x2 y2 z2 a2 b2 c2 d2

v = octonion(x1, y1, z1, a1, b1, c1, d1);
u = octonion(x2, y2, z2, a2, b2, c2, d2);
assume(scalar_product(u, v) == 0); % Declares u perpendicular to v
if ~isAlways(u .* v == -v * u)
    error('Symbolic test fails: perpendicular vectors do not negate on re-ordering.')
end

% Test that a unit pure octonion squares to minus 1. We can do this test
% in more than one way.

if simplify(unit(u).^2) ~= -1
   error('Symbolic test fails: unit pure octonion squared fails to simplify to -1.') 
end

assume(normo(v) == 1);
if ~isAlways(v^2 == -1) % NB We use the matrix square here, but elementwise
                        % square in the previous test, deliberately.
   error('Symbolic test fails: unit pure octonion fails to square to -1.') 
end

% Now try a Moufang identity (also tested in the numeric case) in several
% different ways in order to exercise various functions.

clear variables

syms x_w x_x x_y x_z x_a x_b x_c x_d
syms y_w y_x y_y y_z y_a y_b y_c y_d
syms z_w z_x z_y z_z z_a z_b z_c z_d

x = octonion(x_w, x_x, x_y, x_z, x_a, x_b, x_c, x_d);
y = octonion(y_w, y_x, y_y, y_z, y_a, y_b, y_c, y_d);
z = octonion(z_w, z_x, z_y, z_z, z_a, z_b, z_c, z_d);

check(isAlways((x .* z .* x) .* y == x .* (z .* (x .* y))), ...
      'Moufang identity fails isAlways test for equality');
  
check(isAlways((x .* z .* x) .* y - x .* (z .* (x .* y)) == 0), ...
      'Moufang identity fails isAlways test for zero difference');

check(eval(expand((x .* z .* x) .* y) - ...
           expand(x .* (z .* (x .* y))) == 0), ...
      'Moufang identity fails expand test');
  
disp('Passed')

end

% $Id: test_symbolic.m 1129 2021-09-02 21:03:39Z sangwine $
