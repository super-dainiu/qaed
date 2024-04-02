function test_inv
% Test code for the quaternion inv function.

% Copyright © 2011 Stephen J. Sangwine and Nicolas Le Bihan.
% See the file : Copyright.m for further details.

disp('Testing inverse ...')

T = 1e-10;

% Test scalar inverse on a few real and complex quaternions.

for i=1:10
    q = randq;
    compare(q * inv(q), onesq, T, 'quaternion/inv failed test 1') %#ok<*MINV>
    b = complex(q, randq) .* randn;
    compare(b * inv(b), onesq, T, 'quaternion/inv failed test 2')
end

A = randq(10);
compare(A * inv(A), eyeq(10), T, 'quaternion/inv failed test 3')

B = complex(randq(10), randq(10)) .* randn(10);
compare(B * inv(B), eyeq(10), T, 'quaternion/inv failed test 3')

% Test a case where the matrix has sub-blocks which are zero. This matrix
% has an inverse, but the code prior to April 2022 could not compute it
% because of NaNs (the new code uses adjoint/unadjoint to work around
% this). The result is exact because all the quaternions in A have values
% of zero or one.

A = [ 0, qi, qj, qk; ...
      0,  0, qi, qj; ...
     qk,  0,  0, qi; ...
      0,  0,  0, qk];

B = inv(A);

check(A * B == eyeq(4), 'quaternion/inv fails test 4.')

disp('Passed');

% $Id: test_inv.m 1145 2022-09-07 14:44:23Z sangwine $
