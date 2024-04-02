function test_qr
% Test code for the quaternion qr function.

% Copyright © 2011 Stephen J. Sangwine and Nicolas Le Bihan.
% See the file : Copyright.m for further details.

disp('Testing QR decomposition ...')

% TODO Add check that R is upper triangular, to within a tolerance, by
% using the triu function to extract the upper triangular part, leaving
% only rounding errors, which compare can check for us.

T = 1e-12;

A = randq(10);

[Q, R] = qr(A);
compare(Q * R, A, T,   'quaternion/qr failed test 1')
check(isunitary(Q, T), 'quaternion/qr test 1 has non-unitary Q result')

A = randq(10,13);

[Q, R] = qr(A);
compare(Q * R, A, T,   'quaternion/qr failed test 2')
check(isunitary(Q, T), 'quaternion/qr test 2 has non-unitary Q result')

[Q, R] = qr(A.');
compare(Q * R, A.', T, 'quaternion/qr failed test 3')
check(isunitary(Q, T), 'quaternion/qr test 3 has non-unitary Q result')

% Now repeat the whole lot with complex data, but smaller matrices, keeping
% the same tolerance. No particular reason, but no need to keep the same
% sizes.

% Commented out complex tests while we sort out the problems with
% biquaternion unitary matrices.
disp('Passed')
return

% T = 1e-9; % Relax the requirements a little.
% 
% A = complex(randq(5), randq(5));
% 
% [Q, R] = qr(A);
% compare(Q * R, A, T,   'quaternion/qr failed test 4')
% 
% A = complex(randq(5,7), randq(5,7));
% 
% [Q, R] = qr(A);
% compare(Q * R, A, T,   'quaternion/qr failed test 5')
% 
% [Q, R] = qr(A.');
% compare(Q * R, A.', T, 'quaternion/qr failed test 6')
% 
% disp('Passed');

% $Id: test_qr.m 1131 2021-12-21 21:52:48Z sangwine $
