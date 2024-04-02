function test_octonion_matrix_inverses
% Test code for the octonion matrix inverse functions.

% Copyright © 2021 Stephen J. Sangwine and Nicolas Le Bihan.
% See the file : Copyright.m for further details.

disp('Testing octonion matrix inverses ...');

T = 1e-2; % This is a loose tolerance, the algorithms are known to be inaccurate.

% TODO Provide a more reliable tolerance or a different way to check the
% result that gives fewer false negative results. Ideally scale the
% tolerance according to the matrix size.

for k = 2:5
    A = rando(k); B = complex(A, rando(k)) .* randn(k);
    
    I = octonion(eye(k));

    compare(linv(A) * A, I, T, ['linv failed test k = ' num2str(k)])   
    compare(A * rinv(A), I, T, ['rinv failed test k = ' num2str(k)])

end

disp('Passed');

end

% $Id: test_octonion_matrix_inverses.m 1113 2021-02-01 18:41:09Z sangwine $
