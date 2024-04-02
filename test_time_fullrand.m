addpath('/home/yjshao/MATLAB/qtfm')
tol = eps;
rtol = 1e-12;

order_n = [64, 128, 256, 512, 1024, 2048, 4096, 8192];
for n = order_n
    rng(0)

    fprintf('Testing random matrix of order_n = %d ...\n', n);
    A = randq(n) .* randn(n);
    [~, A] = hessq(A);

    disp('Using aed ...');
    tic
    [Q, T] = aedq(A, tol, true, true);
    toc
    isunitary(Q, rtol); istriu_(T, 0, rtol); 
    compare_(Q' * A * Q, T, rtol, "Q' * A * Q and T does not match!");
    save(sprintf('full_aed_n%d.mat', n), 'A', 'Q', 'T');

    disp('Using iqr ...');
    tic
    [Q, T] = iqrq(A, tol, true, true);
    toc
    isunitary(Q, rtol); istriu_(T, 0, rtol);
    compare_(Q' * A * Q, T, rtol, "Q' * A * Q and T does not match!");
    save(sprintf('full_iqr_n%d.mat', n), 'A', 'Q', 'T');

end

% Helper functions
function y = istridiagonal(A)
% Returns true if A is tridiagonal
y = all(abs(tril(A, -2)) < eps, "all") & all(abs(triu(A, 2)) < eps, "all");
end


function y = istriu_(X, K, E)
% ISTRIU Test function to check that the quaternion matrix A is upper
% triangular.
%
% Y = ISTRIU(X, K, E) returns true if X is upper triangular to within a
% tolerance E, and false otherwise. X is upper triangular if X(i, j) = 0
% for all i > j + K. If K is not specified, it defaults to 0. If E is not
% specified, it defaults to eps.

narginchk(1, 3), nargoutchk(0, 1)

y = true;

if nargin < 2
    K = 0;
end

if nargin < 3
    E = eps;
end

dist = max(max(abs(tril(X, K - 1)))) / sqrt(size(X, 1));

if dist > E
    fprintf('Failure! A is not an upper triangular matrix! abs_err = %.4d\n', dist);
    spy(abs(X) > E);
    y = false;
end
end


function y = isdiag_(X, E)
% ISDIAG Test function to check that the quaternion matrix A is diagonal.
%
% Y = ISDIAG(X, E) returns true if X is diagonal to within a tolerance E,
% and false otherwise. X is diagonal if X(i, j) = 0 for all i ~= j. If E is
% not specified, it defaults to eps.

narginchk(1, 2), nargoutchk(0, 1)

y = true;

if nargin < 2
    E = eps;
end

dist = max(max(abs(X - diag(diag(X))))) / sqrt(size(X, 1));

if dist > E
    fprintf('Failure! A is not a diagonal matrix! abs_err = %.4d\n', dist);
    spy(abs(X) > E);
    y = false;
end
end


function compare_(A, B, T, E)
% COMPARE Test function to check that two quaternion arrays (real or
% complex) are equal, to within a tolerance.
%
% COMPARE(A, B, T, E) returns true if A and B are equal to within a
% tolerance T, and false otherwise. If A and B are not equal, an error
% message E is output. If T is not specified, it defaults to eps. If E is
% not specified, it defaults to ''. This will also work for non-quaternion
% arrays.

narginchk(2, 4), nargoutchk(0, 0)

if nargin < 3
    T = eps;
end

if nargin < 4
    E = '';
end

dist = norm(abs(A - B), 'fro') / norm(abs(A), 'fro');
if dist > T
    disp(dist)
    error(E);
end
end