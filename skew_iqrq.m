function [Q, D] = skew_iqrq(H, rtol, verbose, progress_bar)
% SKEW_IQRQ Implicit Quaternion QR Algorithm for Skew-Symmetric Hessenberg Matrices
%
% [Q, D] = SKEW_IQRQ(H) computes the diagonal form D of the n-by-n quaternion 
% skew-symmetric Hessenberg matrix H, via implicit QR shifts. 
% Q is a unitary quaternion matrix such that:
%
%   Q' * H * Q = D
%
% D is a diagonal matrix.
%
% [Q, D] = SKEW_IQRQ(H, RTOL) uses tolerance RTOL instead of the default eps.
%
% [Q, D] = SKEW_IQRQ(H, RTOL, VERBOSE) prints additional 
% information about the computation if VERBOSE is true.
%
% [Q, D] = SKEW_IQRQ(H, RTOL, VERBOSE, PROGRESS_BAR) displays
% a progress bar during computation if PROGRESS_BAR is true.
%
% Example Usage:
% [Q, D] = skew_iqrq(H); % Default usage
% [Q, D] = skew_iqrq(H, 1e-12, true, true); % With all options specified
%
% References:
% [1] Bunse-Gerstner, et al., A Quaternion QR Algorithm, Numer. Math. 55,
%     83–95 (1989)
% [2] Golub, G. H. and Van Loan, C. F., Matrix Computations, Algorithm 7.5.2,  
%     Johns Hopkins University Press, 1996.

% Input validation
narginchk(1, 4);
nargoutchk(0, 2);

if ~ismatrix(H) 
    error('H must be a square matrix');
end

if ~isa(H, 'quaternion')
    error('H must be a quaternion matrix'); 
end

if ~istridiagonal(H)
    error('H must be a skew-symmetric upper Hessenberg matrix');
end

% Set defaults
if nargin < 4, progress_bar = false; end  
if nargin < 3, verbose = false; end
if nargin < 2, rtol = eps; end

% Helper functions
householder_lapply = @(u, x) x - u * (u' * x);
householder_rapply = @(u, x) x - (x * u) * u';

% Initialize timers and variables
H_mult_time = 0;
Q_time = 0;
House_time = 0;

n = size(H, 1);
ilo = 1;
ihi = n;
Q = eyeq(n);
r = eye(n, 1);
GS = 0;   % steps
Tz = tic;

HH = H;

if progress_bar
    b = waitbar(0, strcat('skew_iqrq, n = ', num2str(n)));
end

while ihi > 1
    while ihi - ilo > 2
        % Compute implicit QR shift using schur2
        x = shift(H(ihi-1:ihi, ihi-1:ihi));
        GS = GS + 1;
        si = -2 * x.w;
        ti = abs(x)^2;
        v = [
            H(ilo, ilo)^2 + H(ilo, ilo+1) * H(ilo+1, ilo) + si * H(ilo, ilo) + ti;
            H(ilo+1, ilo) * H(ilo, ilo) + H(ilo+1, ilo+1) * H(ilo+1, ilo) + si * H(ilo+1, ilo);
            H(ilo+2, ilo+1) * H(ilo+1, ilo)
        ];

        T0 = tic;
        [u, ~] = householder_vector(v, r(1:3));
        House_time = House_time + toc(T0);

        T0 = tic;
        H(ilo:ilo+3, ilo:ilo+2) = householder_rapply(u, H(ilo:ilo+3, ilo:ilo+2));
        H(ilo:ilo+2, ilo:ilo+3) = householder_lapply(u, H(ilo:ilo+2, ilo:ilo+3));
        H_mult_time = H_mult_time + toc(T0);

        T0 = tic;
        Q(:, ilo:ilo+2) = householder_rapply(u, Q(:, ilo:ilo+2));
        Q_time = Q_time + toc(T0);

        for i = ilo : (ihi - 2)
            e  = min(i + 3, ihi);
            ee = min(e + 1, ihi);

            T0 = tic;
            [u, ~] = householder_vector(H(i+1:e, i), r(1:e-i));
            House_time = House_time + toc(T0);

            T0 = tic;
            H(i+1:e, i:ee) = householder_lapply(u, H(i+1:e, i:ee));
            H(i:ee, i+1:e) = householder_rapply(u, H(i:ee, i+1:e));
            H_mult_time = H_mult_time + toc(T0);

            T0 = tic;
            Q(:, i+1:e) = householder_rapply(u, Q(:, i+1:e));
            Q_time = Q_time + toc(T0);

            % Check for deflation using schur1
            if abs(H(i+1, i)) <= rtol * (abs(H(i, i)) + abs(H(i+1, i+1)))
                [U, ~] = schur1(H(i, i));
                Q(:, i) = Q(:, i) * U;
                H(:, i) = H(:, i) * U;
                H(i, :) = U' * H(i, :);
                H(i+1, i) = 0;
                ilo = i + 1;
            end
        end

        % Check for deflation
        while ihi > ilo + 1 && (abs(H(ihi, ihi-1)) <= rtol * (abs(H(ihi, ihi)) + abs(H(ihi-1, ihi-1))))
            [U, ~] = schur1(H(ihi, ihi));
            Q(:, ihi) = Q(:, ihi) * U;
            H(:, ihi) = H(:, ihi) * U;
            H(ihi, :) = U' * H(ihi, :);
            H(ihi, ihi-1) = 0;
            ihi = ihi - 1;
            if progress_bar
                waitbar(1 - ihi / n, b);
            end
        end

        % Ensure matrix remains skew-symmetric and Hessenberg
        H = triu(H, -1);
        H = tril(H, 1);
    end
    [U, H( ilo : ihi , ilo : ihi )] = iqrq(H( ilo : ihi , ilo : ihi ), rtol, false, false);
    Q(  :  , ilo : ihi ) = Q( : , ilo : ihi ) * U;

    ihi = ilo - 1;
    ilo = 1;
end

% Extract diagonal elements to form D
D = diag(diag(H));

if progress_bar
    delete(b);
end

if verbose
    fprintf('Total steps of Skew IQR: %d\n', GS);
    fprintf('Time for H matrix multiplications: %f s\n', H_mult_time);
    fprintf('Time for computing Q matrix: %f s\n', Q_time);
    fprintf('Time for Householder computations: %f s\n', House_time);
    fprintf('Total time: %f s, total deflation: %d\n', toc(Tz), n);
    fprintf('Time without constructing Q: %f s\n', toc(Tz) - Q_time);
end

end

% Helper functions
function y = istridiagonal(A)
% Returns true if A is tridiagonal
y = normq(tril(A, -2)) / normq(A) < eps & normq(triu(A, 2)) / normq(A) < eps;
end

function [Q, T] = schur1(A)
% Standardize a quaternion
% A: 1-by-1 quaternion matrix
[Q, ~] = schur(adjoint(A));
Q = quaternion(real(Q(1)), imag(Q(1)), -real(Q(2)), imag(Q(2)));
T = Q' * A * Q;
end

function [Q, T] = schur2(A)
% Compute an ordered Schur decomposition of A.
% A: 2-by-2 quaternion matrix
[Q, T] = schur(adjoint(A));
[Q, ~] = ordschur(Q, T, groupByValue(abs(diag(T))));
Q = quaternion( ...
    real(Q([1,2], [1,3])), ...
    imag(Q([1,2], [1,3])), ...
    -real(Q([3,4], [1,3])), ...
    imag(Q([3,4], [1,3])) ...
);
T = Q' * A * Q;
end

function group = groupByValue(A)
[~, idx] = sort(A); group = zeros(size(A));
% Assign group numbers
current = 1;
for i = 2:length(A)
    if abs(A(idx(i)) - A(idx(i-1))) < eps
        group(idx(i)) = current;
    else
        current = current + 1;
        group(idx(i)) = current;
    end
end
% Assign group number to the first element
group(idx(1)) = 1;
end

function s = shift(A)
% Selects shifts for the QR algorithm
% A: 2-by-2 matrix
[~, T1] = schur1(A(2, 2));
[~, T2] = schur2(A);
if abs(T2(1, 1) - T1) < abs(T2(2, 2) - T1)
    s = T2(1, 1);
else
    s = T2(2, 2);
end
end

function err = rel_err(A, B)
err = norm(A - B, 'fro') / max(norm(A, 'fro'), norm(B, 'fro'));
end
    