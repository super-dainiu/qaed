function [Q, T] = aedq(H, rtol, verbose, progress_bar)
% AEDQ Eigenvalue and Eigenvector Computation via Aggressive Early Deflation
% for Quaternion Hessenberg Matrices
%
% [Q, T] = AEDQ(H) applies the aggressive early deflation algorithm to the
% n-by-n quaternion Hessenberg matrix H for computing its right eigenvalues
% and corresponding eigenvectors. The algorithm enhances eigenvalue
% computation efficiency and accuracy.
% Q is a unitary quaternion matrix and T is an upper triangular quaternion
% matrix satisfying:
%
%   Q' * H * Q = T
%
% [Q, T] = AEDQ(H, RTOL) uses tolerance RTOL instead of the default eps.
%
% [Q, T] = AEDQ(H, RTOL, VERBOSE) prints additional 
% information about the computation if VERBOSE is true.
%
% [Q, T] = AEDQ(H, RTOL, VERBOSE, PROGRESS_BAR) displays
% a progress bar during computation if PROGRESS_BAR is true.
%
%
% Example Usage:
% [Q, D] = aedq(H); % Default usage
% [Q, D] = aedq(H, 1e-12, true, true); % With all options specified
%
%
% References:
% [1] Bunse-Gerstner, et al., A Quaternion QR Algorithm, Numer. Math. 55,
%     83–95 (1989) 
% [2] Golub, G. H. and Van Loan, C. F., Matrix Computations, Algorithm
%     7.5.2, Johns Hopkins University Press, 1996.
% [3] Braman, K., Byers, R. & Mathias, R. The Multishift QR Algorithm. Part II:
%     Aggressive Early Deflation. Siam J Matrix Anal A 23, 948–973 (2002).
% [4] Granat, R., Kågström, B., Kressner, D., and Shao, M. Algorithm 953:
%     Parallel Library Software for the Multishift QR Algorithm with Aggressive
%     Early Deflation ACM Trans. Math. Softw. 41, 4, Article 29 (2015).
%


% Input validation
narginchk(1, 4), nargoutchk(0, 2);

if ~ismatrix(H) 
    error('A must be a square matrix');
end

if ~isa(H, 'quaternion')
    error('A must be a quaternion matrix'); 
end

if ~isHessenberg(H)
    error('A must be an upper Hessenberg matrix');
end

% Set defaults
if nargin < 4, progress_bar = false; end  
if nargin < 3, verbose = false; end
if nargin < 2, rtol = eps; end

% Helper functions
householder_lapply = @(u, x) x - u * (u' * x);
householder_rapply = @(u, x) x - (x * u) * u';

% Initialize timers
H_mult_time = 0;        % Time for H matrix multiplications
Q_time = 0;             % Time to compute Q
AED_time = 0;           % Time for AED
House_time = 0;         % Time to compute Householder matrices

n    = size(H, 1); ihi = n; 
Q    = eyeq(n);
r    = eye(n, 1);
atol = rtol * norm(H, 'fro');
Tz   = tic;
GS   = 0;   % total shifts
DA   = 0;   % aed deflated

if progress_bar, b = waitbar(0, strcat('aed, n = ', num2str(n))); end

ilo = 1;
HH  = H;
while ihi >= 1
    while (ihi - ilo > 1)
        NS = aed_num_shifts(ihi - ilo + 1);  
        WS = min(aed_win_size(ihi - ilo + 1, NS), ihi - ilo);
        sp = max(ihi - WS, ilo);

        % case 1: aggressive early deflation
        if ihi - ilo + 1 > aed_min_size() && sp > ilo
            T0 = tic;
            [U, H(sp:ihi, sp:ihi), shifts] = aed_step(H(sp:ihi, sp:ihi), rtol);
            H(sp:ihi, ihi+1:n) = U' * H(sp:ihi, ihi+1:n);
            H(1:sp-1,  sp:ihi) = H(1:sp-1,   sp:ihi) * U;
            AED_time = AED_time + toc(T0);

            T0 = tic;
            Q( :   , sp:ihi) = Q( :   , sp:ihi) * U;
            Q_time = Q_time + toc(T0);

            DA = DA + ihi - sp - numel(shifts); ihi = sp + numel(shifts);
        % case 2: wilkinson shift
        else
            shifts = shift(H(ihi-1:ihi, ihi-1:ihi));
        end

        if (ihi - sp + 1) / WS < (1 - nibble())   % if sufficient deflation
            continue
        end

        LS = 0;
        NS = min(numel(shifts), NS);

        % case 3: implicit qr steps
        while (ihi > ilo + 1) && (LS < NS)
            % small-bulge shifts
            GS = GS + 1;
            LS = LS + 1;
            x  = shifts(LS);
            si = - 2 * x.w;
            ti = abs(x) ^ 2;
            v  = [
                H(ilo, ilo) ^ 2 + H(ilo, ilo+1) * H(ilo+1, ilo) + si * H(ilo, ilo) + ti; ...
                H(ilo+1, ilo) * H(ilo, ilo) + H(ilo+1, ilo+1) * H(ilo+1, ilo) + si * H(ilo+1, ilo); ...
                H(ilo+2, ilo+1) * H(ilo+1, ilo);
            ];

            T0 = tic;
            [u, ~] = householder_vector(v, r(1:3));
            House_time = House_time + toc(T0);

            T0 = tic;
            H(ilo:ilo+2, ilo:n) = householder_lapply(u, H(ilo:ilo+2, ilo:n));
            H(1:ihi, ilo:ilo+2) = householder_rapply(u, H(1:ihi, ilo:ilo+2));
            H_mult_time = H_mult_time + toc(T0);

            T0 = tic;
            Q( :   , ilo:ilo+2) = householder_rapply(u, Q( :   , ilo:ilo+2));
            Q_time = Q_time + toc(T0);

            for i = ilo : (ihi - 2)
                e  = min(i + 3, ihi);
                sp = max(i - 1, ilo);

                T0 = tic;
                [u, ~] = householder_vector(H(i+1:e, i), r(1:e-i));
                House_time = House_time + toc(T0);

                T0 = tic;
                H(i+1:e, sp:n ) = householder_lapply(u, H(i+1:e, sp:n ));
                H(1:ihi, i+1:e) = householder_rapply(u, H(1:ihi, i+1:e));
                H(i+2:e, i    ) = 0;  % Zero subdiagonal element
                H_mult_time = H_mult_time + toc(T0);

                T0 = tic;
                Q( :   , i+1:e) = householder_rapply(u, Q( :   , i+1:e));
                Q_time = Q_time + toc(T0);

                if abs(H(i+1, i)) <= atol && abs(H(i+1, i)) <= rtol * (abs(H(i, i)) + abs(H(i+1, i+1)))
                    [U, ~] = schur1(H(i, i));
                    Q( : , i) = Q( : ,  i) * U;
                    H( : , i) = H( : ,  i) * U;
                    H(i,  : ) = U' * H(i,  : );
                    H(i+1, i) = 0; ilo = i+1;
                end
            end

            while abs(H(ihi, ihi-1)) <= rtol * (abs(H(ihi, ihi)) + abs(H(ihi-1, ihi-1)))
                [U, ~] = schur1(H(ihi, ihi));
                Q( : , ihi) = Q( : ,  ihi) * U;
                H( : , ihi) = H( : ,  ihi) * U;
                H(ihi,  : ) = U' * H(ihi,  : );
                H(ihi, ihi-1) = 0; ihi = ihi - 1;
                if progress_bar, waitbar(1 - ihi / n, b); end
                if ihi <= ilo, break; end
            end
        end
    end
    if ihi > ilo
        [U, ~] = schur2(H(ilo:ihi, ilo:ihi));
        Q(   :   , ilo:ihi) = Q(  :   , ilo:ihi) * U;
        H(   :   , ilo:ihi) = H(  :   , ilo:ihi) * U;
        H(ilo:ihi, ilo:  n) = U' * H(ilo:ihi, ilo:n);
        H(ihi, ilo) = 0;
    elseif ihi == ilo
        [U, ~] = schur1(H(ilo, ilo));
        Q( : , ilo) = Q( : , ilo) * U;
        H( : , ilo) = H( : , ilo) * U;
        H(ilo,  : ) = U' * H(ilo,  : );
    end
    ihi = ilo - 1; ilo = 1;
end

T = triu(H);

if verbose
    fprintf('Total step of AED %d\n'                   , GS               );
    fprintf('Time for H matrix multiplications: %f s\n', H_mult_time      );
    fprintf('Time for computing Q matrix: %f s\n'      , Q_time           );
    fprintf('Time for Householder computations: %f s\n', House_time       );
    fprintf('Time for AED: %f s, AED deflation: %d\n'  , AED_time, DA     );
    fprintf('Total time: %f s, total deflation: %d\n'  , toc(Tz), n       );
    fprintf('Time without constructing Q: %f s\n'      , toc(Tz) - Q_time );
end
end


% Helper functions
function y = isHessenberg(A)
% Returns true if A is upper Hessenberg
y = all(abs(tril(A, -2)) < eps); 
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

function [Q, H, shifts] = aed_step(H, rtol)
% AED_STEP Aggressive Early Deflation step
n = size(H, 1); ihi = n;
[Q, H(2:n, 2:n)] = iqrq(H(2:n, 2:n), rtol);
H(1, 2:n) = H(1, 2:n)  * Q;
H(2:n, 1) = Q' * H(2:n, 1);
Q = blkdiag(eyeq(1), Q);

for i = 2 : n
    if abs(H(ihi, 1)) <= rtol * abs(H(ihi, ihi))
        [U, ~] = schur1(H(ihi, ihi)); H(ihi, 1) = 0;
        Q( :   , ihi) = Q( :   , ihi) * U;
        H( :   , ihi) = H( :   , ihi) * U;
        H(ihi,   :  ) = U' * H(ihi,  :  );
        ihi = ihi - 1;
    else
        U = ordschurq(H(2:ihi, 2:ihi));
        Q( :   , 2:ihi) = Q( :   ,  2:ihi)  * U;
        H( :   , 2:ihi) = H( :   ,  2:ihi)  * U;
        H(2:ihi,  :  ) = U'  * H(2:ihi,  :    );
        H(2:ihi, 2:ihi) = triu(H(2:ihi, 2:ihi));
    end
end
shifts = diag(H(2:ihi, 2:ihi));

[U, H(1:ihi, 1:ihi)] = hessq(H(1:ihi, 1:ihi));
H(1:ihi, ihi+1:n) = U' * H(1:ihi, ihi+1:n);
Q( :   ,   1:ihi) = Q( :   ,    1:ihi) * U;
end

function NS = aed_num_shifts(ihi)
if ihi < 30
    NS = 2;
elseif ihi < 60
    NS = 4;
elseif ihi < 150
    NS = 10;
elseif ihi < 590
    NS = max(0, ihi / round(log2(double(ihi))));
elseif ihi < 3000
    NS = 64;
else
    NS = 128;
end

NS = max(2, NS - mod(NS, 2));
end


function WS = aed_win_size(ihi, NS)
if ihi <= 500
    WS = NS;
else
    WS = (3 * NS) / 2;
end

WS = max(4, WS - mod(WS, 2));
end


function MS = aed_min_size()
MS = 12;
end


function N = nibble()
N = 0.14;
end


function err = rel_err(A, B)
err = norm(A - B, 'fro') / max(norm(A, 'fro'), norm(B, 'fro'));
end
