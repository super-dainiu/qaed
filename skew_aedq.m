function [Q, D] = skew_aedq(H, rtol, verbose, progress_bar)
% SKEW_AEDQ Aggressive Early Deflation Algorithm for Skew-Symmetric Quaternion Hessenberg Matrices
%
% [Q, D] = SKEW_AEDQ(H) applies the aggressive early deflation algorithm to the
% n-by-n skew-symmetric quaternion Hessenberg matrix H for computing its eigenvalues
% and corresponding eigenvectors. The algorithm enhances eigenvalue
% computation efficiency and accuracy for skew-symmetric matrices.
% Q is a unitary quaternion matrix and D is a diagonal quaternion
% matrix satisfying:
%
%   Q' * H * Q = D
%
% [Q, D] = SKEW_AEDQ(H, RTOL) uses tolerance RTOL instead of the default eps.
%
% [Q, D] = SKEW_AEDQ(H, RTOL, VERBOSE) prints additional 
% information about the computation if VERBOSE is true.
%
% [Q, D] = SKEW_AEDQ(H, RTOL, VERBOSE, PROGRESS_BAR) displays
% a progress bar during computation if PROGRESS_BAR is true.
%
% Example Usage:
% [Q, D] = skew_aedq(H); % Default usage
% [Q, D] = skew_aedq(H, 1e-12, true, true); % With all options specified

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
AED_time = 0;
House_time = 0;

n = size(H, 1);
ilo = 1;
ihi = n;
Q = eyeq(n);
r = eye(n, 1);
atol = rtol * norm(H, 'fro');
Tz = tic;
GS = 0;   % total shifts
DA = 0;   % aed deflated
HH = H;

if progress_bar
    b = waitbar(0, strcat('skew_aedq, n = ', num2str(n)));
end

while ihi > 1
    while ihi - ilo > 2
        NS = aed_num_shifts(ihi - ilo + 1);  
        WS = aed_win_size(ihi - ilo + 1, NS);
        win = min(WS, ihi - ilo + 1);
        whi = ihi;
        wlo = max(ihi - win + 1, ilo);
        sp = wlo - 1;

        %% case 1: aggressive early deflation
        if ihi - ilo + 1 >= aed_min_size() && sp > ilo && win > 4
            T0 = tic;
            [U, D] = skew_iqrq(H(wlo:whi, wlo:whi), rtol, false, false);
            D = diag(D);
            spike = U' * H(wlo:whi, sp); 
            [~, index] = sort(abs(spike), "descend");
            spike = spike(index);
            H(wlo:whi, wlo:whi) = diag(D(index));

            Q(:, wlo:whi) = Q(:, wlo:whi) * U;
            Q(:, wlo:whi) = Q(:, sp + index);
            Q_time = Q_time + toc(T0);

            for i = win:-1:1
                if abs(spike(i)) <= rtol * abs(D(i))
                    spike(i) = 0;
                    whi = whi - 1;
                end
            end
            H(wlo:ihi, sp) = spike;
            H(sp, wlo:ihi) = -spike';
            shifts = diag(H(wlo:whi, wlo:whi));
            ihi = whi;
            shifts = shifts(end:-1:1);

            % Hessenberg reduction
            T0 = tic;
            [U, H(sp:ihi, sp:ihi)] = hessq(H(sp:ihi, sp:ihi));
            Q(:, sp:ihi) = Q(:, sp:ihi) * U;
            Q_time = Q_time + toc(T0);

            AED_time = AED_time + toc(T0);
            DA = DA + win - numel(shifts);
        else
        %% case 2: Wilkinson shift
            shifts = shift(H(ihi-1:ihi, ihi-1:ihi));
        end

        if (ihi - sp + 1) / win < (1 - nibble())   % if sufficient deflation
            continue
        end

        LS = 0;
        NS = min(numel(shifts), NS);

        %% case 3: implicit qr steps
        while (ihi > ilo + 1) && (LS < NS)
            % small-bulge shifts
            GS = GS + 1;
            LS = LS + 1;
            x = shifts(LS);
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
                e = min(i + 3, ihi);
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
    fprintf('Total steps of AED: %d\n', GS);
    fprintf('Time for H matrix multiplications: %f s\n', H_mult_time);
    fprintf('Time for computing Q matrix: %f s\n', Q_time);
    fprintf('Time for Householder computations: %f s\n', House_time);
    fprintf('Time for AED: %f s, AED deflation: %d\n', AED_time, DA);
    fprintf('Total time: %f s, total deflation: %d\n', toc(Tz), n);
    fprintf('Time without constructing Q: %f s\n', toc(Tz) - Q_time);
end

end

% Helper functions
function y = istridiagonal(A)
% Returns true if A is tridiagonal
y = normq(tril(A, -2)) / normq(A) < eps & normq(triu(A, 2)) / normq(A) < eps;
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
% Compute an ordered Schur decomposition of a 2x2 skew-symmetric matrix A.
[Q, T] = schur(adjoint(A));
[Q, ~] = ordschur(Q, T, groupByValue(abs(diag(T))));
Q = quaternion( ...
    real(Q([1,2], [1,3])), ...
    imag(Q([1,2], [1,3])), ...
    -real(Q([3,4], [1,3])), ...
    imag(Q([3,4], [1,3])) ...
);
T = Q' * A * Q;
% Ensure the result is skew-symmetric
T = (T - T') / 2;
end

function group = groupByValue(A)
[~, idx] = sort(A);
group = zeros(size(A));
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