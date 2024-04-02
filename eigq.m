function [P, D] = eigq(A, rtol, verbose, progress_bar)
% EIGQ Quaternion Eigenvalue and Eigenvector Computation
%
% [P, D] = EIGQ(A) computes the eigenvalues 
% (D) and the corresponding eigenvectors (P) of a quaternion matrix A.
% The function can handle both general and tridiagonal quaternion matrices,
% applying different algorithms based on the matrix size and properties.
%
% [P, D] = EIGQ(A, RTOL) uses tolerance RTOL instead of the default eps.
%
% [P, D] = EIGQ(A, RTOL, VERBOSE) prints additional information about the
% computation if VERBOSE is true.
%
% [P, D] = EIGQ(A, RTOL, VERBOSE, PROGRESS_BAR) displays a progress bar
% during computation if PROGRESS_BAR is true.
%
%
% The function first checks if A can be reduced to a tridiagonal matrix.
% For tridiagonal matrices, it chooses between the `tridiag_qr` and
% `tridiag_aed` algorithms based on the size of A. For non-tridiagonal
% matrices, it performs a preliminary Hessenberg decomposition and then
% applies either the `iqrq` or `aedq` algorithm depending on the matrix
% size. The choice of algorithm optimizes computational efficiency for
% different matrix sizes and types.
%
% If A is not tridiagonal, the function also performs back-substitution to 
% compute the eigenvectors.
%
% Example Usage:
% [P, D] = eigq(A); % Default usage
% [P, D] = eigq(A, 1e-12, true, true); % With all options specified
%
% References:
%
% [1] Bunse-Gerstner, et al., A Quaternion QR Algorithm, Numer. Math. 55,
%     83–95 (1989) 
% [2] Golub, G. H. and Van Loan, C. F., Matrix Computations, Algorithm
%     7.4.2 & Algorithm 7.5.2, Johns Hopkins University Press, 1996.
% [3] Braman, K., Byers, R. & Mathias, R. The Multishift QR Algorithm. Part II:
%     Aggressive Early Deflation. Siam J Matrix Anal A 23, 948–973 (2002).
% [4] Granat, R., Kågström, B., Kressner, D., and Shao, M. Algorithm 953:
%     Parallel Library Software for the Multishift QR Algorithm with Aggressive
%     Early Deflation ACM Trans. Math. Softw. 41, 4, Article 29 (2015).
% [5] Jia, Z., Wei, M., Zhao, M.-X. & Chen, Y. A new real structure-preserving 
%     quaternion QR algorithm. J Comput Appl Math 343, 26–48 (2018).

narginchk(1, 4), nargoutchk(0, 2)

% Set defaults
if nargin < 4, progress_bar = false; end  
if nargin < 3, verbose = false; end
if nargin < 2, rtol = eps; end

[P, H] = hessq(A); n = size(H);

if ~istridiagonal(H)
    if n < 1000
        [Q, T] = iqrq(H, rtol, true, verbose, progress_bar);
    else
        [Q, T] = aedq(H, rtol, true, verbose, progress_bar);
    end

    [P, D] = eigvec(P * Q, T);
else
    if n < 100
        [Q, D] = tridiag_qr(H, rtol, true, verbose, progress_bar);
    else
        [Q, D] = tridiag_aed(H, rtol, true, verbose, progress_bar);
    end

    P = P * Q;
    D = diag(D);
end

if nargout == 1, P = diag(D); end
end


% Helper functions
function y = istridiagonal(A)
% Returns true if A is tridiagonal
y = all(abs(tril(A, -2)) < eps, "all") & all(abs(triu(A, 2)) < eps, "all");
end