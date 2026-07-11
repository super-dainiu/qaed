function [P, D] = eigvec(Q, T)
% EIGVEC - Compute the eigenvectors from Schur decomposition
%
% [P, D] = EIGVEC(Q, T) computes the eigenvectors P and eigenvalues D from
% the Schur decomposition of the matrix A = Q * T * Q'.
%
% Q * T * Q' = P * D * P^{-1}
%
% D is diagonal, T is upper triangular with standard Schur form
% and P is the matrix of eigenvectors.
%
% Reference:
% [1] Golub, G. H., and Van Loan, C. F. (2013). Matrix Computations (4th ed.).
%     Johns Hopkins University Press.

% Dispatch to the C++ core when available (see qaed_accel).
if qaed_accel()
    [P, D] = eigvec_cpp(Q, T);
    return
end

P = Q;
n = size(T, 1);

for i = n:-1:2
    x = sylvesterc_tri(T(1:i-1, 1:i-1), T(i, i), T(1:i-1, i));
    P(1:n, i) = P(1:n, i) - P(1:n, 1:i-1) * x;
    P(1:n, i) = P(1:n, i) / norm(P(1:n, i));
end

D = tril(T);
end

function [P, D] = eigvec_cpp(Q, T)
% EIGVEC_CPP Dispatch EIGVEC to the C++ core (cpp/qaed_mex).
[Qw, Qx, Qy, Qz] = qmat_parts(Q);
[Tw, Tx, Ty, Tz] = qmat_parts(T);
[Pw, Px, Py, Pz, Dw, Dx, Dy, Dz] = qaed_mex('eigvec', Qw, Qx, Qy, Qz, Tw, Tx, Ty, Tz);
P = quaternion(Pw, Px, Py, Pz);
D = quaternion(diag(Dw), diag(Dx), diag(Dy), diag(Dz));
end
