function [Q, T] = schur1(A)
% SCHUR1 Standardize a quaternion: unit U with U'*A*U = w + |v|*i.
% A: 1-by-1 quaternion matrix.
[Q, ~] = schur(adjoint(A));
Q = quaternion(real(Q(1)), imag(Q(1)), -real(Q(2)), imag(Q(2)));
T = Q' * A * Q;
end
