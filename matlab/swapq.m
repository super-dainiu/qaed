function [Q, T] = swapq(T)
% SWAPQ Swap eigenvalues in 2x2 quaternion Schur form
%
% [Q, T] = SWAPQ(T) computes a unitary Q such that Q' * T * Q swaps
% the diagonal elements of the 2x2 upper triangular T.
%
% References:
% [1] Bai, Z., & Demmel, J. W. (1993). On swapping diagonal blocks in 
%     real Schur form. Linear Algebra and its Applications, 186, 75-95.
% 

narginchk(1, 1), nargoutchk(0, 2)

assert(isa(T, 'quaternion')   , 'Not a quaternion matrix.'    );
assert(isequal(size(T), [2 2]), 'Not a quaternion 2x2 block.' );
assert(abs(T(2, 1)) == 0      , 'Not a quaternion Schur form.');

x = sylvesterc(T(1, 1), T(2, 2), -T(1, 2));
z = sqrt(1 + abs(x) ^ 2); c = x / z; s = 1 / z;
Q = [c, -s; s, c'];

if nargout > 1
    T = Q' * T * Q;
    T = triu(T);
end

end