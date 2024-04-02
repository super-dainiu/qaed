function [Q, T] = ordschurq(T)
% ORDSCHURQ Reorder eigenvalues in Schur form 
%
% [Q, T] = ORDSCHURQ(T) moves the (n,n) element of the upper triangular
% matrix T to the (1,1) position by applying a sequence of swaps.
%
% Q is a unitary quaternion matrix such that:
%
%   Q' * T * Q = Reordered T
%
% References:
% [1] D. Kressner, Block algorithms for reordering standard and
%     generalized Schur forms, ACM Transactions on Mathematical Software
%     (TOMS), 32.4 (2006), pp. 521-532.

narginchk(1, 1), nargoutchk(0, 2)

n = size(T, 1); ihi = n; Q = eyeq(n);

while ihi > 1
    ilo = ihi - 1;
    U   = swapq(T(ilo : ihi, ilo : ihi));

    T(ilo : ihi,     :    ) = U' * T(ilo : ihi,     :    );
    T(    :    , ilo : ihi) = T(    :    , ilo : ihi) *  U;
    Q(    :    , ilo : ihi) = Q(    :    , ilo : ihi) *  U;

    ihi = ilo;
end

end