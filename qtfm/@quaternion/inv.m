function R = inv(a)
% INV   Inverse of a quaternion (matrix).
% (Quaternion overloading of standard Matlab function.)
%
% For a single quaternion, this function computes the inverse,
% that is the conjugate divided by the modulus. The result will
% be a NaN if the quaternion has zero modulus. For matrices it
% computes the matrix inverse using an analytical formula based
% on partitioning the matrix into sub-matrices. A better method
% may be substituted in the future.

% Copyright © 2005, 2022 Stephen J. Sangwine and Nicolas Le Bihan.
% See the file : Copyright.m for further details.

% Reference: Wikipedia article on 'Matrix Inversion', 12 July 2005.
%            Wikipedia article on 'Schur Complement', 13 July 2005.
%
% A similar formula is given in:
%
% R. A. Horn and C. R. Johnson, 'Matrix Analysis',
% Cambridge University Press, 1985, §0.7.3, page 18.
%
% and in:
%
% Granino A. Korn and Theresa M. Korn, 'Mathematical Handbook for
% Scientists and Engineers' Dover, 2000. ISBN 978-0-486-41147-7.
% Section 20.3-4 equation (20.3-19)

% TODO There is currently no check for a singular matrix. A result is
% returned which is not a valid inverse.

[r, c] = size(a);

if r ~= c
    error('Matrix must be square.');
end

if r == 1 % c must also be 1, since we have checked that r == c above.
    R = conj(a) ./ normq(a);
elseif r == 2

    % We handle the 2 by 2 case separately since this can be computed
    % using quaternion products, rather than quaternion matrix products.

    A = subsref(a, substruct('()', {1,1})); % A = a(1,1);
    B = subsref(a, substruct('()', {1,2})); % B = a(1,2);
    C = subsref(a, substruct('()', {2,1})); % C = a(2,1);
    D = subsref(a, substruct('()', {2,2})); % D = a(2,2);

    IA = inv(A); % If A is zero, this will cause a NaN and the whole
                 % result will be NaNs.

    if any(isinf(IA), 'all') || any(isnan(IA), 'all')
        R = inv_by_adjoint(a); % Use a matrix isomorphism instead.
        return
    end

    IAB = IA .* B; % These two subexpressions are needed more than
    CIA = C .* IA; % once each, and are computed once for efficiency.

    T  = inv(D - CIA .* B); % This might contain Infs or NaNs.

    if any(isinf(IA), 'all') || any(isnan(IA), 'all')
        R = inv_by_adjoint(a); % Use a matrix isomorphism instead.
        return
    end

    TCIA = T .* CIA; % This is needed twice, so compute it once and reuse.

    R = [[IA + IAB .* TCIA, -IAB .* T];...
         [           -TCIA,         T]];
else

    % For matrices larger than 2 by 2, we partition the matrix into sub
    % blocks of half the height/width of the input matrix. If the input
    % matrix has an even number of rows (and therefore columns,  since it
    % is square),  then the sub-blocks will be exactly half the height
    % (width),  otherwise the left block will have one less element,
    % because we round towards zero. The layout of the sub-blocks is:
    %
    %                 [ A B ]
    %                 [ C D ]
    %
    % the same as in the 2 by 2 case above (except that here they are
    % matrices and not scalars).

    p = fix(r./2); % Calculate half the number of rows/columns rounded down.
    q = p + 1;     % Index of first element of right or bottom block.

    A = subsref(a, substruct('()', {1:p, 1:p})); % A = a(1:p, 1:p);
    B = subsref(a, substruct('()', {1:p, q:c})); % B = a(1:p, q:c);
    C = subsref(a, substruct('()', {q:r, 1:p})); % C = a(q:r, 1:p);
    D = subsref(a, substruct('()', {q:r, q:c})); % D = a(q:r, q:c);

    IA = inv(A); % This might contain Infs or NaNs.

    if any(isinf(IA), 'all') || any(isnan(IA), 'all')
        R = inv_by_adjoint(a); % Use a matrix isomorphism instead.
        return
    end

    IAB = IA * B; %#ok<MINV> % These two subexpressions are needed more than
    CIA = C * IA; %#ok<MINV> % once each, and are computed once for efficiency.

    T  = inv(D - CIA * B);   % This might contain Infs or NaNs.

    if any(isinf(IA), 'all') || any(isnan(IA), 'all')
        R = inv_by_adjoint(a); % Use a matrix isomorphism instead.
        return
    end

    TCIA = T * CIA; %#ok<MINV> % This is needed twice, so compute it once and reuse.

    R = [[IA + IAB * TCIA, -IAB * T];...
         [          -TCIA,        T]]; %#ok<MINV>
end

% NOTE A better method could be to use the LU or QR decomposition. In the
% LU case, L and U are easily inverted since they are triangular. In the QR
% case, Q is unitary (but not in the biquaternion case) and is therefore
% trivial to invert and R is triangular. See Numerical Recipes in C, 2nd
% edition, section 2.2 for a discussion of back-substitution. Having
% obtained the inverses, the inverse of the original matrix can be obtained
% from the product of the inverses: if IA is the inverse of A, then Q * R =
% A, IA = IR * IQ (note the reversed order).

% Note on subscripted references: these do not work inside a class method
% (which this is). See the file 'implementation notes.txt', item 8.

end

function Y = inv_by_adjoint(X)
% Computes the inverse of the real/complex matrix X using the
% adjoint/unadjoint isomorphism pair.
T = adjoint(X, 'real', 'block');
warning('off', 'MATLAB:singularMatrix');
U = inv(T); % This is the MATLAB inv on the real/complex matrix T.
warning('on',  'MATLAB:singularMatrix');
% TODO It might be better to handle the warning suppression below in
% unadjoint by recognising that NaNs and Infs are not accuracy issues.
warning('off', 'QTFM:inaccuracy'); % Suppress warnings about inaccurate
Y = unadjoint(U, 'real', 'block'); % adjoints. We expect that here because
warning('on',  'QTFM:inaccuracy'); % we are computing a difficult inverse.
end

% $Id: inv.m 1166 2022-11-19 20:26:29Z sangwine $
