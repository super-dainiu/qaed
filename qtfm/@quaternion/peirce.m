function [chi, d, a, b] = peirce(q)
% PEIRCE    Peirce decomposition of a quaternion.
%
% Computes a decomposition of a real quaternion into a pair of conjugate
% complex eigenvalues. The eigenvalues permit functions of the quaternion
% argument to be computed using complex functions (for example, powers,
% trigonometric functions).
%
% chi is the 'positive' eigenvalue. The complex conjugate of chi gives the
% 'negative' eigenvalue.
%
% d is the 'positive' idempotent. This is a biquaternion with real scalar
% part and imaginary vector part. d^2 = d (definition of an idempotent).
%
% a and b are the products of the eigenvalues and idempotents such that
% q = a + b (see code for the formula).
%
% The utility of this decomposition is that it satisfies:
%
% f(q) = f(chi) .* d* + f(chi*) .* d, where d* is the complex conjugate of
% d and chi* is the complex conjugate of chi, for a wide range of
% functions, and the RHS may be evaluated using only COMPLEX
% implementations of f(x). Note that because the two terms are conjugates,
% it is possible to compute the RHS as real(f(chi) .* d*) .* 2 with less
% computation (indeed one could go further and avoid computing the
% imaginary part rather than discarding it afterwards).

% Copyright © 2019 Stephen J. Sangwine and Nicolas Le Bihan.
% See the file : Copyright.m for further details.

% Reference:
%
% Roger M. Oba, 'Eigen-Decomposition of Quaternions', Advances in Applied
% Clifford Algebras, Oct 2018, 28(5), p94. doi:10.1007/s00006-018-0911-6

narginchk(1, 1), nargoutchk(0, 4)

if ~isreal(q)
    error('The Peirce decomposition is not defined for biquaternions')
end

vq = vector(q); % We need this twice below, so compute it once here.

% Compute the positive idempotent. This is given in Oba, equations 2.1 and
% 3.1 (the latter defines the unit vector part of q). The negative
% idempotent (d_minus in the paper) is just the complex conjugate of the
% value below, so we do not compute it.

A = abs(vq);

% There is no complex() function in the symbolic toolbox, which is why we
% use explicit multiplication by 1i and addition to construct a complex
% quantity. ones(class(q.x)) is necessary to avoid an error when the
% imaginary part is real.

if nargout > 1
    if isnumeric(q.x)
        d = complex(1, (vq ./ A)) ./ 2;
    elseif isa(q.x, 'sym') || islogical(q.x)
        % There is no complex() function in the symbolic toolbox, which is
        % why we use explicit multiplication by 1i and addition to
        % construct a complex quantity.
        d = (1 + 1i .* (vq ./ A)) ./ 2;
    else
        error(['Unable to handle quaternion with ' class(q.x) ' components'])
    end
end

% Compute the eigenvalues. These are a complex conjugate pair of complex
% numbers, so again we need compute only one of them.

chi = scalar(q) + A .* 1i; % Oba, equation 4.3.

% Now compute the two components of the decomposition, if required (many
% applications of this function would not need them). Oba, Theorem 4.7.

if nargout > 2
    a = chi .* conj(d, 'complex'); 
end

if nargout > 3
    b = conj(chi) .* d;
end

end

% $Id: peirce.m 1113 2021-02-01 18:41:09Z sangwine $
