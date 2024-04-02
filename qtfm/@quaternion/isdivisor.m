function tf = isdivisor(q, tol)
% ISDIVISOR  True where any element of q is a divisor of zero (to within
% the tolerance given (optionally) by the second parameter if q is
% numeric).

% Copyright © 2019, 2022 Stephen J. Sangwine and Nicolas Le Bihan.
% See the file : Copyright.m for further details.

narginchk(1, 2), nargoutchk(0, 1)

% Reference:
%
% Stephen J. Sangwine and Daniel Alfsmann,
% Determination of the Biquaternion Divisors of Zero, Including the
% Idempotents and Nilpotents, Advances in Applied Clifford Algebras, 20,
% (2010), 401–410. DOI 10.1007/s00006-010-0202-3.

% Theorem 1 of the above paper gives the conditions for a quaternion to be
% a divisor of zero. Two quantities have to be zero as computed immediately
% below.

normdiff = abs(normq(real(q)) - normq(imag(q)));
sproduct = abs(scalar_product(real(q), imag(q)));

if isnumeric(x(q))
    if nargin == 1
        tol = eps; % The tolerance was not specified, supply a default.
    end
    tf = normdiff < tol & sproduct < tol;
elseif isa(x(q), 'sym')
    if nargin == 2
        warning('First parameter is symbolic, ignoring tolerance parameter.')
    end
    tf = eval(normdiff == 0 & sproduct == 0);
else
    error('First argument is neither numeric nor symbolic, cannot handle.')
end

% $Id: isdivisor.m 1138 2022-02-13 20:59:56Z sangwine $
