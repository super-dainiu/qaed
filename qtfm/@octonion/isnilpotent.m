function tf = isnilpotent(q, tol)
% ISNILPOTENT  True where any element of q is a nilpotent (to within the
% tolerance given (optionally) by the second parameter, if q is numeric).

% Copyright © 2019, 2022 Stephen J. Sangwine and Nicolas Le Bihan.
% See the file : Copyright.m for further details.

narginchk(1, 2), nargoutchk(0, 1)

% Reference:
%
% Stephen J. Sangwine and Daniel Alfsmann,
% Determination of the Biquaternion Divisors of Zero, Including the
% Idempotents and Nilpotents, Advances in Applied Clifford Algebras, 20,
% (2010), 401–410. DOI 10.1007/s00006-010-0202-3.

% Theorem 5 of the above paper gives the conditions for a biquaternion to
% be a nilpotent. The biquaternion must be pure, and the norm must be zero.
% The same holds for octonions, although not covered in the paper.

if isnumeric(x(q.a))

    if nargin == 1
        tol = 8 .* eps; % The tolerance was not specified, supply a default.
    end

    % We also need to check for a non-zero imaginary part, otherwise we
    % would return true for a value of zero, which is not a nilpotent in
    % the sense needed here.

    if ispure(q)
        tf = abs(normo(imag(q))) > tol & abs(abs(normo(q))) < tol;
    else
        % Nilpotents must be pure, so we check that the scalar part is
        % close to zero, then check the norm of the vector part is also
        % close to zero.

        tf = abs(normo(imag(q))) > tol & ...
             s(q) < tol & abs(abs(normo(v(q)))) < tol;
    end

elseif isa(x(q.a), 'sym')

    if nargin == 2
        warning('First parameter is symbolic, ignoring tolerance parameter.')
    end

    % In the symbolic case, we are able to do an exact check on the norm.

    if isreal(q)
        tf = false(size(q)); % A real quaternion cannot be a nilpotent.
    else
        if ispure(q)
            tf = eval(normo(q) == 0); % The norm must be identically zero,
                                      % even though it is complex.
        else
            % If q is not pure, there could be a zero scalar part, so we
            % check for that, and then check whether the norm is zero as
            % above.
            tf = eval(s(q) == 0) & eval(normo(q) == 0);
        end
    end

else
    error('First argument is neither numeric nor symbolic, cannot handle.')
end

end

% $Id: isnilpotent.m 1138 2022-02-13 20:59:56Z sangwine $
