function Y = sqrt(X)
% SQRT   Square root.
% (Quaternion overloading of standard Matlab function.)

% Copyright © 2005, 2006 Stephen J. Sangwine and Nicolas Le Bihan.
% See the file : Copyright.m for further details.

narginchk(1, 1), nargoutchk(0, 1)

if isreal(X) && ~ispure(X)
    
    % X is a full real quaternion, and we compute the square root of an
    % isomorphic complex number using the standard Matlab square root
    % function, then construct a quaternion with the same axis as the
    % original quaternion.
    
    Y = isoquaternion(sqrt(isocomplex(X)), X);
else
    
    % X is either a complex quaternion or a pure quaternion, and therefore
    % we cannot use the method above, because it is not possible to
    % construct an isomorphic complex number. Therefore we use polar form
    % and halve the argument. Note that the modulus and argument here are
    % complex, so the square root of the modulus is complex.
    
    % Preconstruct a result array of zero quaternions, of the same size as
    % X and with components of the same class.
    
    Y = zerosq(size(X), class(x(X)));
    
    % If the vector part of any element of X is zero, computing its axis
    % will result in an undefined axis warning. We avoid this warning by
    % computing just the square root of the scalar part in these cases.
    % There may be no such cases, or all the elements in X may have
    % undefined axis.
    
    UN = undefined(X);
    DN = ~UN;
        
    % In order to perform the subscripted assignment using the logical
    % array UN (for undefined) we have to use subsasgn and substruct
    % because normal indexing does not work here in a class method.
    % See the file 'Implementation notes.txt', item 8, for more details.

    if nnz(UN) > 0
        % There are some cases of undefined axis.
        U = scalar(subsref(X, substruct('()',  {UN})));
        Y = subsasgn(Y, substruct('()', {UN}), quaternion(sqrt(U)));
    end
    
    if nnz(DN) > 0
        % There are some cases where the axis is defined.
        D = subsref(X, substruct('()', {DN}));    
        Y = subsasgn(Y, substruct('()', {DN}), ...
                    sqrt(abs(D)) .* exp(axis(D) .* angle(D) ./ 2));
    end
end

end

% TODO Consider whether a better algorithm could be devised for the complex
% case based on the Cartesian form. A formula for the complex square root
% is given as eqn 12 on page 17 of:
%
% John H. Mathews,
% 'Basic Complex Variables for Mathematics and Engineering',
% Allyn and Bacon, Boston, 1982. ISBN 0-205-07170-8.

% There is a formula given in the following book, which may be worth
% investigation. It is: Y = (1 + X) ./ sqrt(2 .* (1 + X.w)). Notice this
% requires adding 1 to the scalar part, and then dividing by a real value.
% The formula is valid only for unit quaternions, so an adjustment would be
% needed to handle the general case.
%
% Andrew J. Hanson, 'Visualizing quaternions',
% Morgan Kaufmann, 2006, ISBN 0-12-088400-3, Section F2, p452.
%
% Hanson's formula was tried in March 2021, and looked promising but it
% failed for one of the test cases in the test code. For the record here is
% what was tried. Conclusion - needs more work yet.
%
% % Hanson's formula is valid only for unit X, so we need some extra
% % steps to deal with |X| /= 1.
% 
% N = abs(X);
% Y = X ./ N; % Normalise to a unit quaternion.
% 
% % TODO The formula fails for Y = -1, resulting from a complex quaternion
% % with zero vector part.
% 
% Z = sqrt(N) .* (1 + Y) ./ sqrt(2 .* (1 + scalar(Y))); % scalar avoids [].

% A further possibility is Hero or Heron's method, which is a special case
% of Newton's method. This is iterative, but it seems to converge fast even
% for a biquaternion. It may be better than finding the angle and axis.

% $Id: sqrt.m 1120 2021-03-10 20:04:13Z sangwine $
