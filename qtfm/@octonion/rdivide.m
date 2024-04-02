function r = rdivide(a, b)
% ./  Right array divide.
% (Octonion overloading of standard Matlab function.)

% Copyright © 2011 Stephen J. Sangwine and Nicolas Le Bihan.
% See the file : Copyright.m for further details.

narginchk(2, 2), nargoutchk(0, 1)

% There are three cases to be handled:
%
% 1. Left and right arguments are octonions.
% 2. The left argument is an octonion, the right is not.
% 3. The right argument is an octonion, the left is not.
%
% In fact, cases 1 and 3 can be handled by the same code.
% Case 2 requires different handling.

if isa(b, 'octonion')
    
    % The right argument is an octonion. We can handle this case by
    % forming its elementwise inverse and then multiplying. Of course,
    % if any elements have zero norm, this will result in NaNs.
     
    r = a .* b .^ -1;
    
elseif isnumeric(b) || isa(b, 'sym')
    
    % The right argument is numeric or symbolic. We assume therefore that
    % if we divide components of the left argument by the right argument,
    % that Matlab will do the rest. Obviously if the right argument is
    % zero, there will be a divide by zero error. There may also be errors
    % caused by type mismatch (for example the components of a are double,
    % but b is int8), but Matlab will catch these.

    r = overload(mfilename, a, b);
    
elseif isa(b, 'quaternion')
    error('Cannot divide an octonion by a quaternion')
else
    error(['Unable to divide octonion by a ' class(b)])
end

end

% $Id: rdivide.m 1125 2021-08-15 20:43:38Z sangwine $
