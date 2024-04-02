function a = abs(o)
% ABS Absolute value, or modulus, of an octonion.
% (Octonion overloading of standard Matlab function.)

% Copyright © 2011, 2015, 2021, 2022
%             Stephen J. Sangwine and Nicolas Le Bihan.
% See the file : Copyright.m for further details.

narginchk(1, 1), nargoutchk(0, 1)

isnum = isnumeric(x(o.a));

if isnum && isreal(o)

    % We use here a method based on Cayley-Dickson form and the Matlab
    % function hypot, which avoids overflow for large values. This can only
    % be used for the numeric real case.

    a = hypot(abs(o.a), abs(o.b));

elseif isnum || isa(x(o.a), 'sym') || islogical(x(o.a))

    % If o is numeric at this point, it must be complex, since we handled
    % the real case above.

    a = sqrt(normq(o.a) + normq(o.b));

else
    error(['Unable to handle octonion with ' class(o.x) ' components'])
end

end

% $Id: abs.m 1136 2022-02-11 20:53:24Z sangwine $
