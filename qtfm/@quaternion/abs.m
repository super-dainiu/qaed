function a = abs(q)
% ABS Absolute value, or modulus, of a quaternion.
% (Quaternion overloading of standard Matlab function.)

% Copyright © 2005, 2015, 2020, 2021, 2022
%             Stephen J. Sangwine and Nicolas Le Bihan.
% See the file : Copyright.m for further details.

narginchk(1, 1), nargoutchk(0, 1)

isnum = isnumeric(q.x);

if isnum && isreal(q)
        
        % We use here a method based on Cayley-Dickson form and the Matlab
        % function hypot, which avoids overflow for large values. This can
        % only be used for the numeric real case.
        
        if isempty(q.w)
            a = hypot(abs(q.x), abs(complex(q.y, q.z)));
        else
            a = hypot(abs(complex(q.w, q.x)), abs(complex(q.y, q.z)));
        end


elseif isnum || isa(q.x, 'sym') || islogical(q.x)

    % If q is numeric at this point, it must be complex, since we handled
    % the real case above.

    if isempty(q.w)
        a = sqrt(         q.x.^2 + q.y.^2 + q.z.^2);
    else
        a = sqrt(q.w.^2 + q.x.^2 + q.y.^2 + q.z.^2);
    end
    
else
    error(['Unable to handle quaternion with ' class(q.x) ' components'])
end

end

% $Id: abs.m 1136 2022-02-11 20:53:24Z sangwine $
