function r = unit(o)
% UNIT octonion. Divides an octonion by its own modulus.
% The result is an octonion with unit modulus.

% Copyright © 2011, 2021 Stephen J. Sangwine and Nicolas Le Bihan.
% See the file : Copyright.m for further details.

narginchk(1, 1), nargoutchk(0, 1) 

m = abs(o);

if isa(x(o.a), 'sym')
    % The input argument is symbolic, therefore we can check for equality
    % with zero.
    % TODO Should we not check for zero modulus? We could be dealing with
    % a divisor of zero. Testing for equality with zero won't check this.

    if any(any(abs(m) == 0))
        warning(['At least one element of the argument to the unit '...
                 'function is zero, and divide by zero will occur.']);
    end
    
    r = o ./ m;
    
    return % We must not do the numerical checks below on a symbolic result
    
else
    % Since m could be complex, the warning check below has to
    % use abs again to get a real result for comparison with eps.
    
    if any(any(abs(m) < eps))
        warning('QTFM:information', ...
            ['At least one element has zero modulus, '...
             'and divide by zero will occur.']);
    end
end

r = o ./ m;
 
% Dividing by a small modulus can result in numerical errors such that the
% result does not have unit modulus. This is especially likely with complex
% octonions, so we perform a check here.

% The modulus (and norm) of each element of r should be either 1 or 1 + 0i. 
% In either case, subtracting ones should result in a real or complex
% number with very small real and imaginary parts, which we compare with
% epsilon using an arbitrary scale factor which we choose so that the test
% is not too sensitive.

d = normo(r) - 1;
n = nnz(abs(real(d)) > 1e3.*eps | abs(imag(d)) > 1e3.*eps);

if n > 0
    if n == 1
        warning('QTFM:information', ...
            ['One element of the result of the unit '...
             'function has a modulus which is not accurately 1.']);
    else
        warning('QTFM:information', ...
            [num2str(n), ' elements of the result of the unit '...
             'function have a modulus which is not accurately 1.']);
    end
end

% Discussion note. See the note in the corresponding quaternion function
% for possible issues with accuracy of unit complex octonions.

end

% $Id: unit.m 1127 2021-08-30 19:35:41Z sangwine $

