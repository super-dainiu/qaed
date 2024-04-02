function y = axis(x)
% AXIS  Axis of quaternion.
% If Q = A + mu .* B where A and B are real/complex, and mu is a unit pure
% quaternion, then axis(Q) = mu.

% Copyright © 2005, 2009, 2020 Stephen J. Sangwine and Nicolas Le Bihan.
% See the file : Copyright.m for further details.

narginchk(1, 1), nargoutchk(0, 1)

if isa(x.z, 'sym')

    % The symbolic case is simpler than the numeric, since we don't need to
    % be concerned about divide by zero, which is what complicates the code
    % below for numeric quaternions.
    
    y = unit(v(x));
    return
end

% In the warning check that follows, abs is applied to the normq result in
% order to handle complex quaternions with a complex semi-norm and modulus.
% It has no effect in the real case because the norm is positive and real.

undefined = abs(normq(v(x))) < eps;

if any(undefined(:))
    warning('QTFM:information', ...
            'At least one element may have an undefined axis.');
end

z = v(x);

% z is now the vector part of x, but in some cases this value may be zero.
% These cases will correspond to the entries in undefined (maybe some
% others with very small modulus). We set all these undefined values to a
% safe temporary value, then replace the temporary value with a zero vector
% after evaluating the unit function. In order to perform the subscripted
% assignment using the logical array undefined we have to use subsasgn and
% substruct because normal indexing does not work here in a class method.
% See the file 'implementation notes.txt', item 8, for more details.

z = subsasgn(z, substruct('()', {undefined}), quaternion(1,1,1));
y = unit(z);
y = subsasgn(y, substruct('()', {undefined}), quaternion(0,0,0));

end

% $Id: axis.m 1125 2021-08-15 20:43:38Z sangwine $
