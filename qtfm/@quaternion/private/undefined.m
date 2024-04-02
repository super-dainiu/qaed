function U = undefined(Q)
% UNDEFINED  Detect elements of quaternion array with undefined axis value.

% Copyright © 2021 Stephen J. Sangwine and Nicolas Le Bihan.
% See the file : Copyright.m for further details.

narginchk(1, 1), nargoutchk(0, 1)

% Q must be a quaternion array. The result is a logical array which is true
% in positions where the axis of Q is identically zero (symbolic case), or
% less than eps (numeric case).

M = normq(vector(Q));

if isnumeric(Q.x)
    U = abs(M) < eps; % We have to use abs here in case M is complex.
elseif isa(Q.x, 'sym')

    % Some elements of M may be identically zero and isAlways may be able
    % to prove them zero. However, others may be undecidable, and we
    % suppress the warning message with the 2nd and 3rd parameters.
    
    U = isAlways(M == 0, 'Unknown', 'false');
else
    % TODO It is possible we could deal with logical without difficulty but
    % do we need to?
    error(['Cannot handle argument of class: ', class(Q), ...
           ' in private quaternion function UNDEFINED.'])
end

end

% $Id: undefined.m 1120 2021-03-10 20:04:13Z sangwine $
