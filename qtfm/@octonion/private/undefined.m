function U = undefined(Q)
% UNDEFINED  Detect elements of octonion array with undefined axis value.

% Copyright © 2021 Stephen J. Sangwine and Nicolas Le Bihan.
% See the file : Copyright.m for further details.

narginchk(1, 1), nargoutchk(0, 1)

% Q must be an octonion array. The result is a logical array which is true
% in positions where the axis of Q is identically zero (symbolic case), or
% less than eps (numeric case).

M = normo(vector(Q));

if isnumeric(x(Q.a))
    U = abs(M) < eps; % We have to use abs here in case M is complex.
elseif isa(x(Q.a), 'sym')

    % Some elements of M may be identically zero and isAlways may be able
    % to prove them zero. However, others may be undecidable, and we
    % suppress the warning message with the 2nd and 3rd parameters.
    
    U = isAlways(M == 0, 'Unknown', 'false');
else
    % TODO It is possible we could deal with logical without difficulty but
    % do we need to?
    error(['Cannot handle argument of class: ', class(Q), ...
           ' in private octonion function UNDEFINED.'])
end

end

% $Id: undefined.m 1125 2021-08-15 20:43:38Z sangwine $

% Created automatically from the quaternion
% function of the same name on 15-Aug-2021.