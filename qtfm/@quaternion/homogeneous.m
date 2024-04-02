function h = homogeneous(p)
% HOMOGENEOUS - function to construct or normalise a full quaternion in
% homogeneous coordinates. Given a pure quaternion argument, inserts unity
% scalar part to give a full quaternion result in homogeneous coordinates.
% Given a full quaternion argument, normalises so that the scalar part is
% unity by dividing by the scalar part, and returns the vector part alone.

% Copyright © 2006, 2016, 2017 Stephen J. Sangwine.
% See the file : Copyright.m for further details.

narginchk(1, 1), nargoutchk(0, 1)

if ispure(p)
    
    % Insert scalar part(s) of ones to convert the argument to homogeneous
    % coordinates.
    
    h = ones(size(p), class(x(p))) + p;
else
    
    % Normalise the argument to have unit scalar part(s) and return the
    % vector part only.
                   
    h = v(p) ./ s(p);
    
    if any(isinf(h(:)))
        warning('Some elements of the result have infinite values.')
    end

end

end

% $Id: homogeneous.m 1113 2021-02-01 18:41:09Z sangwine $
