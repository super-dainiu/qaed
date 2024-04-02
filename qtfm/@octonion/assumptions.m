function t = assumptions(q)
% ASSUMPTIONS  List assumptions that apply to components of octonion.
% (Octonion overloading of standard Matlab function.)

% Copyright © 2020 Stephen J. Sangwine and Nicolas Le Bihan.
% See the file : Copyright.m for further details.

narginchk(1, 1), nargoutchk(0, 1)

if ~isa(q, 'octonion') || ~isa(q.x, 'sym')
    error('First argument must be an octonion with symbolic components.')
end

switch nargout
    case 0
        if ~isempty(q.w)
            assumptions(q.w)
        end
        
        assumptions(q.x)
        assumptions(q.y)
        assumptions(q.z)
    case 1
        t = [assumptions(q.x), assumptions(q.y), assumptions(q.z)];
        if ~isempty(q.w)
            % The unique function removes duplicated entries in t. These
            % can occur if an expression applies to coefficients of q, for
            % example, abs(q) == 1.
            t = unique(cat(2, assumptions(q.w), t));
        end
end

end

% $Id: assumptions.m 1122 2021-08-13 19:09:07Z sangwine $

% Created automatically from the quaternion
% function of the same name on 13-Aug-2021.