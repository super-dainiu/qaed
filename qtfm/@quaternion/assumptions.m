function t = assumptions(q)
% ASSUMPTIONS  List assumptions that apply to components of quaternion.
% (Quaternion overloading of standard Matlab function.)

% Copyright © 2020 Stephen J. Sangwine and Nicolas Le Bihan.
% See the file : Copyright.m for further details.

narginchk(1, 1), nargoutchk(0, 1)

if ~isa(q, 'quaternion') || ~isa(q.x, 'sym')
    error('First argument must be a quaternion with symbolic components.')
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

% $Id: assumptions.m 1113 2021-02-01 18:41:09Z sangwine $
