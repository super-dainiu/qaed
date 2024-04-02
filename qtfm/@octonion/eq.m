function tf = eq(l, r)
% ==  Equal.
% (Octonion overloading of standard Matlab function.)

% Copyright © 2011 Stephen J. Sangwine and Nicolas Le Bihan.
% See the file : Copyright.m for further details.

narginchk(2, 2), nargoutchk(0, 1) 

if isa(l, 'octonion') && isa(r, 'octonion')
    
  tf = (l.a == r.a) & (l.b == r.b);
    
else
    % One of the arguments is not an octonion (the other must be, or Matlab
    % would not call this function). The non-octonion argument must be a
    % numeric or symbolic (if we don't impose this restriction it could be
    % anything such as a cell array or string which makes no sense at all
    % to compare with).
    
    if isa(l, 'octonion') && (isa(r, 'numeric') || isa(r, 'sym'))
        tf = (l.a == r) & (l.b == 0);
    elseif (isa(l, 'numeric') || isa(l, 'sym')) && isa(r, 'octonion')
        tf = (l == r.a) & (0 == r.b);
    else
        if isa(r, 'quaternion') || isa(l, 'quaternion')
            error('Cannot compare an octonion with a quaternion')
        else
            if isa(l, 'octonion')
                error(['Cannot compare an octonion with a ', class(r)]);
            else
                error(['Cannot compare a ', class(l), ' with an octonion']);
            end
        end
    end
end

end

% $Id: eq.m 1125 2021-08-15 20:43:38Z sangwine $
