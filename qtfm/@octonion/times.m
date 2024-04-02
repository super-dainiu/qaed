function o = times(l, r)
% .*  Array multiply.
% (Octonion overloading of standard Matlab function.)
 
% Copyright © 2011 Stephen J. Sangwine and Nicolas Le Bihan.
% See the file : Copyright.m for further details.

ol = isa(l, 'octonion');
or = isa(r, 'octonion');

if ol && or
    
    % Both arguments are octonions.
    
    % Use the Cayley-Dickson formula for the octonion product. This is
    % given in Ward (1997), or the Wikipedia page on octonions as:
    % (a,b)(c,d) = [(ac - d*b), (da + bc*)] where the * denotes a quaternion
    % conjugate.

    o = l; % Create a result, avoiding the use of the constructor.
    o.a = l.a .* r.a - conj(r.b) .* l.b;
    o.b = r.b .* l.a + l.b .* conj(r.a);
else

    % One of the arguments is not an octonion. If it is numeric, logical,
    % or symbolic, we can handle it.
    
    if ol && (isa(r, 'numeric') || isa(r, 'logical') || isa(r, 'sym'))
        o = l; % The left operand is an octonion, so use it for the result
               % to avoid calling the constructor.     
        o.a = o.a .* r; o.b = o.b .* r;
    elseif (isa(l, 'numeric') || isa(l, 'logical') || isa(l, 'sym')) && or
        o = r; % The right operand is an octonion, so use it for the
               % result to avoid calling the constructor.  
        o.a = l .* o.a; o.b = l .* o.b;

    else
        if isa(l, 'quaternion') || isa(r, 'quaternion')
            error('Cannot multiply quaternions and octonions')
        else
            if isa(l, 'octonion')
                error(['Cannot multiply an octonion with a ', class(r)]);
            else
                error(['Cannot multiply a ', class(l), ' with an octonion']);
            end
        end
    end
end

% if isa(x(o.a), 'sym')
%    o = elide(o); 
% end

end

% $Id: times.m 1127 2021-08-30 19:35:41Z sangwine $
