function tf = isempty(q)
% ISEMPTY True for empty matrix.
% (Octonion overloading of standard Matlab function.)

% Copyright © 2011 Stephen J. Sangwine and Nicolas Le Bihan.
% See the file : Copyright.m for further details.
     
tf = isempty(q.a);
if tf
    assert(isempty(q.b));
end

end

% $Id: isempty.m 1113 2021-02-01 18:41:09Z sangwine $
