function b = subsref(a, ss)
% SUBSREF Subscripted reference.
% (Octonion overloading of standard Matlab function.)

% Copyright © 2011 Stephen J. Sangwine and Nicolas Le Bihan.
% See the file : Copyright.m for further details.

if length(ss) ~= 1
    error('Only one level of subscripting is currently supported for octonions.');
end

switch ss.type
case '()'
    if length(ss) ~= 1
        error('Multiple levels of subscripting are not supported for octonions.')
    end

    b = overload(mfilename, a, ss);
case '{}'
    error('Cell array indexing is not valid for octonions.')
case '.'
    error('Structure indexing is not valid for octonions.')
otherwise
    error('subsref received an invalid subscripting type.')
end

end

% $Id: subsref.m 1113 2021-02-01 18:41:09Z sangwine $
