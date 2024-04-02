function m = mean(X, dim)
% MEAN   Average or mean value.
% (Octonion overloading of standard Matlab function.)

% Copyright © 2015 Stephen J. Sangwine and Nicolas Le Bihan.
% See the file : Copyright.m for further details.

narginchk(1, 2), nargoutchk(0, 1) 

if nargin == 1
    m = overload(mfilename, X);
else
    if ~isnumeric(dim)
        error('Dimension argument must be numeric');
    end
    
    if ~isscalar(dim) || ~ismember(dim, 1:ndims(X))
        error(['Dimension argument must be a positive'...
               ' integer scalar within indexing range.']);
    end

    m = overload(mfilename, X, dim);
end

end

% $Id: mean.m 1113 2021-02-01 18:41:09Z sangwine $
