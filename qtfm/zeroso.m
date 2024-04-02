function E = zeroso(varargin)
% ZEROSO   N-by-N octonion matrix of zeros. Takes the same parameters as
% the Matlab function ZEROS (q.v.).

% Copyright © 2011 Stephen J. Sangwine and Nicolas Le Bihan.
% See the file : Copyright.m for further details.

E = octonion(zeros(varargin{:}));

end

% $Id: zeroso.m 1113 2021-02-01 18:41:09Z sangwine $
