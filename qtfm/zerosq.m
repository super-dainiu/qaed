function E = zerosq(varargin)
% ZEROSQ   N-by-N quaternion matrix of zeros. Takes the same parameters as
% the Matlab function ZEROS (q.v.).

% Copyright © 2008 Stephen J. Sangwine and Nicolas Le Bihan.
% See the file : Copyright.m for further details.

E = quaternion(zeros(varargin{:}));

end

% $Id: zerosq.m 1113 2021-02-01 18:41:09Z sangwine $
