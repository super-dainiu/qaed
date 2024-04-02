function varargout = arrayfun(varargin)
% ARRAYFUN Apply a function to each element of an array
% (Octonion overloading of standard Matlab function.)

% Copyright © 2020 Stephen J. Sangwine and Nicolas Le Bihan.
% See the file : Copyright.m for further details.

narginchk(2, 2), nargoutchk(0, 1)

F = varargin{1}; % We support only the calling syntax arrayfun(func, A) at
A = varargin{2}; % present.

if ~isscalar(F) || ~isa(F, 'function_handle')
    error('First parameter must be a function handle.')
end

% The second parameter must be an octonion, otherwise this function would
% not have been called. So we don't check this.

% The method below may not be the best way to implement arrayfun but it is
% better than no implementation at all. The MATLAB function won't work for
% octonions. The method below uses linear indexing: make a vector out of
% whatever array is provided, apply the function F one by one to the
% elements and then reshape the result back to the dimensions and size of
% the input argument.

B = subsref(A, substruct('()', {':'})); % B = A(:);

for j = 1:length(B)
    r = substruct('()', {j});
    B = subsasgn(B, r, F(subsref(B, r)));  % B(j) = F(B(j));
end

varargout{1} = reshape(B, size(A));

end

% $Id: arrayfun.m 1113 2021-02-01 18:41:09Z sangwine $
