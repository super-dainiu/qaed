function Y = iqfft2ds(X, Axes)
% IQFFT2DS Inverse Fast Quaternion 2D Two-Sided Fourier transform.

% Copyright © 2022 Stephen J. Sangwine and Nicolas Le Bihan.
% See the file : Copyright.m for further details.

narginchk(2, 2), nargoutchk(0, 1)

if ~iscell(Axes)
    error('The second parameter must be a cell array.')
end
      
if length(Axes) ~= 2
    error('The second parameter must have 2 elements.')
end

% We compute the inverse here by negating the two elements of Axes. This
% trick doesn't apply a scale factor, so we have to do that here.

Y = qfft2ds(X, cellfun(@uminus, Axes, 'UniformOutput', false)) ./ numel(X);

end

% References: see the function QFFT2DS which computes the matching forward
% transform.

% $Id: iqfft2ds.m 1167 2022-11-20 12:31:11Z sangwine $
