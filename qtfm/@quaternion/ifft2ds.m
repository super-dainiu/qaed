function Y = ifft2ds(X)
% IFFT2DS Inverse Fast Quaternion Two-Sided Two-Dimensional Fourier transform.

% Copyright © 2022 Stephen J. Sangwine and Nicolas Le Bihan.
% See the file : Copyright.m for further details.

narginchk(1, 1), nargoutchk(0, 1) 

Y = iqfft2ds(X, {qj, qk});

end

% References: see the forward transform function FFT2DS.

% $Id: ifft2ds.m 1162 2022-11-15 21:34:09Z sangwine $
