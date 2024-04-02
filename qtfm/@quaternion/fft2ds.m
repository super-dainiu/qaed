function Y = fft2ds(X)
% FFT2DS Fast Quaternion Two-Sided Two-Dimensional Fourier transform.

% Copyright © 2022 Stephen J. Sangwine and Nicolas Le Bihan.
% See the file : Copyright.m for further details.

narginchk(1, 1), nargoutchk(0, 1) 

Y = qfft2ds(X, {qj, qk});

end

% References: the transform computed by this function is that of Ell (1992,
% 1993) in discrete form as published by Sangwine (1996).
%
% T. A. Ell, 'Hypercomplex Spectral Transformations'. PhD thesis,
% University of Minnesota, June 1992.
%
% T. A. Ell, 'Quaternion-Fourier transforms for analysis of 2-dimensional
% linear time-invariant partial-differential systems'. In Proceedings of
% the 32nd Conference on Decision and Control, pages 1830–1841, San
% Antonio, Texas, USA, 15–17 December 1993. IEEE Control Systems Society.
%
% S. J. Sangwine, 'Fourier transforms of colour images using quaternion, or
% hypercomplex, numbers', Electronics Letters, 32(21):1979–1980, 10 October
% 1996. DOI: 10.1049/el:19961331.

% $Id: fft2ds.m 1162 2022-11-15 21:34:09Z sangwine $
