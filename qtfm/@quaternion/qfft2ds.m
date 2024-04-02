function Y = qfft2ds(X, A)
% QFFT2DS Fast Quaternion 2D Two-Sided Fourier transform.
%
% This function calculates a fast quaternion 2D two-sided separable Fourier
% transform of the two dimensional quaternion array X. The second parameter
% specifies the left and right transform axes (the direction in 3-space of
% the vector part of the left and right hypercomplex exponentials). It must
% be a cell array of two pure quaternions (real or complex), not
% necessarily unit quaternions (the code here uses unit versions of the
% supplied values). The transform computed here is a generalisation of the
% 2D QFFT of Todd Ell (1992) (see below for references).

% Copyright © 2022 Stephen J. Sangwine and Nicolas Le Bihan.
% See the file : Copyright.m for further details.

narginchk(2, 2), nargoutchk(0, 1)

if ~iscell(A)
    error('The second parameter must be a cell array.')
end
      
if length(A) ~= 2
    error('The second parameter must have 2 elements.')
end

L = A{1};
R = A{2};
    
if ~isscalar(L) || ~isscalar(R)
    error('The transform axes must be scalars (not matrices or vectors).');
end

if ~isa(L, 'quaternion') || ~isa(R, 'quaternion')   
    error('The transform axes must be quaternions.')
end

if ~isempty(L.w) || ~isempty(R.w)     
    error('The transform axes must be pure quaternions.')
end

f = unit(L); % Now we know that L and R are pure quaternions, we can safely
g = unit(R); % normalise them, for use below.

[XP, XM] = opd(X, f, g); % Decompose the input array X using the orthogonal
                         % plane decomposition.

% We use Theorem 4.2 in arxiv:1306.2157 (see below). The negation of the
% first index value requires that we perform an FFTFLIP on the first
% decomposed component before computing the FFT2. This reverses/negates the
% row indices (i.e. it flips the array up and down the columns).

Y = qfft2(fftflip(XP), g, 'R') + qfft2(XM, g, 'R');

end

% References
% 
% The function which is computed above is a simple generalization of the 2D
% transform of Todd Ell using generalized roots of -1.  The function FFT2DS
% computes Todd's original transform by providing j and k as the roots of
% -1 as in Todd's thesis:
%
% T. A. Ell, Hypercomplex spectral transformations, Ph.D. thesis,
% University of Minnesota, 1992.
%
% Ell, T. A. and Le Bihan, N. and Sangwine, S. J., 'Quaternion Fourier
% Transforms for Signal and Image Processing', ISTE-Wiley,
% ISBN 978-1-84821-478-1, May 2014. doi:10.1002/9781118930908. See §3.2.
%
% The method by which the transform is here implemented (an orthogonal
% plane decomposition followed by two one-sided transforms) is given in:
%
% Eckhard Hitzer and Stephen J. Sangwine, The Orthogonal 2D Planes Split of
% Quaternions and Steerable Quaternion Fourier Transformations,
% arxiv:1306.2157, 10 June 2013. Lemma 3.4 is crucial to the algorithm.
%
% The following also sets out relationships between one and two-sided QFFTs
% and thus provides a route to implementation as used here.
%
% Min-Hung Yeh, 'Relationships Among Various 2-D Quaternion Fourier
% Transforms' IEEE Signal Processing Letters, 15, 669-672, (2008).
% DOI: 10.1109/LSP.2008.2002714.

% $Id: qfft2ds.m 1163 2022-11-16 12:41:20Z sangwine $
