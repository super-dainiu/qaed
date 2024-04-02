function Y = qdft2ds(X, A)
% QDFT2DS Quaternion 2D Two-Sided Discrete Fourier Transform.
%
% This function calculates a quaternion 2D two-sided separable discrete
% Fourier transform of the two dimensional quaternion array X. The second
% parameter specifies the left and right transform axes (the direction in
% 3-space of the vector part of the left and right hypercomplex
% exponentials). It must be a cell array of two pure quaternions (real or
% complex), not necessarily unit quaternions (the code here uses unit
% versions of the supplied values). The transform computed here is a
% generalisation of the 2D QDFT of Todd Ell (1992) (see below for
% references). This is not a fast implementation - it is provided as a
% check on the FFT version.

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

% The code here directly follows equation 4 in the 1996 Electronics Letters
% paper cited below, except that j and k are replaced by f and g, which are
% arbitrary unit pure quaternions, and there is no scale factor of 1/√(MN).

[M, N] = size(X);

Y = zerosq(M, N);

for u = 1:M % Loop indices are 1 greater than in the equation in the paper
    for v = 1:N % in order to index correctly into the MATLAB arrays. This
                % means they must be decremented in the expressions for
                % the exponentials.
        S = substruct('()', {u , v}); % See below inside the inner loops.
        for m = 1:M
            for n = 1:N
                % To index into X and Y here we need substruct/subsref
                % because this is a class method.
                % See Implementation_notes.txt, section on Indexing in
                % Class Methods.
                % Y(u,v) = Y(u,v) + ...
                Y = subsasgn(Y, S, subsref(Y, S) + ...
                         exp(-f .* 2 .* pi .* (m - 1) .* (u - 1) ./ M) ...
                         ... % .* X(m, n) .* is implemented on the next line.
                         .* subsref(X, substruct('()', {m, n})) .* ...
                         exp(-g .* 2 .* pi .* (n - 1) .* (v - 1) ./ N) ...
                         );
            end
        end
    end
end

end

% References
% 
% The function which is computed above is a simple generalization of the 2D
% transform of Todd Ell using generalized roots of -1.
%
% T. A. Ell, Hypercomplex spectral transformations, Ph.D. thesis,
% University of Minnesota, 1992.
%
% Sangwine, S. J., 'Fourier transforms of colour images using quaternion,
% or hypercomplex, numbers , Electronics Letters, 32(21), October 10 1996,
% 1979–80. DOI: 10.1049/el:19961331.

% $Id: qdft2ds.m 1167 2022-11-20 12:31:11Z sangwine $
