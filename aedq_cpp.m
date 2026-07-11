function [Q, T] = aedq_cpp(H, rtol, verbose, progress_bar, alpha) %#ok<INUSD>
% AEDQ_CPP Drop-in replacement for AEDQ backed by the C++ core (qaed_mex).
% Same signature as AEDQ; VERBOSE and PROGRESS_BAR are accepted but ignored.
narginchk(1, 5), nargoutchk(0, 2)
if nargin < 5, alpha = 0.25; end
if nargin < 2, rtol = eps; end

if ~ismatrix(H), error('A must be a square matrix'); end
if ~isa(H, 'quaternion'), error('A must be a quaternion matrix'); end
if ~all(abs(tril(H, -2)) < eps, 'all'), error('A must be an upper Hessenberg matrix'); end

[w, x, y, z] = qmat_parts(H);
[Qw, Qx, Qy, Qz, Tw, Tx, Ty, Tz] = qaed_mex('aed', w, x, y, z, rtol, alpha);
Q = quaternion(Qw, Qx, Qy, Qz);
T = quaternion(Tw, Tx, Ty, Tz);
end
