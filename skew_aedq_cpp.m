function [Q, D] = skew_aedq_cpp(H, rtol, verbose, progress_bar) %#ok<INUSD>
% SKEW_AEDQ_CPP Drop-in replacement for SKEW_AEDQ backed by the C++ core.
% Same signature as SKEW_AEDQ; VERBOSE and PROGRESS_BAR are ignored.
narginchk(1, 4), nargoutchk(0, 2)
if nargin < 2, rtol = eps; end

[w, x, y, z] = qmat_parts(H);
[Qw, Qx, Qy, Qz, Dw, Dx, Dy, Dz] = qaed_mex('skew_aed', w, x, y, z, rtol);
Q = quaternion(Qw, Qx, Qy, Qz);
D = quaternion(diag(Dw), diag(Dx), diag(Dy), diag(Dz));
end
