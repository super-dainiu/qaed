function [Q, D] = skew_iqrq_cpp(H, rtol, verbose, progress_bar) %#ok<INUSD>
% SKEW_IQRQ_CPP Drop-in replacement for SKEW_IQRQ backed by the C++ core.
% Same signature as SKEW_IQRQ; VERBOSE and PROGRESS_BAR are ignored.
narginchk(1, 4), nargoutchk(0, 2)
if nargin < 2, rtol = eps; end

[w, x, y, z] = qmat_parts(H);
[Qw, Qx, Qy, Qz, Dw, Dx, Dy, Dz] = qaed_mex('skew_iqr', w, x, y, z, rtol);
Q = quaternion(Qw, Qx, Qy, Qz);
D = quaternion(diag(Dw), diag(Dx), diag(Dy), diag(Dz));
end
