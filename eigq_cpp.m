function [P, D] = eigq_cpp(A, rtol, verbose, progress_bar) %#ok<INUSD>
% EIGQ_CPP Drop-in replacement for EIGQ backed by the C++ core (qaed_mex).
% Same signature as EIGQ; VERBOSE and PROGRESS_BAR are accepted but ignored.
% Note: uses the general (Hessenberg + QR) path for all inputs; the
% skew/tridiagonal specialization is not yet ported.
narginchk(1, 4), nargoutchk(0, 2)
if nargin < 2, rtol = eps; end

[w, x, y, z] = qmat_parts(A);
[Pw, Px, Py, Pz, Dw, Dx, Dy, Dz] = qaed_mex('eig', w, x, y, z, rtol);
P = quaternion(Pw, Px, Py, Pz);
D = quaternion(diag(Dw), diag(Dx), diag(Dy), diag(Dz));

if nargout == 1, P = diag(D); end
end
