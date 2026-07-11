function [P, D] = eigvec_cpp(Q, T)
% EIGVEC_CPP Drop-in replacement for EIGVEC backed by the C++ core (qaed_mex).
% Same signature and semantics as EIGVEC.
narginchk(2, 2), nargoutchk(0, 2)

[Qw, Qx, Qy, Qz] = qmat_parts(Q);
[Tw, Tx, Ty, Tz] = qmat_parts(T);
[Pw, Px, Py, Pz, Dw, Dx, Dy, Dz] = qaed_mex('eigvec', Qw, Qx, Qy, Qz, Tw, Tx, Ty, Tz);
P = quaternion(Pw, Px, Py, Pz);
D = quaternion(diag(Dw), diag(Dx), diag(Dy), diag(Dz));
end
