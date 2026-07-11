function [Q, H] = hessq_cpp(A)
% HESSQ_CPP Drop-in replacement for HESSQ backed by the C++ core (qaed_mex).
% Same signature and semantics as HESSQ.
narginchk(1, 1), nargoutchk(0, 2)

% Mirror HESSQ: skip the reduction entirely if A is already Hessenberg.
if all(abs(tril(A, -2)) < eps, 'all')
    Q = eyeq(size(A, 1));
    H = triu(A, -1);
    if nargout < 2, Q = H; end
    return
end

[w, x, y, z] = qmat_parts(A);
[Qw, Qx, Qy, Qz, Hw, Hx, Hy, Hz] = qaed_mex('hess', w, x, y, z);
Q = quaternion(Qw, Qx, Qy, Qz);
H = quaternion(Hw, Hx, Hy, Hz);

if nargout < 2
    Q = H;
end
end
