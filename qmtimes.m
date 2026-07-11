function C = qmtimes(A, B)
% QMTIMES Quaternion matrix product via the C++ quaternion BLAS (qaed_mex).
%
% C = QMTIMES(A, B) computes A*B for qtfm quaternion matrices using the
% native quaternion GEMM kernel (no mapping to complex/real BLAS).
% Falls back to qtfm's mtimes when qaed_mex is not built.
narginchk(2, 2)

if exist('qaed_mex', 'file') ~= 3
    C = A * B;
    return
end

[Aw, Ax, Ay, Az] = qmat_parts(A);
[Bw, Bx, By, Bz] = qmat_parts(B);
[Cw, Cx, Cy, Cz] = qaed_mex('gemm', Aw, Ax, Ay, Az, Bw, Bx, By, Bz);
C = quaternion(Cw, Cx, Cy, Cz);
end
