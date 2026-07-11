// qaed_mex.cpp - MEX gateway to the C++ qaed core.
//
// Usage (from the thin MATLAB wrappers hessq_cpp / iqrq_cpp / aedq_cpp /
// eigq_cpp; not intended to be called directly):
//
//   [Qw,Qx,Qy,Qz, Tw,Tx,Ty,Tz] = qaed_mex(mode, Hw,Hx,Hy,Hz, rtol, alpha)
//
//   mode: 'hess' -> Q, H      (Hessenberg reduction)
//         'iqr'  -> Q, T      (implicit QR Schur form)
//         'aed'  -> Q, T      (QR with aggressive early deflation)
//         'eig'  -> P, D      (eigenvectors; D returned as n-by-1 diagonal)
//
// Build (inside MATLAB, in the cpp/ directory):
//   mex -R2018a -DQAED_ILP64 CXXFLAGS='$CXXFLAGS -O3 -std=c++17 -fopenmp' ...
//       LDFLAGS='$LDFLAGS -fopenmp' qaed_mex.cpp -lmwlapack -lmwblas
#include "qaed.hpp"
#include "mex.h"
#include <cstring>
#include <string>

using namespace qaed;

static QMat unpack(const mxArray* w, const mxArray* x, const mxArray* y,
                   const mxArray* z) {
    const mwSize m = mxGetM(x), n = mxGetN(x);
    QMat A(static_cast<int>(m), static_cast<int>(n));
    const size_t sz = static_cast<size_t>(m) * n;
    auto fill = [&](const mxArray* a, double Quat::*p) {
        if (a == nullptr || mxIsEmpty(a)) return;  // pure quaternion: w = []
        if (mxGetM(a) != m || mxGetN(a) != n || !mxIsDouble(a) || mxIsComplex(a))
            mexErrMsgIdAndTxt("qaed:input", "component arrays must be real double and same size");
        const double* d = mxGetPr(a);
        for (size_t t = 0; t < sz; ++t) A.data()[t].*p = d[t];
    };
    fill(w, &Quat::w); fill(x, &Quat::x); fill(y, &Quat::y); fill(z, &Quat::z);
    return A;
}

static void pack(mxArray* plhs[], int base, const QMat& A) {
    const mwSize m = A.rows(), n = A.cols();
    const size_t sz = static_cast<size_t>(m) * n;
    double Quat::* comps[4] = {&Quat::w, &Quat::x, &Quat::y, &Quat::z};
    for (int c = 0; c < 4; ++c) {
        plhs[base + c] = mxCreateDoubleMatrix(m, n, mxREAL);
        double* d = mxGetPr(plhs[base + c]);
        for (size_t t = 0; t < sz; ++t) d[t] = A.data()[t].*comps[c];
    }
}

void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]) {
    if (nrhs < 5 || !mxIsChar(prhs[0]))
        mexErrMsgIdAndTxt("qaed:usage",
            "usage: [Qw..Qz,Tw..Tz] = qaed_mex(mode, Hw,Hx,Hy,Hz [,rtol,alpha])");
    char modebuf[16];
    mxGetString(prhs[0], modebuf, sizeof modebuf);
    std::string mode(modebuf);

    if (mode == "gemm") {
        // [Cw..Cz] = qaed_mex('gemm', Aw..Az, Bw..Bz) - quaternion BLAS GEMM
        if (nrhs < 9)
            mexErrMsgIdAndTxt("qaed:usage", "gemm mode needs 8 component arrays");
        const QMat A = unpack(prhs[1], prhs[2], prhs[3], prhs[4]);
        const QMat B = unpack(prhs[5], prhs[6], prhs[7], prhs[8]);
        if (A.cols() != B.rows())
            mexErrMsgIdAndTxt("qaed:input", "inner dimensions must agree");
        pack(plhs, 0, qgemm(A, B));
        return;
    }

    if (mode == "eigvec") {
        // [Pw..Pz, Dw..Dz] = qaed_mex('eigvec', Qw..Qz, Tw..Tz)
        if (nrhs < 9)
            mexErrMsgIdAndTxt("qaed:usage", "eigvec mode needs 8 component arrays");
        const QMat Q = unpack(prhs[1], prhs[2], prhs[3], prhs[4]);
        const QMat T = unpack(prhs[5], prhs[6], prhs[7], prhs[8]);
        if (Q.rows() != T.rows() || T.rows() != T.cols())
            mexErrMsgIdAndTxt("qaed:input", "Q and T must be square and same size");
        try {
            QMat P;
            std::vector<Quat> D;
            eigvec(Q, T, P, D);
            pack(plhs, 0, P);
            if (nlhs > 4) {
                QMat Dm(static_cast<int>(D.size()), 1);
                for (int i = 1; i <= Dm.rows(); ++i) Dm(i, 1) = D[i - 1];
                pack(plhs, 4, Dm);
            }
        } catch (const std::exception& e) {
            mexErrMsgIdAndTxt("qaed:internal", "%s", e.what());
        }
        return;
    }

    const QMat H = unpack(prhs[1], prhs[2], prhs[3], prhs[4]);
    if (H.rows() != H.cols())
        mexErrMsgIdAndTxt("qaed:input", "matrix must be square");
    double rtol  = (nrhs > 5 && !mxIsEmpty(prhs[5])) ? mxGetScalar(prhs[5]) : 2.220446049250313e-16;
    double alpha = (nrhs > 6 && !mxIsEmpty(prhs[6])) ? mxGetScalar(prhs[6]) : 0.25;
    if (nlhs > 8) mexErrMsgIdAndTxt("qaed:usage", "at most 8 outputs");

    try {
        QMat Q, T;
        if (mode == "skew_iqr" || mode == "skew_aed") {
            // [Qw..Qz, Dw..Dz] = qaed_mex('skew_iqr'|'skew_aed', Hw..Hz, rtol)
            std::vector<Quat> D;
            if (mode == "skew_iqr") skew_iqrq(H, rtol, Q, D);
            else                    skew_aedq(H, rtol, Q, D);
            T = QMat(static_cast<int>(D.size()), 1);
            for (int i = 1; i <= T.rows(); ++i) T(i, 1) = D[i - 1];
            pack(plhs, 0, Q);
            if (nlhs > 4) pack(plhs, 4, T);
            return;
        }
        if (mode == "hess") {
            hessq(H, Q, T);
        } else if (mode == "iqr") {
            iqrq(H, rtol, Q, T);
        } else if (mode == "aed") {
            aedq(H, rtol, Q, T, alpha);
        } else if (mode == "eig") {
            std::vector<Quat> D;
            eigq(H, rtol, Q, D);
            T = QMat(static_cast<int>(D.size()), 1);
            for (int i = 1; i <= T.rows(); ++i) T(i, 1) = D[i - 1];
        } else {
            mexErrMsgIdAndTxt("qaed:mode", "mode must be hess|iqr|aed|eig");
        }
        pack(plhs, 0, Q);
        if (nlhs > 4) pack(plhs, 4, T);
    } catch (const std::exception& e) {
        mexErrMsgIdAndTxt("qaed:internal", "%s", e.what());
    }
}
