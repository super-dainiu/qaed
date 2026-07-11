// qmat.hpp - Dense column-major quaternion matrix + kernels.
//
// All algorithm-facing helpers use 1-BASED INCLUSIVE index ranges so the code
// mirrors the MATLAB reference implementation line by line.
#pragma once
#include "quat.hpp"
#include <vector>
#include <cassert>
#include <complex>
#include <stdexcept>

// Inside a MEX file we link MATLAB's own BLAS/LAPACK (libmwblas/libmwlapack),
// which use 64-bit integers; standalone builds link LP64 FlexiBLAS.
#ifdef QAED_ILP64
using blas_int = long long;
#else
using blas_int = int;
#endif

extern "C" {
void zgemm_(const char* transa, const char* transb, const blas_int* m, const blas_int* n,
            const blas_int* k, const std::complex<double>* alpha,
            const std::complex<double>* a, const blas_int* lda,
            const std::complex<double>* b, const blas_int* ldb,
            const std::complex<double>* beta, std::complex<double>* c,
            const blas_int* ldc);
void zgeev_(const char* jobvl, const char* jobvr, const blas_int* n,
            std::complex<double>* a, const blas_int* lda, std::complex<double>* w,
            std::complex<double>* vl, const blas_int* ldvl, std::complex<double>* vr,
            const blas_int* ldvr, std::complex<double>* work, const blas_int* lwork,
            double* rwork, blas_int* info);
}

namespace qaed {

class QMat {
public:
    QMat() : m_(0), n_(0) {}
    QMat(int m, int n) : m_(m), n_(n), d_(static_cast<size_t>(m) * n) {}

    static QMat eye(int n) {
        QMat I(n, n);
        for (int i = 1; i <= n; ++i) I(i, i) = Quat::one();
        return I;
    }

    int rows() const { return m_; }
    int cols() const { return n_; }

    // 1-based access
    Quat& operator()(int i, int j) { return d_[static_cast<size_t>(j - 1) * m_ + (i - 1)]; }
    const Quat& operator()(int i, int j) const { return d_[static_cast<size_t>(j - 1) * m_ + (i - 1)]; }

    Quat* data() { return d_.data(); }
    const Quat* data() const { return d_.data(); }

    double norm_fro() const {
        double s = 0;
        for (const Quat& q : d_) s += q.norm2();
        return std::sqrt(s);
    }

    QMat block(int r1, int r2, int c1, int c2) const {  // inclusive
        QMat B(r2 - r1 + 1, c2 - c1 + 1);
        for (int j = c1; j <= c2; ++j)
            for (int i = r1; i <= r2; ++i)
                B(i - r1 + 1, j - c1 + 1) = (*this)(i, j);
        return B;
    }
    void set_block(int r1, int c1, const QMat& B) {
        for (int j = 1; j <= B.cols(); ++j)
            for (int i = 1; i <= B.rows(); ++i)
                (*this)(r1 + i - 1, c1 + j - 1) = B(i, j);
    }

    QMat ctranspose() const {
        QMat B(n_, m_);
        for (int j = 1; j <= n_; ++j)
            for (int i = 1; i <= m_; ++i)
                B(j, i) = (*this)(i, j).conj();
        return B;
    }

    void triu(int k = 0) {  // keep upper triangle (diag k), zero the rest
        for (int j = 1; j <= n_; ++j)
            for (int i = 1; i <= m_; ++i)
                if (j - i < k) (*this)(i, j) = Quat::zero();
    }
    void tril(int k = 0) {
        for (int j = 1; j <= n_; ++j)
            for (int i = 1; i <= m_; ++i)
                if (j - i > k) (*this)(i, j) = Quat::zero();
    }

private:
    int m_, n_;
    std::vector<Quat> d_;
};

// ---------------------------------------------------------------------------
// Quaternion GEMM: C = A * B via 4 complex GEMMs on the split A = A1 + A2*j.
// (A1+A2 j)(B1+B2 j) = (A1 B1 - A2 conj(B2)) + (A1 B2 + A2 conj(B1)) j
// ---------------------------------------------------------------------------
inline QMat qmul(const QMat& A, const QMat& B) {
    const blas_int m = A.rows(), k = A.cols(), n = B.cols();
    assert(B.rows() == static_cast<int>(k));
    using cd = std::complex<double>;
    std::vector<cd> A1(static_cast<size_t>(m) * k), A2(A1.size());
    std::vector<cd> B1(static_cast<size_t>(k) * n), B2(B1.size()), B2c(B1.size()), B1c(B1.size());
    std::vector<cd> C1(static_cast<size_t>(m) * n), C2(C1.size());

    const Quat* a = A.data();
    for (size_t t = 0; t < A1.size(); ++t) { A1[t] = a[t].c1(); A2[t] = a[t].c2(); }
    const Quat* b = B.data();
    for (size_t t = 0; t < B1.size(); ++t) {
        B1[t] = b[t].c1(); B2[t] = b[t].c2();
        B1c[t] = std::conj(B1[t]); B2c[t] = std::conj(B2[t]);
    }

    const cd one(1, 0), zero(0, 0), neg(-1, 0);
    const char N = 'N';
    // C1 = A1*B1 - A2*conj(B2)
    zgemm_(&N, &N, &m, &n, &k, &one, A1.data(), &m, B1.data(), &k, &zero, C1.data(), &m);
    zgemm_(&N, &N, &m, &n, &k, &neg, A2.data(), &m, B2c.data(), &k, &one, C1.data(), &m);
    // C2 = A1*B2 + A2*conj(B1)
    zgemm_(&N, &N, &m, &n, &k, &one, A1.data(), &m, B2.data(), &k, &zero, C2.data(), &m);
    zgemm_(&N, &N, &m, &n, &k, &one, A2.data(), &m, B1c.data(), &k, &one, C2.data(), &m);

    QMat C(static_cast<int>(m), static_cast<int>(n));
    Quat* c = C.data();
    for (size_t t = 0; t < C1.size(); ++t) c[t] = from_c(C1[t], C2[t]);
    return C;
}

// In-place block updates: X(r1:r2, c1:c2) = X(r1:r2, c1:c2) * U  and  U' * X.
// Small U is applied with direct scalar loops (zgemm overhead dominates there).
inline void right_mul(QMat& X, int r1, int r2, int c1, int c2, const QMat& U) {
    if (r2 < r1 || c2 < c1) return;
    const int k = U.rows();
    if (k <= 16) {
        std::vector<Quat> t(k);
        for (int i = r1; i <= r2; ++i) {
            for (int j = 1; j <= k; ++j) {
                Quat s = Quat::zero();
                for (int l = 1; l <= k; ++l) s += X(i, c1 + l - 1) * U(l, j);
                t[j - 1] = s;
            }
            for (int j = 1; j <= k; ++j) X(i, c1 + j - 1) = t[j - 1];
        }
        return;
    }
    X.set_block(r1, c1, qmul(X.block(r1, r2, c1, c2), U));
}
inline void left_mul_ct(QMat& X, int r1, int r2, int c1, int c2, const QMat& U) {
    if (r2 < r1 || c2 < c1) return;
    const int k = U.rows();
    if (k <= 16) {
        std::vector<Quat> t(k);
        for (int j = c1; j <= c2; ++j) {
            for (int i = 1; i <= k; ++i) {
                Quat s = Quat::zero();
                for (int l = 1; l <= k; ++l) s += U(l, i).conj() * X(r1 + l - 1, j);
                t[i - 1] = s;
            }
            for (int i = 1; i <= k; ++i) X(r1 + i - 1, j) = t[i - 1];
        }
        return;
    }
    X.set_block(r1, c1, qmul(U.ctranspose(), X.block(r1, r2, c1, c2)));
}

// Scale column j (rows all) by unit quaternion u on the right; row i by conj(u) on the left.
inline void col_scale(QMat& X, int j, const Quat& u) {
    for (int i = 1; i <= X.rows(); ++i) X(i, j) = X(i, j) * u;
}
inline void row_scale_ct(QMat& X, int i, const Quat& u) {
    for (int j = 1; j <= X.cols(); ++j) X(i, j) = u.conj() * X(i, j);
}

// ---------------------------------------------------------------------------
// Householder: qtfm-compatible householder_vector with v = e1 (column case).
// Returns u (norm sqrt(2)) and zeta such that (I - u u') a = e1 * (zeta*alpha).
// ---------------------------------------------------------------------------
inline void householder_vector(const std::vector<Quat>& a, std::vector<Quat>& u, Quat& zeta) {
    const size_t n = a.size();
    u.assign(n, Quat::zero());
    double alpha2 = 0;
    for (const Quat& q : a) alpha2 += q.norm2();
    double alpha = std::sqrt(alpha2);
    if (alpha == 0.0) { zeta = Quat::one(); return; }

    Quat romega = a[0];           // a.' * e1 (no conjugation)
    double r = romega.abs();
    zeta = (r != 0.0) ? -(romega / r) : Quat::one();
    double mu = std::sqrt(alpha * (alpha + r));

    u[0] = (a[0] - zeta * alpha) / mu;
    for (size_t i = 1; i < n; ++i) u[i] = a[i] / mu;
}

// x = x - u * (u' * x) on X(r1:r2, c1:c2); u has length r2-r1+1.
inline void householder_lapply(QMat& X, const std::vector<Quat>& u, int r1, int r2, int c1, int c2) {
    const int len = r2 - r1 + 1;
    assert(static_cast<int>(u.size()) == len);
#pragma omp parallel for schedule(static) if ((c2 - c1) > 64)
    for (int j = c1; j <= c2; ++j) {
        Quat s = Quat::zero();
        for (int i = 0; i < len; ++i) s += u[i].conj() * X(r1 + i, j);
        for (int i = 0; i < len; ++i) X(r1 + i, j) -= u[i] * s;
    }
}

// x = x - (x * u) * u' on X(r1:r2, c1:c2); u has length c2-c1+1.
inline void householder_rapply(QMat& X, const std::vector<Quat>& u, int r1, int r2, int c1, int c2) {
    const int len = c2 - c1 + 1;
    assert(static_cast<int>(u.size()) == len);
    const int nr = r2 - r1 + 1;
    std::vector<Quat> v(nr, Quat::zero());
    // v = X * u  (column-major friendly: accumulate column by column)
    for (int j = 0; j < len; ++j) {
        const Quat& uj = u[j];
        if (uj.norm2() == 0.0) continue;
        for (int i = 0; i < nr; ++i) v[i] += X(r1 + i, c1 + j) * uj;
    }
    // X -= v * u'
#pragma omp parallel for schedule(static) if (len > 64)
    for (int j = 0; j < len; ++j) {
        Quat ujc = u[j].conj();
        if (ujc.norm2() == 0.0) continue;
        for (int i = 0; i < nr; ++i) X(r1 + i, c1 + j) -= v[i] * ujc;
    }
}

// ---------------------------------------------------------------------------
// schur2: ordered Schur decomposition of a 2x2 quaternion matrix.
// Q' * A * Q = T upper triangular, diagonal standardized (complex, Im >= 0),
// ordered by decreasing |lambda| (as MATLAB ordschur+groupByValue).
// ---------------------------------------------------------------------------
Quat swapq(const QMat& T, QMat& Qout);  // fwd decl (defined below)

inline void schur2(const QMat& A, QMat& Q, QMat& T) {
    assert(A.rows() == 2 && A.cols() == 2);
    using cd = std::complex<double>;
    // 4x4 complex adjoint M = [[A1, A2], [-conj(A2), conj(A1)]]
    cd M[16];  // column-major 4x4
    auto put = [&](int i, int j, cd v) { M[j * 4 + i] = v; };  // 0-based here
    for (int j = 0; j < 2; ++j)
        for (int i = 0; i < 2; ++i) {
            const Quat& q = A(i + 1, j + 1);
            put(i,     j,     q.c1());
            put(i,     j + 2, q.c2());
            put(i + 2, j,     -std::conj(q.c2()));
            put(i + 2, j + 2, std::conj(q.c1()));
        }
    blas_int n4 = 4, lwork = 64, info = 0;
    cd w[4], vr[16], work[64];
    double rwork[8];
    const char Nc = 'N', V = 'V';
    zgeev_(&Nc, &V, &n4, M, &n4, w, nullptr, &n4, vr, &n4, work, &lwork, rwork, &info);
    if (info != 0) throw std::runtime_error("zgeev failed in schur2");

    // Pick the eigenvalue with the largest modulus among those with Im >= 0
    // (eigenvalues come in conjugate pairs).
    int best = -1;
    double bestabs = -1;
    for (int t = 0; t < 4; ++t) {
        if (w[t].imag() < 0) continue;
        if (std::abs(w[t]) > bestabs) { bestabs = std::abs(w[t]); best = t; }
    }
    if (best < 0) {  // numerical safety: all Im < 0, take largest modulus
        for (int t = 0; t < 4; ++t)
            if (std::abs(w[t]) > bestabs) { bestabs = std::abs(w[t]); best = t; }
    }

    // Quaternion eigenvector v = top - conj(bottom)*j, A v = v lambda.
    std::vector<Quat> v(2);
    for (int i = 0; i < 2; ++i)
        v[i] = from_c(vr[best * 4 + i], -std::conj(vr[best * 4 + i + 2]));
    double vn = std::sqrt(v[0].norm2() + v[1].norm2());
    v[0] = v[0] / vn; v[1] = v[1] / vn;

    // Unitary Q with first column = v: Q = (I - u u') * diag(zeta, 1).
    std::vector<Quat> u;
    Quat zeta;
    householder_vector(v, u, zeta);
    Q = QMat::eye(2);
    householder_rapply(Q, u, 1, 2, 1, 2);   // Q = I - u u' (Hermitian)
    Q(1, 1) = Q(1, 1) * zeta;
    Q(2, 1) = Q(2, 1) * zeta;

    T = qmul(qmul(Q.ctranspose(), A), Q);
    T(2, 1) = Quat::zero();

    // Standardize the diagonal.
    Quat U1, t1, U2, t2;
    schur1(T(1, 1), U1, t1);
    T(1, 1) = t1;
    T(1, 2) = U1.conj() * T(1, 2);
    col_scale(Q, 1, U1);
    schur1(T(2, 2), U2, t2);
    T(2, 2) = t2;
    T(1, 2) = T(1, 2) * U2;
    col_scale(Q, 2, U2);

    // Order by decreasing modulus (MATLAB ordschur puts the largest first).
    if (T(2, 2).abs() > T(1, 1).abs() * (1.0 + 1e-15)) {
        QMat U;
        swapq(T, U);
        T = qmul(qmul(U.ctranspose(), T), U);
        T.triu();
        Q = qmul(Q, U);
    }
}

// swapq: unitary U such that U' * T * U swaps the diagonal of upper
// triangular 2x2 T (with complex standardized diagonal).
inline Quat swapq(const QMat& T, QMat& U) {
    Quat x = sylvesterc(T(1, 1), T(2, 2), -T(1, 2));
    double zn = std::sqrt(1.0 + x.norm2());
    Quat c = x / zn;
    double s = 1.0 / zn;
    U = QMat(2, 2);
    U(1, 1) = c;              U(1, 2) = Quat(-s);
    U(2, 1) = Quat(s);        U(2, 2) = c.conj();
    return x;
}

// shift: Wilkinson-like shift from trailing 2x2 block (mirrors MATLAB shift()).
inline Quat shift2(const QMat& A) {
    Quat u1, T1;
    schur1(A(2, 2), u1, T1);
    QMat Q2, T2;
    schur2(A, Q2, T2);
    if ((T2(1, 1) - T1).abs() < (T2(2, 2) - T1).abs())
        return T2(1, 1);
    return T2(2, 2);
}

// ordschurq: move T(n,n) to T(1,1) by a sequence of adjacent swaps.
// Returns unitary Q; T is updated in place.
inline QMat ordschurq(QMat& T) {
    const int n = T.rows();
    QMat Q = QMat::eye(n);
    for (int ihi = n; ihi > 1; --ihi) {
        int ilo = ihi - 1;
        QMat U;
        swapq(T.block(ilo, ihi, ilo, ihi), U);
        left_mul_ct(T, ilo, ihi, 1, T.cols(), U);
        right_mul(T, 1, T.rows(), ilo, ihi, U);
        right_mul(Q, 1, n, ilo, ihi, U);
    }
    return Q;
}

} // namespace qaed
