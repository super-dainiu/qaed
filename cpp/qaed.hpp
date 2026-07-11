// qaed.hpp - C++ port of the qaed algorithms (hessq / iqrq / aedq / eigvec /
// eigq). Faithful, line-by-line translation of the MATLAB reference using
// 1-based inclusive indexing throughout.
#pragma once
#include "qmat.hpp"
#include <algorithm>
#include <cstdio>

namespace qaed {

// ---------------------------------------------------------------------------
// hessq: Q' * A * Q = H upper Hessenberg.
// ---------------------------------------------------------------------------
inline void hessq(QMat A, QMat& Q, QMat& H) {
    const int n = A.rows();
    Q = QMat::eye(n);
    std::vector<Quat> a, u;
    Quat zeta;
    for (int i = 1; i <= n - 2; ++i) {
        a.assign(n - i, Quat::zero());
        for (int k = i + 1; k <= n; ++k) a[k - i - 1] = A(k, i);
        householder_vector(a, u, zeta);
        householder_lapply(A, u, i + 1, n, i, n);
        householder_rapply(A, u, 1, n, i + 1, n);
        householder_rapply(Q, u, 1, n, i + 1, n);
    }
    A.triu(-1);
    H = std::move(A);
}

// Standardize diagonal entry (i,i): H(i,i) <- U' H(i,i) U, update row/col of H
// and column i of Q. Mirrors the schur1-deflation idiom in the MATLAB code.
inline void standardize_diag(QMat& H, QMat& Q, int i) {
    Quat U, t;
    schur1(H(i, i), U, t);
    col_scale(Q, i, U);
    col_scale(H, i, U);
    row_scale_ct(H, i, U);
}

// ---------------------------------------------------------------------------
// iqrq: implicit quaternion QR on a Hessenberg matrix. Q' * H0 * Q = T.
// ---------------------------------------------------------------------------
struct IqrStats { long steps = 0; };

inline void iqrq(QMat H, double rtol, QMat& Q, QMat& T, IqrStats* stats = nullptr) {
    const int n = H.rows();
    int ilo = 1, ihi = n;
    const double atol = rtol * H.norm_fro();
    Q = QMat::eye(n);
    std::vector<Quat> a, u;
    Quat zeta;
    long GS = 0;

    while (ihi >= 1) {
        while (ihi - ilo > 1) {
            Quat x = shift2(H.block(ihi - 1, ihi, ihi - 1, ihi));
            ++GS;
            double si = -2.0 * x.w;
            double ti = x.norm2();
            a.assign(3, Quat::zero());
            a[0] = H(ilo, ilo) * H(ilo, ilo) + H(ilo, ilo + 1) * H(ilo + 1, ilo)
                 + si * H(ilo, ilo) + Quat(ti);
            a[1] = H(ilo + 1, ilo) * H(ilo, ilo) + H(ilo + 1, ilo + 1) * H(ilo + 1, ilo)
                 + si * H(ilo + 1, ilo);
            a[2] = H(ilo + 2, ilo + 1) * H(ilo + 1, ilo);

            householder_vector(a, u, zeta);
            householder_lapply(H, u, ilo, ilo + 2, ilo, n);
            householder_rapply(H, u, 1, ihi, ilo, ilo + 2);
            householder_rapply(Q, u, 1, n, ilo, ilo + 2);

            for (int i = ilo; i <= ihi - 2; ++i) {
                int e = std::min(i + 3, ihi);
                int sp = std::max(i - 1, ilo);

                a.assign(e - i, Quat::zero());
                for (int k = i + 1; k <= e; ++k) a[k - i - 1] = H(k, i);
                householder_vector(a, u, zeta);

                householder_lapply(H, u, i + 1, e, sp, n);
                householder_rapply(H, u, 1, ihi, i + 1, e);
                for (int k = i + 2; k <= e; ++k) H(k, i) = Quat::zero();
                householder_rapply(Q, u, 1, n, i + 1, e);

                double sub = H(i + 1, i).abs();
                if (sub <= atol && sub <= rtol * (H(i, i).abs() + H(i + 1, i + 1).abs())) {
                    standardize_diag(H, Q, i);
                    H(i + 1, i) = Quat::zero();
                    ilo = i + 1;
                }
            }

            while (H(ihi, ihi - 1).abs() <= rtol * (H(ihi, ihi).abs() + H(ihi - 1, ihi - 1).abs())) {
                standardize_diag(H, Q, ihi);
                H(ihi, ihi - 1) = Quat::zero();
                --ihi;
                if (ihi <= ilo) break;
            }
        }
        if (ihi > ilo) {
            QMat U, T2;
            schur2(H.block(ilo, ihi, ilo, ihi), U, T2);
            right_mul(Q, 1, n, ilo, ihi, U);
            right_mul(H, 1, n, ilo, ihi, U);
            left_mul_ct(H, ilo, ihi, ilo, n, U);
            H(ihi, ilo) = Quat::zero();
        } else if (ihi == ilo) {
            standardize_diag(H, Q, ilo);
        }
        ihi = ilo - 1;
        ilo = 1;
    }

    H.triu();
    T = std::move(H);
    if (stats) stats->steps = GS;
}

// ---------------------------------------------------------------------------
// aedq: implicit QR with aggressive early deflation. Q' * H0 * Q = T.
// ---------------------------------------------------------------------------
namespace detail {
inline int aed_num_shifts(int ihi, double alpha) {
    return std::max(2, static_cast<int>(std::lround(alpha * ihi)));
}
inline int aed_win_size(int ihi, int NS) {
    int WS = (ihi <= 500) ? NS : (3 * NS) / 2;
    return std::max(4, WS - WS % 2);
}
constexpr int    aed_min_size = 12;
constexpr double nibble       = 0.14;

// AED step on the window H (whole matrix passed is the window itself).
// Returns unitary Q (window-sized) and the shift list; H is updated in place.
inline void aed_step(QMat& H, double rtol, QMat& Q, std::vector<Quat>& shifts) {
    const int n = H.rows();
    int ihi = n;

    {   // [Q1, H(2:n,2:n)] = iqrq(H(2:n,2:n))
        QMat Q1, T1;
        iqrq(H.block(2, n, 2, n), rtol, Q1, T1);
        H.set_block(2, 2, T1);
        right_mul(H, 1, 1, 2, n, Q1);      // H(1,2:n) = H(1,2:n) * Q1
        left_mul_ct(H, 2, n, 1, 1, Q1);    // H(2:n,1) = Q1' * H(2:n,1)
        Q = QMat::eye(n);
        Q.set_block(2, 2, Q1);
    }

    for (int i = 2; i <= n; ++i) {
        if (H(ihi, 1).abs() <= rtol * H(ihi, ihi).abs()) {
            H(ihi, 1) = Quat::zero();
            standardize_diag(H, Q, ihi);
            --ihi;
        } else {
            QMat Tw = H.block(2, ihi, 2, ihi);
            QMat U = ordschurq(Tw);
            right_mul(Q, 1, n, 2, ihi, U);
            right_mul(H, 1, n, 2, ihi, U);
            left_mul_ct(H, 2, ihi, 1, n, U);
            Tw = H.block(2, ihi, 2, ihi);
            Tw.triu();
            H.set_block(2, 2, Tw);
        }
    }

    shifts.clear();
    for (int i = 2; i <= ihi; ++i) shifts.push_back(H(i, i));

    // Rebuild Hessenberg form of the spike block.
    QMat U, Hh;
    hessq(H.block(1, ihi, 1, ihi), U, Hh);
    H.set_block(1, 1, Hh);
    left_mul_ct(H, 1, ihi, ihi + 1, n, U);
    right_mul(Q, 1, n, 1, ihi, U);
}
} // namespace detail

struct AedStats { long steps = 0; long aed_deflated = 0; };

inline void aedq(QMat H, double rtol, QMat& Q, QMat& T, double alpha = 0.25,
                 AedStats* stats = nullptr) {
    const int n = H.rows();
    int ilo = 1, ihi = n;
    const double atol = rtol * H.norm_fro();
    Q = QMat::eye(n);
    std::vector<Quat> a, u, shifts;
    Quat zeta;
    long GS = 0, DA = 0;

    while (ihi >= 1) {
        while (ihi - ilo > 1) {
            int NS = detail::aed_num_shifts(ihi - ilo + 1, alpha);
            int WS = std::min(detail::aed_win_size(ihi - ilo + 1, NS), ihi - ilo);
            int sp = std::max(ihi - WS, ilo);

            if (ihi - ilo + 1 > detail::aed_min_size && sp > ilo) {
                QMat Hw = H.block(sp, ihi, sp, ihi);
                QMat U;
                detail::aed_step(Hw, rtol, U, shifts);
                H.set_block(sp, sp, Hw);
                left_mul_ct(H, sp, ihi, ihi + 1, n, U);
                right_mul(H, 1, sp - 1, sp, ihi, U);
                right_mul(Q, 1, n, sp, ihi, U);
                DA += ihi - sp - static_cast<long>(shifts.size());
                ihi = sp + static_cast<int>(shifts.size());
            } else {
                shifts.assign(1, shift2(H.block(ihi - 1, ihi, ihi - 1, ihi)));
            }

            if (static_cast<double>(ihi - sp + 1) / WS < (1.0 - detail::nibble))
                continue;  // sufficient deflation

            int LS = 0;
            NS = std::min<int>(static_cast<int>(shifts.size()), NS);

            while (ihi > ilo + 1 && LS < NS) {
                ++GS;
                ++LS;
                Quat x = shifts[LS - 1];
                double si = -2.0 * x.w;
                double ti = x.norm2();
                a.assign(3, Quat::zero());
                a[0] = H(ilo, ilo) * H(ilo, ilo) + H(ilo, ilo + 1) * H(ilo + 1, ilo)
                     + si * H(ilo, ilo) + Quat(ti);
                a[1] = H(ilo + 1, ilo) * H(ilo, ilo) + H(ilo + 1, ilo + 1) * H(ilo + 1, ilo)
                     + si * H(ilo + 1, ilo);
                a[2] = H(ilo + 2, ilo + 1) * H(ilo + 1, ilo);

                householder_vector(a, u, zeta);
                householder_lapply(H, u, ilo, ilo + 2, ilo, n);
                householder_rapply(H, u, 1, ihi, ilo, ilo + 2);
                householder_rapply(Q, u, 1, n, ilo, ilo + 2);

                for (int i = ilo; i <= ihi - 2; ++i) {
                    int e = std::min(i + 3, ihi);
                    int s2 = std::max(i - 1, ilo);

                    a.assign(e - i, Quat::zero());
                    for (int k = i + 1; k <= e; ++k) a[k - i - 1] = H(k, i);
                    householder_vector(a, u, zeta);

                    householder_lapply(H, u, i + 1, e, s2, n);
                    householder_rapply(H, u, 1, ihi, i + 1, e);
                    for (int k = i + 2; k <= e; ++k) H(k, i) = Quat::zero();
                    householder_rapply(Q, u, 1, n, i + 1, e);

                    double sub = H(i + 1, i).abs();
                    if (sub <= atol && sub <= rtol * (H(i, i).abs() + H(i + 1, i + 1).abs())) {
                        standardize_diag(H, Q, i);
                        H(i + 1, i) = Quat::zero();
                        ilo = i + 1;
                    }
                }

                while (H(ihi, ihi - 1).abs() <= rtol * (H(ihi, ihi).abs() + H(ihi - 1, ihi - 1).abs())) {
                    standardize_diag(H, Q, ihi);
                    H(ihi, ihi - 1) = Quat::zero();
                    --ihi;
                    if (ihi <= ilo) break;
                }
            }
        }
        if (ihi > ilo) {
            QMat U, T2;
            schur2(H.block(ilo, ihi, ilo, ihi), U, T2);
            right_mul(Q, 1, n, ilo, ihi, U);
            right_mul(H, 1, n, ilo, ihi, U);
            left_mul_ct(H, ilo, ihi, ilo, n, U);
            H(ihi, ilo) = Quat::zero();
        } else if (ihi == ilo) {
            standardize_diag(H, Q, ilo);
        }
        ihi = ilo - 1;
        ilo = 1;
    }

    H.triu();
    T = std::move(H);
    if (stats) { stats->steps = GS; stats->aed_deflated = DA; }
}

// ---------------------------------------------------------------------------
// sylvesterc_tri: solve T X - X b = c with T upper triangular (complex diag).
// ---------------------------------------------------------------------------
inline std::vector<Quat> sylvesterc_tri(const QMat& T, const Quat& b, std::vector<Quat> c) {
    const int n = T.rows();
    std::vector<Quat> X(n, Quat::zero());
    for (int i = n; i >= 1; --i) {
        X[i - 1] = sylvesterc(T(i, i), b, c[i - 1]);
        for (int k = 1; k <= i - 1; ++k) c[k - 1] -= T(k, i) * X[i - 1];
    }
    return X;
}

// ---------------------------------------------------------------------------
// eigvec: eigenvectors from the Schur decomposition Q * T * Q'.
// ---------------------------------------------------------------------------
inline void eigvec(const QMat& Q, const QMat& T, QMat& P, std::vector<Quat>& D) {
    const int n = T.rows();
    P = Q;
    std::vector<Quat> c;
    for (int i = n; i >= 2; --i) {
        c.assign(i - 1, Quat::zero());
        for (int k = 1; k <= i - 1; ++k) c[k - 1] = T(k, i);
        std::vector<Quat> x = sylvesterc_tri(T.block(1, i - 1, 1, i - 1), T(i, i), c);
        // P(:,i) -= P(:,1:i-1) * x ; then normalize
#pragma omp parallel for schedule(static) if (n > 128)
        for (int r = 1; r <= n; ++r) {
            Quat s = Quat::zero();
            for (int k = 1; k <= i - 1; ++k) s += P(r, k) * x[k - 1];
            P(r, i) -= s;
        }
        double nn = 0;
        for (int r = 1; r <= n; ++r) nn += P(r, i).norm2();
        nn = std::sqrt(nn);
        for (int r = 1; r <= n; ++r) P(r, i) = P(r, i) / nn;
    }
    D.assign(n, Quat::zero());
    for (int i = 1; i <= n; ++i) D[i - 1] = T(i, i);
}

// ---------------------------------------------------------------------------
// eigq: full eigen-decomposition A = P * diag(D) * P^{-1}.
// (General path; the skew/tridiagonal specialization falls back to it.)
// ---------------------------------------------------------------------------
inline void eigq(const QMat& A, double rtol, QMat& P, std::vector<Quat>& D) {
    QMat Qh, H;
    hessq(A, Qh, H);
    QMat Q, T;
    if (H.rows() < 1000) iqrq(std::move(H), rtol, Q, T);
    else                 aedq(std::move(H), rtol, Q, T);
    eigvec(qmul(Qh, Q), T, P, D);
}

} // namespace qaed
