// qaed.hpp - C++ port of the qaed algorithms (hessq / iqrq / aedq / eigvec /
// eigq). Faithful, line-by-line translation of the MATLAB reference using
// 1-based inclusive indexing throughout.
#pragma once
#include "qmat.hpp"
#include <algorithm>
#include <chrono>
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

// Cumulative timers matching the MATLAB verbose breakdown (time spent
// updating the Schur vectors Q, and time spent inside the AED branch).
struct QaedTimers {
    double q_time = 0, aed_time = 0;
    void reset() { q_time = aed_time = 0; }
};
inline QaedTimers& qaed_timers() { static QaedTimers t; return t; }
inline double qaed_now() {
    return std::chrono::duration<double>(
        std::chrono::steady_clock::now().time_since_epoch()).count();
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
// qr_sweep: one implicit-shift QR sweep on the active block [ilo, ihi] with
// shift class x (standardized complex eigenvalue). Mirrors the reference
// sweep (bulge introduction + chase with immediate updates); shared by iqrq
// and aedq. Updates ilo in place on mid-sweep deflations.
// ---------------------------------------------------------------------------
inline void qr_sweep(QMat& H, QMat& Q, const Quat& x, int& ilo, int ihi,
                     double atol, double rtol) {
    const int n = H.cols();
    const double si = -2.0 * x.w;
    const double ti = x.norm2();
    std::vector<Quat> a, u;
    Quat zeta;

    a.assign(3, Quat::zero());
    a[0] = H(ilo, ilo) * H(ilo, ilo) + H(ilo, ilo + 1) * H(ilo + 1, ilo)
         + si * H(ilo, ilo) + Quat(ti);
    a[1] = H(ilo + 1, ilo) * H(ilo, ilo) + H(ilo + 1, ilo + 1) * H(ilo + 1, ilo)
         + si * H(ilo + 1, ilo);
    a[2] = H(ilo + 2, ilo + 1) * H(ilo + 1, ilo);

    householder_vector(a, u, zeta);
    householder_lapply(H, u, ilo, ilo + 2, ilo, n);
    householder_rapply(H, u, 1, ihi, ilo, ilo + 2);
    double t0 = qaed_now();
    householder_rapply(Q, u, 1, Q.rows(), ilo, ilo + 2);
    qaed_timers().q_time += qaed_now() - t0;

    for (int i = ilo; i <= ihi - 2; ++i) {
        const int e = std::min(i + 3, ihi);
        const int sp = std::max(i - 1, ilo);

        a.assign(e - i, Quat::zero());
        for (int k = i + 1; k <= e; ++k) a[k - i - 1] = H(k, i);
        householder_vector(a, u, zeta);

        householder_lapply(H, u, i + 1, e, sp, n);
        householder_rapply(H, u, 1, ihi, i + 1, e);
        for (int k = i + 2; k <= e; ++k) H(k, i) = Quat::zero();
        t0 = qaed_now();
        householder_rapply(Q, u, 1, Q.rows(), i + 1, e);
        qaed_timers().q_time += qaed_now() - t0;

        double sub = H(i + 1, i).abs();
        if (sub <= atol && sub <= rtol * (H(i, i).abs() + H(i + 1, i + 1).abs())) {
            standardize_diag(H, Q, i);
            H(i + 1, i) = Quat::zero();
            ilo = i + 1;
        }
    }
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
    long GS = 0;

    while (ihi >= 1) {
        while (ihi - ilo > 1) {
            Quat x = shift2(H.block(ihi - 1, ihi, ihi - 1, ihi));
            ++GS;
            qr_sweep(H, Q, x, ilo, ihi, atol, rtol);

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
    std::vector<Quat> shifts;
    long GS = 0, DA = 0;

    while (ihi >= 1) {
        while (ihi - ilo > 1) {
            int NS = detail::aed_num_shifts(ihi - ilo + 1, alpha);
            int WS = std::min(detail::aed_win_size(ihi - ilo + 1, NS), ihi - ilo);
            int sp = std::max(ihi - WS, ilo);

            if (ihi - ilo + 1 > detail::aed_min_size && sp > ilo) {
                double t0 = qaed_now();
                QMat Hw = H.block(sp, ihi, sp, ihi);
                QMat U;
                detail::aed_step(Hw, rtol, U, shifts);
                H.set_block(sp, sp, Hw);
                left_mul_ct(H, sp, ihi, ihi + 1, n, U);
                right_mul(H, 1, sp - 1, sp, ihi, U);
                qaed_timers().aed_time += qaed_now() - t0;
                t0 = qaed_now();
                right_mul(Q, 1, n, sp, ihi, U);
                qaed_timers().q_time += qaed_now() - t0;
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
                qr_sweep(H, Q, shifts[LS - 1], ilo, ihi, atol, rtol);

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
// Skew-Hermitian tridiagonal specialization (H' = -H). Faithful port of
// skew_iqrq.m / skew_aedq.m: band-limited sweeps, deflation criteria and the
// AED-with-sorted-spike step, including the fixes of 2026-07 (spike/eigen-
// value pairing, trailing-only deflation).
// ---------------------------------------------------------------------------

// One band-limited implicit-shift sweep on the tridiagonal block [ilo, ihi].
inline void skew_qr_sweep(QMat& H, QMat& Q, const Quat& x, int& ilo, int ihi,
                          double rtol) {
    const double si = -2.0 * x.w;
    const double ti = x.norm2();
    std::vector<Quat> a, u;
    Quat zeta;

    a.assign(3, Quat::zero());
    a[0] = H(ilo, ilo) * H(ilo, ilo) + H(ilo, ilo + 1) * H(ilo + 1, ilo)
         + si * H(ilo, ilo) + Quat(ti);
    a[1] = H(ilo + 1, ilo) * H(ilo, ilo) + H(ilo + 1, ilo + 1) * H(ilo + 1, ilo)
         + si * H(ilo + 1, ilo);
    a[2] = H(ilo + 2, ilo + 1) * H(ilo + 1, ilo);
    householder_vector(a, u, zeta);

    householder_rapply(H, u, ilo, ilo + 3, ilo, ilo + 2);
    householder_lapply(H, u, ilo, ilo + 2, ilo, ilo + 3);
    double t0 = qaed_now();
    householder_rapply(Q, u, 1, Q.rows(), ilo, ilo + 2);
    qaed_timers().q_time += qaed_now() - t0;

    for (int i = ilo; i <= ihi - 2; ++i) {
        const int e  = std::min(i + 3, ihi);
        const int ee = std::min(e + 1, ihi);

        a.assign(e - i, Quat::zero());
        for (int k = i + 1; k <= e; ++k) a[k - i - 1] = H(k, i);
        householder_vector(a, u, zeta);

        householder_lapply(H, u, i + 1, e, i, ee);
        householder_rapply(H, u, i, ee, i + 1, e);
        // Trim the wake (kept banded by skew symmetry up to rounding).
        for (int k = i + 2; k <= e; ++k) { H(k, i) = Quat::zero(); H(i, k) = Quat::zero(); }
        t0 = qaed_now();
        householder_rapply(Q, u, 1, Q.rows(), i + 1, e);
        qaed_timers().q_time += qaed_now() - t0;

        if (H(i + 1, i).abs() <= rtol * (H(i, i).abs() + H(i + 1, i + 1).abs())) {
            standardize_diag(H, Q, i);
            H(i + 1, i) = Quat::zero();
            H(i, i + 1) = Quat::zero();
            ilo = i + 1;
        }
    }
    // Trim rounding fill in the trailing corner.
    for (int c = std::max(1, ihi - 4); c <= ihi; ++c)
        for (int r = std::max(1, ihi - 4); r <= ihi; ++r)
            if (r - c > 1 || c - r > 1) { H(r, c) = Quat::zero(); }
}

// Trailing deflation loop shared by skew_iqrq / skew_aedq.
inline void skew_deflate_trailing(QMat& H, QMat& Q, int ilo, int& ihi, double rtol) {
    while (ihi > ilo + 1 &&
           H(ihi, ihi - 1).abs() <= rtol * (H(ihi, ihi).abs() + H(ihi - 1, ihi - 1).abs())) {
        standardize_diag(H, Q, ihi);
        H(ihi, ihi - 1) = Quat::zero();
        H(ihi - 1, ihi) = Quat::zero();
        --ihi;
    }
}

struct SkewStats { long steps = 0; long aed_deflated = 0; };

// skew_iqrq: Q' * H0 * Q = D (diagonal, standardized) for skew-Hermitian
// tridiagonal H0.
inline void skew_iqrq(QMat H, double rtol, QMat& Q, std::vector<Quat>& D,
                      SkewStats* stats = nullptr) {
    const int n = H.rows();
    int ilo = 1, ihi = n;
    Q = QMat::eye(n);
    long GS = 0;

    while (ihi > 1) {
        while (ihi - ilo > 2) {
            Quat x = shift2(H.block(ihi - 1, ihi, ihi - 1, ihi));
            ++GS;
            skew_qr_sweep(H, Q, x, ilo, ihi, rtol);
            skew_deflate_trailing(H, Q, ilo, ihi, rtol);
        }
        {   // Finish the remaining (<= 3x3) block with the general solver.
            QMat U, T;
            iqrq(H.block(ilo, ihi, ilo, ihi), rtol, U, T);
            H.set_block(ilo, ilo, T);
            double t0 = qaed_now();
            right_mul(Q, 1, n, ilo, ihi, U);
            qaed_timers().q_time += qaed_now() - t0;
        }
        ihi = ilo - 1;
        ilo = 1;
    }

    D.assign(n, Quat::zero());
    for (int i = 1; i <= n; ++i) D[i - 1] = H(i, i);
    if (stats) stats->steps = GS;
}

namespace detail {
inline int skew_num_shifts(int m) {
    double NS;
    if (m < 30) NS = 2;
    else if (m < 60) NS = 4;
    else if (m < 150) NS = 10;
    else if (m < 590) NS = static_cast<double>(m) / std::lround(std::log2(static_cast<double>(m)));
    else if (m < 3000) NS = 64;
    else NS = 128;
    int k = static_cast<int>(NS);
    return std::max(2, k - k % 2);
}
} // namespace detail

// skew_aedq: skew_iqrq with aggressive early deflation.
inline void skew_aedq(QMat H, double rtol, QMat& Q, std::vector<Quat>& D,
                      SkewStats* stats = nullptr) {
    const int n = H.rows();
    int ilo = 1, ihi = n;
    Q = QMat::eye(n);
    long GS = 0, DA = 0;
    constexpr int MIN_SIZE = 12;
    constexpr double NIBBLE = 0.14;
    std::vector<Quat> shifts;

    while (ihi > 1) {
        while (ihi - ilo > 2) {
            int NS = detail::skew_num_shifts(ihi - ilo + 1);
            int WS = (ihi - ilo + 1 <= 500) ? NS : (3 * NS) / 2;
            WS = std::max(4, WS - WS % 2);
            const int win = std::min(WS, ihi - ilo + 1);
            int whi = ihi;
            const int wlo = std::max(ihi - win + 1, ilo);
            const int sp = wlo - 1;

            if (ihi - ilo + 1 >= MIN_SIZE && sp > ilo && win > 4) {
                // AED: Schur form of the window, sort the spike, deflate its
                // negligible trailing entries.
                double t0 = qaed_now();
                QMat U;
                std::vector<Quat> Dw;
                skew_iqrq(H.block(wlo, whi, wlo, whi), rtol, U, Dw);
                QMat spike = qmul(U.ctranspose(), H.block(wlo, whi, sp, sp));

                std::vector<int> idx(win);
                for (int i = 0; i < win; ++i) idx[i] = i;
                std::sort(idx.begin(), idx.end(), [&](int a2, int b2) {
                    return spike(a2 + 1, 1).abs() > spike(b2 + 1, 1).abs();
                });
                qaed_timers().aed_time += qaed_now() - t0;

                t0 = qaed_now();
                right_mul(Q, 1, n, wlo, whi, U);
                QMat Qw = Q.block(1, n, wlo, whi);   // permute window columns
                for (int j = 0; j < win; ++j)
                    for (int r = 1; r <= n; ++r)
                        Q(r, wlo + j) = Qw(r, idx[j] + 1);
                qaed_timers().q_time += qaed_now() - t0;

                t0 = qaed_now();
                std::vector<Quat> spikeS(win), Ds(win);
                for (int j = 0; j < win; ++j) {
                    spikeS[j] = spike(idx[j] + 1, 1);
                    Ds[j] = Dw[idx[j]];
                }
                // Window block <- diag(Ds); spike column/row rewritten below.
                for (int r = wlo; r <= whi; ++r)
                    for (int c = wlo; c <= whi; ++c)
                        H(r, c) = (r == c) ? Ds[r - wlo] : Quat::zero();

                for (int j = win - 1; j >= 0; --j) {
                    if (spikeS[j].abs() <= rtol * Ds[j].abs()) {
                        spikeS[j] = Quat::zero();
                        --whi;
                    } else {
                        break;
                    }
                }
                for (int j = 0; j < win; ++j) {
                    H(wlo + j, sp) = spikeS[j];
                    H(sp, wlo + j) = -spikeS[j].conj();
                }
                shifts.clear();
                for (int r = whi; r >= wlo; --r) shifts.push_back(H(r, r));
                DA += ihi - whi;
                ihi = whi;

                // Restore the tridiagonal form of the arrow block.
                QMat U2, Hh;
                hessq(H.block(sp, ihi, sp, ihi), U2, Hh);
                Hh.tril(1);   // skew: rounding fill above the band
                H.set_block(sp, sp, Hh);
                qaed_timers().aed_time += qaed_now() - t0;

                t0 = qaed_now();
                right_mul(Q, 1, n, sp, ihi, U2);
                qaed_timers().q_time += qaed_now() - t0;
            } else {
                shifts.assign(1, shift2(H.block(ihi - 1, ihi, ihi - 1, ihi)));
            }

            if (static_cast<double>(ihi - sp + 1) / win < (1.0 - NIBBLE))
                continue;

            int LS = 0;
            NS = std::min<int>(static_cast<int>(shifts.size()), NS);
            while (ihi > ilo + 1 && LS < NS) {
                ++GS;
                ++LS;
                skew_qr_sweep(H, Q, shifts[LS - 1], ilo, ihi, rtol);
            }
            skew_deflate_trailing(H, Q, ilo, ihi, rtol);
        }
        {
            QMat U, T;
            iqrq(H.block(ilo, ihi, ilo, ihi), rtol, U, T);
            H.set_block(ilo, ilo, T);
            double t0 = qaed_now();
            right_mul(Q, 1, n, ilo, ihi, U);
            qaed_timers().q_time += qaed_now() - t0;
        }
        ihi = ilo - 1;
        ilo = 1;
    }

    D.assign(n, Quat::zero());
    for (int i = 1; i <= n; ++i) D[i - 1] = H(i, i);
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
        // Rescale to avoid overflow when eigenvalues are clustered (the
        // equation is linear, so scaling X and c together preserves the
        // solution direction; callers normalize the result anyway).
        double ax = X[i - 1].abs();
        if (ax > 1e100) {
            double s = 1.0 / ax;
            for (int k = i - 1; k < n; ++k) X[k] = X[k] * s;
            for (int k = 0; k < i - 1; ++k) c[k] = c[k] * s;
        }
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
