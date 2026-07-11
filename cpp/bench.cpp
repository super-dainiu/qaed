// bench.cpp - benchmark + self-verification driver for the C++ qaed port.
//
// Usage:
//   ./qaed_bench --alg aed|iqr|eig|hess --n 256 [--seed 0] [--rtol 2.2e-16]
//                [--alpha 0.25]
//   ./qaed_bench --selftest
#include "qaed.hpp"
#include <chrono>
#include <cstring>
#include <random>
#include <string>

using namespace qaed;
using clk = std::chrono::steady_clock;

static double tictoc(clk::time_point t0) {
    return std::chrono::duration<double>(clk::now() - t0).count();
}

// A = randq(n) .* randn(n): uniform unit quaternion scaled by N(0,1).
static QMat rand_qmat(int n, unsigned seed) {
    std::mt19937_64 gen(seed);
    std::normal_distribution<double> N(0.0, 1.0);
    QMat A(n, n);
    for (int j = 1; j <= n; ++j)
        for (int i = 1; i <= n; ++i) {
            Quat q(N(gen), N(gen), N(gen), N(gen));
            double a = q.abs();
            if (a == 0) q = Quat::one(); else q = q / a;
            A(i, j) = q * N(gen);
        }
    return A;
}

// skewrand: A = randq(n).*rand(n); A = A - A'; reduce to tridiagonal form.
static QMat rand_skew_tridiag(int n, unsigned seed, double* t_hess = nullptr) {
    std::mt19937_64 gen(seed);
    std::normal_distribution<double> N(0.0, 1.0);
    std::uniform_real_distribution<double> U(0.0, 1.0);
    QMat A(n, n);
    for (int j = 1; j <= n; ++j)
        for (int i = 1; i <= n; ++i) {
            Quat q(N(gen), N(gen), N(gen), N(gen));
            double a = q.abs();
            if (a == 0) q = Quat::one(); else q = q / a;
            A(i, j) = q * U(gen);
        }
    for (int j = 1; j <= n; ++j)          // A <- A - A'
        for (int i = 1; i < j; ++i) {
            Quat v = A(i, j) - A(j, i).conj();
            A(i, j) = v;
            A(j, i) = -v.conj();
        }
    for (int i = 1; i <= n; ++i) {
        Quat d = A(i, i);
        A(i, i) = (d - d.conj()) * 0.5;   // pure vector part
    }
    QMat Qh, H;
    auto t0 = clk::now();
    hessq(A, Qh, H);
    if (t_hess) *t_hess = tictoc(t0);
    H.tril(1);                            // skew-Hermitian Hessenberg = tridiagonal
    return H;
}

// ---- verification helpers -------------------------------------------------
static double unitarity_err(const QMat& Q) {
    QMat E = qmul(Q.ctranspose(), Q);
    for (int i = 1; i <= E.rows(); ++i) E(i, i) -= Quat::one();
    return E.norm_fro() / std::sqrt(static_cast<double>(Q.rows()));
}

static double schur_residual(const QMat& A, const QMat& Q, const QMat& T) {
    QMat R = qmul(qmul(Q.ctranspose(), A), Q);
    for (int j = 1; j <= R.cols(); ++j)
        for (int i = 1; i <= R.rows(); ++i)
            R(i, j) -= T(i, j);
    return R.norm_fro() / A.norm_fro();
}

static double lower_mass(const QMat& T, int k) {  // ||tril(T,k)|| / ||T||
    double s = 0;
    for (int j = 1; j <= T.cols(); ++j)
        for (int i = 1; i <= T.rows(); ++i)
            if (j - i <= k) s += T(i, j).norm2();
    return std::sqrt(s) / T.norm_fro();
}

static double eig_residual(const QMat& A, const QMat& P, const std::vector<Quat>& D) {
    // || A*P - P*diag(D) || / (||A|| * ||P||)
    QMat AP = qmul(A, P);
    for (int j = 1; j <= P.cols(); ++j)
        for (int i = 1; i <= P.rows(); ++i)
            AP(i, j) -= P(i, j) * D[j - 1];
    return AP.norm_fro() / (A.norm_fro() * P.norm_fro());
}

// ---- self tests -----------------------------------------------------------
static int selftest() {
    int fail = 0;
    auto check = [&](const char* name, double err, double tol) {
        bool ok = err < tol;
        printf("  %-28s err = %.3e  %s\n", name, err, ok ? "ok" : "FAIL");
        if (!ok) ++fail;
    };

    {   // schur1
        std::mt19937_64 g(1);
        std::normal_distribution<double> N;
        double worst = 0;
        for (int t = 0; t < 1000; ++t) {
            Quat q(N(g), N(g), N(g), N(g)), u, s;
            schur1(q, u, s);
            worst = std::max(worst, std::abs(u.abs() - 1.0));
            Quat r = u.conj() * q * u;
            worst = std::max(worst, (r - s).abs());
            worst = std::max(worst, std::abs(r.y) + std::abs(r.z));
            if (s.x < -1e-14) worst = 1;
        }
        check("schur1 (1000 random)", worst, 1e-12);
    }
    {   // sylvesterc: a x - x b = c
        std::mt19937_64 g(2);
        std::normal_distribution<double> N;
        double worst = 0;
        for (int t = 0; t < 1000; ++t) {
            Quat a(N(g), N(g), 0, 0), b(N(g), N(g), 0, 0), c(N(g), N(g), N(g), N(g));
            Quat x = sylvesterc(a, b, c);
            worst = std::max(worst, (a * x - x * b - c).abs() / c.abs());
        }
        check("sylvesterc (1000 random)", worst, 1e-10);
    }
    {   // schur2
        std::mt19937_64 g(3);
        std::normal_distribution<double> N;
        double worst = 0;
        for (int t = 0; t < 500; ++t) {
            QMat A(2, 2);
            for (int j = 1; j <= 2; ++j)
                for (int i = 1; i <= 2; ++i) A(i, j) = Quat(N(g), N(g), N(g), N(g));
            QMat Q, T;
            schur2(A, Q, T);
            worst = std::max(worst, unitarity_err(Q));
            worst = std::max(worst, schur_residual(A, Q, T));
            worst = std::max(worst, T(2, 1).abs());
            if (T(1, 1).abs() + 1e-12 < T(2, 2).abs()) worst = 1;   // ordering
            if (std::abs(T(1, 1).y) + std::abs(T(1, 1).z) > 1e-12) worst = 1;
        }
        check("schur2 (500 random)", worst, 1e-10);
    }
    {   // swapq via ordschurq on 2x2
        std::mt19937_64 g(4);
        std::normal_distribution<double> N;
        double worst = 0;
        for (int t = 0; t < 200; ++t) {
            QMat A(2, 2);
            for (int j = 1; j <= 2; ++j)
                for (int i = 1; i <= 2; ++i) A(i, j) = Quat(N(g), N(g), N(g), N(g));
            QMat Q, T;
            schur2(A, Q, T);
            QMat T0 = T;
            QMat U = ordschurq(T);
            worst = std::max(worst, unitarity_err(U));
            worst = std::max(worst, schur_residual(T0, U, T));
            worst = std::max(worst, T(2, 1).abs() / T0.norm_fro());
        }
        check("swapq/ordschurq (200)", worst, 1e-10);
    }
    {   // hessq
        QMat A = rand_qmat(60, 5), Q, H;
        hessq(A, Q, H);
        check("hessq unitarity", unitarity_err(Q), 1e-12);
        check("hessq residual", schur_residual(A, Q, H), 1e-12);
        check("hessq lower mass", lower_mass(H, -2), 1e-13);
    }
    for (int n : {24, 60, 150}) {  // iqrq + aedq
        QMat A = rand_qmat(n, 100 + n), Qh, H;
        hessq(A, Qh, H);
        QMat Q, T;
        iqrq(H, 2.2e-16, Q, T);
        char buf[64];
        snprintf(buf, sizeof buf, "iqrq n=%d resid", n);
        check(buf, schur_residual(H, Q, T), 1e-11);
        snprintf(buf, sizeof buf, "iqrq n=%d unitary", n);
        check(buf, unitarity_err(Q), 1e-11);

        QMat Q2, T2;
        aedq(H, 2.2e-16, Q2, T2);
        snprintf(buf, sizeof buf, "aedq n=%d resid", n);
        check(buf, schur_residual(H, Q2, T2), 1e-11);
        snprintf(buf, sizeof buf, "aedq n=%d unitary", n);
        check(buf, unitarity_err(Q2), 1e-11);
    }
    {   // eigq end-to-end
        QMat A = rand_qmat(80, 42), P;
        std::vector<Quat> D;
        eigq(A, 2.2e-16, P, D);
        check("eigq n=80 residual", eig_residual(A, P, D), 1e-10);
    }
    for (int n : {24, 60, 150}) {   // skew tridiagonal specialization
        QMat H = rand_skew_tridiag(n, 300 + n);
        QMat Q1, Q2;
        std::vector<Quat> D1, D2;
        skew_iqrq(H, 2.2e-16, Q1, D1);
        skew_aedq(H, 2.2e-16, Q2, D2);
        char buf[64];
        auto diag_resid = [&](const QMat& Q, const std::vector<Quat>& D) {
            QMat R = qmul(qmul(Q.ctranspose(), H), Q);
            for (int i = 1; i <= n; ++i) R(i, i) -= D[i - 1];
            return R.norm_fro() / H.norm_fro();
        };
        snprintf(buf, sizeof buf, "skew_iqrq n=%d resid", n);
        check(buf, diag_resid(Q1, D1), 1e-11);
        snprintf(buf, sizeof buf, "skew_aedq n=%d resid", n);
        check(buf, diag_resid(Q2, D2), 1e-11);
        snprintf(buf, sizeof buf, "skew unitary n=%d", n);
        check(buf, std::max(unitarity_err(Q1), unitarity_err(Q2)), 1e-11);
        // eigenvalue multisets agree (standardized |v| values)
        std::vector<double> a1, a2;
        for (int i = 0; i < n; ++i) { a1.push_back(D1[i].abs()); a2.push_back(D2[i].abs()); }
        std::sort(a1.begin(), a1.end());
        std::sort(a2.begin(), a2.end());
        double worst = 0, scale = a1.back();
        for (int i = 0; i < n; ++i) worst = std::max(worst, std::abs(a1[i] - a2[i]));
        snprintf(buf, sizeof buf, "skew eig match n=%d", n);
        check(buf, worst / scale, 1e-11);
    }
    printf(fail ? "SELFTEST: %d FAILURES\n" : "SELFTEST: all passed\n", fail);
    return fail ? 1 : 0;
}

int main(int argc, char** argv) {
    std::string alg = "aed", type = "full";
    int n = 256;
    unsigned seed = 0;
    double rtol = 2.220446049250313e-16, alpha = 0.25;
    bool do_selftest = false;

    for (int i = 1; i < argc; ++i) {
        auto next = [&]() -> const char* { return (i + 1 < argc) ? argv[++i] : ""; };
        if (!strcmp(argv[i], "--alg")) alg = next();
        else if (!strcmp(argv[i], "--type")) type = next();
        else if (!strcmp(argv[i], "--n")) n = atoi(next());
        else if (!strcmp(argv[i], "--seed")) seed = static_cast<unsigned>(atoi(next()));
        else if (!strcmp(argv[i], "--rtol")) rtol = atof(next());
        else if (!strcmp(argv[i], "--alpha")) alpha = atof(next());
        else if (!strcmp(argv[i], "--selftest")) do_selftest = true;
        else { fprintf(stderr, "unknown arg: %s\n", argv[i]); return 2; }
    }
    if (do_selftest) return selftest();

    printf("alg = %s, type = %s, n = %d, rtol = %.3e\n", alg.c_str(), type.c_str(), n, rtol);

    if (alg == "skew_iqr" || alg == "skew_aed") {
        double t_hess = 0;
        QMat H = rand_skew_tridiag(n, seed, &t_hess);
        printf("tridiagonalization time: %.3f s\n", t_hess);
        QMat Q;
        std::vector<Quat> D;
        SkewStats st;
        qaed_timers().reset();
        auto t0 = clk::now();
        if (alg == "skew_iqr") skew_iqrq(H, rtol, Q, D, &st);
        else                   skew_aedq(H, rtol, Q, D, &st);
        double dt = tictoc(t0);
        printf("%s time: %.3f s (QR steps: %ld, AED deflated: %ld)\n",
               alg.c_str(), dt, st.steps, st.aed_deflated);
        printf("time to construct Q: %.3f s, time for AED: %.3f s\n",
               qaed_timers().q_time, qaed_timers().aed_time);
        // e1, e2, e3 (for a normal matrix the Schur vectors are eigenvectors)
        QMat R = qmul(qmul(Q.ctranspose(), H), Q);
        for (int i = 1; i <= n; ++i) R(i, i) -= D[i - 1];
        double dnorm = 0;
        for (const Quat& d : D) dnorm += d.norm2();
        dnorm = std::sqrt(dnorm);
        QMat HQ = qmul(H, Q);
        for (int j = 1; j <= n; ++j)
            for (int i = 1; i <= n; ++i)
                HQ(i, j) -= Q(i, j) * D[j - 1];
        printf("unitarity ||Q'Q - I||/sqrt(n)      = %.3e\n", unitarity_err(Q));
        printf("residual  ||Q'HQ - D||/||H||       = %.3e\n", R.norm_fro() / H.norm_fro());
        printf("eigvec    ||HQ - QD||/((|H|+|D|)|Q|) = %.3e\n",
               HQ.norm_fro() / ((H.norm_fro() + dnorm) * Q.norm_fro()));
        return 0;
    }

    QMat A = rand_qmat(n, seed);
    if (type == "hess")  // hessrand: directly generated upper Hessenberg
        A.triu(-1);

    if (alg == "eig") {
        QMat P;
        std::vector<Quat> D;
        auto t0 = clk::now();
        eigq(A, rtol, P, D);
        double dt = tictoc(t0);
        printf("eigq time: %.3f s\n", dt);
        printf("||A*P - P*D|| / (||A||*||P||) = %.3e\n", eig_residual(A, P, D));
        return 0;
    }

    QMat Qh, H;
    auto th = clk::now();
    hessq(A, Qh, H);
    double dth = tictoc(th);
    printf("hessq time: %.3f s (unitarity %.3e, residual %.3e)\n",
           dth, unitarity_err(Qh), schur_residual(A, Qh, H));
    if (alg == "hess") return 0;

    QMat Q, T;
    double dt = 0;
    qaed_timers().reset();
    if (alg == "iqr") {
        IqrStats st;
        auto t0 = clk::now();
        iqrq(H, rtol, Q, T, &st);
        dt = tictoc(t0);
        printf("iqrq time: %.3f s (QR steps: %ld)\n", dt, st.steps);
        printf("time to construct Q: %.3f s\n", qaed_timers().q_time);
    } else if (alg == "aed") {
        AedStats st;
        auto t0 = clk::now();
        aedq(H, rtol, Q, T, alpha, &st);
        dt = tictoc(t0);
        printf("aedq time: %.3f s (QR steps: %ld, AED deflated: %ld)\n",
               dt, st.steps, st.aed_deflated);
        printf("time to construct Q: %.3f s, time for AED: %.3f s\n",
               qaed_timers().q_time, qaed_timers().aed_time);
    } else {
        fprintf(stderr, "unknown alg %s\n", alg.c_str());
        return 2;
    }

    printf("unitarity ||Q'Q - I||/sqrt(n)      = %.3e\n", unitarity_err(Q));
    printf("residual  ||Q'HQ - T||/||H||       = %.3e\n", schur_residual(H, Q, T));
    printf("lower mass ||tril(T,-1)||/||T||    = %.3e\n", lower_mass(T, -1));

    return 0;
}
