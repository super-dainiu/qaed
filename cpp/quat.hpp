// quat.hpp - Quaternion scalar type for qaed C++ port.
#pragma once
#include <cmath>
#include <complex>
#include <cstdio>

namespace qaed {

struct Quat {
    double w{0}, x{0}, y{0}, z{0};

    Quat() = default;
    Quat(double w_, double x_, double y_, double z_) : w(w_), x(x_), y(y_), z(z_) {}
    explicit Quat(double s) : w(s) {}

    static Quat zero() { return Quat(); }
    static Quat one()  { return Quat(1, 0, 0, 0); }

    Quat operator+(const Quat& b) const { return {w + b.w, x + b.x, y + b.y, z + b.z}; }
    Quat operator-(const Quat& b) const { return {w - b.w, x - b.x, y - b.y, z - b.z}; }
    Quat operator-() const { return {-w, -x, -y, -z}; }

    Quat& operator+=(const Quat& b) { w += b.w; x += b.x; y += b.y; z += b.z; return *this; }
    Quat& operator-=(const Quat& b) { w -= b.w; x -= b.x; y -= b.y; z -= b.z; return *this; }

    Quat operator*(const Quat& b) const {
        return {
            w * b.w - x * b.x - y * b.y - z * b.z,
            w * b.x + x * b.w + y * b.z - z * b.y,
            w * b.y - x * b.z + y * b.w + z * b.x,
            w * b.z + x * b.y - y * b.x + z * b.w
        };
    }
    Quat& operator*=(const Quat& b) { *this = (*this) * b; return *this; }

    Quat operator*(double s) const { return {w * s, x * s, y * s, z * s}; }
    Quat operator/(double s) const { return {w / s, x / s, y / s, z / s}; }
    friend Quat operator*(double s, const Quat& q) { return q * s; }
    Quat operator+(double s) const { return {w + s, x, y, z}; }
    friend Quat operator+(double s, const Quat& q) { return q + s; }

    Quat conj() const { return {w, -x, -y, -z}; }
    double norm2() const { return w * w + x * x + y * y + z * z; }
    double abs() const { return std::sqrt(norm2()); }

    Quat inv() const { double n = norm2(); return {w / n, -x / n, -y / n, -z / n}; }
    Quat operator/(const Quat& b) const { return (*this) * b.inv(); }

    // Complex split q = c1 + c2 * j with c1 = w + x*i, c2 = y + z*i.
    std::complex<double> c1() const { return {w, x}; }
    std::complex<double> c2() const { return {y, z}; }
};

inline Quat conj(const Quat& q) { return q.conj(); }
inline double abs(const Quat& q) { return q.abs(); }
inline Quat from_c(std::complex<double> c1, std::complex<double> c2) {
    return {c1.real(), c1.imag(), c2.real(), c2.imag()};
}

// schur1: standardize a quaternion. Returns unit quaternion u such that
// conj(u) * q * u = w + |v| * i  (complex, nonnegative imaginary part).
inline void schur1(const Quat& q, Quat& u, Quat& t) {
    double a = std::sqrt(q.x * q.x + q.y * q.y + q.z * q.z);
    if (a == 0.0) { u = Quat::one(); t = q; return; }
    double xh = q.x / a, yh = q.y / a, zh = q.z / a;
    if (xh >= 0.0) {
        // u = (1 - vhat*i) / sqrt(2*(1+xh)) = ((1+xh) - zh*j + yh*k)/mu
        double mu = std::sqrt(2.0 * (1.0 + xh));
        u = Quat(1.0 + xh, 0.0, -zh, yh) / mu;
    } else {
        // u = ((1 + vhat*i)/sqrt(2*(1-xh))) * j ; 1 + vhat*i = (1-xh) + zh*j - yh*k
        double mu = std::sqrt(2.0 * (1.0 - xh));
        u = (Quat(1.0 - xh, 0.0, zh, -yh) / mu) * Quat(0, 0, 1, 0);
    }
    t = Quat(q.w, a, 0.0, 0.0);
}

// sylvesterc: solve a*x - x*b = c where a, b are complex (as quaternions with
// zero j,k parts) and c is a quaternion.
inline Quat sylvesterc(const Quat& a, const Quat& b, const Quat& c) {
    std::complex<double> a1(a.w, a.x);
    std::complex<double> b1(b.w, b.x), b2(b.w, -b.x);
    std::complex<double> c1(c.w, c.x), c2(c.y, c.z);
    std::complex<double> x1 = c1 / (a1 - b1);
    std::complex<double> x2 = c2 / (a1 - b2);
    return from_c(x1, x2);
}

} // namespace qaed
