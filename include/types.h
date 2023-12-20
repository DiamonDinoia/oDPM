//
// Created by mbarbone on 7/28/21.
//

#ifndef ODPM_CPU_INCLUDE_TYPES_H_
#define ODPM_CPU_INCLUDE_TYPES_H_

#include <absl/types/span.h>

#include <cstring>
#include <experimental/mdspan>
#include <iostream>

#include "elementary_functions.h"
#include "utils.h"

namespace opmc {

#ifdef SINGLE
using real_type = float;
#else
using real_type = double;
#endif

using unsigned_type = unsigned int;

template <typename T, typename V>
struct Pair {
    T first;
    V second;
};

struct Point2D {
    real_type x, y;
};

// data that describes a single point in the table
struct Point4D {
    // the discretized stochastic variable value
    real_type x;
    // the corresponding p.d.f. value (not necessarily normalised over fXdata)
    real_type y;
    // the correspnding alias probability (not necessarily normalised over fXdata)
    real_type alias_w;
    // the corresponding alias index
    int alias_idx;

    friend std::ostream &operator<<(std::ostream &os, const Point4D &point) {
        os << "fXdata: " << point.x << " fYdata: " << point.y << " fAliasW: " << point.alias_w
           << " fAliasIndx: " << point.alias_idx;
        return os;
    }
};

struct Point5D {
    real_type x, y, b, c, d;

    friend std::ostream &operator<<(std::ostream &os, const Point5D &point);
};

template <typename T>
using span1d = absl::Span<T>;

template <typename T, auto size_0 = std::experimental::dynamic_extent, auto size_1 = std::experimental::dynamic_extent>
using span2d = std::experimental::mdspan<T, std::experimental::extents<std::size_t, size_0, size_1>>;

template <typename T = real_type>
class DoubleVector {
   public:
    T x;
    T y;

    constexpr DoubleVector() noexcept : x(0), y(0){};

    constexpr DoubleVector(T x, T y) noexcept : x(x), y(y){};

    constexpr DoubleVector(const DoubleVector &other) noexcept = default;

    constexpr DoubleVector(DoubleVector &&other) noexcept      = default;

    //
    DoubleVector<T> &operator=(const DoubleVector<T> &other) noexcept = default;

    DoubleVector<T> &operator=(DoubleVector<T> &&other) noexcept      = default;

    friend std::ostream &operator<<(std::ostream &os, const DoubleVector &vector) {
        os << "x: " << vector.x << " y: " << vector.y;
        return os;
    }
};

template <typename T = real_type>
class ThreeVector {
   public:
    T x;
    T y;
    T z;

    ODPM_INLINE constexpr ThreeVector() noexcept : x(0), y(0), z(0){};

    ODPM_INLINE constexpr ThreeVector(T x, T y, T z) noexcept : x(std::move(x)), y(std::move(y)), z(std::move(z)) {}

    template <typename S1, typename S2, typename S3>
    ODPM_INLINE constexpr ThreeVector(S1 &&x, S2 &&y, S3 &&z) noexcept
        : x(std::forward<T>(x)), y(std::forward<T>(y)), z(std::forward<T>(z)) {}

    ODPM_INLINE constexpr ThreeVector(const ThreeVector<T> &other) noexcept = default;

    ODPM_INLINE constexpr ThreeVector(ThreeVector<T> &&other) noexcept      = default;

    template <typename V>
    ODPM_INLINE constexpr friend ThreeVector<T> &operator+=(ThreeVector<T> &vector,
                                                            const ThreeVector<V> &other) noexcept {
        vector.x += other.x;
        vector.y += other.y;
        vector.z += other.z;
        return vector;
    }

    template <typename V>
    ODPM_INLINE constexpr friend ThreeVector<T> &operator+=(ThreeVector<T> &vector, const V &other) noexcept {
        constexpr ThreeVector<V> b{other, other, other};
        vector += b;
        return vector;
    }

    template <typename V>
    ODPM_INLINE constexpr friend ThreeVector<T> operator+(ThreeVector<T> vector, const ThreeVector<V> &other) noexcept {
        vector += other;
        return vector;
    }

    template <typename V>
    ODPM_INLINE constexpr friend ThreeVector<T> operator+(ThreeVector<T> vector, const V &other) noexcept {
        vector += other;
        return vector;
    }

    template <typename V>
    ODPM_INLINE constexpr friend ThreeVector<T> &operator-=(ThreeVector<T> &vector,
                                                            const ThreeVector<V> &other) noexcept {
        vector.x -= other.x;
        vector.y -= other.y;
        vector.z -= other.z;
        return vector;
    }

    template <typename V>
    ODPM_INLINE constexpr friend ThreeVector<T> &operator-=(ThreeVector<T> &vector, const V &other) noexcept {
        vector -= ThreeVector<V>{other, other, other};
        return vector;
    }

    template <typename V>
    ODPM_INLINE constexpr friend ThreeVector<T> operator-(ThreeVector<T> vector, const ThreeVector<V> &other) noexcept {
        vector -= other;
        return vector;
    }

    template <typename V>
    ODPM_INLINE constexpr friend ThreeVector<T> operator-(ThreeVector<T> vector, const V &other) noexcept {
        vector -= other;
        return vector;
    }

    ODPM_INLINE constexpr ThreeVector<T> operator-() noexcept {
        x = -x;
        y = -y;
        z = -z;
        return *this;
    }

    template <typename V>
    ODPM_INLINE constexpr friend ThreeVector<T> &operator*=(ThreeVector<T> &vector,
                                                            const ThreeVector<V> &other) noexcept {
        vector.x *= other.x;
        vector.y *= other.y;
        vector.z *= other.z;
        return vector;
    }

    template <typename V>
    ODPM_INLINE constexpr friend ThreeVector<T> &operator*=(ThreeVector<T> &vector, const V &other) noexcept {
        const ThreeVector<V> b{other, other, other};
        vector *= b;
        return vector;
    }

    template <typename V>
    ODPM_INLINE constexpr friend ThreeVector<T> operator*(ThreeVector<T> vector, const ThreeVector<V> &other) noexcept {
        vector *= other;
        return vector;
    }

    template <typename V>
    ODPM_INLINE constexpr friend ThreeVector<T> operator*(ThreeVector<T> vector, const V &other) noexcept {
        vector *= other;
        return vector;
    }

    template <typename V>
    ODPM_INLINE constexpr friend ThreeVector<T> &operator/=(ThreeVector<T> &vector,
                                                            const ThreeVector<V> &other) noexcept {
        vector.x /= other.x;
        vector.y /= other.y;
        vector.z /= other.z;
        return vector;
    }

    template <typename V>
    ODPM_INLINE constexpr friend ThreeVector<T> &operator/=(ThreeVector<T> &vector, const V &other) noexcept {
        vector /= ThreeVector<V>{other, other, other};
        return vector;
    }

    template <typename V>
    ODPM_INLINE constexpr friend ThreeVector<T> operator/(ThreeVector<T> vector, const ThreeVector<V> &other) noexcept {
        vector /= other;
        return vector;
    }

    template <typename V>
    ODPM_INLINE constexpr friend ThreeVector<T> operator/(ThreeVector<T> vector, const V &other) noexcept {
        vector /= other;
        return vector;
    }

    template <typename V>
    ODPM_INLINE constexpr friend ThreeVector<T> &operator%=(ThreeVector<T> &vector,
                                                            const ThreeVector<V> &other) noexcept {
        vector.x %= other.x;
        vector.y %= other.y;
        vector.z %= other.z;
        return vector;
    }

    template <typename V>
    ODPM_INLINE constexpr friend ThreeVector<T> &operator%=(ThreeVector<T> &vector, const V &other) noexcept {
        vector %= ThreeVector<V>{other, other, other};
        return vector;
    }

    template <typename V>
    ODPM_INLINE constexpr friend ThreeVector<T> operator%(ThreeVector<T> vector, const ThreeVector<V> &other) noexcept {
        vector %= other;
        return vector;
    }

    template <typename V>
    ODPM_INLINE constexpr friend ThreeVector<T> operator%(ThreeVector<T> vector, const V &other) noexcept {
        vector %= other;
        return vector;
    }

    ODPM_INLINE constexpr friend bool operator==(ThreeVector<T> vector, const ThreeVector<T> &other) noexcept {
        return vector.x == other.x && vector.y == other.y && vector.z == other.z;
    }

    ODPM_INLINE constexpr friend bool operator!=(ThreeVector<T> vector, const ThreeVector<T> &other) noexcept {
        return vector.x != other.x || vector.y != other.y || vector.z != other.z;
    }

    /* needed to implement and std::map<ThreeVector, T> */
    template <typename V>
    ODPM_INLINE constexpr friend bool operator<(ThreeVector<T> vector, const ThreeVector<V> &other) noexcept {
        return vector.x < other.x && vector.y < other.y && vector.z < other.z;
    }

    template <typename V>
    ODPM_INLINE constexpr friend bool operator>=(ThreeVector<T> vector, const ThreeVector<V> &other) noexcept {
        return vector.x >= other.x && vector.y >= other.y && vector.z >= other.z;
    }

    template <typename V>
    ODPM_INLINE constexpr friend bool operator>(ThreeVector<T> vector, const ThreeVector<V> &other) noexcept {
        return vector.x > other.x && vector.y > other.y && vector.z > other.z;
    }

    ODPM_INLINE ThreeVector<T> &operator=(const ThreeVector<T> &other)     = default;

    ODPM_INLINE ThreeVector<T> &operator=(ThreeVector<T> &&other) noexcept = default;

    ODPM_INLINE constexpr friend std::ostream &operator<<(std::ostream &os, const ThreeVector<T> &vector) noexcept {
        return os << "[ " << vector.x << " " << vector.y << " " << vector.z << " ]";
    }

    template <typename V>
    ODPM_INLINE ODPM_PURE constexpr operator ThreeVector<V>() const noexcept {
        return {static_cast<V>(x), static_cast<V>(y), static_cast<V>(z)};
    }

    ODPM_INLINE constexpr ThreeVector<T> &floor() noexcept {
        x = std::floor(x);
        y = std::floor(y);
        z = std::floor(z);
        return *this;
    }

    ODPM_INLINE ODPM_PURE constexpr T min() const noexcept { return math::min(math::min(x, y), z); }

    ODPM_INLINE constexpr void rotate(T cosTheta, T phi) noexcept {
        assert(-1 <= cosTheta && cosTheta <= 1);

        const T rho2      = x * x + y * y;
        const T sinPhi    = math::sin(phi);
        const T cosPhi    = math::cos(phi);
        const T sinTheta2 = 1. - cosTheta * cosTheta;

        if (rho2 > 0.) {
            T sthrho = math::sqrt(sinTheta2 / rho2);
            T urho   = x * sthrho;
            T vrho   = y * sthrho;
            x        = x * cosTheta - vrho * sinPhi + z * urho * cosPhi;
            y        = y * cosTheta + urho * sinPhi + z * vrho * cosPhi;
            z        = z * cosTheta - rho2 * sthrho * cosPhi;
        } else {
            T sinTheta = math::sqrt(sinTheta2);
            y          = sinTheta * sinPhi;
            if (z > 0.) {
                x = sinTheta * cosPhi;
                z = cosTheta;
            } else {
                x = -sinTheta * cosPhi;
                z = -cosTheta;
            }
        }
    }

    template <typename V>
    ODPM_INLINE ODPM_PURE constexpr T dot(const ThreeVector<V> &other) const noexcept {
        return x * other.x + y * other.y + z * other.z;
    }
};

}  // namespace opmc

#endif  // ODPM_CPU_INCLUDE_TYPES_H_
