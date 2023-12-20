//
// Created by mbarbone on 11/30/21.
//

#ifndef ODPM_CPU_INCLUDE_ELEMENTARY_FUNCTIONS_H_
#define ODPM_CPU_INCLUDE_ELEMENTARY_FUNCTIONS_H_

#include <bit>
#include <cmath>

#include "utils.h"

namespace opmc::math {

namespace {
/* compute natural logarithm,
 * source https://github.com/ekmett/approximate/blob/dc1ee7cef58a6b31661edde6ef4a532d6fc41b18/cbits/fast.c#L127
 */
ODPM_INLINE ODPM_PURE constexpr float approximate_logf(const float a) noexcept {
    return (std::bit_cast<int32_t>(a) - 1064866805) * 8.262958405176314e-8f; /* 1 / 12102203.0; */
}

/* compute natural logarithm, maximum error 0.85089 ulps
 * source https://stackoverflow.com/questions/39821367/very-fast-approximate-logarithm-natural-log-function-in-c
 */
ODPM_INLINE ODPM_PURE constexpr float faster_logf(const float a) noexcept {
    const int32_t e = (std::bit_cast<int32_t>(a) - 0x3f2aaaab) & 0xff800000;
    const float m   = std::bit_cast<float>(std::bit_cast<int32_t>(a) - e);
    const float i   = static_cast<float>(e) * 1.19209290e-7f;  // 0x1.0p-23
    /* m in [2/3, 4/3] */
    const float f = m - 1.0f;
    const float s = f * f;
    /* Compute log1p(f) for f in [-1/3, 1/3] */
    float r       = std::fma(0.230836749f, f, -0.279208571f);  // 0x1.d8c0f0p-3, -0x1.1de8dap-2
    const float t = std::fma(0.331826031f, f, -0.498910338f);  // 0x1.53ca34p-2, -0x1.fee25ap-2
    r             = std::fma(r, s, t);
    r             = std::fma(r, s, f);
    r             = std::fma(i, 0.693147182f, r);  // 0x1.62e430p-1 // log(2)
    return r;
}

ODPM_INLINE ODPM_PURE constexpr float fast_logf(float a) noexcept {
    float i = 0.0f;
    if (a < 1.175494351e-38f) {  // 0x1.0p-126
        a = a * 8388608.0f;      // 0x1.0p+23
        i = -23.0f;
    }
    float e = (std::bit_cast<int32_t>(a) - std::bit_cast<int32_t>(0.666666667f)) & 0xff800000;
    float m = std::bit_cast<float>(std::bit_cast<int32_t>(a) - e);
    i       = std::fmaf((float)e, 1.19209290e-7f, i);  // 0x1.0p-23
    /* m in [2/3, 4/3] */
    m       = m - 1.0f;
    float s = m * m;
    /* Compute log1p(m) for m in [-1/3, 1/3] */
    float r = -0.130310059f;                   // -0x1.0ae000p-3
    float t = 0.140869141f;                    //  0x1.208000p-3
    r       = std::fmaf(r, s, -0.121483512f);  // -0x1.f198b2p-4
    t       = std::fmaf(t, s, 0.139814854f);   //  0x1.1e5740p-3
    r       = std::fmaf(r, s, -0.166846126f);  // -0x1.55b36cp-3
    t       = std::fmaf(t, s, 0.200120345f);   //  0x1.99d8b2p-3
    r       = std::fmaf(r, s, -0.249996200f);  // -0x1.fffe02p-3
    r       = std::fmaf(t, m, r);
    r       = std::fmaf(r, m, 0.333331972f);   //  0x1.5554fap-2
    r       = std::fmaf(r, m, -0.500000000f);  // -0x1.000000p-1
    r       = std::fmaf(r, s, m);
    r       = std::fmaf(i, 0.693147182f, r);  //  0x1.62e430p-1 // log(2)
    if (!((a > 0.0f) && (a < INFINITY))) {
        r = a + a;                              // silence NaNs if necessary
        if (a < 0.0f) r = INFINITY - INFINITY;  //  NaN
        if (a == 0.0f) r = -INFINITY;
    }
    return r;
}
}  // namespace

template <typename T>
ODPM_INLINE ODPM_PURE constexpr T sin(T x) noexcept {
    return std::sin(x);
}

template <typename T>
ODPM_INLINE ODPM_PURE constexpr T cos(T x) noexcept {
    return std::cos(x);
}

template <typename T>
ODPM_INLINE ODPM_PURE constexpr T log(T x) noexcept {
    return std::log(x);
}

template <typename T>
ODPM_INLINE ODPM_PURE constexpr T log1p(T x) noexcept {
    return std::log1p(x);
}

template <typename T>
ODPM_INLINE ODPM_PURE constexpr T exp(T x) noexcept {
    return std::exp(x);
}

template <typename T>
ODPM_INLINE ODPM_PURE constexpr T sqrt(T x) noexcept {
    return std::sqrt(x);
}

template <typename T, typename V>
ODPM_INLINE ODPM_PURE constexpr T pow(T x, V y) noexcept {
    return std::pow(x, y);
}

template <typename T, typename V>
ODPM_INLINE ODPM_PURE constexpr T min(T x, V y) noexcept {
    return x < y ? x : y;
}

template <typename T, typename V>
ODPM_INLINE ODPM_PURE constexpr T max(T x, V y) noexcept {
    return x > y ? x : y;
}

template <typename T>
    requires std::is_floating_point_v<T>
ODPM_INLINE ODPM_PURE constexpr bool signbit(T x) noexcept {
    return __builtin_signbit(x);
}

template <typename T>
ODPM_INLINE ODPM_PURE constexpr T abs(T x) noexcept {
    return x > 0 ? x : -x;
    ;
}

}  // namespace opmc::math

#endif  // ODPM_CPU_INCLUDE_ELEMENTARY_FUNCTIONS_H_
