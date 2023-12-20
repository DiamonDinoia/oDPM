//
// Created by mbarbone on 12/13/21.
//

#ifndef ODPM_CPU_INCLUDE_UTILS_H_
#define ODPM_CPU_INCLUDE_UTILS_H_

#ifndef __CUDACC__
#define ODPM_UNLIKELY(expr) __builtin_expect(!!(expr), 0)
#define ODPM_LIKELY(expr) __builtin_expect(!!(expr), 1)
#else
#define ODPM_UNLIKELY(expr) expr
#define ODPM_LIKELY(expr) expr
#endif

#define ODPM_INLINE inline __attribute__((always_inline))
#define ODPM_NOINLINE __attribute__((noinline))
#define ODPM_CONST __attribute__((const))
#define ODPM_CONSTEXPR constexpr
#define ODPM_PURE __attribute__((pure))
#define ODPM_UNREACHABLE __builtin_unreachable()
#define ODPM_ASSUME(expr)                     \
    do {                                      \
        if (!(expr)) __builtin_unreachable(); \
    } while (0)
#endif  // ODPM_CPU_INCLUDE_UTILS_H_
