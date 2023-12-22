//
// Created by mbarbone on 7/28/21.
//

#ifndef ODPM_CPU_INCLUDE_RANDOM_H_
#define ODPM_CPU_INCLUDE_RANDOM_H_

#include <mixmax/mixmax.h>

#include <random>

#include "types.h"

namespace opmc {

class Random {
   public:
    explicit Random(unsigned long seed) : m_rng(seed) {}

    explicit Random(unsigned long seed, unsigned long stream) : m_rng(seed, stream) {}

    explicit Random(unsigned int seed0, unsigned int seed1, unsigned int seed2, unsigned int seed3)
        : m_rng(seed0, seed1, seed2, seed3) {}

    Random(const Random &)            = default;

    Random(Random &&)                 = default;

    Random &operator=(const Random &) = default;

    Random &operator=(Random &&)      = default;

    ODPM_INLINE constexpr real_type getUniform() noexcept { return m_rng.Uniform(); }

    ODPM_INLINE constexpr real_type getExponential() noexcept { return -math::faster_logf(m_rng.Uniform()); }

   private:
    MIXMAX::MixMaxRng17 m_rng;
};

}  // namespace opmc

#endif  // ODPM_CPU_INCLUDE_RANDOM_H_
