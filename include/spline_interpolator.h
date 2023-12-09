//
// Copyright M. Novak: 2021
// Modified by M.Barbone 2021
//

#ifndef ODPM_CPU_INCLUDE_SPLINE_INTERPOLATOR_H_
#define ODPM_CPU_INCLUDE_SPLINE_INTERPOLATOR_H_

#include <cassert>
#include <ostream>
#include <vector>

#include "constants.h"
#include "types.h"

namespace opmc {

class SplineInterpolator {
   public:
    SplineInterpolator(std::vector<Point5D> data) : m_data(std::move(data)) {
        m_span = span1d<Point5D>(m_data.data(), m_data.size());
        initialize();
    }

    SplineInterpolator(const SplineInterpolator &) = default;

    SplineInterpolator(SplineInterpolator &&)      = default;

    virtual ~SplineInterpolator()                  = default;

    ODPM_INLINE ODPM_PURE constexpr real_type eval(real_type xval) const noexcept {
        const auto logxval = math::approximate_logf(static_cast<float>(xval));
        auto ilow          = static_cast<unsigned>(((logxval - m_log_min_x) * m_inv_log_delta));
        // ilow is computed using an approximation of the logarithm which is much faster than math::log()
        // However it might lead to slightly wrong results that can be fixed with the if statement below.
        ilow += m_span[ilow + 1].x <= xval;
        ilow -= m_span[ilow].x > xval & ilow > 0;
        ilow -= m_span[ilow].x > xval & ilow > 0;

        assert(ilow <= unsigned(m_num_data - 2));
        ODPM_ASSUME(ilow <= m_num_data - 2);
        xval -= m_span[ilow].x;
        return m_span[ilow].y + xval * (m_span[ilow].b + xval * (m_span[ilow].c + xval * m_span[ilow].d));
    }

    void initialize();

   protected:
    span1d<Point5D> m_span;

   private:
    std::vector<Point5D> m_data;
    // Number of data stored
    std::size_t m_num_data;
    // The logarithm of the minimum x-value
    float m_log_min_x;
    // Inverse delta log x i.e. 1./ln[x_{i+1}/x_i]
    float m_inv_log_delta;
};

class MaterialSplineInterpolator : public SplineInterpolator {
   public:
    MaterialSplineInterpolator(std::vector<Point5D> data)
        : SplineInterpolator(std::move(data)),
          m_min_input{m_span[0].x},
          m_min_value{eval(m_min_input)},
          m_max_input{m_span.back().x - constants::k_EMaxThreshold},
          m_max_value{eval(m_max_input)} {}

    MaterialSplineInterpolator(const MaterialSplineInterpolator &) = default;

    MaterialSplineInterpolator(MaterialSplineInterpolator &&)      = default;

    virtual ~MaterialSplineInterpolator()                          = default;
    //
    ODPM_INLINE ODPM_CONST constexpr real_type sample(real_type ekin) const noexcept {
        // make sure that E_min <= ekin < E_max
        if (ekin <= m_min_input) {
            return m_min_value;
        }
        if (ekin >= m_max_input) {
            return m_max_value;
        }
        return eval(ekin);
    }

    real_type getMinInput() const noexcept { return m_min_input; }

    real_type getMaxInput() const noexcept { return m_max_input; }

   private:
    const real_type m_min_input, m_min_value, m_max_input, m_max_value;
};

}  // namespace opmc
#endif  // ODPM_CPU_INCLUDE_SPLINE_INTERPOLATOR_H_
