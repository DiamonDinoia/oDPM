//
// Copyright M. Novak: 2021
// Modified by M.Barbone 2021
//

#ifndef ODPM_CPU_SRC_LINEARINTERPOLATOR_H_
#define ODPM_CPU_SRC_LINEARINTERPOLATOR_H_

#include <cassert>

#include "types.h"

#ifdef GPU
#include "gpu_array.h"
#endif

namespace opmc {
class LinearInterpolator {
   public:
    explicit LinearInterpolator(std::vector<Point2D> data)
        : m_min_x{data[0].x},
          m_max_x{data[data.size() - 2].x},
          m_inv_delta{1. / (data[1].x - data[0].x)},
          m_data(data),
          m_span(m_data.data(), m_data.size()) {}

    LinearInterpolator(const LinearInterpolator &)     = default;

    LinearInterpolator(LinearInterpolator &&) noexcept = default;

    virtual ~LinearInterpolator()                      = default;

    ODPM_INLINE ODPM_PURE constexpr real_type eval(const real_type xval) const noexcept {
        auto ilow = static_cast<unsigned int>((xval - m_min_x) * m_inv_delta);
        assert(ilow <= m_span.size() - 2);
        return (m_data[ilow + 1].y - m_data[ilow].y) * (xval - m_data[ilow].x) / (m_data[ilow + 1].x - m_data[ilow].x) +
               m_data[ilow].y;
    }

   private:
    // Minimum x-value
    const real_type m_min_x;
    // Maximum x-value
    const real_type m_max_x;
    // Inverse delta x i.e. 1./[x_{i+1}-x_i]
    const real_type m_inv_delta;
    // interpolation points
    // data -- using span because vector are not constexpr in cpp < 20
    // this also simplifies GPU support
    const std::vector<Point2D> m_data;

    const span1d<const Point2D> m_span;
};

}  // namespace opmc

#endif  // ODPM_CPU_SRC_LINEARINTERPOLATOR_H_
