//
// Copyright M. Novak: 2021
// Modified by M.Barbone 2021
//

#ifndef ODPM_CPU_INCLUDE_WALKERS_APPROXIMATOR_H_
#define ODPM_CPU_INCLUDE_WALKERS_APPROXIMATOR_H_

#include <absl/container/inlined_vector.h>

#include <ostream>

#include "constants.h"
#include "types.h"

namespace opmc {

class WalkersApproximator {
   public:
    explicit WalkersApproximator(std::size_t numPrimaryEnergies, std::size_t numSecondaryEnergies,
                                 real_type minPrimaryEnergy, real_type invLogDeltaPrimaryEnergy,
                                 std::vector<Point4D> data)
        : m_num_primary_energies(numPrimaryEnergies),
          m_num_sampling_secondary_energies(numSecondaryEnergies),
          m_min_primary_energy(minPrimaryEnergy),
          m_log_min_primary_energy(math::log(m_min_primary_energy)),
          m_inv_log_delta_primary_energy(invLogDeltaPrimaryEnergy),
          m_data(data),
          m_tables(m_data.data(), m_num_primary_energies, m_num_sampling_secondary_energies) {}

    WalkersApproximator(const WalkersApproximator &)     = default;

    WalkersApproximator(WalkersApproximator &&) noexcept = default;

    virtual ~WalkersApproximator()                       = default;

   protected:
    ODPM_INLINE ODPM_PURE constexpr auto get_index(const auto x, const auto y) const noexcept {
        return x * m_num_sampling_secondary_energies + y;
    }

    // produce sample from the represented distribution u
    constexpr real_type sample(const unsigned ergyindx, const real_type rndm1, const real_type rndm2) const {
        // get the lower index of the bin by using the alias part
        const real_type rest = rndm1 * (m_num_sampling_secondary_energies - 1);
        int indxl            = static_cast<int>(rest);
        assert(indxl >= 0 && indxl < m_num_sampling_secondary_energies - 1);
        if (m_data[get_index(ergyindx, indxl)].alias_w < rest - indxl) {
            indxl = m_data[get_index(ergyindx, indxl)].alias_idx;
            assert(indxl >= 0 && indxl <= m_num_sampling_secondary_energies - 1);
        }
        // sample value within the selected bin by using linear aprox. of the p.d.f.
        real_type xval   = m_data[get_index(ergyindx, indxl)].x;
        real_type xdelta = m_data[get_index(ergyindx, indxl + 1)].x - xval;
        if (m_data[get_index(ergyindx, indxl)].y > 0.0) {
            real_type dum = (m_data[get_index(ergyindx, indxl + 1)].y - m_data[get_index(ergyindx, indxl)].y) /
                            m_data[get_index(ergyindx, indxl)].y;
            if (math::abs(dum) > 0.1) {
                return xval - xdelta / dum * (1.0 - math::sqrt(1.0 + rndm2 * dum * (dum + 2.0)));
            } else {  // use second order Taylor around dum = 0.0
                return xval + rndm2 * xdelta * (1.0 - 0.5 * dum * (rndm2 - 1.0) * (1.0 + dum * rndm2));
            }
        }
        return xval + xdelta * math::sqrt(rndm2);
    }

    // The number of primary energies
    const std::size_t m_num_primary_energies;
    // The size of the individual sampling tables
    const std::size_t m_num_sampling_secondary_energies;
    // Minimum primary energy from which sampling tables are built
    const real_type m_min_primary_energy;
    // The logarithm of the minimum primary energy.
    const real_type m_log_min_primary_energy;
    // Inverse delta log kinetic energy (i.e. bin size in log-scale).
    const real_type m_inv_log_delta_primary_energy;
    // data -- using span because vector is not 2D
    const std::vector<Point4D> m_data;
    // reference to the sampling points
    const span2d<const Point4D> m_tables;
};

}  // namespace opmc
#endif  // ODPM_CPU_INCLUDE_WALKERS_APPROXIMATOR_H_
