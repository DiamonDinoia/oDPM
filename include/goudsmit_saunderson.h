//
// Copyright M. Novak: 2021
// Modified by M.Barbone 2021
//

#ifndef ODPM_CPU_SRC_GOUDSMITSAUNDERSON_H_
#define ODPM_CPU_SRC_GOUDSMITSAUNDERSON_H_

#include "types.h"

#ifdef GPU
#include "gpu_array.h"
#endif

class ElectronData;

namespace opmc {

// sampling tables to provide rejection free sampling of the angular deflection
// due to multiple Coulomb scattering in the reference material along the
// maximally allowed MSC step length that is a function of the energy and
// determined by a (sigmoid-like) function of the `slow`, `shigh` and `ecross`
// parameters
class GoudsmitSaunderson {
   public:
    struct OnePoint {
        real_type var_u;    // the transfomred variable `u` that corresponds this cumulative value
        real_type param_a;  // interpolation parameters for the approximate inversion of the cumulative
        real_type param_b;
    };

    GoudsmitSaunderson(int numPrimaryEnergies, int numCumData, real_type minPrimaryEnergy,
                       std::vector<real_type> primaryEnergyGrid, real_type logMinPrimaryEnergy,
                       real_type invLogDeltaPrimaryEnergy, real_type invDeltaCumulative,
                       std::vector<real_type> transformParams, std::vector<OnePoint> tables)
        : m_num_primary_energies(numPrimaryEnergies),
          m_sampling_table_size(numCumData),
          m_min_primary_energy(minPrimaryEnergy),
          m_primary_energy_grid_data(std::move(primaryEnergyGrid)),
          m_log_min_primary_energy(logMinPrimaryEnergy),
          m_inv_log_delta_primary_energy(invLogDeltaPrimaryEnergy),
          m_inv_delta_cumulative(invDeltaCumulative),
          m_transform_params_data(std::move(transformParams)),
          m_data(std::move(tables)) {
        m_primary_energy_grid = span1d<real_type>(m_primary_energy_grid_data.data(), m_num_primary_energies);
        m_transform_params    = span1d<real_type>(m_transform_params_data.data(), m_num_primary_energies);
        m_tables              = span2d<OnePoint>(m_data.data(), m_num_primary_energies, m_sampling_table_size);
    }

    GoudsmitSaunderson(const ElectronData &electron_data, const int id);

    GoudsmitSaunderson(const GoudsmitSaunderson &) = default;

    GoudsmitSaunderson(GoudsmitSaunderson &&)      = default;

    virtual ~GoudsmitSaunderson()                  = default;

    // The number of discrete primary energies at which sampling tables are built.
    int m_num_primary_energies;

    int m_sampling_table_size;

    // Minimum primary energy from which sampling tables are built
    real_type m_min_primary_energy;

    std::vector<real_type> m_primary_energy_grid_data;
    // The electron kinetic energy grid
    span1d<real_type> m_primary_energy_grid;  // [#fNumPrimaryEnergies]
    // The logarithm of the minimum primary energy.
    real_type m_log_min_primary_energy;
    // Inverse delta log kinetic energy (i.e. bin size in log-scale).
    real_type m_inv_log_delta_primary_energy;
    // Delta cumulative value: 1/(fSamplingTableSize-1)
    real_type m_inv_delta_cumulative;

    std::vector<real_type> m_transform_params_data;
    // transformation paraneter `a` used to generate the corresponding smooth pdf
    span1d<real_type> m_transform_params;

    std::vector<OnePoint> m_data;
    // The sampling table over the discrete primary energy grid for this material.
    span2d<OnePoint> m_tables;

    inline constexpr real_type angularDeflection(real_type eprim, real_type rndm1, real_type rndm2) const noexcept {
        // determine electron energy lower grid point and sample if that or one above is used now
        real_type lpenergy = math::log(eprim);
        real_type phigher  = (lpenergy - m_log_min_primary_energy) * m_inv_log_delta_primary_energy;
        int penergyindx    = static_cast<int>(phigher);
        // keep the lower index of the energy bin
        const int ielow = penergyindx;
        phigher -= penergyindx;
        penergyindx += rndm1 < phigher;
        // should always be fine if electron-cut < eprim < E_max but make sure
        //  penergyindx      = std::min(fNumPrimaryEnergies-1, penergyindx);
        // sample the transformed variable \xi
        // lower index of the (common) discrete cumulative bin and the residual fraction
        const int indxl       = static_cast<int>(rndm2 / m_inv_delta_cumulative);
        const real_type resid = rndm2 - indxl * m_inv_delta_cumulative;
        // compute value `u` by using ratin based numerical inversion
        const real_type parA = m_tables(penergyindx, indxl).param_a;
        const real_type parB = m_tables(penergyindx, indxl).param_b;
        const real_type u0   = m_tables(penergyindx, indxl).var_u;
        const real_type u1   = m_tables(penergyindx, indxl + 1).var_u;
        const real_type dum1 = (1.f + parA + parB) * m_inv_delta_cumulative * resid;
        const real_type dum2 = m_inv_delta_cumulative * m_inv_delta_cumulative + parA * m_inv_delta_cumulative * resid +
                               parB * resid * resid;
        const real_type theU = u0 + dum1 / dum2 * (u1 - u0);
        // transform back the sampled `u` to `mu(u)` using the transformation parameter `a`
        // mu(u) = 1 - 2au/[1-u+a] as given by Eq.(34)
        // interpolate (linearly) the transformation parameter to E
        const real_type a0        = m_transform_params[ielow];
        const real_type a1        = m_transform_params[ielow + 1];
        const real_type e0        = m_primary_energy_grid[ielow];
        const real_type e1        = m_primary_energy_grid[ielow + 1];
        const real_type parTransf = (a1 - a0) / (e1 - e0) * (eprim - e0) + a0;
        return 1.f - 2.f * parTransf * theU / (1.f - theU + parTransf);
    }
};
}  // namespace opmc
#endif  // ODPM_CPU_SRC_GOUDSMITSAUNDERSON_H_
