//
// Copyright M. Novak: 2021
// Modified by M.Barbone 2021
//


#ifndef ODPM_CPU_INCLUDE_ELECTRON_TABLES_H_
#define ODPM_CPU_INCLUDE_ELECTRON_TABLES_H_

#include "spline_interpolator.h"
#include "walkers_approximator.h"

#ifdef GPU
#include "gpu_array.h"
#endif

class ElectronData;

namespace opmc {

// sampling tables to provide rejection free sampling of the energy that is
// transferred to the secondary gamma in electron bremsstrahlung according to
// the Seltzer-Berger DCS
    class SeltzerBerger : public WalkersApproximator {
    public:
        SeltzerBerger(const ElectronData &data, const int id);
//    const std::string datafile() const final { return "/brem_SBDtrData.dat"; }
        ODPM_INLINE ODPM_PURE
        constexpr real_type energyTransfer(real_type energy, real_type rndm1, real_type rndm2,
                                           real_type rndm3) const {
            // determine the primary electron energy lower grid point and sample if that or one above is used now
            real_type lpenergy = math::log(energy);
            real_type phigher = (lpenergy - m_log_min_primary_energy) * m_inv_log_delta_primary_energy;
            auto penergyindx = static_cast<unsigned>(phigher);
            phigher -= penergyindx;
            penergyindx += (rndm1 < phigher);
            // should always be fine if gamma-cut < eprim < E_max but make sure
            //  penergyindx       = std::min(fNumPrimaryEnergies-2, penergyindx);
            // sample the transformed variable
            const real_type xi = sample(penergyindx, rndm2, rndm3);
            // transform it back to kappa then to gamma energy (fMinPrimaryEnergy = gcut
            // and fLogMinPrimaryEnergy = log(gcut) so log(kappac) = log(gcut/eprim) =
            // = fLogMinPrimaryEnergy - lpenergy = -(lpenergy - fLogMinPrimaryEnergy) that
            // is computed above but keep it visible here)
            const real_type kappac = m_min_primary_energy / energy;
            const real_type kappa = kappac * math::exp(-xi * (m_log_min_primary_energy - lpenergy));
            return kappa * energy;
        }
    };

// sampling tables to provide rejection free sampling of the energy transferred
// to the secondary electron in Moller interaction
    class MollerEnergyTransfer : public WalkersApproximator {
    public:
        MollerEnergyTransfer(size_t numPrimaryEnergies, size_t numSecondaryEnergies, real_type minPrimaryEnergy,
                             real_type invLogDeltaPrimaryEnergy, std::vector<Point4D> data) : WalkersApproximator(
                numPrimaryEnergies, numSecondaryEnergies, minPrimaryEnergy, invLogDeltaPrimaryEnergy, data) {}

        MollerEnergyTransfer(const MollerEnergyTransfer &) = default;

        MollerEnergyTransfer(MollerEnergyTransfer &&) noexcept = default;

//    const std::string datafile() const final { return "/ioni_MollerDtrData.dat"; }
        ODPM_INLINE ODPM_PURE
        constexpr real_type sampleEnergyTransfer(real_type energy, real_type rndm1, real_type rndm2,
                                                 real_type rndm3) const noexcept {
            // determine the primary electron energy lower grid point and sample if that or one above is used now
            real_type lpenergy = math::log(energy);
            real_type phigher = (lpenergy - m_log_min_primary_energy) * m_inv_log_delta_primary_energy;
            auto penergyindx = static_cast<unsigned>(phigher);
            phigher -= penergyindx;
            penergyindx += rndm1 < phigher;
            // should always be fine if gamma-cut < eprim < E_max but make sure
            //  penergyindx       = std::min(fNumPrimaryEnergies-2, penergyindx);
            // sample the transformed variable xi=[kappa-ln(T_cut/T_0)]/[ln(T_cut/T_0)-ln(T_max/T_0)]
            // where kappa = ln(eps) with eps = T/T_0
            // so xi= [ln(T/T_0)-ln(T_cut/T_0)]/[ln(T_cut/T_0)-ln(T_max/T_0)] that is in [0,1]
            const real_type xi = sample(penergyindx, rndm2, rndm3);
            // transform it back to kappa then to gamma energy (fMinPrimaryEnergy = gcut
            // and fLogMinPrimaryEnergy = log(gcut) so log(kappac) = log(gcut/eprim) =
            // = fLogMinPrimaryEnergy - lpenergy = -(lpenergy - fLogMinPrimaryEnergy) that
            // is computed above but keep it visible here)
            const real_type dum1 = lpenergy - m_log_min_primary_energy;
            // return with the sampled kinetic energy transfered to the electron
            // (0.5*fMinPrimaryEnergy is the electron production cut)
            return math::exp(xi * dum1) * 0.5 * m_min_primary_energy;
        }
    };

    class ODPMElectronData {
    public:

        explicit ODPMElectronData(const ElectronData &data);

        ODPMElectronData(const ODPMElectronData &) = default;

        ODPMElectronData(ODPMElectronData &&) noexcept = default;

        virtual ~ODPMElectronData() = default;

        // samples value of energy transferred to the secondary electron in Moller
        // interaction at the given `eprim` primary electron energy. The three
        // additional input arguments are uniformly random values on [0,1].
        //
        // NOTE: it is assumed that: 2 x electron-cut < eprim < E_max (also for e+)
        ODPM_INLINE ODPM_PURE
        constexpr real_type mollerEnergyTransfer(real_type energy, real_type rndm1, real_type rndm2,
                                                 real_type rndm3) const noexcept {
            assert(constants::k_ElectronCut * 2 < energy);
            assert(constants::k_MaxEkin >= energy);
            ODPM_ASSUME(constants::k_ElectronCut * 2 < energy);
            ODPM_ASSUME(constants::k_MaxEkin >= energy);
            return m_mollerEnergyTransfer.sampleEnergyTransfer(energy, rndm1, rndm2, rndm3);
        }
        // sampling tables to provide rejection free sampling of the energy transferred
        // to the secondary electron in Moller interaction
        MollerEnergyTransfer m_mollerEnergyTransfer;

    };

}  // namespace opmc
#endif  // ODPM_CPU_INCLUDE_ELECTRON_TABLES_H_
