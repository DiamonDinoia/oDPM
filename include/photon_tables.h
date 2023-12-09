//
// Copyright M. Novak: 2021
// Modified by M.Barbone 2021
//

#ifndef ODPM_CPU_SRC_PHOTON_TABLES_H_
#define ODPM_CPU_SRC_PHOTON_TABLES_H_

#include "linear_interpolator.h"
#include "walkers_approximator.h"

#ifdef GPU
#include "gpu_array.h"
#endif

class PhotonData;

namespace opmc {

class PhotonLinearInterpolator : public LinearInterpolator {
   public:
    explicit PhotonLinearInterpolator(std::vector<Point2D> data)
        : LinearInterpolator(data),
          m_min_input{data.front().x},
          m_min_value{eval(m_min_input)},
          m_max_input{data.back().x - constants::k_EMaxThreshold},
          m_max_value{eval(m_max_input)} {};

    PhotonLinearInterpolator(const PhotonLinearInterpolator &)     = default;

    PhotonLinearInterpolator(PhotonLinearInterpolator &&) noexcept = default;

    virtual ~PhotonLinearInterpolator()                            = default;

    ODPM_INLINE ODPM_PURE constexpr real_type sample(real_type ekin) const noexcept {
        // make sure that E_min <= ekin < E_max
        if (ODPM_UNLIKELY(ekin <= m_min_input)) {
            return m_min_value;
        }
        if (ODPM_UNLIKELY(ekin >= m_max_input)) {
            return m_max_value;
        }
        return eval(ekin);
    }

    real_type GetMinInput() const { return m_min_input; }

    real_type GetMaxInput() const { return m_max_input; }

   private:
    const real_type m_min_input, m_min_value, m_max_input, m_max_value;
};

// sampling tables to provide rejection free sampling of the (reduced) energy
// transferred to the secondary electron in Compton scattering according to
// the Klein-Nishina DCS
class KleinNishina : public WalkersApproximator {
   public:
    explicit KleinNishina(size_t samplingTableSize, size_t numPrimaryEnergies, real_type minPrimaryEnergy,
                          real_type invLogDeltaPrimaryEnergy, std::vector<Point4D> data)
        : WalkersApproximator(samplingTableSize, numPrimaryEnergies, minPrimaryEnergy, invLogDeltaPrimaryEnergy, data) {
    }

    KleinNishina(const KleinNishina &)     = default;

    KleinNishina(KleinNishina &&) noexcept = default;

    virtual ~KleinNishina()                = default;

    constexpr real_type sample_energy_transfer(const real_type egamma, const real_type rndm1, const real_type rndm2,
                                               const real_type rndm3) const {
        // determine the primary electron energy lower grid point and sample if that or one above is used now
        real_type lpenergy = math::log(egamma);
        real_type phigher  = (lpenergy - m_log_min_primary_energy) * m_inv_log_delta_primary_energy;
        int penergyindx    = static_cast<int>(phigher);
        phigher -= penergyindx;
        penergyindx += rndm1 < phigher;
        // should always be fine if gamma-cut < egamma < E_max but make sure
        //  penergyindx      = std::min(fNumPrimaryEnergies-1, penergyindx);
        // sample the transformed variable xi=[\alpha-ln(ep)]/\alpha (where \alpha=ln(1/(1+2\kappa)))
        // that is in [0,1] when ep is in [ep_min=1/(1+2\kappa),ep_max=1] (that limits comes from energy and momentum
        // conservation in case of scattering on free electron at rest).
        // where ep = E_1/E_0 and kappa = E_0/(mc^2)
        real_type xi = sample(penergyindx, rndm2, rndm3);
        // transform it back to eps = E_1/E_0
        // \epsion(\xi) = \exp[ \alpha(1-\xi) ] = \exp [\ln(1+2\kappa)(\xi-1)]
        real_type kappa = egamma * k_INVEMC2;
        return math::exp(math::log(1. + 2. * kappa) * (xi - 1.));  // eps = E_1/E
    }

   private:
    constexpr static real_type k_INVEMC2 = 1. / constants::k_EMC2;
};

class ODPMPhotonData {
   public:
    ODPMPhotonData(const PhotonData &data);

    ODPMPhotonData(const ODPMPhotonData &)     = default;

    ODPMPhotonData(ODPMPhotonData &&) noexcept = default;

    virtual ~ODPMPhotonData()                  = default;

    // the global maximum inverse MFP in [1/mm] units
    constexpr real_type iMFPMaxPhoton(real_type ekin) const { return m_IMFPMaxPhoton.sample(ekin); }

    // samples (reduced) energy transfer in Compton scattering according to the
    // Klein-Nishina DCS at the given `egamma` primary gamma energy. The three
    // additional input arguments are uniformly random values on [0,1].
    //
    // NOTE: it is assumed that: gamma-cut < egamma < E_max
    constexpr real_type kleinNishina(real_type ekin, real_type rndm1, real_type rndm2, real_type rndm3) const {
        assert(constants::k_GammaCut < ekin);
        assert(constants::k_MaxEkin >= ekin);
        return m_KleinNishina.sample_energy_transfer(ekin, rndm1, rndm2, rndm3);
    }

   protected:
    // global maximum of the total IMFP over the entire geometry (i.e. across the
    // different materials) at each individual discrete photon energy in [1/mm]
    // units.
    PhotonLinearInterpolator m_IMFPMaxPhoton;
    // sampling tables to provide rejection free sampling of the (reduced) energy
    // transferred to the secondary electron in Compton scattering according to
    // the Klein-Nishina DCS
    KleinNishina m_KleinNishina;
};
}  // namespace opmc
#endif  // ODPM_CPU_SRC_PHOTON_TABLES_H_
