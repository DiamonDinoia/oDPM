//
// Copyright M. Novak: 2021
// Modified by M.Barbone 2021
//

#ifndef ODPM_MATERIALS_H
#define ODPM_MATERIALS_H

#include "constants.h"
#include "electron_tables.h"
#include "goudsmit_saunderson.h"
#include "particles.h"
#include "photon_tables.h"
#include "spline_interpolator.h"
#include "types.h"
#include "walkers_approximator.h"

class ElectronData;
class PhotonData;

namespace opmc {

class Material {
   public:
    Material(const ElectronData &electronData, const PhotonData &photonData, const std::string &name, const int id);

    virtual ~Material()            = default;

    Material(const Material &)     = default;

    Material(Material &&) noexcept = default;

    // the inverse Tr1-MFP in [1/mm] [cm3/g] scalled units
    ODPM_INLINE ODPM_PURE constexpr real_type iTr1MFPElasticElectron(real_type ekin) const noexcept {
        return m_ITr1MFPElastic.sample(ekin);
    }
    // the maximum allowed scattering strength (K_1(E)) in a single MSC step in the
    // reference material approximated as `S_max(E)/tr1-mfp(E)` (no units) where
    // S_max(E) is a sigmoid-like function of the maximum allowed MSC step length
    // determined by the `slow` and `shigh` parameters with the smooth transition
    // around `ecross`.
    ODPM_INLINE ODPM_PURE constexpr real_type maxScatteringStrengthElectron(real_type ekin) const noexcept {
        return m_MaxScatStrength.sample(ekin);
    }

    // the inverse MFP in [1/mm] [cm3/g] scaled units
    ODPM_INLINE ODPM_PURE constexpr real_type iMFPMollerElectron(real_type ekin) const noexcept {
        // FIXME: std::max should not be here if the table is set up properly
        return math::max(m_IMFPMoller.sample(ekin), constants::k_lowInvMFPVacuum);
    }
    // the inverse MFP in [1/mm] [cm3/g] scaled units
    ODPM_INLINE ODPM_PURE constexpr real_type iMFPBremElectron(real_type ekin) const noexcept {
        return m_IMFPBrem.sample(ekin);
    }
    // the stopping power in [MeV/mm] [cm3/g] scaled units
    ODPM_INLINE ODPM_PURE constexpr real_type stoppingPowerElectron(real_type ekin) const noexcept {
        return m_StoppingPower.sample(ekin);
    }

    // samples values of energy transferred to the emitted gamma according to the
    // Seltzer-Berger DCS at the given `eprim` primary electron energy and material
    // specified by its `imat` index. The three additional input arguments are
    // uniformly random values on [0,1].
    //
    // NOTE: it is assumed that: gamma-cut < eprim < E_max
    ODPM_INLINE ODPM_PURE constexpr real_type seltzerBergerElectron(real_type ekin, real_type rng1, real_type rng2,
                                                                    real_type rng3) const noexcept {
        assert(constants::k_GammaCut < ekin);
        assert(constants::k_MaxEkin >= ekin);
        ODPM_ASSUME(constants::k_GammaCut < ekin);
        ODPM_ASSUME(constants::k_MaxEkin >= ekin);
        return m_SeltzerBerger.energyTransfer(ekin, rng1, rng2, rng3);
    }
    // samples cosine of the angular deflection at the given `ekin` electron energy
    // (travelling the predefined path length in the predefined referecne material).
    // The two additional input arguments are uniformly random values on [0,1].
    //
    // NOTE: it is assumed that the `eprim` electron energy: electron-cut < eprim <E_max
    ODPM_INLINE ODPM_PURE constexpr real_type angularDeflectionElectron(real_type ekin, real_type rng1,
                                                                        real_type rng2) const {
        assert(constants::k_ElectronCut < ekin);
        assert(constants::k_MaxEkin >= ekin);
        ODPM_ASSUME(constants::k_ElectronCut < ekin);
        ODPM_ASSUME(constants::k_MaxEkin >= ekin);
        return m_GoudsmitSaunderson.angularDeflection(ekin, rng1, rng2);
    }
    // the inverse IMFP in [1/mm] [cm3/g] scalled units
    ODPM_INLINE ODPM_PURE constexpr real_type iMFPTotalPhoton(real_type ekin) const noexcept {
        return m_IMFPTotalPhoton.sample(ekin);
    }
    // the inverse IMFP in [1/mm] [cm3/g] scalled units
    ODPM_INLINE ODPM_PURE constexpr real_type iMFPComptonPhoton(real_type ekin) const noexcept {
        return m_IMFPComptonPhoton.sample(ekin);
    }
    // the inverse IMFP in [1/mm] [cm3/g] scalled units
    ODPM_INLINE ODPM_PURE constexpr real_type iMFPPairProdPhoton(real_type ekin) const noexcept {
        return m_IMFPPairProdPhoton.sample(ekin);
    }
    // density [g/cm3]
    ODPM_INLINE ODPM_CONST constexpr real_type mass_density() const noexcept { return m_density; }
    // [Z/w]/[Z/w]_ref
    ODPM_INLINE ODPM_CONST constexpr real_type mollerIMFPScaling() const noexcept { return m_MollerIMFPScaling; }

    constexpr unsigned short id() const noexcept { return k_id; };

    const std::string &getName() const;

   protected:
    const unsigned short k_id;
    const std::string name;

    // MollerIMFPScaling
    const real_type m_MollerIMFPScaling;
    const real_type m_density; //[g/cm3]

    // Electron data

    // inverse first transprt mean free path for elastic interaction only for
    // the reference material and scalled (divided) by the material density in [g/cm3]
    // So it will have the proper [1/mm] units only after scaling back (multiplying)
    // by a given material density in [g/cm3] units.
    MaterialSplineInterpolator m_ITr1MFPElastic;
    // maximum scattering strength (K_1(E)) in MSC step allowed in the reference
    // material that is approximated as K_1(E) ~ S_max(E)/tr1-mfp(E) with S_max(E)
    // being the max-MSC-step length determined by the `slow`, `shigh` and `ecross`
    // parameters.
    // Note, data are loaded only for the reference material and not scalled by
    // the density (i.e. it has no units).
    MaterialSplineInterpolator m_MaxScatStrength;
    // restriced IMFP (macroscopic cross section) for Moller interaction only for
    // the reference material and scalled (divided) by the material density in [g/cm3]
    // So it will have the proper [1/mm] units only after scaling back (multiplying)
    // by a given material density in [g/cm3] units.
    MaterialSplineInterpolator m_IMFPMoller;
    // restriced IMFP (macroscopic cross section) for Brem. interaction for each of
    // the individual material and scalled (divided) by the material density in [g/cm3]
    // So it will have the proper [1/mm] units only after scaling back (multiplying)
    // by a given material density in [g/cm3] units.
    MaterialSplineInterpolator m_IMFPBrem;
    // restriced total (radiative plus collisonal) stopping power for each of the
    // individual material and scalled (divided) by the material density in [g/cm3]
    // So it will have the proper [MeV/mm] units only after scaling back (multiplying)
    // by a given material density in [g/cm3] units.
    MaterialSplineInterpolator m_StoppingPower;
    // sampling tables to provide rejection free sampling of the energy that is
    // transferred to the secondary gamma in electron bremsstrahlung according to
    // the Seltzer-Berger DCS
    SeltzerBerger m_SeltzerBerger;
    // sampling tables to provide rejection free sampling of the angular deflection
    // due to multiple Coulomb scattering in the reference material along the
    // maximally allowed MSC step length that is a function of the energy and
    // determined by a (sigmoid-like) function of the `slow`, `shigh` and `ecross`
    // parameters
    GoudsmitSaunderson m_GoudsmitSaunderson;
    // Photon data

    // total (sum of Compton, Pair-production and Photoelectric) inverse MFP for
    // each of the individual material and scalled (divided) by the material density
    // in [g/cm3]. So it will have the proper [1/mm] units only after scaling back
    // (multiplying) by a given (voxel) material density in [g/cm3] units.
    PhotonLinearInterpolator m_IMFPTotalPhoton;
    // IMFP for Compton scattering for each of the individual material and scalled
    // (divided) by the material density in [g/cm3]. So it will have the proper
    // [1/mm] units only after scaling back (multiplying) by a given (voxel) material
    // density in [g/cm3] units.
    PhotonLinearInterpolator m_IMFPComptonPhoton;
    // IMFP for Pair-production for each of the individual material and scalled
    // (divided) by the material density in [g/cm3]. So it will have the proper
    // [1/mm] units only after scaling back (multiplying) by a given (voxel) material
    // density in [g/cm3] units.
    PhotonLinearInterpolator m_IMFPPairProdPhoton;
};
}  // namespace opmc
#endif  // ODPM_MATERIALS_H
