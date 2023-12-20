//
// Copyright M. Novak: 2021
// Modified by M.Barbone 2021
//

#ifndef ODPM_CPU_INCLUDE_PHYSICS_H_
#define ODPM_CPU_INCLUDE_PHYSICS_H_

#include "geometry.h"
#include "materials.h"
#include "particle_queue.h"
#include "particles.h"
#include "random.h"
#include "types.h"

namespace opmc {

// Possible interaction events saved as an enum instead of string or integer for clarity and debug reasons.
enum Event { NOTHING = 0, BREMSSTRAHLUNG, MOLLER, HINGE, OUT_OF_BOUNDARY, ABSORBED, MICROSTEP };

/**
 * Stateful wrapper for the physics.
 * It has the functionality to simulate particles {Electrons, Photons, Positrons} until they run out of energy ot escape
 * the geometry.
 */
template <typename T, typename V>
class Physics {
   public:
    explicit Physics(Random &rng, T &geometry, V &stack, const Material &reference_material,
                     const ODPMPhotonData &photon_data, const ODPMElectronData &electron_data)
        : r_rng(rng),
          r_geometry(geometry),
          r_stack(stack),
          r_reference_material(reference_material),
          k_photonData(photon_data),
          k_electronData(electron_data){};

    /**
     * Samples cos(theta) according to the G&S distribution.
     * Uses interpolated data for bw and the q surface and this
     * latter quantity to perform a rejection procedure.
     * NOTE: that angular deflection is applied beween the pre- and post-step point
     * while it needs to be computed at the e- energy at the pre-step point.
     * The pre-step point kinetic energy is `initialEnergy` while the current kinetic
     * energy is `particle.kineticEnergy`.
     * @param material material of the voxel
     * @param particle particle to simulate
     */
    ODPM_INLINE void hingeInteraction(Electron &particle, real_type initialEnergy) noexcept {
        const real_type polarAngle =
            r_reference_material.angularDeflectionElectron(initialEnergy, r_rng.getUniform(), r_rng.getUniform());
        // smaple \phi: uniform in [0,2Pi] <== spherical symmetry of the scattering potential
        // rotate new direction from the scattering to the lab frame
        particle.direction.rotate(polarAngle, samplePhi());
    }

    /**
     * Creates a new secondary electron from a Moller interaction
     * @param primaryElectron electron originating the event
     * @return the secondary electron
     * It is assumed that track.fEkin > 2*electron-cut!
     * (Interaction is not possible otherwise)
     */
    ODPM_INLINE Electron mollerInteraction(Electron &primaryElectron) noexcept {
        const auto primary_energy        = primaryElectron.energy;
        const real_type secondary_energy = k_electronData.mollerEnergyTransfer(
            primaryElectron.energy, r_rng.getUniform(), r_rng.getUniform(), r_rng.getUniform());
        Electron secondElectron{primaryElectron.position, primaryElectron.direction, secondary_energy};

        const real_type cosTheta = math::sqrt(secondary_energy * (primary_energy + constants::k_2EMC2) /
                                              (primary_energy * (secondary_energy + constants::k_2EMC2)));

        const real_type phi      = samplePhi();

        secondElectron.direction.rotate(cosTheta, phi);
        // decrease primary energy: DMP do not deflect the primary
        primaryElectron.energy -= secondary_energy;
        return secondElectron;
    }

    /**
     * Creates a new secondary bremsstrahlung photon
     * @param material material of the voxel
     * @param electron the electron originating  the event
     * @return the photon if generated
     * It is assumed that track.fEkin > gamma-cut!
     * (Interaction is not possible otherwise)
     */
    ODPM_INLINE Photon bremsstrahlungInteraction(const Material &material, Electron &electron) noexcept {
        // sample energy transferred to the emitted gamma photon
        const real_type energyLoss =
            material.seltzerBergerElectron(electron.energy, r_rng.getUniform(), r_rng.getUniform(), r_rng.getUniform());
        // compute emission direction (rough approximation in DPM by the mean)
        // and no deflection of the primary e-
        const real_type dum0      = constants::k_HalfSqrt2EMC2 / (electron.energy + constants::k_EMC2);
        const real_type cos_theta = 1.0 - dum0 * dum0;
        const real_type phi       = samplePhi();
        auto direction            = electron.direction;
        direction.rotate(math::cos(cos_theta), phi);
        // decrease the primary energy:
        electron.energy -= energyLoss;
        // return the secondary particle
        return {electron.position, direction, energyLoss};
    }

    /**
     * Pair photon production
     * @param positron
     * @return
     */
    ODPM_INLINE Pair<Photon, Photon> annihilation(const Positron &positron) noexcept {
        // isotropic direction
        const real_type cost = 1.0 - 2.0 * r_rng.getUniform();
        const real_type sint = math::sqrt(1.0 - cost * cost);
        const real_type phi  = samplePhi();
        const real_type rx   = sint * math::cos(phi);
        const real_type ry   = sint * math::sin(phi);
        const real_type rz   = cost;
        return {{positron.position, {rx, ry, rz}, constants::k_EMC2},
                {positron.position, {-rx, -ry, -rz}, constants::k_EMC2}};
    }

    /**
     * Creates a secondary electron from a Compton interaction
     * @param primaryParticle the photon originating the scattering event
     * @return the electron generated
     */
    ODPM_INLINE Electron comptonScattering(Photon &primaryParticle) noexcept {
        const real_type theEps = k_photonData.kleinNishina(primaryParticle.energy, r_rng.getUniform(),
                                                           r_rng.getUniform(), r_rng.getUniform());
        const real_type kappa  = primaryParticle.energy * constants::k_InvEMC2;
        const real_type phCost = 1.0 - (1.0 - theEps) / (theEps * kappa);  // 1- (1-cost)
        const real_type phEner = theEps * primaryParticle.energy;
        const real_type phPhi  = samplePhi();
        const real_type elEner = primaryParticle.energy - phEner;
        Electron electron{primaryParticle.position, primaryParticle.direction, elEner};
        // This is needed only if the particle will be simulated, otherwise just
        // score the energy
        if (elEner >= constants::k_ElectronCut) {
            const real_type e0     = primaryParticle.energy * constants::k_InvEMC2;
            const real_type elCost = (1.0 + e0) * math::sqrt((1.0 - theEps) / (e0 * (2.0 + e0 * (1.0 - theEps))));
            const real_type phi    = phPhi + constants::k_PI;
            electron.direction.rotate(elCost, phi);
        }
        // This is needed only if the particle will be simulated, otherwise just
        // score the energy
        primaryParticle.energy = phEner;
        if (phEner >= constants::k_GammaCut) {
            primaryParticle.direction.rotate(phCost, phPhi);
        }
        // Given this implementation the main simulation loop should check the
        // energy and score it if necessary
        return electron;
    }

    /** Performs one tracking iteration it is possible that no interaction happens and the particle is:
     * 1. advanced
     * 2. Continuous energy loos is accounted
     */
    ODPM_INLINE Event electronTracking(Electron &electron, ParticleState &state) noexcept {
        Event event{NOTHING};
        const auto position = r_geometry[electron.position];
        // the geometry returns an optional. If the particle escaped the optional is empty
        if (ODPM_UNLIKELY(!position)) {
            electron.energy = 0;
            return OUT_OF_BOUNDARY;
        }

        // compute the distance to boundary: this will be the current value of the maximal step length
        // init the current step lenght to that of the distance to boundary: we might
        // or might not go that far (depending on what the 3 physics tells) but for
        // sure that is the maximum step length becasue 'boundary crossing interaction'
        // happens after travelling that far. So set what happen to that (i.e. = 0.)
        real_type stepLength = r_geometry.distanceToBoundary(electron);
        // if the stewleng is < than 10e4 (predefined value) then this is a microstep, the particle is advanced and the
        // tracking starts over
        if (stepLength < constants::k_minStepLength) {
            return MICROSTEP;
        }
        // The Moller mfp energy dependence is not considered in DPM so we do the same:
        // we will keep using the Moller mfp evaluated at the initial energy for the reference
        // material. Only the material scaling will be applied below:
        // \lam' = \lam_ref(E') [Z\rho/A]_ref [A/(Z\rho)]_actual (for 1/\lam' to compute delta #mfp)
        Voxel voxel = *position;

        if (electron.energy < constants::k_ElectronCut) {
            voxel.addDose(electron.energy);
            electron.energy = 0;
            return ABSORBED;
        }

        // the material scaling factor for the Moller inverse-mf: [A/(Z\rho/)]_ref [(Z\rho)/A]_actual
        // or more exactly its [A/Z)]_ref [(Z)/A]_actual part
        const real_type scalMolMFP = voxel.material.mollerIMFPScaling();
        // WE ASSUME HERE NOW THAT EACH VOXEL IS A CLEAR MATERIAL SO WE WILL
        // USE theVoxelMaterialDensity = theVoxelBaseMaterialDensity. HOWEVER,
        // THIS MIGHT BE CHANGED ANYTIME WHEN THE GEOMETRY CAN PROVIDE A SUITABLE
        // VOXEL MATERIAL DENSITY.
        //
        // NOTE: density must be in g/cm3 units !!!!
        const real_type massDensity = voxel.material.mass_density();
        //
        // Here we compute the decrese of the #mfp/#tr1-mfp for the 3 interaction with the current,
        // maximum step length (i.e. the distance to boundary): #mfp' = #mfp - ds/mfp' or - ds/tr1_mfp
        // for MSC (where ' indicates values at the end point)
        //
        // compute the mid-point energy along this step by assuming:
        // - constant dEdx along the step, i.e. dEdx=dEdx(E_0) and dE = s dEdx --> E_mid = E_0 - 0.5 s dEdx
        // - the step equal to the current one, i.e. `stepLength` (dist. to boundary)
        // the restricted stopping power for this material: for the referecne material and scalled with the current
        // density
        real_type theHalfDEDX = -voxel.material.stoppingPowerElectron(electron.energy) * massDensity * 0.5;
        // make sure that do not go below the minim e- energy
        real_type midStepE = math::max(std::fma(stepLength, theHalfDEDX, electron.energy), constants::k_ElectronCut);
        // elastic: #tr1-mfp' = #tr1-mfp - ds/tr1-mfp' so the change in #tr1-mfp is ds/tr1-mfp' and
        //          1/mfp' is computed here
        real_type delNumTr1MFP = voxel.material.iTr1MFPElasticElectron(midStepE) * massDensity;
        // moller: see above the details
        const real_type delNumMollerMFP = state.inv_moller_mfp * scalMolMFP * massDensity;
        // brem: #mfp' = #mfp - ds/mfp' with mfp = brem_mfp so the change in #mfp is ds/mfp' and
        //       1/mfp' is computed here
        const real_type delNumBremMFP = voxel.material.iMFPBremElectron(midStepE) * massDensity;
        //
        //
        // Now we will see how far actually we go by trying to decrese each of the 3 #mfp/#tr1-mfp
        // by the changes in the number of mfp/tr1-mfp computed above as `delNum` :
        // - if we could manage to decrese all the 3 #mfp/tr1-mfp such that they are > 0
        //   then actually we reached the boundary: we cross and perform an other step.
        // - if any of the 3 #mfp/#tr1-mfp goes down to zero, then the lowest (i.e. shortest
        //   path to) will be considered to happen: the given inetraction need to beinvoked
        // In all cases, the energy loss along the given step needs to be computed!
        //
        // In DMP, the #mfp are tried to decresed with the current step length (that
        // is always the current shortest) as #mfp: n' = n - ds/mfp' then corrected
        // back if goes below zero, etc...
        // We compute the step length ds = n x mfp' (or n x tr1-mfp') i.e. that
        // would bring the given number #mfp/#tr1-mfp down to 0 (since n'=n-ds/mfp' or ds/tr1-mfp).
        // If this step is the shortest then we take this as current step length and we determine
        // the overal shortest and the corresponding interaction will happen (including boundary crossing).
        //
        double stepBrem    = state.num_brem_mfp / delNumBremMFP;
        double stepElastic = state.num_tr_1_mfp / delNumTr1MFP;
        double stepMoller  = state.num_moller_mfp / delNumMollerMFP;

        if (stepBrem < stepLength) {
            // discrete bremsstrahlung (might) happen before reaching the boundary:
            stepLength = stepBrem;
            event      = BREMSSTRAHLUNG;
        }

        if (stepMoller < stepLength) {
            // discrete Moller (might) happen even before bremsstrahlung:
            stepLength = stepMoller;
            event      = MOLLER;
        }

        if (stepElastic < stepLength) {
            // elastic interaction happens:
            // - either the hinge: sample and apply deflection and update numElMFP to the
            //                     remaining part
            // - end of 2nd step : nothing to do just resample the #mfp left since all
            //                     has been eaten up by the step length travelled so far
            // Before anything, refine the computation of the mfp (i.e. the first transport
            // mean free path in case of elastic) regarding its energy dependence.
            // NOTE: that the 1/mfp values were computed at the mid-point energy (brem
            //       and elastic since Moller is assumed to be constant), assuming that
            //       the geometry step will be taken and the dEdx is constant along this
            //       step (i.e. no energy dependence).
            //       Here we know that actually not the geometry step, but the stepElastic
            //       is taken since that is the shortest. So we recompute the mid-step-point
            //       energy according to the step length of stepElastic and re-evaluate
            //       the 1./mfp i.e. 1/tr1mfp at this energy value
            //        delNumTr1MFP = elData.ITr1MFPElastic.GetITr1MFPPerDensity(electron.kineticEnergy, materialID) *
            //        massDensity;
            delNumTr1MFP = voxel.material.iTr1MFPElasticElectron(electron.energy) * massDensity;
            //        assert(delNumTr1MFP == new_delNumTr1MFP);
            stepElastic  = state.num_tr_1_mfp / delNumTr1MFP;
            midStepE     = math::max(std::fma(stepElastic, theHalfDEDX, electron.energy), constants::k_ElectronCut);
            delNumTr1MFP = voxel.material.iTr1MFPElasticElectron(midStepE) * massDensity;
            // don't let longer than the original in order to make sure that it
            // is still the minimum of all step lengths
            stepElastic = state.num_tr_1_mfp / delNumTr1MFP;
            // don't let longer than the original in order to make sure that it is still the
            // minimum of all step lengths
            if (ODPM_LIKELY(stepElastic < stepLength)) {
                stepLength = stepElastic;
                event      = HINGE;
            }
        }
        state.num_brem_mfp   = (event == BREMSSTRAHLUNG) ? 0 : state.num_brem_mfp - stepLength * delNumBremMFP;
        state.num_moller_mfp = (event == MOLLER) ? 0 : state.num_moller_mfp - stepLength * delNumMollerMFP;
        state.num_tr_1_mfp   = (event == HINGE) ? 0 : state.num_tr_1_mfp - stepLength * delNumTr1MFP;
        //
        // Compte the (sub-treshold, i.e. along step) energy loss:
        // - first the mid-step energy using the final value of the step lenght and the
        //   pre-step point dEdx (assumed to be constant along the step).
        midStepE = math::max(std::fma(stepLength, theHalfDEDX, electron.energy), constants::k_ElectronCut);
        // - then the dEdx at this energy
        // - then the energy loss along the step using the mid-step dEdx (as
        // constant)
        //   and the final energy
        real_type energyLoss  = stepLength * voxel.material.stoppingPowerElectron(midStepE) * massDensity;
        real_type finalEnergy = electron.energy - energyLoss;
        // check if energy dropped below tracking cut, i.e. below seconday e- production threshold
        // NOTE: HERE THERE IS A SUB-THRESHOLD TRACKING CONDITION IN DPM BUT WE WILL NEED TO SEE THAT !!!
        // ALSO: IF THE SELECTED STEP LENGHT BRINGS THE EKIN BELOW THRESHOLD WE DON'T TRY TO FIND THE
        //       STEP LENGTH (a smaller than selected) THAT ACTUALLY BRINGS EKIN EXACTLY TO THE THRESHOLD.
        //       SO THE TRACK LENGTH IS NOT PRECISE, SINCE WE KNOW THAT THE e- WAS STOPPED IN THIS VOLUME/BOX
        //       BUT WE DON'T COMPUTE THE EXCT POSITION
        if (finalEnergy < constants::k_ElectronCut) {
            energyLoss  = electron.energy;
            finalEnergy = 0;
            event       = ABSORBED;
        }
        // Update particle position, track length etc.
        electron.position += electron.direction * stepLength;
        electron.energy = finalEnergy;
        // Score the continuous energy loss before going back to perform the discrete interaction
        voxel.addDose(energyLoss);

        return event;
    }

    /** Keeps tracking an e- till one of the following condition is reached:
     * 1. a discrete brem interaction take place
     * 2. a discrete ioni interaction take place
     * 3. MSC hinge or end of MSC step take place
     * 4. the electron energy drops below zero so stops
     */
    ODPM_INLINE Event continuousElectronTracking(Electron &electron, ParticleState &state) noexcept {
        Event event;
        while (event = electronTracking(electron, state), event == NOTHING)
            ;
        return event;
    }
    /**
     * Initialize MFPs for the particle and saves it in State
     */
    ODPM_INLINE constexpr ParticleState initializeTracking(real_type kineticEnergy) noexcept {
        ParticleState state;
        //
        // Track e-/e+ ortherwise:
        //
        // WE ASSUME HERE NOW THAT EACH VOXEL IS A CLEAR MATERIAL SO WE WILL
        // USE theVoxelMaterialDensity = theVoxelBaseMaterialDensity. HOWEVER,
        // THIS MIGH BE CHANGED ANYTIME WHEN THE GEOMETRY CAN PROVIDE A SUITABLE
        // VOXEL MATERIAL DENSITY.
        //
        // NOTE: density must be in g/cm3 units. (cannot be vacuum at this point)
        //          double theVoxelMatDensity = geom.GetVoxelMaterialDensity(theVoxelMatIndx);
        //
        // this will be used to alter between an hinge and remaining sub MSC-step
        state.is_msc_hinge = true;
        // this will be used to keep track of the pre-step point energy that we need
        // when we reach the msc hinge point (init. is not important cause will be set below)
        state.initial_energy = kineticEnergy;
        //
        //
        // Compute the initial number of mfp left till the different interactions
        // which is -ln(R) = s/mfp or U(R) = s/tr1-mfp for the MSC hinge interaction:
        //
        // 1. elastic interaction with msc: s/tr1mfp and sample the hinge position as well
        //    NOTE: data are used for the reference material and not scaled by the density
        //          so numTr1MFP ~ s_max(E)/tr1-mfp(E) [no units]
        state.num_tr_1_mfp = r_reference_material.maxScatteringStrengthElectron(kineticEnergy);
        // travell this #tr1-mfp after the MSC-hinge took place
        state.num_tr_1_mfp_0 = r_rng.getUniform() * state.num_tr_1_mfp;
        // travel this #tr1-mfp till the MSC-hinge
        state.num_tr_1_mfp = state.num_tr_1_mfp - state.num_tr_1_mfp_0;
        //
        // 2. Moller:
        // NOTE: as in DPM, the mfp for Moller will be assumed to have no energy dependence!!!
        //       So the update of the number-of-mfp-left will contain only the material
        //       dependence related scaling (that's again approximate but DPM does this).
        //       This material dependent update(scaling) relies on the \lam ~ [A/(Z\rho)]
        //       dependence: \lam' = \lam_ref(E') [Z\rho/A]_ref [A/(Z\rho)]_actual (such
        //       that in case of Moller (in DPM) even E'=E_0 in each step as dicussed above).
        state.num_moller_mfp = r_rng.getExponential();
        // again, the reference material Moller IMFP value is used
        state.inv_moller_mfp = r_reference_material.iMFPMollerElectron(kineticEnergy);
        //
        // 3. bremsstrahlung:
        // NOTE: IMFP for bremsstrahlung is important so it's interpolated by using values for
        //       the actual material and kinetic energy in the `KeepTracking` code. Here we
        //       sample a `number-of-interaction-left` value, i.e. -ln(R) only (=s/imfp(E,mat))
        state.num_brem_mfp = r_rng.getExponential();
        return state;
    }

    /** Keeps tracking a photon till one of the following condition is reached:
     * 1. pair-production take place
     * 2. Compton scattering take place
     * 3. Photoelectric absorption take place
     * Unlike in case of electrons, this function also performs the interactions
     * themselves since the photon interactions are very simple in a DPM like simulation
     */
    ODPM_INLINE void simulate(Photon &photon) noexcept {
        while (ODPM_LIKELY(photon.energy > 0.0)) {
            simulateStep(photon);
        }
    }

    /** Performs only one step of photon simulation */
    ODPM_INLINE constexpr void simulateStep(Photon &photon) noexcept {
        // get the global max-macroscopic cross section and use it for samppling the
        // the length till the next photon intercation (that includes delta interaction
        // as well)
        const real_type globalMaxMFP = 1.0 / k_photonData.iMFPMaxPhoton(photon.energy);
        const real_type stepLength   = globalMaxMFP * r_rng.getExponential();
        // Update particle position, track length etc.
        photon.position += photon.direction * stepLength;
        // determine current voxel index
        const auto voxel_optional = r_geometry[photon.position];
        // out of bounds
        if (ODPM_UNLIKELY(!voxel_optional)) {
            photon.energy = 0;
            return;
        }
        auto voxel = *voxel_optional;
        /*
         * this is present in the original simulation, but it is not necessary.
         * Moreover, I believe that it is also wrong to update the photon position to the boundary
         * in case of micro steps.
            const auto distance = geometry.distanceToBoundary(photon);
            if (distance < distanceLowCut) {
                photon.position += photon.direction * distance;
                voxel = geometry[photon.position];
                if (!voxel) {
                    photon.kineticEnergy = 0;
                    return;
                }
            }
        */
        const real_type massDensity = voxel.material.mass_density();
        const real_type totalIMFP   = voxel.material.iMFPTotalPhoton(photon.energy) * massDensity;
        // P(no-inetaction) = 1.0-mxsecTotal/mxsecGlobalMax
        const real_type r1 = r_rng.getUniform();
        real_type theProb  = 1.0 - totalIMFP * globalMaxMFP;
        if (r1 < theProb) {
            return;  // with the same globalMaxMFP since the energy did not change !!!
        }
        //
        // otherwise: check which interaction happened P(i) = mxsec-i/mxsecTotal
        // compute cumulated probability of adding Compton prob
        theProb += voxel.material.iMFPComptonPhoton(photon.energy) * massDensity * globalMaxMFP;
        if (r1 < theProb) {
            // Compton interaction: Klein-Nishina like
            // the photon scattering angle and post-interafctin energy fraction
            auto electron = comptonScattering(photon);
            if (electron.energy < constants::k_ElectronCut) {
                voxel.addDose(electron.energy);
            } else {
                addParticle(std::move(electron));
            }
            if (photon.energy < constants::k_GammaCut) {
                voxel.addDose(photon.energy);
                photon.energy = 0;
            }
            return;
        }

        // compute cumulated probability of adding Pair-production prob
        theProb += voxel.material.iMFPPairProdPhoton(photon.energy) * massDensity * globalMaxMFP;
        if (r1 < theProb && photon.energy > 2.0 * constants::k_EMC2) {
            // Pair-production interaction:
            const real_type sumEkin = photon.energy - 2.0 * constants::k_EMC2;

            // simple uniform share of the energy between the e- and e+ going to the
            // same direction as the original photon.
            // no difference between the e- and e+ transport till the end:
            // - when the e+ stops, 2 photons are emitted
            // we will assume that e1 is the e+
            const real_type e1 = r_rng.getUniform() * sumEkin;
            const real_type e2 = sumEkin - e1;
            // insert the e- and e+ only if their energy is above the tracking cut
            // the e-
            if (e2 < constants::k_ElectronCut) {
                voxel.addDose(e2);
            } else {
                addParticle(Electron{photon.position, photon.direction, e2});
            }
            // the e+
            if (e1 < constants::k_ElectronCut) {
                voxel.addDose(e1);
                auto [p1, p2] = annihilation({photon.position, photon.direction, e1});
                addParticle(std::move(p1));
                addParticle(std::move(p1));
            } else {
                addParticle(Positron{photon.position, photon.direction, e1});
            }
            // kill the primary photon
            photon.energy = 0;
            return;
        }
        // if we are here then Photoelectric effect happens that absorbs the photon:
        // - score the current photon energy and stopp the photon
        voxel.addDose(photon.energy);
        photon.energy = 0;
    }

    ODPM_INLINE constexpr void simulate(Electron &electron) noexcept {
        auto state = initializeTracking(electron.energy);
        while (electron.energy > 0) {
            const auto event = continuousElectronTracking(electron, state);
            simulateStep(electron, state, event);
        }
    }

    /** this function is executed after tracking, it handles the interaction event, which can be:
     * 1. a discrete brem interaction
     * 2. a discrete ioni interaction
     * 3. MSC hinge or end of MSC step
     * 4. a MICROSTEP
     * When the event is handled
     */
    ODPM_INLINE constexpr void simulateStep(Electron &electron, ParticleState &state, Event event) noexcept {
        switch (event) {
            // discrete bremsstrahlung interaction should be sampled:
            //     - sample energy transfer to the photon (if any)
            case BREMSSTRAHLUNG: {
                auto voxel = *r_geometry[electron.position];
                // perform bremsstrahlung interaction but only if E0 > gcut
                if (electron.energy > constants::k_GammaCut) {
                    addParticle(bremsstrahlungInteraction(voxel.material, electron));
                }
                // check if the post-interaction electron energy dropped
                // below the tracking cut and stop tracking if yes
                if (electron.energy < constants::k_ElectronCut) {
                    // deposit the whole energy and stop tracking
                    voxel.addDose(electron.energy);
                    electron.energy = 0;
                    return;
                }
                // if the primary electron (or e+) survived, i.e. if we are here
                // then, re-sample the #mfp to travel till the next brem event
                // and break
                state.num_brem_mfp = r_rng.getExponential();
                break;
            }
                // (2) discrete Moller interaction is sampled:
                // NOTE: no energy dependence is considered in case of Moller in DPM so
                //       the numMollerMFP is evaluated only at this point and assumed to be
                //       constant along the entire `KeepTracking` part (only material scaling is applied).
                //       Furthermore, the kinetic energies of both the post interaction primary
                //       and secondly electrons are guarantied to be above the secondary electron
                //       production threshold so no need to check if their energy dropped below after
                //       the interaction
                //       Furthermore, note that Moller interaction is independent of Z
            case MOLLER: {
                // perform ionisation (Moller) intraction but only if E0 > 2cut
                if (electron.energy > constants::k_ElectronCut * 2) {
                    addParticle(mollerInteraction(electron));
                }
                // Resample #mfp left and interpolate the IMFP since the energy has been changed.
                // Again, the reference material Moller IMFP value is used
                state.num_moller_mfp = r_rng.getExponential();
                state.inv_moller_mfp = r_reference_material.iMFPMollerElectron(electron.energy);
                break;
            }
                // msc interaction happened: either hinge or just or end of an MSC step
            case HINGE: {
                if (state.is_msc_hinge) {
                    // Sample angular deflection from GS distr. and apply it
                    // -----------------------------------------------------
                    hingeInteraction(electron, state.initial_energy);
                    // -----------------------------------------------------
                    // set the #tr1-mfp left to the remaining, i.e. after hinge part
                    state.num_tr_1_mfp = state.num_tr_1_mfp_0;
                    // the end point is the next msc stop and not the hinge
                    state.is_msc_hinge = false;
                } else {
                    // end point so resample #tr1-mfp left and the hinge point
                    state.initial_energy = electron.energy;
                    // again, the reference material K_1(E) is used
                    state.num_tr_1_mfp = r_reference_material.maxScatteringStrengthElectron(state.initial_energy);
                    // travell this #tr1-mfp after the MSC-hinge took place
                    state.num_tr_1_mfp_0 = r_rng.getUniform() * state.num_tr_1_mfp;
                    // travell this #tr1-mfp till the MSC-hinge
                    state.num_tr_1_mfp -= state.num_tr_1_mfp_0;
                    // hinge will be the next msc stop
                    state.is_msc_hinge = true;
                }
                break;
            }
                // The electron is advanced to the boundary of the next voxel
            case MICROSTEP: {
                auto min_distance = constants::k_minStepLength;
                electron.position += electron.direction * min_distance;
                break;
            }
                // the kinetic energy dropped below the tracking cut
                // Fall through since in these cases only the energy needs to be set to 0 or do nothing and return.
            case ABSORBED:
            case OUT_OF_BOUNDARY:
                electron.energy = 0;
            case NOTHING:
                return;
        }
    }

    ODPM_INLINE constexpr void simulate(Positron &positron) noexcept {
        simulate(static_cast<Electron &>(positron));
        // perform annihilation (case of e+)
        auto [p1, p2] = annihilation(positron);
        addParticle(std::move(p1));
        addParticle(std::move(p2));
    }

    // Auxiliary functions to add particles to the queue.
    template <typename Q>
    ODPM_INLINE constexpr void addParticle(Q &&particle) noexcept {
        r_stack.push(std::forward<Q>(particle));
    }

   private:
    Random &r_rng;
    T &r_geometry;
    V &r_stack;
    const Material &r_reference_material;
    const ODPMPhotonData &k_photonData;
    const ODPMElectronData &k_electronData;

    ODPM_INLINE constexpr real_type samplePhi() noexcept { return r_rng.getUniform() * 2 * M_PI; }
};

// Nice printing function for Events to have nicer logs
std::ostream &operator<<(std::ostream &out, Event value);

}  // namespace opmc

#endif  // ODPM_CPU_INCLUDE_PHYSICS_H_