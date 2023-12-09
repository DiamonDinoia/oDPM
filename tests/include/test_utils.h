//
// Created by mbarbone on 5/26/22.
//

#ifndef ODPM_CPU_TESTS_INCLUDE_TEST_UTILS_H_
#define ODPM_CPU_TESTS_INCLUDE_TEST_UTILS_H_

#include <Track.hh>
#include <ostream>

#include "materials.h"
#include "photon_tables.h"

namespace opmc {
    class TestMaterial : public Material {
    public:
        TestMaterial(const ElectronData &electronData, const PhotonData &photonData,
                     const std::string &name,
                     const int id) : Material(electronData, photonData, name, id) {}

        const auto &GetITr1MFPElastic() { return m_ITr1MFPElastic; }

        const auto &GetMaxScatStrength() { return m_MaxScatStrength; }

        const auto &GetIMFPMoller() { return m_IMFPMoller; }

        const auto &GetIMFPBrem() { return m_IMFPBrem; }

        const auto &GetStoppingPower() { return m_StoppingPower; }

        const auto &GetIMFPTotalPhoton() { return m_IMFPTotalPhoton; }

        const auto &GetIMFPComptonPhoton() { return m_IMFPComptonPhoton; }

        const auto &GetIMFPPairProdPhoton() { return m_IMFPPairProdPhoton; }
    };

    class TestPhotonData : public ODPMPhotonData {
    public:
        TestPhotonData(const ODPMPhotonData &data) : ODPMPhotonData(data) {}

        const auto &GetIMFPMaxPhoton() { return m_IMFPMaxPhoton; }

        const auto &GetKleinNishina() { return m_KleinNishina; }
    };

    class TestElectronData : public ODPMElectronData {
    public:
        TestElectronData(const ODPMElectronData &data) : ODPMElectronData(data) {}

        const auto &GetMollerEnergyTransfer() { return m_mollerEnergyTransfer; }
    };

    struct ParticleInterface {
        real_type pos_x, pos_y, pos_z, dir_x, dir_y, dir_z, energy;

        constexpr void fromParticle(const Particle &particle) {
            pos_x = particle.position.x;
            pos_y = particle.position.y;
            pos_z = particle.position.z;
            dir_x = particle.direction.x;
            dir_y = particle.direction.y;
            dir_z = particle.direction.z;
            energy = particle.energy;
        }

        constexpr void toParticle(Particle &particle) const {
            particle.position.x = pos_x;
            particle.position.y = pos_y;
            particle.position.z = pos_z;
            particle.direction.x = dir_x;
            particle.direction.y = dir_y;
            particle.direction.z = dir_z;
            particle.energy = energy;
        }

        constexpr void toTrack(Track &primaryTrack) const {
            primaryTrack.fType = -1;
            primaryTrack.fEkin = energy;
            primaryTrack.fPosition[0] = pos_x;
            primaryTrack.fPosition[1] = pos_y;
            primaryTrack.fPosition[2] = pos_z;
            primaryTrack.fDirection[0] = dir_x;
            primaryTrack.fDirection[1] = dir_y;
            primaryTrack.fDirection[2] = dir_z;
        }

        friend std::ostream &operator<<(std::ostream &os, const ParticleInterface &an_interface) {
            os << "pos_x: " << an_interface.pos_x << " pos_y: " << an_interface.pos_y << " pos_z: "
               << an_interface.pos_z
               << " dir_x: " << an_interface.dir_x << " dir_y: " << an_interface.dir_y << " dir_z: "
               << an_interface.dir_z
               << " energy: " << an_interface.energy;
            return os;
        }
    };

}  // namespace opmc

#endif  // ODPM_CPU_TESTS_INCLUDE_TEST_UTILS_H_
