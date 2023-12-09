//
// Created by mbarbone on 8/13/21.
//

#include <Geom.hh>
#include <SimDPMLike.hh>
#include <SimElectronData.hh>
#include <SimMaterialData.hh>
#include <Track.hh>
#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>

#include "constants.h"
#include "geometry.h"
#include "physics.h"
#include "random.h"
#include "initialise_tables.h"

using namespace opmc;
using Catch::Approx;

static const std::vector<std::string> test_materials{
        "G4_WATER",
};


TEST_CASE("Electron Tracking") {
    auto seed = std::random_device()();
    INFO("SEED " << seed);
    Random rng{seed};
    std::default_random_engine generator(rng.getUniform());

    static const DPMTables tables{test_materials};
    static const Material reference_material{tables.getElectronData(), tables.getPhotonData(),
                                tables.getMaterialName(0), 0};
    static const PhotonData &photon_data = tables.getPhotonData();
    static const ElectronData &electron_data = tables.getElectronData();

    SimMaterialData refMData;
    SimElectronData refElData;
    refMData.Load(OUTPUT_DATA_DIR);
    refElData.Load(OUTPUT_DATA_DIR);

    unsigned long event_distribution[7] = {0};

    const auto tests = 16384;
    const auto voxel_size = 2.;
    const auto geom_size = 200;
    const auto epsilon = 10e-5;

    ThreeVector<real_type> resolution(voxel_size, voxel_size, voxel_size);
    Geom geom(voxel_size, &refMData);
    HalfDistanceVoxelCube voxel_map{{geom_size, geom_size, geom_size * 10}, resolution, tables.populateMaterials()};
    DummyQueue queue;

    std::uniform_real_distribution<real_type> energy_distribution(opmc::constants::k_ElectronCut,
                                                                  opmc::constants::k_MaxEkin);
    std::uniform_real_distribution<real_type> position_distribution(-voxel_size / 2, voxel_size / 2);
    std::uniform_real_distribution<real_type> direction_distribution(-1, 1);

    Physics physics{rng, voxel_map, queue, reference_material, photon_data, electron_data};

    auto getEnergy = [&generator, &energy_distribution] { return energy_distribution(generator); };
    auto getPosition = [&generator, &position_distribution] { return position_distribution(generator); };
    auto getDirection = [&generator, &direction_distribution] { return direction_distribution(generator); };

    SECTION("EVENT") {

        for (int i = 0; i < tests; ++i) {
            INFO("TEST " << i);

            const auto energy = getEnergy();
            const real_type position[3] = {getPosition(), getPosition(), getPosition()};
            const real_type direction[3] = {getDirection(), getDirection(), getDirection()};
            ThreeVector<real_type> position_vector{position[0], position[1], position[2]};
            ThreeVector<real_type> direction_vector{direction[0], direction[1], direction[2]};

            INFO("Position " << position_vector);
            INFO("Direction " << direction_vector);
            INFO("Energy " << energy);

            Electron electron{position_vector, direction_vector, energy};
            ParticleState state = physics.initializeTracking(energy);

            real_type numTr1MFP = state.num_tr_1_mfp;
            real_type numMollerMFP = state.num_moller_mfp;
            real_type invMollerMFP = state.inv_moller_mfp;
            real_type numBremMFP = state.num_brem_mfp;

            INFO("Initial state " << numTr1MFP << " " << numMollerMFP << " " << invMollerMFP << " " << numBremMFP);

            auto result = physics.continuousElectronTracking(electron, state);
            event_distribution[result] += 1;

            INFO("DISTRIBUTION " << event_distribution[0] << " " << event_distribution[1] << " "
                                 << event_distribution[2] << " " << event_distribution[3] << " "
                                 << event_distribution[4] << " " << event_distribution[5] << " "
                                 << event_distribution[6]);

            REQUIRE(result != opmc::NOTHING);

            Track primaryTrack{};
            primaryTrack.fEkin = energy;
            primaryTrack.fType = -1;
            primaryTrack.fDirection[0] = direction[0];
            primaryTrack.fDirection[1] = direction[1];
            primaryTrack.fDirection[2] = direction[2];
            primaryTrack.fPosition[0] = position[0];
            primaryTrack.fPosition[1] = position[1];
            primaryTrack.fPosition[2] = position[2];

            // These are extra events not contained in the reference implementation
            if (result == opmc::MICROSTEP || result == opmc::ABSORBED) {
                continue;
            }

            Event reference = static_cast<Event>(KeepTrackingElectron(refElData, refMData, geom, primaryTrack,
                                                                      numTr1MFP, numMollerMFP,
                                                                      invMollerMFP, numBremMFP));

            if (result == opmc::OUT_OF_BOUNDARY || reference == opmc::OUT_OF_BOUNDARY) {
                continue;
            }

            REQUIRE(result == reference);

            REQUIRE(electron.energy == Approx(primaryTrack.fEkin));

            REQUIRE(electron.position.x == Approx(primaryTrack.fPosition[0]).epsilon(epsilon));
            REQUIRE(electron.position.y == Approx(primaryTrack.fPosition[1]).epsilon(epsilon));
            REQUIRE(electron.position.z == Approx(primaryTrack.fPosition[2]).epsilon(epsilon));

            REQUIRE(electron.direction.x == Approx(primaryTrack.fDirection[0]).epsilon(epsilon));
            REQUIRE(electron.direction.y == Approx(primaryTrack.fDirection[1]).epsilon(epsilon));
            REQUIRE(electron.direction.z == Approx(primaryTrack.fDirection[2]).epsilon(epsilon));
        }
    }

    std::cout << "DISTRIBUTION " << event_distribution[0] << " " << event_distribution[1] << " "
              << event_distribution[2] << " " << event_distribution[3] << " " << event_distribution[4] << " "
              << event_distribution[5] << " " << event_distribution[6] << std::endl;
}
