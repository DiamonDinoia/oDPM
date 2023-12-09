//
// Created by mbarbone on 10/1/21.
//

#define CATCH_CONFIG_MAIN

#include <Geom.hh>
#include <SimDPMLike.hh>
#include <SimMaterialData.hh>
#include <SimPhotonData.hh>
#include <Track.hh>
#include <catch2/catch_test_macros.hpp>
#include <iomanip>

#include "beams.h"
#include "constants.h"
#include "ksTest.h"
#include "physics.h"
#include "random.h"
#include "initialise_tables.h"

using namespace std;
using namespace opmc;

static const std::vector<std::string> test_materials{
        "G4_WATER",
};

DPMTables tables{test_materials};
Material reference_material{tables.getElectronData(), tables.getPhotonData(), tables.getMaterialName(0),
                            0};
PhotonData photon_data{tables.getPhotonData()};
ElectronData electron_data{tables.getElectronData()};
SimMaterialData refMData;
SimPhotonData refPhData;
bool init = true;

TEST_CASE("Photon simulation") {
    const auto voxel_size = 2;
    const auto geom_size = 200;
    const auto histories = 10e6;
    const auto alpha = 5.0E-9;
    const auto seed = std::random_device()();
    const auto treshold = 0.00;
    INFO("SEED " << seed);
    Random rng{seed};

    if (init) {
        refMData.Load(OUTPUT_DATA_DIR);
        refPhData.Load(OUTPUT_DATA_DIR);
        init = false;
    }
    for (int test_energy = 1; test_energy < 21; test_energy++) {
        SECTION("Energy " + std::to_string(test_energy) + "Mev") {
            const real_type primaryEnergy = test_energy;
            Geom geom(voxel_size, &refMData);
            Track primaryTrack;
            for (int i = 0; i < histories; ++i) {
                primaryTrack.fEkin = primaryEnergy;  // Primary kinetic energy in [MeV]
                primaryTrack.fType = 0;              // e-(-1) or photon(0)

                primaryTrack.fDirection[0] = 0.0;  // initial direction is [0,0,1]
                primaryTrack.fDirection[1] = 0.0;
                primaryTrack.fDirection[2] = 1.0;

                primaryTrack.fPosition[0] = 0.0;  // initial position is [0,0, theRZ0]
                primaryTrack.fPosition[1] = 0.0;
                primaryTrack.fPosition[2] = -0.5 * voxel_size;
                KeepTrackingPhoton(refPhData, refMData, geom, primaryTrack);
            }

            auto reference = geom.histogram();
            reference.erase(
                    std::remove_if(reference.begin(), reference.end(),
                                   [treshold](real_type x) { return x <= treshold; }),
                    reference.end());


            real_type norm = 1. / *max_element(reference.begin(), reference.end());
            for (auto &item: reference) {
                item *= norm;
            }


            ThreeVector<real_type> resolution(voxel_size, voxel_size, voxel_size);
            HalfDistanceVoxelCube voxel_map{{geom_size, geom_size, reference.size()}, resolution,
                                            tables.populateMaterials()};
            voxel_map.setDefaultMaterial(0);

            DummyQueue queue;
            Physics physics{rng, voxel_map, queue, reference_material, photon_data, electron_data};

            const ThreeVector<real_type> initialPosition{0, 0, -0.5 * voxel_size};
            const ThreeVector<real_type> initialDirection{0, 0, 1};

            PhotonPencilBeam beam{initialPosition, initialDirection, primaryEnergy};

            for (int i = 0; i < histories; ++i) {
                auto particle = beam.generateParticle(rng);
                physics.simulate(particle);
            }

            auto result = voxel_map.doseDepthDistribution();
            std::cout << std::setprecision(10) << std::endl;
            reference.pop_back();
            result.pop_back();

            for (auto &item: result) {
                std::cout << item << " ";
            }
            std::cout << std::endl;

            for (auto &item: reference) {
                std::cout << item << " ";
            }
            std::cout << std::endl;
            REQUIRE(opmc::KsTest(reference, result) >= alpha);
        }
    }
}