//
// Created by mbarbone on 5/30/22.
//
#define CATCH_CONFIG_MAIN
#include <Geom.hh>
#include <SimDPMLike.hh>
#include <SimMaterialData.hh>
#include <SimPhotonData.hh>
#include <Track.hh>
#include <catch2/catch.hpp>

#include "beams.h"
#include "constants.h"
#include "ksTest.h"
#include "physics.h"
#include "random.h"

namespace opmc {

constexpr auto voxel_size        = 2;
constexpr auto geom_size         = 200;
constexpr auto histories         = 1000192 * 4;
constexpr auto alpha             = 5.0E-8;
static constexpr auto block_size = 256;
static constexpr auto grid_size  = div_rounding_up(histories, block_size);

static ODPM_KERNEL void photonKernel(HalfDistanceVoxelCube* voxel_map, DummyQueue* queue, Water* material,
                                     PhotonData* photon_data, ElectronData* electron_data,
                                     const real_type primaryEnergy, const int gpu_histories, const unsigned seed) {
    Random rng{seed};
    Physics physics{rng, *voxel_map, *queue, *material, *photon_data, *electron_data};
    ThreeVector<real_type> initialPosition{0, 0, -0.5 * voxel_size};
    ThreeVector<real_type> initialDirection{0, 0, 1};
    PhotonPencilBeam beam{initialPosition, initialDirection, primaryEnergy};
    for (int i = 0; i < gpu_histories; ++i) {
        auto particle = beam.generateParticle(rng);
        physics.simulate(particle);
    }
}

TEST_CASE("GPU Photon simulation") {
    const unsigned seed = std::random_device()();
    INFO("SEED " << seed);
    SimMaterialData refMData;
    refMData.Load(DATA_DIR);
    SimPhotonData refPhData;
    refPhData.Load(DATA_DIR);
    auto device = cuda::device::current::get();
    ThreeVector<real_type> resolution(voxel_size, voxel_size, voxel_size);
    Water reference_material{DATA_DIR};
    reference_material.initializeGPU(device);
    auto gpu_material = cuda::make_unique(device, reference_material);

    PhotonData photon_data{DATA_DIR};
    photon_data.initializeGPU(device);
    auto gpu_photon_data = cuda::make_unique(device, photon_data);

    ElectronData electron_data{DATA_DIR};
    electron_data.initializeGPU(device);
    auto gpu_electron_data = cuda::make_unique(device, electron_data);

    DummyQueue queue;
    auto gpu_queue = cuda::make_unique(device, queue);

    for (int test_energy = 1; test_energy < 21; test_energy++) {
        SECTION("Energy " + std::to_string(test_energy) + "Mev") {
            const real_type primaryEnergy = test_energy;
            Geom geom(voxel_size, &refMData);
            Track primaryTrack;
            for (int i = 0; i < histories; ++i) {
                primaryTrack.fEkin         = primaryEnergy;  // Primary kinetic energy in [MeV]
                primaryTrack.fType         = 0;              // e-(-1) or photon(0)

                primaryTrack.fDirection[0] = 0.0;  // initial direction is [0,0,1]
                primaryTrack.fDirection[1] = 0.0;
                primaryTrack.fDirection[2] = 1.0;

                primaryTrack.fPosition[0]  = 0.0;  // initial position is [0,0, theRZ0]
                primaryTrack.fPosition[1]  = 0.0;
                primaryTrack.fPosition[2]  = -0.5 * voxel_size;
                KeepTrackingPhoton(refPhData, refMData, geom, primaryTrack);
            }

            auto reference = geom.histogram();
            while (!reference.empty() && reference[reference.size() - 1] == 0) {
                reference.pop_back();
            }
            reference.pop_back();

            real_type norm = 1. / *max_element(reference.begin(), reference.end());
            for (auto& item : reference) {
                item *= norm;
            }

            HalfDistanceVoxelCube voxel_map{{geom_size, geom_size, reference.size()}, resolution};
            voxel_map.setDefaultMaterial(gpu_material.get());
            voxel_map.initializeGPU(device);
            auto gpu_voxel_map = cuda::make_unique(device, voxel_map);

            cuda::launch(photonKernel, cuda::make_launch_config(grid_size, block_size), gpu_voxel_map.get(),
                         gpu_queue.get(), gpu_material.get(), gpu_photon_data.get(), gpu_electron_data.get(),
                         primaryEnergy, 1, seed);

            voxel_map.fromGPU();
            const auto result = voxel_map.doseDepthDistribution();
            for (auto& item : result) {
                std::cout << item << " ";
            }
            std::cout << std::endl;

            for (auto& item : reference) {
                std::cout << item << " ";
            }
            std::cout << std::endl;
            REQUIRE(opmc::KsTest(reference, result) >= alpha);
        }
        std::cout << std::endl;
    }
}
}  // namespace opmc