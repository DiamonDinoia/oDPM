//
// Created by mbarbone on 8/13/21.
//
#define CATCH_CONFIG_MAIN

#include <Geom.hh>
#include <SimDPMLike.hh>
#include <SimElectronData.hh>
#include <SimMaterialData.hh>
#include <Track.hh>
#include <catch2/catch.hpp>

#include "constants.h"
#include "geometry.h"
#include "physics.h"
#include "random.h"
#include "test_utils.h"

using namespace std;
using namespace opmc;

const auto tests                 = 16384;
const auto voxel_size            = 2;
const auto geom_size             = 200;
const auto epsilon               = 10e-5;

static constexpr auto block_size = 256;
static constexpr auto grid_size  = div_rounding_up(tests, block_size);

static ODPM_KERNEL void trackingKernel(HalfDistanceVoxelCube* voxel_map, DummyQueue* queue, Water* material,
                                       PhotonData* photon_data, ElectronData* electron_data,
                                       ParticleInterface* electrons, ParticleState* particle_state, Event* events) {
    const auto i = blockIdx.x * blockDim.x + threadIdx.x;
    Random rng{0};
    Physics physics{rng, *voxel_map, *queue, *material, *photon_data, *electron_data};
    Electron electron;
    electrons[i].toParticle(electron);
    physics.initializeTracking(electron.energy, particle_state[i]);
    ParticleState new_state = particle_state[i];
    events[i]               = physics.continuousElectronTracking(electron, new_state);
    electrons[i].fromParticle(electron);
}

TEST_CASE("GPU Electron Tracking") {
    auto seed = std::random_device()();
    INFO("SEED " << seed);
    Random rng{seed};
    std::default_random_engine generator(rng.getUniform());
    SimMaterialData refMData;
    refMData.Load(DATA_DIR);
    SimElectronData refElData;
    refElData.Load(DATA_DIR);
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

    Geom geom(voxel_size, &refMData);

    HalfDistanceVoxelCube voxel_map{{geom_size, geom_size, geom_size * 10}, resolution};
    voxel_map.setDefaultMaterial(gpu_material.get());
    voxel_map.initializeGPU(device);
    auto gpu_voxel_map = cuda::make_unique(device, voxel_map);

    std::uniform_real_distribution<real_type> energy_distribution(opmc::constants::k_MinEkin,
                                                                  opmc::constants::k_MaxEkin);
    std::uniform_real_distribution<real_type> position_distribution(-voxel_size / 2, voxel_size / 2);
    std::uniform_real_distribution<real_type> direction_distribution(-1, 1);

    auto getEnergy    = [&generator, &energy_distribution] { return energy_distribution(generator); };
    auto getPosition  = [&generator, &position_distribution] { return position_distribution(generator); };
    auto getDirection = [&generator, &direction_distribution] { return direction_distribution(generator); };

    std::vector<Electron> initial_electrons(tests);
    std::vector<ParticleState> initial_state(tests);
    std::vector<ParticleInterface> gpu_particles(tests);
    for (int i = 0; i < tests; ++i) {
        const auto energy            = getEnergy();
        const real_type position[3]  = {getPosition(), getPosition(), getPosition()};
        const real_type direction[3] = {getDirection(), getDirection(), getDirection()};
        ThreeVector<real_type> position_vector{position[0], position[1], position[2]};
        ThreeVector<real_type> direction_vector{direction[0], direction[1], direction[2]};
        initial_electrons[i] = {position_vector, direction_vector, energy};
        gpu_particles[i].fromParticle(initial_electrons[i]);
    }

    cuda::gpu_array<ParticleInterface> gpu_electrons(cuda::device::current::get(), gpu_particles);
    cuda::gpu_array<ParticleState> gpu_state(cuda::device::current::get(), initial_state);
    cuda::gpu_array<Event> gpu_events(cuda::device::current::get(), tests);

    cuda::launch(trackingKernel, cuda::make_launch_config(grid_size, block_size), gpu_voxel_map.get(), gpu_queue.get(),
                 gpu_material.get(), gpu_photon_data.get(), gpu_electron_data.get(), gpu_electrons.data(),
                 gpu_state.data(), gpu_events.data());
    const auto final_electrons = gpu_electrons.to_vector();
    initial_state              = gpu_state.to_vector();
    std::vector<Event> events  = gpu_events.to_vector();

    SECTION("EVENT") {
        Track primaryTrack;
        ParticleInterface primary_electron;
        Electron final_electron;
        for (int i = 0; i < tests; ++i) {
            final_electrons[i].toParticle(final_electron);
            INFO("TEST " << i)  //
            INFO("Initial state " << initial_state[i])
            INFO("Event" << events[i]);
            INFO("Initial Electron " << initial_electrons[i])
            INFO("Final Electron " << final_electron)
            REQUIRE(events[i] != opmc::NOTHING);
            if (events[i] == opmc::MICROSTEP || events[i] == opmc::ABSORBED || events[i] == opmc::OUT_OF_BOUNDARY) {
                continue;
            }
            ParticleState state = initial_state[i];
            primary_electron.fromParticle(initial_electrons[i]);
            primary_electron.toTrack(primaryTrack);
            const auto reference = KeepTrackingElectron(refElData, refMData, geom, primaryTrack, state.num_tr_1_mfp,
                                                        state.num_moller_mfp, state.inv_moller_mfp, state.num_brem_mfp);
            REQUIRE(events[i] == reference);
            REQUIRE(final_electron.energy == Approx(primaryTrack.fEkin));
            REQUIRE(final_electron.position.x == Approx(primaryTrack.fPosition[0]).epsilon(epsilon));
            REQUIRE(final_electron.position.y == Approx(primaryTrack.fPosition[1]).epsilon(epsilon));
            REQUIRE(final_electron.position.z == Approx(primaryTrack.fPosition[2]).epsilon(epsilon));
            REQUIRE(final_electron.direction.x == Approx(primaryTrack.fDirection[0]).epsilon(epsilon));
            REQUIRE(final_electron.direction.y == Approx(primaryTrack.fDirection[1]).epsilon(epsilon));
            REQUIRE(final_electron.direction.z == Approx(primaryTrack.fDirection[2]).epsilon(epsilon));
        }
    }
}
