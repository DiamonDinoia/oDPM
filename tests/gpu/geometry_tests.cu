//
// Created by mbarbone on 5/27/22.
//
//
// Created by mbarbone on 7/29/21.
//

#define CATCH_CONFIG_MAIN
#include <catch2/catch.hpp>

#include "Geom.hh"
#include "geometry.h"
#include "gpu_array.h"
#include "gpu_unique_ptr.h"
#include "gpu_utils.h"
#include "test_utils.h"

namespace opmc {

static constexpr auto size       = 2;
static constexpr auto geom_size  = 50;
static constexpr auto tests      = 16384;

static constexpr auto block_size = 256;
static constexpr auto grid_size  = div_rounding_up(tests, block_size);

static ODPM_KERNEL void distanceKernel(ParticleInterface* inputElectrons, HalfDistanceVoxelCube* voxel_map,
                                       real_type* distance) {
    const auto i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i < tests) {
        Particle electron;
        inputElectrons[i].toParticle(electron);
        distance[i] = voxel_map->distanceToBoundary(electron);
        inputElectrons[i].fromParticle(electron);
    }
}

TEST_CASE("GPU GEOMETRY") {
    const auto seed = std::random_device()();
    INFO("SEED " << seed)
    std::default_random_engine rng(seed);
    SimMaterialData sim;
    sim.Load(DATA_DIR);
    const auto energy = 0;
    Geom geom(size, &sim);
    ThreeVector<real_type> resolution(size, size, size);
    TestHalfDistanceVoxelCube voxel_map{{geom_size, geom_size, geom_size}, resolution};
    auto device = cuda::device::current::get();

    std::uniform_real_distribution<real_type> dist(0, size);
    std::uniform_real_distribution<real_type> direction_dist(-1, 1);

    std::uniform_int_distribution<int> position_dist(-geom_size / 2, geom_size / 2);

    std::vector<Electron> initial_electrons(tests), final_electrons(tests);
    std::vector<real_type> reference_distance(tests);
    std::vector<std::array<int, 3>> reference_iVoxel(tests);
    std::vector<ParticleInterface> gpu_electrons(tests);

    for (int i = 0; i < tests; ++i) {
        INFO("TEST " << i)
        real_type position[3]  = {dist(rng) + position_dist(rng), dist(rng) + position_dist(rng),
                                  dist(rng) + position_dist(rng) + geom_size};
        real_type direction[3] = {direction_dist(rng), direction_dist(rng), direction_dist(rng)};
        ThreeVector<real_type> position_vector{position[0], position[1], position[2]};
        ThreeVector<real_type> direction_vector{direction[0], direction[1], direction[2]};
        initial_electrons[i] = {position_vector, direction_vector, energy};
        gpu_electrons[i].fromParticle(initial_electrons[i]);
        reference_distance[i] = geom.DistanceToBoundary(position, direction, reference_iVoxel[i].data());
    }

    // GPU computation

    // Copy data to GPU
    voxel_map.initializeGPU(device);
    auto gpu_voxel_map = cuda::make_unique(device, voxel_map);
    cuda::gpu_array<real_type> gpu_distance(cuda::device::current::get(), tests);
    cuda::gpu_array<ParticleInterface> gpu_particles(cuda::device::current::get(), gpu_electrons);
    // execute the computation
    cuda::launch(distanceKernel, cuda::make_launch_config(grid_size, block_size), gpu_particles.data(),
                 gpu_voxel_map.get(), gpu_distance.data());
    // copy data back
    auto result_distance = gpu_distance.to_vector();
    gpu_electrons        = gpu_particles.to_vector();
    for (int i = 0; i < tests; ++i) {
        gpu_electrons[i].toParticle(final_electrons[i]);
    }

    SECTION("DISTANCE TO BOUNDARY") {
        for (int i = 0; i < tests; ++i) {
            INFO("TEST " << i)
            INFO("electron " << initial_electrons[i])
            const real_type reference = reference_distance[i];
            const real_type result    = result_distance[i];
            if (result <= constants::k_minStepLength || reference <= constants::k_minStepLength) {
                // Microsteps are handled differently in the two codebases
                continue;
            }
            REQUIRE(reference == Approx(result));
            REQUIRE(reference_iVoxel[i][0] + geom_size / 2 == voxel_map.getPosition(final_electrons[i].position).x);
            REQUIRE(reference_iVoxel[i][1] + geom_size / 2 == voxel_map.getPosition(final_electrons[i].position).y);
            REQUIRE(reference_iVoxel[i][2] == voxel_map.getPosition(final_electrons[i].position).z);
        }
    }
}
}  // namespace opmc
