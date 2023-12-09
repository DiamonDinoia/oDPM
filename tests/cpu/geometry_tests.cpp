//
// Created by mbarbone on 7/29/21.
//

#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>
#include <random>

#include "Geom.hh"
#include "geometry.h"
#include "initialise_tables.h"

using namespace opmc;
using Catch::Approx;

TEST_CASE("GEOMETRY") {
    const auto size = 2;
    const auto geom_size = 50;
    const auto seed = std::random_device()();
    INFO("SEED " << seed);
    std::default_random_engine rng(seed);
    const auto tests = 16384;
    SimMaterialData sim;
    sim.Load(OUTPUT_DATA_DIR);
    const auto energy = 0;
    Geom geom(size, &sim);
    static const std::vector<std::string> test_materials{
            "G4_WATER",
    };
    ThreeVector<real_type> resolution(size, size, size);
    DPMTables tables{test_materials};
    TestHalfDistanceVoxelCube voxel_map{{geom_size, geom_size, geom_size}, resolution, tables.populateMaterials()};
    std::uniform_real_distribution<real_type> dist(0, size);
    std::uniform_real_distribution<real_type> direction_dist(-1, 1);

    std::uniform_int_distribution<int> position_dist(-geom_size / 2, geom_size / 2);
    int iVoxel[3] = {0, 0, 0};

    SECTION("DISTANCE TO BOUNDARY") {
        for (int i = 0; i < tests; ++i) {
            INFO("TEST " << i);
            real_type position[3] = {dist(rng) + position_dist(rng), dist(rng) + position_dist(rng),
                                     dist(rng) + position_dist(rng) + geom_size};
            real_type direction[3] = {direction_dist(rng), direction_dist(rng), direction_dist(rng)};
            ThreeVector<real_type> position_vector{position[0], position[1], position[2]};

            ThreeVector<real_type> direction_vector{direction[0], direction[1], direction[2]};
            INFO("Position " << position_vector);
            INFO("Direction" << direction_vector);
            const Electron electron{position_vector, direction_vector, energy};
            const real_type reference = geom.DistanceToBoundary(position, direction, iVoxel);
            const real_type result = voxel_map.distanceToBoundary(electron);
            if (result < 1.0E-3) {
                continue;
            }
            REQUIRE(reference == Approx(result));
            REQUIRE(iVoxel[0] + geom_size / 2 == voxel_map.getPosition(position_vector).x);
            REQUIRE(iVoxel[1] + geom_size / 2 == voxel_map.getPosition(position_vector).y);
            REQUIRE(iVoxel[2] == voxel_map.getPosition(position_vector).z);
        }
    }
}
