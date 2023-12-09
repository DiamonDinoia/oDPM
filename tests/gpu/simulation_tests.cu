//
// Created by mbarbone on 10/26/21.
//
#define CATCH_CONFIG_ENABLE_BENCHMARKING
#define CATCH_CONFIG_MAIN
#include <catch2/catch.hpp>
#include <chrono>

#include "Geom.hh"
#include "SimDPMLike.hh"
#include "SimElectronData.hh"
#include "SimMaterialData.hh"
#include "SimPhotonData.hh"
#include "gpu_run.h"
#include "ksTest.h"

using namespace std;
using namespace opmc;

const auto voxel_size = 2;
const auto geom_size  = 200;
const auto alpha      = 5.0E-9;
const auto treshold   = 0.00;

template <typename T, typename V>
void testSimulation(const real_type primaryEnergy, const unsigned long histories, SimMaterialData &materialData,
                    SimElectronData &electronData, SimPhotonData &photonData) {
    const ThreeVector<real_type> resolution(voxel_size, voxel_size, voxel_size);
    constexpr real_type theRZ0   = std::is_same<V, HalfDistanceVoxelCube>::value ? -0.5 * voxel_size : 0;
    constexpr real_type initialX = std::is_same<V, HalfDistanceVoxelCube>::value ? 0 : geom_size * voxel_size * .5;
    constexpr real_type initialY = std::is_same<V, HalfDistanceVoxelCube>::value ? 0 : geom_size * voxel_size * .5;

    auto geom      = Simulate(histories, primaryEnergy, std::is_same<T, Electron>::value, voxel_size, materialData,
                              electronData, photonData, 0);

    auto reference = geom.histogram();

    real_type norm = 1. / *max_element(reference.begin(), reference.end());
    for (auto &item : reference) {
        item *= norm;
    }

    reference.erase(std::remove_if(reference.begin(), reference.end(), [](real_type x) { return x <= treshold; }),
                    reference.end());

    ThreeVector<unsigned long> geometrySize{geom_size, geom_size, reference.size()};

    const ThreeVector<real_type> initialPosition{initialX, initialY, theRZ0};
    const ThreeVector<real_type> initialDirection{0, 0, 1};

    PencilBeam<T> beam{initialPosition, initialDirection, primaryEnergy};

    const auto seed0 = std::random_device()();
    const auto seed1 = std::random_device()();
    std::cout << "SEED SIMULATION " << seed0 << " " << seed1 << std::endl;

    GPURun<PencilBeam<T>, V> run{histories, geometrySize, resolution, beam, seed0, seed1, DATA_DIR};

    INFO("" + std::string(typeid(V).name()) + " energy " + std::to_string(primaryEnergy));
    std::cout << "" + std::string(typeid(V).name()) + " energy " + std::to_string(primaryEnergy) << std::endl;
    //
    std::vector<real_type> result;
    const auto start     = std::chrono::steady_clock::now();
    const auto &voxelMap = run.simulate();
    const auto end       = std::chrono::steady_clock::now();
    result               = voxelMap.doseDepthDistribution();
    //
    std::chrono::duration<double, std::milli> timeRequired = (end - start);
    INFO("Monte Carlo Simulation required " << timeRequired.count() << " milliseconds");
    std::cout << "Monte Carlo Simulation required " << timeRequired.count() << " milliseconds" << std::endl;

    for (auto &item : result) {
        std::cout << item << " ";
    }
    std::cout << std::endl;

    for (auto &item : reference) {
        std::cout << item << " ";
    }
    std::cout << std::endl;

    REQUIRE(opmc::KsTest(reference, result) >= alpha);
}

TEST_CASE("GPU Complete simulation") {
    SimMaterialData materialData{};
    SimElectronData electronData{};
    SimPhotonData photonData{};
    materialData.Load(DATA_DIR);
    electronData.Load(DATA_DIR);
    photonData.Load(DATA_DIR);
    const unsigned long el_histories = 10e4;
    const unsigned long ph_histories = 10e4;
    for (int i = 1; i < 21; ++i) {
        SECTION("Electron pencil beam " + std::to_string(i) + "Mev") {
            std::cout << "Electron-------------------------------------------------------------------------"
                      << std::endl;
            testSimulation<Electron, VoxelCube>(i, el_histories, materialData, electronData, photonData);
        }
        SECTION("Electron pencil beam " + std::to_string(i) + "Mev") {
            std::cout << "Electron-------------------------------------------------------------------------"
                      << std::endl;
            testSimulation<Electron, HalfDistanceVoxelCube>(i, el_histories, materialData, electronData, photonData);
        }
    }
    for (int i = 1; i < 21; ++i) {
        SECTION("Photon pencil beam " + std::to_string(i) + "Mev") {
            std::cout << "Photon----------------------------------------------------------------------------"
                      << std::endl;
            testSimulation<Photon, VoxelCube>(i, ph_histories, materialData, electronData, photonData);
        }
        SECTION("Photon pencil beam " + std::to_string(i) + "Mev") {
            std::cout << "Photon----------------------------------------------------------------------------"
                      << std::endl;
            testSimulation<Photon, HalfDistanceVoxelCube>(i, ph_histories, materialData, electronData, photonData);
        }
    }
}