
#include <array>

#include "DPMG4Run.h"
#include "beams.h"
#include "run.h"

using namespace opmc;

int main(int argc, char **argv) {
    std::array<G4long, 2> seeds{1, 2};

    opmc::G4Run g4Run(seeds, G4DATASET_DIR);
    //    g4Run.initializeDicom("/home/mbarbone/repos/photonMonteCarlo/lung_ct/1-", 142);
    g4Run.initializeWater(32, 32, 32, 1.);
    //    g4Run.initializeGeometry();
    g4Run.Run("e-", {1.5, 1.5, 1000}, {0, 0, 1}, 6 * CLHEP::MeV, static_cast<int>(1));
    auto result                = g4Run.Run("e-", {0.5, 0.5, 0.001}, {0, 0, 1}, 6 * CLHEP::MeV, static_cast<int>(10e4));
    auto center                = g4Run.center();
    auto doseDepthDistribution = result.doseDepthDistribution();
    //    auto doseDepthDistribution2 = result2.doseDepthDistribution();
    auto sum1 = std::accumulate(doseDepthDistribution.begin(), doseDepthDistribution.end(), 0.);
    //    auto sum2 = std::accumulate(doseDepthDistribution2.begin(), doseDepthDistribution2.end(), 0.);
    //    for (auto i = 0; i < doseDepthDistribution.size(); ++i) {
    //        std::cout << doseDepthDistribution[i] << " " << doseDepthDistribution2[i] << std::endl;
    //    }
    //     PRINT THE SUM of th two distributions
    //    std::cout << sum1 << " " << sum2 << std::endl;
    //    std::cout << sum2 - sum1 << std::endl;
    // print the percentage difference
    //    std::cout << (sum2 - sum1) / sum1 << std::endl;

    auto voxelMap   = g4Run.convertGeometry();
    const auto seed = 42;

    auto DPMRrun    = Run<HalfDistanceVoxelCube>(seed);

    ThreeVector<real_type> initialPosition{0, 0, -static_cast<double>(voxelMap.resolution().z) * 0.49};
    constexpr ThreeVector<real_type> initialDirection{0, 0, 1};
    constexpr real_type initialEnergy = 6;
    auto beam                         = ElectronPencilBeam{initialPosition, initialDirection, initialEnergy};

    const auto tables                 = g4Run.getTables();
    auto referenceMaterial            = tables->referenceMaterial();
    auto electronData                 = tables->getElectronData();
    auto photonData                   = tables->getPhotonData();
    DPMRrun.simulate(voxelMap, beam, referenceMaterial, electronData, photonData, 10e4);

    auto dpmDoseDepthDistribution = voxelMap.doseDepthDistribution();

    for (auto i = 0; i < doseDepthDistribution.size(); ++i) {
        std::cout << doseDepthDistribution[i] << " " << dpmDoseDepthDistribution[i] << std::endl;
    }
    //    auto gamma = validateDistributions(result.doseDistribution(), voxelMap.doseDistribution());
    //    std::cout << "gamma " << gamma << std::endl;
    return 0;
}