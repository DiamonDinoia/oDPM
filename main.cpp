
#include <G4RunManager.hh>
#include <Shielding.hh>
#include "G4RunManagerFactory.hh"
#include "g4_utils.h"
#include "run.h"
#include "beams.h"

using namespace opmc;

static G4String inpfile = "Data.dat";


int main(int argc, char **argv) {

    initializeEnvironment(G4DATASET_DIR);
    create_water_dat(32, 32, 32, 1., "WATER");
    ODPMDicomFileMgr theFileMgr;
    auto runManager = std::unique_ptr<G4RunManager>(G4RunManagerFactory::CreateRunManager());
    auto *theGeometry = new OpmcDicomDetectorConstruction();
    runManager->SetUserInitialization(theGeometry);
    auto* phys = new Shielding();
    runManager->SetUserInitialization(phys);
    theFileMgr.Convert(inpfile);
    std::unique_ptr<opmc::DPMTables<G4String>> tables;

    auto voxelMap = convertGeometryToVoxelCube(*theGeometry, theFileMgr, tables);

    const auto seed = 42;

    auto DPMRrun = Run<HalfDistanceVoxelCube>(seed);

    ThreeVector<real_type> initialPosition{0, 0, -static_cast<double>(voxelMap.resolution().z) * 0.49};
    constexpr ThreeVector<real_type> initialDirection{0, 0, 1};
    constexpr real_type initialEnergy = 6;
    auto beam = ElectronPencilBeam{initialPosition, initialDirection, initialEnergy};
    auto referenceMaterial = tables->referenceMaterial();
    auto electronData = tables->getElectronData();
    auto photonData = tables->getPhotonData();
    DPMRrun.simulate(voxelMap, beam, referenceMaterial, electronData, photonData, 10);

    auto dpmDoseDepthDistribution = voxelMap.doseDepthDistribution();

    return 0;
}