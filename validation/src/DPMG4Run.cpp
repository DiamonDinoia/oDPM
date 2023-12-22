
#include "DPMG4Run.h"

#include <FTFP_BERT.hh>
#include <G4EmStandardPhysicsSS.hh>
#include <G4EmStandardPhysics_option4.hh>
#include <G4GeometryManager.hh>
#include <G4LogicalVolumeStore.hh>
#include <G4PhysicalVolumeStore.hh>
#include <G4RunManager.hh>
#include <G4RunManagerFactory.hh>
#include <G4SolidStore.hh>
#include <Shielding.hh>

#include "DicomRun.hh"
#include "g4_utils.h"

namespace opmc {

G4Run::G4Run(const std::array<G4long, 2> &seeds, const std::string &gDataDirEnvVar, const std::string &physListName) {
    initializeEnvironment(gDataDirEnvVar);
    CLHEP::HepRandom::setTheEngine(new CLHEP::MixMaxRng);
    CLHEP::HepRandom::setTheSeeds(seeds.data());
    runManager          = std::unique_ptr<G4RunManager>(G4RunManagerFactory::CreateRunManager());
    const auto nThreads = static_cast<G4int>(std::thread::hardware_concurrency());
    runManager->SetNumberOfThreads(nThreads);
    if (physListName == "Shielding") {
        phys = new Shielding();
    } else if (physListName == "FTFP_BERT") {
        phys = new FTFP_BERT();
        phys->ReplacePhysics(new G4EmStandardPhysicsSS());
    } else if (physListName == "FTFP_BERT_Fast") {
        phys = new FTFP_BERT();
        phys->ReplacePhysics(new G4EmStandardPhysics_option4());
    } else {
        G4Exception("G4Run", "001", FatalException, "Physics List not found!!!");
    }
    runManager->SetUserInitialization(phys);
    theBeam = new OpmcDicomActionInitialization(nThreads);
}

static G4String inpfile = "Data.dat";

void G4Run::initializeDicom(const std::string &filename, const int nFiles) {
    create_dicom_dat(filename, nFiles);
    initializeGeometry();
}

void G4Run::initializeWater(const std::size_t dimx, const std::size_t dimy, const std::size_t dimz,
                            const real_type voxel_size) {
    create_water_dat(dimx, dimy, dimz, voxel_size, "WATER");
    initializeGeometry();
}

void G4Run::initializeGeometry() {
    runManager->GeometryHasBeenModified();
    theGeometry = new OpmcDicomDetectorConstruction();
    theFileMgr.Convert(inpfile);
    // Initialisation of physics, geometry, primary particles ...
    runManager->SetUserInitialization(theGeometry);
    runManager->SetUserInitialization(theBeam);
    ThreeVector<real_type> voxelSize{theGeometry->GetVoxelHalfX(), theGeometry->GetVoxelHalfY(),
                                     theGeometry->GetVoxelHalfZ()};
    std::cout << " Voxel size in initializeGeometry: " << voxelSize << std::endl;
    phys->SetCutValue(.1, "gamma");
    phys->SetCutValue(.1, "e-");
    phys->SetCutValue(.1, "e+");
    phys->SetCutValue(.1, "proton");
    runManager->Initialize();
}

VoxelCube G4Run::Run(const std::string &particle, const ThreeVector<double> &pos, const ThreeVector<double> &dir,
                     const double energy, const std::int32_t histories) {
    theBeam->SetBeam(particle, {pos.x, pos.y, pos.z}, {dir.x, dir.y, dir.z}, energy);
    runManager->BeamOn(histories);
    auto &run       = *dynamic_cast<const DicomRun *>(runManager->GetNonConstCurrentRun());
    auto &runAction = *dynamic_cast<const ODPMDicomRunAction *>(G4RunManager::GetRunManager()->GetUserRunAction());
    const auto id   = runAction.GetDicomRun()->GetRunID();
    std::cout << "Run ID = " << id << std::endl;
    const auto totalVoxels = theGeometry->GetTotalVoxels();
    if (totalVoxels == 0) {
        G4Exception("DicomRunAction", "001", JustWarning, "Number of Voxels is zero!!!");
        return VoxelCube{};
    }
    std::vector<double> doses(totalVoxels, 0);
    doses.reserve(totalVoxels);
    const auto hitsMaps = run.GetNumberOfHitsMap();
    std::cout << "hitsMaps = " << hitsMaps << std::endl;
    for (auto i = 0; i < hitsMaps; ++i) {
        auto &doseDeposit = *run.GetHitsMap(i)->GetMap();
        for (auto [index, dose] : doseDeposit) {
            if (index >= totalVoxels) {
                G4Exception("DicomRunAction", "001", FatalException, "Dose Voxel Index exceeds number of Voxels!!!");
            }
            if (index < totalVoxels) {
                doses[index] += *dose / CLHEP::gray;
            }
        }
    }
    ThreeVector<int> dim{theGeometry->GetNoVoxelsX(), theGeometry->GetNoVoxelsY(), theGeometry->GetNoVoxelsZ()};
    ThreeVector<real_type> voxelSize{theGeometry->GetVoxelHalfX() * 2, theGeometry->GetVoxelHalfY() * 2,
                                     theGeometry->GetVoxelHalfZ() * 2};
    ThreeVector<real_type> origin{theGeometry->GetMinX(), theGeometry->GetMinY(), theGeometry->GetMinZ()};
    ThreeVector<real_type> end{theGeometry->GetMaxX(), theGeometry->GetMaxY(), theGeometry->GetMaxZ()};

    std::cout << "dim = " << dim << std::endl;
    std::cout << "voxelSize = " << voxelSize << std::endl;
    std::cout << "origin = " << origin << std::endl;
    std::cout << "end = " << end << std::endl;
    std::cout << "Z = " << theGeometry->GetNoVoxelsZ() << std::endl;

    VoxelCube G4doses{dim, voxelSize, std::vector<Material>{}};
    std::copy(doses.begin(), doses.end(), G4doses.origin_dose());

    return G4doses;
}

VoxelCube G4Run::convertGeometry() { return convertGeometryToVoxelCube(*theGeometry, theFileMgr, tables); }

const DPMTables<G4String> *G4Run::getTables() const { return tables.get(); }
ThreeVector<double> G4Run::center() const { return getCenter(*theGeometry); }
ThreeVector<double> G4Run::resolution() const { return getResolution(*theGeometry); }

}  // namespace opmc