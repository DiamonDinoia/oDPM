
#include <G4RunManagerFactory.hh>
#include <G4RunManager.hh>
#include <Shielding.hh>
#include <G4GeometryManager.hh>
#include <G4PhysicalVolumeStore.hh>
#include <G4LogicalVolumeStore.hh>
#include <G4SolidStore.hh>

#include "DPMG4Run.h"
#include "DicomRun.hh"
#include "g4_utils.h"


namespace opmc {

    G4Run::G4Run(const std::array<G4long, 2> &seeds, const std::string &gDataDirEnvVar) {
        initializeEnvironment(gDataDirEnvVar);
        CLHEP::HepRandom::setTheEngine(new CLHEP::MixMaxRng);
        CLHEP::HepRandom::setTheSeeds(seeds.data());
        runManager = std::unique_ptr<G4RunManager>(G4RunManagerFactory::CreateRunManager());
        runManager->SetNumberOfThreads(static_cast<G4int>(std::thread::hardware_concurrency()));
        phys = new Shielding();
        runManager->SetUserInitialization(phys);
    }

    void G4Run::initializeDicom(const std::string &filename, const int nFiles) {
        create_dicom_dat(filename, nFiles);
    }

    void G4Run::initializeWater(const std::size_t dimx, const std::size_t dimy, const std::size_t dimz,
                                const real_type voxel_size) {
        create_water_dat(dimx, dimy, dimz, voxel_size, "WATER");
    }

    static G4String inpfile = "Data.dat";


    void G4Run::initializeGeometry() {
        theGeometry = new OpmcDicomDetectorConstruction();
        theFileMgr.Convert(inpfile);
        // Initialisation of physics, geometry, primary particles ...
        runManager->SetUserInitialization(theGeometry);
        ThreeVector<real_type> voxelSize{theGeometry->GetVoxelHalfX(), theGeometry->GetVoxelHalfY(),
                                         theGeometry->GetVoxelHalfZ()};
        phys->SetDefaultCutValue(voxelSize.min());
    }

    HalfDistanceVoxelCube G4Run::Run(const std::string &particle, const ThreeVector<double> &pos,
                                     const ThreeVector<double> &dir, const double energy,
                                     const std::int32_t histories) {
        runManager->SetUserInitialization(
                new OpmcDicomActionInitialization(particle, {pos.x, pos.y, pos.z}, {dir.x, dir.y, dir.z}, energy));
        runManager->Initialize();
        runManager->BeamOn(histories);

        auto &run = *dynamic_cast<const DicomRun *>(runManager->GetNonConstCurrentRun());
        auto &runAction = *dynamic_cast<const ODPMDicomRunAction *>(G4RunManager::GetRunManager()->GetUserRunAction());
        const auto id = runAction.GetDicomRun()->GetRunID();
        std::cout << "Run ID = " << id << std::endl;
        const auto totalVoxels = theGeometry->GetTotalVoxels();
        std::vector<double> doses(totalVoxels, 0);
        doses.reserve(totalVoxels);
        const auto hitsMaps = run.GetNumberOfHitsMap();
        for (auto i = 0; i < hitsMaps; ++i) {
            auto &doseDeposit = *run.GetHitsMap(i)->GetMap();
            for (auto [index, dose]: doseDeposit) {
                if (index >= totalVoxels) {
                    G4Exception("DicomRunAction", "001", FatalException,
                                "Dose Voxel Index exceeds number of Voxels!!!");
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

        HalfDistanceVoxelCube G4doses{dim, voxelSize, std::vector<Material>{}};
        std::copy(doses.begin(), doses.end(), G4doses.origin_dose());

        return G4doses;
    }

    HalfDistanceVoxelCube G4Run::convertGeometry() {
        return convertGeometryToVoxelCube(*theGeometry, theFileMgr, tables);
    }

    const DPMTables<G4String>* G4Run::getTables() const {
        return tables.get();
    }

}