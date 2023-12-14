//
// Created by mbarbone on 12/3/23.
//

#ifndef ODPM_DPMG4RUN_H
#define ODPM_DPMG4RUN_H

#include <G4RunManager.hh>
#include <array>

#include "g4_utils.h"
#include "geometry.h"
#include "initialise_tables.h"

class G4VModularPhysicsList;

namespace opmc {

class G4Run {
   public:
    G4Run(const std::array<G4long, 2> &seeds, const std::string &gDataDirEnvVar,
          const std::string &physListName = "Shielding");

    VoxelCube Run(const std::string &particle, const ThreeVector<double> &position,
                              const ThreeVector<double> &direction, const double energy, const std::int32_t histories);

    void initializeDicom(const std::string &filename, const int nFiles);

    void initializeWater(const std::size_t dimx, const std::size_t dimy, const std::size_t dimz,
                         const real_type voxel_size);

    void initializeGeometry();

    VoxelCube convertGeometry();

    const DPMTables<G4String> *getTables() const;

    ThreeVector<double> center() const;

   private:
    std::unique_ptr<G4RunManager> runManager;
    OpmcDicomDetectorConstruction *theGeometry;
    ODPMDicomFileMgr theFileMgr;
    G4VModularPhysicsList *phys;
    std::unique_ptr<DPMTables<G4String>> tables;
    std::vector<G4String> theMaterials;
};

}  // namespace opmc
#endif  // ODPM_DPMG4RUN_H
