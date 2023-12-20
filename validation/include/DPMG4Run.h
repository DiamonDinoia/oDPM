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
                  const ThreeVector<double> &direction, double energy, std::int32_t histories);

    void initializeDicom(const std::string &filename, int nFiles);

    void initializeWater(std::size_t dimx, std::size_t dimy, std::size_t dimz, real_type voxel_size);

    VoxelCube convertGeometry();

    const DPMTables<G4String> *getTables() const;

    ThreeVector<double> center() const;

   private:
    void initializeGeometry();

    std::unique_ptr<G4RunManager> runManager;
    OpmcDicomDetectorConstruction *theGeometry;
    ODPMDicomFileMgr theFileMgr;
    G4VModularPhysicsList *phys;
    std::unique_ptr<DPMTables<G4String>> tables;
    std::vector<G4String> theMaterials;
    OpmcDicomActionInitialization *theBeam;
};

}  // namespace opmc
#endif  // ODPM_DPMG4RUN_H
