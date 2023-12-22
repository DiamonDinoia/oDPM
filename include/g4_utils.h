#ifndef ODPM_G4_UTILS_H
#define ODPM_G4_UTILS_H

#include <DicomActionInitialization.hh>
#include <DicomEventAction.hh>
#include <DicomFileMgr.hh>
#include <DicomNestedParamDetectorConstruction.hh>
#include <DicomRunAction.hh>
#include <G4ParticleGun.hh>
#include <G4ParticleTable.hh>
#include <G4VUserPrimaryGeneratorAction.hh>
#include <ostream>

#include "geometry.h"
#include "types.h"

namespace opmc {
template <typename T>
class DPMTables;

void create_dicom_dat(const std::string &filename, const int nFiles);

void create_water_dat(const std::size_t dimx, const std::size_t dimy, const std::size_t dimz,
                      const real_type voxel_size, const std::string &filename = "WATER");

void initializeEnvironment(const std::string &gDataDirEnvVar, bool overwrite = true);

// DicomRegularDetectorConstruction
// DicomNestedParamDetectorConstruction
// it can be switched to DicomNestedParamDetectorConstruction
// if there are issues
class OpmcDicomDetectorConstruction : public DicomNestedParamDetectorConstruction {
   public:
    using DicomNestedParamDetectorConstruction::DicomNestedParamDetectorConstruction;

    G4int GetNoVoxelsX() const;
    ;

    G4int GetNoVoxelsY() const;
    ;

    G4int GetNoVoxelsZ() const;
    ;

    G4double GetVoxelHalfX() const;
    ;

    G4double GetVoxelHalfY() const;
    ;

    G4double GetVoxelHalfZ() const;
    ;

    G4double GetMinX() const;
    ;

    G4double GetMinY() const;
    ;

    G4double GetMinZ() const;
    ;

    G4double GetMaxX() const;
    ;

    G4double GetMaxY() const;
    ;

    G4double GetMaxZ() const;
    ;

    G4ThreeVector GetTranslation() const;
    ;
};

class ODPMDicomRunAction : public DicomRunAction {
   public:
    using DicomRunAction::DicomRunAction;

    const auto &GetSDName() const { return fSDName; }
};

class ODPMDicomPrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction {
   public:
    ODPMDicomPrimaryGeneratorAction();

    void GeneratePrimaries(G4Event *anEvent) override;

    void SetBeam(const G4String &particleName, G4ThreeVector position, G4ThreeVector direction, const double energy);

   private:
    G4ParticleGun fParticleGun;
    G4ParticleDefinition *fParticleDefinition;
    G4ThreeVector fPosition;
    G4ThreeVector fDirection;
    double fEnergy;
};

class OpmcDicomActionInitialization : public DicomActionInitialization {
   public:
    explicit OpmcDicomActionInitialization(int nThreads);

    void Build() const override;

    void BuildForMaster() const override;

    void SetBeam(const G4String &particleName, G4ThreeVector position, G4ThreeVector direction, const double energy);

   private:
    std::vector<ODPMDicomPrimaryGeneratorAction *> fPrimaryGeneratorActions;
};

class ODPMDicomFileMgr {
   public:
    ODPMDicomFileMgr();

    void CheckNColumns(std::vector<G4String> wl, size_t vsizeTh);

    G4int GetCompression() const;

    G4String GetFileOutName() const;

    std::vector<G4String> getMaterials() const;

    std::vector<const DicomFileCT *> getCTs() const;

    void SetVoxels(const G4String &fVoxelsX, const G4String &fVoxelsY, const G4String &fVoxelsZ, const G4String &fDim);

    void DumpWaterToTextFile() const;

    void Convert(const G4String &fileName);

   private:
    DicomFileMgr &instance;
    G4int fNoVoxelsX;
    G4int fNoVoxelsY;
    G4int fNoVoxelsZ;
    G4double fDimVox;
};

VoxelCube convertGeometryToVoxelCube(OpmcDicomDetectorConstruction &theGeometry, ODPMDicomFileMgr &theFileMgr,
                                     std::unique_ptr<DPMTables<G4String>> &tables);

struct GammaResults {
    double maxRelativeError;
    double averageRelativeError;
    double gammaPassingRate;
    double minGamma;

    friend std::ostream &operator<<(std::ostream &os, const GammaResults &results);
};

GammaResults validateDistributions(DoseDistribution<> reference_distribution, DoseDistribution<> opmc_distribution);

ThreeVector<real_type> getCenter(const OpmcDicomDetectorConstruction &theGeometry);

ThreeVector<real_type> getResolution(const OpmcDicomDetectorConstruction &theGeometry);

}  // namespace opmc

#endif  // ODPM_G4_UTILS_H
