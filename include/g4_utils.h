#ifndef ODPM_G4_UTILS_H
#define ODPM_G4_UTILS_H


#include <DicomNestedParamDetectorConstruction.hh>
#include <DicomRunAction.hh>
#include <G4VUserPrimaryGeneratorAction.hh>
#include <G4ParticleTable.hh>
#include <G4ParticleGun.hh>
#include <DicomActionInitialization.hh>
#include <DicomEventAction.hh>
#include <DicomFileMgr.hh>
#include <ostream>

#include "types.h"
#include "geometry.h"

namespace opmc {
    template<typename T>
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

        G4int GetNoVoxelsX() const { return fNoVoxelsX; };

        G4int GetNoVoxelsY() const { return fNoVoxelsY; };

        G4int GetNoVoxelsZ() const { return fNoVoxelsZ; };

        G4double GetVoxelHalfX() const { return fVoxelHalfDimX; };

        G4double GetVoxelHalfY() const { return fVoxelHalfDimY; };

        G4double GetVoxelHalfZ() const { return fVoxelHalfDimZ; };

        G4double GetMinX() const { return fMinX; };

        G4double GetMinY() const { return fMinY; };

        G4double GetMinZ() const { return fMinZ; };

        G4double GetMaxX() const { return fMaxX; };

        G4double GetMaxY() const { return fMaxY; };

        G4double GetMaxZ() const { return fMaxZ; };
    };

    class ODPMDicomRunAction : public DicomRunAction {
    public:
        using DicomRunAction::DicomRunAction;

        const auto &GetSDName() const { return fSDName; }
    };

    class ODPMDicomPrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction {
    public:
        ODPMDicomPrimaryGeneratorAction(const G4String &particleName, G4ThreeVector position, G4ThreeVector direction,
                                        double energy);

        virtual void GeneratePrimaries(G4Event *anEvent);

    private:
        G4ParticleGun fParticleGun;
        G4ParticleDefinition *fParticleDefinition;
        const G4ThreeVector fPosition;
        const G4ThreeVector fDirection;
        const double energy;
    };

    class OpmcDicomActionInitialization : public DicomActionInitialization {
    public:
        using DicomActionInitialization::DicomActionInitialization;

        OpmcDicomActionInitialization(G4String particleName, G4ThreeVector position, G4ThreeVector direction,
                                      const double energy);

        void Build() const override;

        void BuildForMaster() const override { SetUserAction(new ODPMDicomRunAction); }

    private:
        const G4String particleName;
        const G4ThreeVector fPosition;
        const G4ThreeVector fDirection;
        const double energy;
    };

    class ODPMDicomFileMgr {
    public:
        ODPMDicomFileMgr();

        void CheckNColumns(std::vector<G4String> wl, size_t vsizeTh);

        G4int GetCompression() const;

        G4String GetFileOutName() const;

        std::vector<G4String> getMaterials() const;

        std::vector<const DicomFileCT *> getCTs() const;

        void SetVoxels(G4String fVoxelsX, G4String fVoxelsY, G4String fVoxelsZ, G4String fDim);

        void DumpWaterToTextFile();

        void Convert(G4String fileName);

    private:
        DicomFileMgr &instance;
        G4int fNoVoxelsX;
        G4int fNoVoxelsY;
        G4int fNoVoxelsZ;
        G4double fDimVox;
    };

    VoxelCube convertGeometryToVoxelCube(OpmcDicomDetectorConstruction &theGeometry,
                                                     ODPMDicomFileMgr &theFileMgr,
                                                     std::unique_ptr<DPMTables<G4String>> &tables);

    struct GammaResults {
        double maxRelativeError;
        double averageRelativeError;
        double gammaPassingRate;
        double minGamma;

        friend std::ostream &operator<<(std::ostream &os, const GammaResults &results);
    };

    GammaResults validateDistributions(DoseDistribution<> reference_distribution,
                                       DoseDistribution<> opmc_distribution);

    ThreeVector<real_type> getCenter(const OpmcDicomDetectorConstruction &theGeometry);
}


#endif //ODPM_G4_UTILS_H
