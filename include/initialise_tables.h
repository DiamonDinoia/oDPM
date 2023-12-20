//
// Created by mbarbone on 6/5/23.
//

#ifndef ODPM_INITIALISE_TABLES_H
#define ODPM_INITIALISE_TABLES_H

#include <ElectronData.hh>
#include <G4String.hh>
#include <PhotonData.hh>

#include "materials.h"

class G4Material;

namespace opmc {

template <typename T>
std::vector<G4String> saveMaterialNames(const std::vector<T> &materials);

template <typename T>
class DPMTables {
   public:
    DPMTables(const std::vector<T> &materials, const std::string &gInputDataDir = INPUT_DATA_DIR,
              const std::string &gOutputDataDir = OUTPUT_DATA_DIR, const std::string &gDataDirEnvVar = G4DATASET_DIR)
        : materialNames(saveMaterialNames(materials)), m_materials(initialiseMaterials(materials)) {
        ElectronData electronData(gInputDataDir, static_cast<int>(m_materials.size()),
                                  opmc::constants::k_ElectronCut * .99, opmc::constants::k_MaxEkin,
                                  opmc::constants::k_NumEkinForElectron);
        PhotonData photonData(static_cast<int>(m_materials.size()), opmc::constants::k_GammaCut * .99,
                              opmc::constants::k_MaxEkin, opmc::constants::k_NumEkinForElectron);
        initialiseTables(electronData, photonData, gOutputDataDir, gDataDirEnvVar);
        m_electronData = std::make_unique<ODPMElectronData>(electronData);
        m_photonData   = std::make_unique<ODPMPhotonData>(photonData);
        populateMaterials(electronData, photonData);
    }

    constexpr ~DPMTables() = default;

    const std::vector<G4Material *> &getMaterials() const;

    const std::string getMaterialName(const int id) const;

    const std::vector<G4String> materialNames;

    const std::vector<Material> &DPMmaterials() const;

    const Material &referenceMaterial() const;

    const ODPMElectronData &getElectronData() const;

    const ODPMPhotonData &getPhotonData() const;

   private:
    std::vector<G4Material *> m_materials;
    std::vector<Material> m_DPMmaterials;

    std::unique_ptr<ODPMElectronData> m_electronData;
    std::unique_ptr<ODPMPhotonData> m_photonData;

    void initialiseTables(ElectronData &electronData, PhotonData &photonData, const std::string &gOutputDataDir,
                          const std::string &gDataDirEnvVar);

    std::vector<G4Material *> initialiseMaterials(const std::vector<T> &materials);

    void populateMaterials(ElectronData &electronData, PhotonData &photonData);
};

template <typename T>
const Material &DPMTables<T>::referenceMaterial() const {
    return m_DPMmaterials.front();
}

}  // namespace opmc
#endif  // ODPM_INITIALISE_TABLES_H
