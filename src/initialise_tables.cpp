//
// Created by mbarbone on 6/5/23.
//

#include "initialise_tables.h"

#include <ElectronData.hh>
#include <G4GeometryManager.hh>
#include <G4LogicalVolumeStore.hh>
#include <G4NistManager.hh>
#include <G4PhysicalVolumeStore.hh>
#include <G4RunManager.hh>
#include <G4Setup.hh>
#include <G4SolidStore.hh>
#include <G4StateManager.hh>
#include <InitElectronData.hh>
#include <InitPhotonData.hh>
#include <PhotonData.hh>

#include "constants.h"

namespace opmc {

static void initialiseEnvironment(const std::string &gDataDirEnvVar) {
    setenv("G4ENSDFSTATEDATA", (gDataDirEnvVar + "/G4ENSDFSTATE2.3").c_str(), 0);
    std::cout << "Setting environtment variable G4ENSDFSTATEDATA = " << getenv("G4ENSDFSTATEDATA") << std::endl;
    setenv("G4LEDATA", (gDataDirEnvVar + "/G4EMLOW8.2").c_str(), 0);
    std::cout << "Setting environtment variable G4LEDATA = " << getenv("G4LEDATA") << std::endl;
}

static void resetGeant4() {
    G4GeometryManager::GetInstance()->OpenGeometry();
    G4PhysicalVolumeStore::GetInstance()->Clean();
    G4LogicalVolumeStore::GetInstance()->Clean();
    G4SolidStore::GetInstance()->Clean();
    G4Material::GetMaterialTable()->clear();
}

template <typename T>
std::vector<G4Material *> DPMTables<T>::initialiseMaterials(const std::vector<T> &input_materials) {
    std::vector<G4Material *> materials;
    if constexpr (std::is_same_v<T, G4String> || std::is_same_v<T, std::string>) {
        for (auto &materialName : input_materials) {
            materials.push_back(G4NistManager::Instance()->FindOrBuildMaterial(materialName));
        }
    } else {
        materials = input_materials;
    }
    // iterate over materials and if the first is not G4_WATER put it in front
    // of the vector
    bool found = false;
    for (auto it = materials.begin(); it != materials.end(); ++it) {
        if ((*it)->GetName() == "G4_WATER") {
            found = true;
            std::swap(*it, *materials.begin());
            break;
        }
    }
    // first element must be water
    if (!found) {
        materials.insert(materials.begin(), G4NistManager::Instance()->FindOrBuildMaterial("G4_WATER"));
    }
    return materials;
}

template <typename T>
void DPMTables<T>::initialiseTables(ElectronData &electronData, PhotonData &photonData,
                                    const std::string &gOutputDataDir, const std::string &gDataDirEnvVar) {
    initialiseEnvironment(gDataDirEnvVar);

    resetGeant4();
    FakeG4Setup(m_materials, opmc::constants::k_ElectronCut, opmc::constants::k_GammaCut);

    std::cout << "\n === Data are generated now (don't care about the "
              << "\n     Geant4 warning: 'WARNING from G4PenelopeIonisationModel...')\n"
              << std::endl;

    InitElectronData(electronData, opmc::constants::k_ElectronCut, opmc::constants::k_GammaCut,
                     opmc::constants::k_MscStepShigh, opmc::constants::k_MscStepSlow, opmc::constants::k_MscEcross);

    InitPhotonData(photonData, opmc::constants::k_ElectronCut, opmc::constants::k_GammaCut);

    std::cout << "\n === All data have been generated and will be written now "
              << "\n     into files under the directory: " << gOutputDataDir << "\n"
              << std::endl;

    electronData.WriteData(gOutputDataDir);
    photonData.WriteData(gOutputDataDir);
}

template <typename T>
const std::vector<G4Material *> &DPMTables<T>::getMaterials() const {
    return m_materials;
}

template <typename T>
const std::string DPMTables<T>::getMaterialName(const int id) const {
    return m_materials[id]->GetName();
}

template <typename T>
void DPMTables<T>::populateMaterials(ElectronData &electronData, PhotonData &photonData) {
    auto i = 0;
    for (const auto &material : m_materials) {
        m_DPMmaterials.emplace_back(electronData, photonData, material->GetName(), i++);
    }
}

template <typename T>
std::vector<G4String> saveMaterialNames(const std::vector<T> &materials) {
    std::vector<G4String> materialNames;
    if constexpr (std::is_same_v<T, G4String> || std::is_same_v<T, std::string>) {
        for (auto &material : materials) {
            materialNames.emplace_back(material);
        }
    } else {
        for (auto &material : materials) {
            materialNames.emplace_back(material->GetName());
        }
    }
    return materialNames;
}

template <typename T>
const ODPMElectronData &DPMTables<T>::getElectronData() const {
    return *m_electronData;
}

template <typename T>
const ODPMPhotonData &DPMTables<T>::getPhotonData() const {
    return *m_photonData;
}

template <typename T>
const std::vector<Material> &DPMTables<T>::DPMmaterials() const {
    return m_DPMmaterials;
}

template std::vector<G4String> saveMaterialNames(const std::vector<G4String> &materials);

template std::vector<G4String> saveMaterialNames(const std::vector<G4Material *> &materials);

template class DPMTables<G4Material *>;

template class DPMTables<G4String>;

template class DPMTables<std::string>;
}  // namespace opmc