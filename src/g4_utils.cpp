

#include "g4_utils.h"

#include <dispatcher.h>

#include <DicomFileCT.hh>
#include <DicomPrimaryGeneratorAction.hh>
#include <DicomRun.hh>
#include <G4ParticleGun.hh>
#include <G4RunManagerFactory.hh>
#include <G4Timer.hh>
#include <G4UImanager.hh>
#include <G4VisExecutive.hh>
#include <G4tgrFileIn.hh>
#include <Shielding.hh>
#include <filesystem>
#include <globals.hh>
#include <iostream>

#include "initialise_tables.h"

// #include <DicomRegularDetectorConstruction.hh>
namespace opmc {

void create_dicom_dat(const std::string &filename, const int nFiles) {
    namespace fs = std::filesystem;
    fs::copy(std::string(DICOM_DATA_DIR) + "/ColourMap.dat", "./ColourMap.dat", fs::copy_options::update_existing);
    fs::copy(std::string(DICOM_DATA_DIR) + "/CT2Density.dat", "./CT2Density.dat", fs::copy_options::update_existing);
    //    fs::copy(std::string(DICOM_DATA_DIR) + "/run.mac", "./run.mac", fs::copy_options::update_existing);
    // open DICOM_DATA_DAT read only
    std::ifstream dicom_dat(std::string(DICOM_DATA_DIR) + "/Data.dat.new");
    if (!dicom_dat.is_open()) {
        throw std::runtime_error(std::string("Error opening") + DICOM_DATA_DIR + "/Data.dat");
    }
    std::ofstream dicom_dat_out("Data.dat");
    if (!dicom_dat_out.is_open()) {
        throw std::runtime_error("Error opening Data.dat");
    }
    std::string line;
    std::getline(dicom_dat, line);
    dicom_dat_out << ":COMPRESSION 1 " << std::endl;
    while (std::getline(dicom_dat, line)) {
        if (line.find("FILE") == std::string::npos) {
            dicom_dat_out << line << std::endl;
        }
    }
    for (int i = 0; i < nFiles; ++i) {
        dicom_dat_out << ":FILE " << filename << std::setw(3) << std::setfill('0') << i + 1 << ".dcm" << std::endl;
    }

    dicom_dat_out << ":FILE_OUT " << filename << ".g4dcm";
}

void create_water_dat(const std::size_t dimx, const std::size_t dimy, const std::size_t dimz,
                      const real_type voxel_size, const std::string &filename) {
    namespace fs = std::filesystem;
    fs::copy(std::string(DATA_DIR) + "/ColourMap.dat", "./ColourMap.dat", fs::copy_options::update_existing);
    fs::copy(std::string(DICOM_DATA_DIR) + "/CT2Density.dat", "./CT2Density.dat", fs::copy_options::update_existing);
    std::ifstream dicom_dat(std::string(DICOM_DATA_DIR) + "/Data.dat.new");
    if (!dicom_dat.is_open()) {
        throw std::runtime_error(std::string("Error opening") + DICOM_DATA_DIR + "/Data.dat");
    }
    std::ofstream dicom_dat_out("Data.dat");
    if (!dicom_dat_out.is_open()) {
        throw std::runtime_error("Error opening Data.dat");
    }
    std::string line;
    std::getline(dicom_dat, line);
    dicom_dat_out << ":COMPRESSION 1 " << std::endl;
    while (std::getline(dicom_dat, line)) {
        if (line.find("FILE") == std::string::npos) {
            dicom_dat_out << line << std::endl;
        }
    }
    dicom_dat_out << ":FILE " << filename << std::endl;
    dicom_dat_out << ":#VOXELS " << dimx << ' ' << dimy << ' ' << dimz << ' ' << voxel_size << std::endl;
    dicom_dat_out << ":FILE_OUT " << filename << ".g4dcm";
}

void initializeEnvironment(const std::string &gDataDirEnvVar, bool overwrite) {
    setenv("G4PARTICLEHPDATA", (gDataDirEnvVar + "/G4TENDL").c_str(), overwrite);
    std::cout << "Setting environment variable G4PARTICLEHPDATA = " << getenv("G4PARTICLEHPDATA") << std::endl;
    setenv("G4NEUTRONHPDATA", (gDataDirEnvVar + "/G4NDL").c_str(), overwrite);
    std::cout << "Setting environment variable G4NEUTRONHPDATA = " << getenv("G4NEUTRONHPDATA") << std::endl;
    setenv("G4LEDATA", (gDataDirEnvVar + "/G4EMLOW").c_str(), overwrite);
    std::cout << "Setting environment variable G4LEDATA = " << getenv("G4LEDATA") << std::endl;
    setenv("G4ENSDFSTATEDATA", (gDataDirEnvVar + "/G4ENSDFSTATE").c_str(), overwrite);
    std::cout << "Setting environment variable G4ENSDFSTATEDATA = " << getenv("G4ENSDFSTATEDATA") << std::endl;
    setenv("G4LEVELGAMMADATA", (gDataDirEnvVar + "/G4PhotonEvaporation").c_str(), overwrite);
    std::cout << "Setting environment variable G4LEVELGAMMADATA = " << getenv("G4LEVELGAMMADATA") << std::endl;
    setenv("G4RADIOACTIVEDATA", (gDataDirEnvVar + "/G4RadioactiveDecay").c_str(), overwrite);
    std::cout << "Setting environment variable G4RADIOACTIVEDATA = " << getenv("G4RADIOACTIVEDATA") << std::endl;
    setenv("G4PARTICLEXSDATA", (gDataDirEnvVar + "/G4PARTICLEXS").c_str(), overwrite);
    std::cout << "Setting environment variable G4PARTICLEXSDATA = " << getenv("G4PARTICLEXSDATA") << std::endl;
    setenv("DICOM_NESTED_PARAM", "1", overwrite);
    std::cout << "Setting environment variable DICOM_NESTED_PARAM = " << getenv("DICOM_NESTED_PARAM") << std::endl;
}

void ODPMDicomPrimaryGeneratorAction::GeneratePrimaries(G4Event *anEvent) {
    fParticleGun.SetParticleDefinition(fParticleDefinition);
    fParticleGun.SetParticleEnergy(fEnergy);
    fParticleGun.SetParticlePosition(fPosition);
    fParticleGun.SetParticleMomentumDirection(fDirection);
    fParticleGun.GeneratePrimaryVertex(anEvent);
}
ODPMDicomPrimaryGeneratorAction::ODPMDicomPrimaryGeneratorAction()
    : fParticleGun(1), fParticleDefinition(nullptr), fPosition(), fDirection(), fEnergy() {}

void ODPMDicomPrimaryGeneratorAction::SetBeam(const G4String &particleName, G4ThreeVector position,
                                              G4ThreeVector direction, const double energy) {
    fParticleDefinition = G4ParticleTable::GetParticleTable()->FindParticle(particleName);
    fPosition           = std::move(position);
    fDirection          = std::move(direction);
    fEnergy             = energy;
}

void OpmcDicomActionInitialization::Build() const {
    SetUserAction(fPrimaryGeneratorActions[G4Threading::G4GetThreadId()]);
    SetUserAction(new ODPMDicomRunAction);
    SetUserAction(new DicomEventAction);
}

void OpmcDicomActionInitialization::SetBeam(const G4String &particleName, G4ThreeVector position,
                                            G4ThreeVector direction, const double energy) {
    for (auto &pga : fPrimaryGeneratorActions) {
        pga->SetBeam(particleName, std::move(position), std::move(direction), energy);
    }
}
OpmcDicomActionInitialization::OpmcDicomActionInitialization(const int nThreads) : fPrimaryGeneratorActions() {
    for (auto i = 0; i < nThreads; ++i) {
        fPrimaryGeneratorActions.emplace_back(new ODPMDicomPrimaryGeneratorAction);
    }
}

ODPMDicomFileMgr::ODPMDicomFileMgr()
    : instance(*DicomFileMgr::GetInstance()), fNoVoxelsX(), fNoVoxelsY(), fNoVoxelsZ(), fDimVox() {}

void ODPMDicomFileMgr::CheckNColumns(std::vector<G4String> wl, size_t vsizeTh) { instance.CheckNColumns(wl, vsizeTh); }

G4int ODPMDicomFileMgr::GetCompression() const { return instance.GetCompression(); }

G4String ODPMDicomFileMgr::GetFileOutName() const { return instance.GetFileOutName(); }

std::vector<G4String> ODPMDicomFileMgr::getMaterials() const { return instance.getMaterials(); }

std::vector<const DicomFileCT *> ODPMDicomFileMgr::getCTs() const { return instance.getCTs(); }

void ODPMDicomFileMgr::SetVoxels(const G4String &fVoxelsX, const G4String &fVoxelsY, const G4String &fVoxelsZ,
                                 const G4String &fDim) {
    fNoVoxelsX = G4UIcommand::ConvertToInt(fVoxelsX);
    fNoVoxelsY = G4UIcommand::ConvertToInt(fVoxelsY);
    fNoVoxelsZ = G4UIcommand::ConvertToInt(fVoxelsZ);
    fDimVox    = G4UIcommand::ConvertToDouble(fDim);
}

void ODPMDicomFileMgr::DumpWaterToTextFile() const {
    G4cout << " DicomFileMgr::Dumping WATER To Text File " << G4endl;  // GDEB

    std::ofstream fout(GetFileOutName());

    G4double fMaxX = 0.5 * fNoVoxelsX * fDimVox;
    G4double fMinX = -0.5 * fNoVoxelsX * fDimVox;
    G4double fMaxY = 0.5 * fNoVoxelsY * fDimVox;
    G4double fMinY = -0.5 * fNoVoxelsY * fDimVox;
    G4double fMaxZ = fNoVoxelsZ * fDimVox;
    G4double fMinZ = 0.0;
    fout << "1" << std::endl;
    fout << "0 \"G4_WATER\"" << std::endl;
    fout << fNoVoxelsX / GetCompression() << " " << fNoVoxelsY / GetCompression() << " " << fNoVoxelsZ << std::endl;
    fout << fMinX << " " << fMaxX << std::endl;
    fout << fMinY << " " << fMaxY << std::endl;
    fout << fMinZ << " " << fMaxZ << std::endl;

    // material values for each voxel:
    for (G4int Zn = 0; Zn < fNoVoxelsZ; ++Zn) {
        for (G4int Yn = 0; Yn < fNoVoxelsY / GetCompression(); ++Yn) {
            for (G4int Xn = 0; Xn < fNoVoxelsX / GetCompression(); ++Xn) {
                fout << "0 ";
            }
            fout << G4endl;
        }
    }

    // density values for each voxel:
    //  0.9982 g/cc or 1.0 g/cm3?
    for (G4int Zn = 0; Zn < fNoVoxelsZ; ++Zn) {
        for (G4int Yn = 0; Yn < fNoVoxelsY / GetCompression(); ++Yn) {
            for (G4int Xn = 0; Xn < fNoVoxelsX / GetCompression(); ++Xn) {
                fout << "1.000 ";
            }
            fout << G4endl;
        }
    }
}

void EraseFileLine(const std::string &path, const std::string &eraseLine) {
    std::string line;
    std::ifstream fin;

    fin.open(path);
    // contents of path must be copied to a temp file then
    // renamed back to the path file
    std::ofstream temp;
    temp.open("temp.txt");

    while (getline(fin, line)) {
        // write all lines to temp other than the line marked for erasing
        if (line.find(eraseLine)) {
            temp << line << std::endl;
        }
    }

    temp.close();
    fin.close();

    // required conversion for remove and rename functions
    const char *p = path.c_str();
    remove(p);
    rename("temp.txt", p);
}

void ODPMDicomFileMgr::Convert(const G4String &fileName) {
    G4tgrFileIn fin = G4tgrFileIn::GetInstance(fileName);
    std::vector<G4String> wl;
    auto waterphantom = false;
    for (int ii = 0;; ii++) {
        if (!fin.GetWordsInLine(wl)) break;
        if (wl[0] == ":FILE") {
            if (wl[1] == "WATER") {
                waterphantom = true;
            }
        } else if (wl[0] == ":#VOXELS") {
            CheckNColumns(wl, 5);
            SetVoxels(wl[1], wl[2], wl[3], wl[4]);
        }
    }
    fin.Close();
    if (waterphantom) {
        EraseFileLine(fileName, ":FILE");
        EraseFileLine(fileName, ":#VOXELS");
    }
    G4tgrFileIn::GetInstance(fileName).OpenNewFile(fileName);
    instance.Convert(fileName);
    //@@@@@@ Process files
    if (waterphantom) {
        DumpWaterToTextFile();
    }
}

VoxelCube convertGeometryToVoxelCube(OpmcDicomDetectorConstruction &theGeometry, ODPMDicomFileMgr &theFileMgr,
                                     std::unique_ptr<DPMTables<G4String>> &tables) {
    const auto &materials = theFileMgr.getMaterials();
    const auto cts        = theFileMgr.getCTs();

    // merge the content of all cts voxels into one flat array
    std::vector<std::uint8_t> voxels(theGeometry.GetTotalVoxels(), -1);
    for (size_t i = 0; i < cts.size(); ++i) {
        auto IDS = cts[i]->GetMaterislIDs();
        std::copy(IDS.begin(), IDS.end(), voxels.begin() + static_cast<std::int64_t>(i * IDS.size()));
    }
    if (cts.empty()) {
        auto water = std::find_if(materials.begin(), materials.end(), [](const auto &m) { return m == "G4_WATER"; });
        if (water == materials.end()) {
            throw std::runtime_error("No water found in the materials vector");
        }
        std::fill(voxels.begin(), voxels.end(), std::distance(materials.begin(), water));
    }

    if (!tables || tables->materialNames != materials) {
        tables = std::make_unique<DPMTables<G4String>>(materials);
    }

    ThreeVector<int> dim{theGeometry.GetNoVoxelsX(), theGeometry.GetNoVoxelsY(), theGeometry.GetNoVoxelsZ()};
    ThreeVector<real_type> voxelSize{theGeometry.GetVoxelHalfX() * 2, theGeometry.GetVoxelHalfY() * 2,
                                     theGeometry.GetVoxelHalfZ() * 2};
    ThreeVector<real_type> origin{theGeometry.GetMinX(), theGeometry.GetMinY(), theGeometry.GetMinZ()};
    ThreeVector<real_type> end{theGeometry.GetMaxX(), theGeometry.GetMaxY(), theGeometry.GetMaxZ()};

    std::cout << "dim = " << dim << std::endl;
    std::cout << "voxelSize = " << voxelSize << std::endl;
    std::cout << "origin = " << origin << std::endl;
    std::cout << "end = " << end << std::endl;
    std::cout << "Z = " << theGeometry.GetNoVoxelsZ() << std::endl;

    const auto &new_materials = tables->getMaterials();

    std::unordered_map<std::size_t, std::uint64_t> materialMap{};

    // populate the materialMap so that it maps the material index in the cts
    // file to the material index in the new_materials vector
    for (size_t i = 0; i < materials.size(); ++i) {
        auto it = std::find_if(new_materials.begin(), new_materials.end(),
                               [&](const auto &material) { return material->GetName() == materials[i]; });
        if (it != new_materials.end()) {
            materialMap[i] = std::distance(new_materials.begin(), it);
            std::cout << "Material " << (*it)->GetName() << " found in the new_materials vector ";
            std::cout << "Mapping [" << i << "]->[" << materialMap[i] << "]" << std::endl;
        } else {
            std::cout << "Material " << materials[i] << " not found in the new_materials vector" << std::endl;
        }
    }

    // update the voxels with the new material indices
    for (auto &v : voxels) {
        v = materialMap[v];
    }
    VoxelCube G4doses{dim, voxelSize, tables->DPMmaterials()};
    std::copy(voxels.begin(), voxels.end(), G4doses.origin_material());
    return G4doses;
}

GammaResults validateDistributions(DoseDistribution<> reference_distribution, DoseDistribution<> opmc_distribution) {
    reference_distribution.filter(.20);
    opmc_distribution.filter(.20);
    const auto [dimx, dimy, dimz]    = static_cast<ThreeVector<int>>(reference_distribution.dimensions());
    const auto [sizex, sizey, sizez] = reference_distribution.voxelSize();
    auto gamma = calculateGamma(0, 3, reference_distribution.data(), opmc_distribution.data(), sizex * .5, sizex, dimx,
                                sizey * .5, sizey, dimy, sizez * .5, sizez, dimz, sizex * .5, sizex, dimx, sizey * .5,
                                sizey, dimy, sizez * .5, sizez, dimz, 1, 1, true, -1, 100);
    std::experimental::mdspan gamma_span(gamma, dimx, dimy, dimz);
    // iterate over the 3d mdspan and find the index of minimum gamma
    auto min_gamma = std::numeric_limits<double>::max();
    ThreeVector<std::int64_t> index;
    auto max_gamma          = std::numeric_limits<double>::min();
    auto sum                = 0.;
    auto nVoxels            = 0;
    auto relative_error     = 0.;
    auto max_relative_error = 0.;
    ThreeVector<std::int64_t> max_index;

    for (auto z = 0; z < dimz; ++z) {
        for (auto y = 0; y < dimy; ++y) {
            for (auto x = 0; x < dimx; ++x) {
                if (reference_distribution[x, y, z] == 0 && opmc_distribution[x, y, z] == 0) {
                    continue;
                }
                if (gamma_span[x, y, z] < min_gamma) {
                    min_gamma = gamma_span[x, y, z];
                    index     = {x, y, z};
                }
                if (gamma_span[x, y, z] > max_gamma) {
                    max_gamma = gamma_span[x, y, z];
                    max_index = {x, y, z};
                }
                const auto current_relative_error =
                    std::abs(reference_distribution[x, y, z] - opmc_distribution[x, y, z]) /
                    std::max(reference_distribution[x, y, z], opmc_distribution[x, y, z]);
                relative_error += current_relative_error;
                max_relative_error = std::max(max_relative_error, current_relative_error);
                sum += gamma_span[x, y, z];
                nVoxels++;
            }
        }
    }

    std::cout << "Gamma pass rate = " << sum / nVoxels << std::endl;
    return {relative_error / nVoxels, max_relative_error, sum / nVoxels, min_gamma};
}

std::ostream &operator<<(std::ostream &os, const GammaResults &results) {
    os << "maxRelativeError: " << results.maxRelativeError << " averageRelativeError: " << results.averageRelativeError
       << " gammaPassingRate: " << results.gammaPassingRate << " minGamma: " << results.minGamma;
    return os;
}

ThreeVector<real_type> getCenter(const OpmcDicomDetectorConstruction &theGeometry) {
    ThreeVector<real_type> dim{theGeometry.GetNoVoxelsX(), theGeometry.GetNoVoxelsY(), theGeometry.GetNoVoxelsZ()};
    ThreeVector<real_type> voxelSize{theGeometry.GetVoxelHalfX() * 2, theGeometry.GetVoxelHalfY() * 2,
                                     theGeometry.GetVoxelHalfZ() * 2};
    ThreeVector<real_type> origin{theGeometry.GetMinX(), theGeometry.GetMinY(), theGeometry.GetMinZ()};
    ThreeVector<real_type> end{theGeometry.GetMaxX(), theGeometry.GetMaxY(), theGeometry.GetMaxZ()};

    std::cout << "dim = " << dim << std::endl;
    std::cout << "voxelSize = " << voxelSize << std::endl;
    std::cout << "origin = " << origin << std::endl;
    std::cout << "end = " << end << std::endl;
    std::cout << "Z = " << theGeometry.GetNoVoxelsZ() << std::endl;

    auto center = dim - end;
    center.z    = std::abs(origin.z);
    std::cout << "center = " << center << std::endl;
    return center;
}
}  // namespace opmc