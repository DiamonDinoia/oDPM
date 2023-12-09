//
// Created by mbarbone on 4/12/22.
//

#include "materials.h"

#include <sstream>

#include "constants.h"


#include <ElectronData.hh>
#include <PhotonData.hh>


namespace opmc {


    static auto getITr1MFPElastic(const ElectronData &data, const int id) {
        auto iTr1MFPElastic = data.GetITr1MFPElasticPerMaterial(id);
        std::vector<Point5D> iTr1MFPElasticSpline;
        for (const auto &[x, y]: iTr1MFPElastic) {
            iTr1MFPElasticSpline.emplace_back(x, y, 0, 0, 0);
        }
        return iTr1MFPElasticSpline;
    }

    static auto getMaxScatteringStrength(const ElectronData &data, const int id) {
        auto maxScatteringStrength = data.GetMaxScatStrengthPerMaterial(id);
        std::vector<Point5D> maxScatteringStrengthSpline;
        for (const auto &[x, y]: maxScatteringStrength) {
            maxScatteringStrengthSpline.emplace_back(x, y, 0, 0, 0);
        }
        return maxScatteringStrengthSpline;
    }

    static auto getIMFPMoller(const ElectronData &data, const int id) {
        auto imfpMoller = data.GetIMFPMollerPerMaterial(id);
        std::vector<Point5D> imfpMollerSpline;
        for (const auto &[x, y]: imfpMoller) {
            imfpMollerSpline.emplace_back(x, y, 0, 0, 0);
        }
        return imfpMollerSpline;
    }

    static auto getIMFPBrem(const ElectronData &data, const int id) {
        auto imfpBrem = data.GetIMFPBremPerMaterial(id);
        std::vector<Point5D> imfpBremSpline;
        for (const auto &[x, y]: imfpBrem) {
            imfpBremSpline.emplace_back(x, y, 0, 0, 0);
        }
        return imfpBremSpline;
    }

    static auto getStoppingPower(const ElectronData &data, const int id) {
        auto stoppingPower = data.GetStoppingPowerPerMaterial(id);
        std::vector<Point5D> stoppingPowerSpline;
        for (const auto &[x, y]: stoppingPower) {
            stoppingPowerSpline.emplace_back(x, y, 0, 0, 0);
        }
        return stoppingPowerSpline;
    }


    static auto getIMFPTotalMaxPhoton(const PhotonData &data, const int id) {
        std::vector<Point2D> imfpTotal;
        for (const auto &[x, y]: data.GetIMFPTotal(id)) {
            imfpTotal.emplace_back(x, y);
        }
        return imfpTotal;
    }

    static auto getIMFPComptonMaxPhoton(const PhotonData &data, const int id) {
        std::vector<Point2D> imfpCompton;
        for (const auto &[x, y]: data.GetIMFPCompton(id)) {
            imfpCompton.emplace_back(x, y);
        }
        return imfpCompton;
    }

    static auto getIMFPPairProdPhoton(const PhotonData &data, const int id) {
        std::vector<Point2D> imfpPhotoelectric;
        for (const auto &[x, y]: data.GetIMFPPairProd(id)) {
            imfpPhotoelectric.emplace_back(x, y);
        }
        return imfpPhotoelectric;
    }

    Material::Material(const ElectronData &electronData, const PhotonData &photonData,
                       const std::string &name,
                       const int id) : k_id(id),
                                       name(name),
                                       m_MollerIMFPScaling(electronData.GetMollerIMFPScaling(id)),
          m_density(electronData.GetMaterialDensity(id)),
                                       m_ITr1MFPElastic(getITr1MFPElastic(electronData, id)),
                                       m_MaxScatStrength(getMaxScatteringStrength(electronData, id)),
                                       m_IMFPMoller(getIMFPMoller(electronData, id)),
                                       m_IMFPBrem(getIMFPBrem(electronData, id)),
                                       m_StoppingPower(getStoppingPower(electronData, id)),
                                       m_SeltzerBerger(electronData, id),
                                       m_GoudsmitSaunderson(electronData, id),
                                       m_IMFPTotalPhoton(getIMFPTotalMaxPhoton(photonData, id)),
                                       m_IMFPComptonPhoton(getIMFPComptonMaxPhoton(photonData, id)),
                                       m_IMFPPairProdPhoton(getIMFPPairProdPhoton(photonData, id)) {}


    const std::string &Material::getName() const {
        return name;
    }

}  // namespace opmc
