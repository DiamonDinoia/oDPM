//
// Created by mbarbone on 5/11/22.
//

#include "photon_tables.h"

#include <PhotonData.hh>

namespace opmc {

static auto getIMFPMax(const PhotonData &data) {
    std::vector<Point2D> iMFPMax;
    for (const auto &[x, y] : data.GetIMFPTotalMax()) {
        iMFPMax.emplace_back(x, y);
    }
    return iMFPMax;
}

static auto getKleinNishinaNumPrimaryEnergies(const PhotonData &data) {
    return data.getFknTables().getFNumSamplingPrimaryEnergies();
}

static auto getKleinNishinaNumSecondaryEnergies(const PhotonData &data) {
    return data.getFknTables().getFNumSamplingSecondaryEnergies();
}

static auto getKleinNishinaMinPrimaryEnergy(const PhotonData &data) {
    return data.getFknTables().getFMinPrimaryEnergy();
}

static auto getKleinNishinaInvLogDeltaPrimaryEnergy(const PhotonData &data) {
    return data.getFknTables().getFInvLogDeltaPrimaryEnergy();
}

static auto getKleinNishinaData(const PhotonData &data) {
    std::vector<Point4D> kleinNishina;
    const auto &dataKleinNishina = data.getFknTables().getFTables();
    for (auto i = 0; i < getKleinNishinaNumPrimaryEnergies(data); ++i) {
        for (int j = 0; j < getKleinNishinaNumSecondaryEnergies(data); ++j) {
            const auto x   = dataKleinNishina[i]->fXdata[j];
            const auto y   = dataKleinNishina[i]->fYdata[j];
            const auto w   = dataKleinNishina[i]->fAliasW[j];
            const auto ind = dataKleinNishina[i]->fAliasIndx[j];
            kleinNishina.emplace_back(x, y, w, ind);
        }
    }
    return kleinNishina;
}

ODPMPhotonData::ODPMPhotonData(const PhotonData &data)
    : m_IMFPMaxPhoton(getIMFPMax(data)),
      m_KleinNishina(getKleinNishinaNumPrimaryEnergies(data), getKleinNishinaNumSecondaryEnergies(data),
                     getKleinNishinaMinPrimaryEnergy(data), getKleinNishinaInvLogDeltaPrimaryEnergy(data),
                     getKleinNishinaData(data)) {}

}  // namespace opmc
