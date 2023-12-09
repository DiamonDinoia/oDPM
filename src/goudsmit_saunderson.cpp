//
// Created by mbarbone on 5/11/22.
//

#include "goudsmit_saunderson.h"

#include <ElectronData.hh>

namespace opmc {

static auto getNumEkin(const ElectronData &electronData, const int id) {
    return electronData.GetGSDtrData(id).fNumEKin;
}

static auto getNumCumData(const ElectronData &electronData, const int id) {
    return electronData.GetGSDtrData(id).fNumCumData;
}

static auto getMinEkin(const ElectronData &electronData, const int id) {
    return electronData.GetGSDtrData(id).fEKinGrid[0];
}

static auto getEkingGrid(const ElectronData &electronData, const int id) {
    return electronData.GetGSDtrData(id).fEKinGrid;
}

static auto getLogEKinMin(const ElectronData &electronData, const int id) {
    return math::log(electronData.GetGSDtrData(id).fEKinGrid[0]);
}

static auto getInvLofDelta(const ElectronData &electronData, const int id) {
    return 1. / (math::log(electronData.GetGSDtrData(id).fEKinGrid[1]) -
                 math::log(electronData.GetGSDtrData(id).fEKinGrid[0]));
}

static auto getInvDelaCum(const ElectronData &electronData, const int id) {
    return 1. / (getNumCumData(electronData, id) - 1);
}

static auto getTransformParams(const ElectronData &electronData, const int id) {
    std::vector<double> data;
    for (const auto &onePoint : electronData.GetGSDtrData(id).fGSDtrData) {
        data.emplace_back(onePoint.fTransformParam);
    }
    return data;
}

static auto getGSDtrData(const ElectronData &electronData, const int id) {
    std::vector<GoudsmitSaunderson::OnePoint> data;
    for (const auto &onePoint : electronData.GetGSDtrData(id).fGSDtrData) {
        for (const auto &[u, a, b] : onePoint.fCumData) {
            data.emplace_back(u, a, b);
        }
    }
    return data;
}

GoudsmitSaunderson::GoudsmitSaunderson(const ElectronData &electronData, const int id)
    : GoudsmitSaunderson::GoudsmitSaunderson(
          getNumEkin(electronData, id), getNumCumData(electronData, id), getMinEkin(electronData, id),
          getEkingGrid(electronData, id), getLogEKinMin(electronData, id), getInvLofDelta(electronData, id),
          getInvDelaCum(electronData, id), getTransformParams(electronData, id), getGSDtrData(electronData, id)) {}

}  // namespace opmc