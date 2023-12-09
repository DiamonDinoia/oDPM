//
// Created by mbarbone on 5/17/22.
//

#include "electron_tables.h"

#include <ElectronData.hh>

namespace opmc {

static auto getMollerNumPrimaryEnergies(const ElectronData &data) {
    return data.GetMollerTables().getFNumSamplingPrimaryEnergies();
}

static auto getMollerNumSecondaryEnergies(const ElectronData &data) {
    return data.GetMollerTables().getFNumSamplingSecondaryEnergies();
}

static auto getMollerMinPrimaryEnergy(const ElectronData &data) {
    return data.GetMollerTables().getFMinPrimaryEnergy();
}

static auto getMollerInvLogDeltaPrimaryEnergy(const ElectronData &data) {
    return data.GetMollerTables().getFInvLogDeltaPrimaryEnergy();
}

static auto getMollerData(const ElectronData &data) {
    const auto primaryEnergies   = getMollerNumPrimaryEnergies(data);
    const auto secondaryEnergies = getMollerNumSecondaryEnergies(data);
    std::vector<Point4D> moller;
    const auto &dataMoller = data.GetMollerTables().getFTables();
    for (auto i = 0; i < primaryEnergies; ++i) {
        for (int j = 0; j < secondaryEnergies; ++j) {
            const auto x   = dataMoller[i]->fXdata[j];
            const auto y   = dataMoller[i]->fYdata[j];
            const auto w   = dataMoller[i]->fAliasW[j];
            const auto ind = dataMoller[i]->fAliasIndx[j];
            assert(ind < secondaryEnergies);
            moller.emplace_back(x, y, w, ind);
        }
    }
    return moller;
}

static auto getSeltzerBergerNumPrimaryEnergies(const ElectronData &data) {
    return data.getFsbTables().getFNumSamplingElecEnergies();
}

static auto getSeltzerBergerNumSecondaryEnergies(const ElectronData &data) {
    return data.getFsbTables().getFNumSamplingPhotEnergies();
}

static auto getSeltzerBergerMinPrimaryEnergy(const ElectronData &data) {
    return data.getFsbTables().getFMinElecEnergy();
}

static auto getSeltzerBergerInvLogDeltaPrimaryEnergy(const ElectronData &data) {
    return data.getFsbTables().getFElEnIlDelta();
}

static auto getSeltzerBergerData(const ElectronData &data, const int id) {
    std::vector<Point4D> seltzerBerger;
    const auto &dataSeltzerBerger = data.getFsbTables().getFTablesPerMaterial(id);
    for (auto i = 0; i < getSeltzerBergerNumPrimaryEnergies(data); ++i) {
        for (int j = 0; j < getSeltzerBergerNumSecondaryEnergies(data); ++j) {
            const auto x   = dataSeltzerBerger.fTheTables[i]->fXdata[j];
            const auto y   = dataSeltzerBerger.fTheTables[i]->fYdata[j];
            const auto w   = dataSeltzerBerger.fTheTables[i]->fAliasW[j];
            const auto ind = dataSeltzerBerger.fTheTables[i]->fAliasIndx[j];
            seltzerBerger.emplace_back(x, y, w, ind);
        }
    }
    return seltzerBerger;
}

SeltzerBerger::SeltzerBerger(const ElectronData &data, const int id)
    : WalkersApproximator(getSeltzerBergerNumPrimaryEnergies(data), getSeltzerBergerNumSecondaryEnergies(data),
                          getSeltzerBergerMinPrimaryEnergy(data), getSeltzerBergerInvLogDeltaPrimaryEnergy(data),
                          getSeltzerBergerData(data, id)) {}

ODPMElectronData::ODPMElectronData(const ElectronData &data)
    : m_mollerEnergyTransfer(getMollerNumPrimaryEnergies(data), getMollerNumSecondaryEnergies(data),
                             getMollerMinPrimaryEnergy(data), getMollerInvLogDeltaPrimaryEnergy(data),
                             getMollerData(data)) {}

}  // namespace opmc