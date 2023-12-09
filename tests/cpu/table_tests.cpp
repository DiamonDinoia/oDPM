//
// Created by mbarbone on 5/4/22.
//

#include "test_utils.h"
#include "initialise_tables.h"

#include <SimElectronData.hh>
#include <SimGSTables.hh>
#include <SimIMFPBrem.hh>
#include <SimIMFPMaxPhoton.hh>
#include <SimIMFPMoller.hh>
#include <SimIMFPPhoton.hh>
#include <SimITr1MFPElastic.hh>
#include <SimKNTables.hh>
#include <SimMaterialData.hh>
#include <SimMaxScatStrength.hh>
#include <SimMollerTables.hh>
#include <SimPhotonData.hh>
#include <SimSBTables.hh>
#include <SimStoppingPower.hh>
#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>
#include <functional>

#include <random>

using Catch::Approx;

namespace opmc {

    static const auto tests = 16384 << 4;

    static void test_table(auto &reference, auto &spline,
                           auto &input_generator) {
        for (int i = 0; i < tests; ++i) {
            const auto energy = input_generator();
            INFO("Test " << i);
            INFO("Energy " << energy);
            // epsilon required otherwise moller tables fail due to rounding error
            REQUIRE(reference(energy) == Approx(spline(energy)).epsilon(1e-4));
        }
    }

    static const std::vector<std::string> test_materials{
            "G4_WATER",
    };
    const auto seed = std::random_device()();
    std::default_random_engine rng{seed};
    DPMTables tables{test_materials};
    TestMaterial material{tables.getElectronData(), tables.getPhotonData(), tables.getMaterialName(0), 0};
    TestPhotonData photonData{tables.getPhotonData()};
    TestElectronData electronData{tables.getElectronData()};
    SimElectronData elData{};
    SimPhotonData phData{};
    SimMaterialData mData{};

    bool init = false;

    TEST_CASE("INTERPOLATION TESTS", "[tables]" ) {
        if (!init) {
            INFO("SEED " << seed);
            elData.Load(OUTPUT_DATA_DIR, 0);
            phData.Load(OUTPUT_DATA_DIR, 0);
            mData.Load(OUTPUT_DATA_DIR, 0);
            init = true;
        }

        SECTION("ITr1MFPElastic") {
            std::uniform_real_distribution<real_type> distribution{material.GetITr1MFPElastic().getMinInput(),
                                                                   material.GetITr1MFPElastic().getMaxInput()};
            auto input_generator = [&distribution]() -> real_type { return distribution(rng); };
            auto reference = [](real_type x) -> real_type {
                return elData.GetITr1MFPElastic()->GetITr1MFPPerDensity(x, material.id());
            };
            auto spline = [](real_type x) -> real_type {
                return material.iTr1MFPElasticElectron(x);
            };
            test_table(reference, spline, input_generator);
        }SECTION("MaxScatStrength") {
            std::uniform_real_distribution<real_type> distribution{material.GetMaxScatStrength().getMinInput(),
                                                                   material.GetMaxScatStrength().getMaxInput()};
            std::function<real_type()> input_generator = [&distribution]() -> real_type {
                return distribution(rng);
            };
            auto reference = [](real_type x) -> real_type {
                return elData.GetMaxScatStrength()->GetMaxScatStrength(x);
            };
            auto spline = [](real_type x) -> real_type {
                return material.maxScatteringStrengthElectron(x);
            };
            test_table(reference, spline, input_generator);
        }SECTION("IMFPMoller") {
            std::uniform_real_distribution<real_type> distribution{material.GetIMFPMoller().getMinInput(),
                                                                   material.GetIMFPMoller().getMaxInput()};
            std::function<real_type()> input_generator = [&distribution]() -> real_type {
                return distribution(rng);
            };
            auto reference = [](real_type x) -> real_type {
                return elData.GetIMFPMoller()->GetIMFPPerDensity(x);
            };
            auto spline = [](real_type x) -> real_type {
                return material.iMFPMollerElectron(x);
            };
            test_table(reference, spline, input_generator);
        }SECTION("IMFPBrem") {
            std::uniform_real_distribution<real_type> distribution{material.GetIMFPBrem().getMinInput(),
                                                                   material.GetIMFPBrem().getMaxInput()};
            std::function<real_type()> input_generator = [&distribution]() -> real_type {
                return distribution(rng);
            };
            auto reference = [](real_type x) -> real_type {
                return elData.GetIMFPBrem()->GetIMFPPerDensity(x, material.id());
            };
            auto spline = [](real_type x) -> real_type {
                return material.iMFPBremElectron(x);
            };
            test_table(reference, spline, input_generator);
        }SECTION("StoppingPower") {
            std::uniform_real_distribution<real_type> distribution{material.GetStoppingPower().getMinInput(),
                                                                   material.GetStoppingPower().getMaxInput()};
            std::function<real_type()> input_generator = [&distribution]() -> real_type {
                return distribution(rng);
            };
            auto reference = [](real_type x) -> real_type {
                return elData.GetDEDX()->GetDEDXPerDensity(x, material.id());
            };
            auto spline = [](real_type x) -> real_type {
                return material.stoppingPowerElectron(x);
            };
            test_table(reference, spline, input_generator);
        }SECTION("SeltzerBerger") {
            std::uniform_real_distribution<real_type> energy_distribution{constants::k_GammaCut, constants::k_MaxEkin};
            std::uniform_real_distribution<real_type> uniform_distribution{};
            for (int i = 0; i < tests; ++i) {
                const auto energy = energy_distribution(rng);
                const auto rng1 = uniform_distribution(rng);
                const auto rng2 = uniform_distribution(rng);
                const auto rng3 = uniform_distribution(rng);
                INFO("Energy " << energy);
                INFO("rng1 " << rng1);
                INFO("rng2 " << rng2);
                INFO("rng3 " << rng3);
                INFO("test: " << i);
                const auto reference =
                        elData.GetTheSBTables()->SampleEnergyTransfer(energy, material.id(), rng1, rng2, rng3);
                const auto result = material.seltzerBergerElectron(energy, rng1, rng2, rng3);
                REQUIRE(reference == Approx(result));
            }
        }SECTION("MollerEnergyTransfer") {
            std::uniform_real_distribution<real_type> energy_distribution{constants::k_ElectronCut * 2,
                                                                          constants::k_MaxEkin};
            std::uniform_real_distribution<real_type> uniform_distribution{};
            for (int i = 0; i < tests; ++i) {
                const auto energy = energy_distribution(rng);
                const auto rng1 = uniform_distribution(rng);
                const auto rng2 = uniform_distribution(rng);
                const auto rng3 = uniform_distribution(rng);
                INFO("Energy " << energy);
                INFO("rng1 " << rng1);
                INFO("rng2 " << rng2);
                INFO("rng3 " << rng3);
                INFO("test: " << i);
                const auto reference = elData.GetTheMollerTables()->SampleEnergyTransfer(energy, rng1, rng2, rng3);
                const auto result = electronData.mollerEnergyTransfer(energy, rng1, rng2, rng3);
                REQUIRE(reference == Approx(result));
            }
        }SECTION("GoudsmitSaunderson") {
            std::uniform_real_distribution<real_type> energy_distribution{constants::k_ElectronCut,
                                                                          constants::k_MaxEkin};
            std::uniform_real_distribution<real_type> uniform_distribution{};
            for (int i = 0; i < tests; ++i) {
                const auto energy = energy_distribution(rng);
                const auto rng1 = uniform_distribution(rng);
                const auto rng2 = uniform_distribution(rng);
                INFO("Energy " << energy);
                INFO("rng1 " << rng1);
                INFO("rng2 " << rng2);
                INFO("test: " << i);
                const auto reference = elData.GetTheGSTables()->SampleAngularDeflection(energy, rng1, rng2);
                const auto result = material.angularDeflectionElectron(energy, rng1, rng2);
                REQUIRE(reference == Approx(result));
            }
        }SECTION("IMFPTotalPhoton") {
            std::uniform_real_distribution<real_type> distribution{material.GetIMFPTotalPhoton().GetMinInput(),
                                                                   material.GetIMFPTotalPhoton().GetMaxInput()};
            auto input_generator = [&distribution]() -> real_type {
                return distribution(rng);
            };
            auto reference = [](real_type x) -> real_type {
                return phData.GetIMFPTotal()->GetIMFPPerDensity(x, material.id());
            };
            auto interpolation = [](real_type x) -> real_type {
                return material.iMFPTotalPhoton(x);
            };
            test_table(reference, interpolation, input_generator);
        }SECTION("IMFPComptonPhoton") {
            std::uniform_real_distribution<real_type> distribution{material.GetIMFPComptonPhoton().GetMinInput(),
                                                                   material.GetIMFPComptonPhoton().GetMaxInput()};
            auto input_generator = [&distribution]() -> real_type {
                return distribution(rng);
            };
            auto reference = [](real_type x) -> real_type {
                return phData.GetIMFPCompton()->GetIMFPPerDensity(x, material.id());
            };
            auto interpolation = [](real_type x) -> real_type {
                return material.iMFPComptonPhoton(x);
            };
            test_table(reference, interpolation, input_generator);
        }SECTION("IMFPPairProdPhoton") {
            std::uniform_real_distribution<real_type> distribution{material.GetIMFPPairProdPhoton().GetMinInput(),
                                                                   material.GetIMFPPairProdPhoton().GetMaxInput()};
            auto input_generator = [&distribution]() -> real_type {
                return distribution(rng);
            };
            auto reference = [](real_type x) -> real_type {
                return phData.GetIMFPPairProd()->GetIMFPPerDensity(x, material.id());
            };
            auto interpolation = [](real_type x) -> real_type {
                return material.iMFPPairProdPhoton(x);
            };
            test_table(reference, interpolation, input_generator);
        }SECTION("IMFPMaxPhoton") {
            std::uniform_real_distribution<real_type> distribution{photonData.GetIMFPMaxPhoton().GetMinInput(),
                                                                   photonData.GetIMFPMaxPhoton().GetMaxInput()};
            auto input_generator = [&distribution]() -> real_type {
                return distribution(rng);
            };
            auto reference = [](real_type x) -> real_type {
                return phData.GetIMFPTotalMax()->GetIMFP(x);
            };
            auto interpolation = [](real_type x) -> real_type {
                return photonData.iMFPMaxPhoton(x);
            };
            test_table(reference, interpolation, input_generator);
        }
        SECTION("KleinNishina") {
            std::uniform_real_distribution<real_type> energy_distribution{constants::k_GammaCut, constants::k_MaxEkin};
            std::uniform_real_distribution<real_type> uniform_distribution{};
            for (int i = 0; i < tests; ++i) {
                const auto energy = energy_distribution(rng);
                assert(energy >= constants::k_GammaCut);
                const auto rng1 = uniform_distribution(rng);
                const auto rng2 = uniform_distribution(rng);
                const auto rng3 = uniform_distribution(rng);
                INFO("Energy " << energy);
                INFO("rng1 " << rng1);
                INFO("rng2 " << rng2);
                INFO("rng3 " << rng3);
                INFO("test: " << i);
                const auto reference = phData.GetTheKNTables()->SampleEnergyTransfer(energy, rng1, rng2, rng3);
                const auto result = photonData.kleinNishina(energy, rng1, rng2, rng3);
                REQUIRE(reference == Approx(result));
            }
        }
    }

}  // namespace opmc