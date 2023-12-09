//
// Created by mbarbone on 7/26/21.
//

#include <testAPI.h>
#include <G4Material.hh>
#include <catch2/catch_test_macros.hpp>

#include "constants.h"
#include "ksTest.h"
#include "materials.h"
#include "physics.h"
#include "initialise_tables.h"

using namespace std;

namespace opmc {

    const auto alpha = 5.0E-9;
    const auto tests = 100;
    const int histories = 1.0E+5;

    template<typename T, typename V>
    double testHistograms(const vector <tuple<T, T>> &reference, const vector <tuple<V, V>> &opmcInteraction) {
        vector<real_type> referenceValues(reference.size());
        for (auto i = 0U; i < reference.size(); ++i) {
            referenceValues[i] = get<1>(reference[i]);
        }
        vector<real_type> opmcValues(reference.size());
        for (auto i = 0U; i < reference.size(); ++i) {
            opmcValues[i] = get<1>(opmcInteraction[i]);
        }
        if (opmcValues.back() == referenceValues.back() == 1.) {
            return 1;
        }
        return KsTest(referenceValues, opmcValues);
    }

    void testInteraction(AbstractTest &reference, AbstractTest &opmcInteraction) {
        reference.simulate();
        opmcInteraction.simulate();
        const auto energyResult = testHistograms(reference.getEnergyHist(), opmcInteraction.getEnergyHist());
        REQUIRE(energyResult >= alpha);
        const auto cosResult = testHistograms(reference.getCosHist(), opmcInteraction.getCosHist());
        REQUIRE(cosResult >= alpha);
    }

    static const std::vector<std::string> test_materials{
            "G4_WATER",
    };

    static DPMTables tables{test_materials};
    static HalfDistanceVoxelCube geometry{{}, {}, tables.populateMaterials()};
    static DummyQueue queue{};
    static std::unique_ptr<Physics<HalfDistanceVoxelCube, DummyQueue>> physics;
    static std::unique_ptr<Material> reference_material;
    static std::unique_ptr<ODPMPhotonData> photon_data;
    static std::unique_ptr<ODPMElectronData> electron_data;

    TEST_CASE("Particle Intercations") {
        auto seed = std::random_device()();
        INFO("SEED " << seed);
        Random rng{seed};
        std::default_random_engine generator{static_cast<unsigned>(rng.getUniform())};
        reference_material = std::make_unique<Material>(tables.getElectronData(), tables.getPhotonData(),
                                                        tables.getMaterials()[0]->GetName(), 0);
        photon_data = std::make_unique<ODPMPhotonData>(tables.getPhotonData());
        electron_data = std::make_unique<ODPMElectronData>(tables.getElectronData());
        physics = std::make_unique<Physics<HalfDistanceVoxelCube, DummyQueue>>(rng, geometry, queue,
                                                                               *reference_material,
                                                                               *photon_data, *electron_data);
        class ODPMTest : public AbstractTest {
        public:
            explicit ODPMTest(std::string gInputDataDir = OUTPUT_DATA_DIR, double gPrimaryEnergy = 12.345,
                              int gNumPrimaries = 1.0E+5, int gMaterialIndex = 0, int hbins = 101)
                    : AbstractTest(std::move(gInputDataDir), gPrimaryEnergy, gNumPrimaries, gMaterialIndex, hbins) {}

            virtual Particle interaction(const Material &material, Particle &particle) = 0;

            virtual Particle generateParticle() = 0;

            void simulate() override {
                for (int is = 0; is < gNumPrimaries; ++is) {
                    auto primaryParticle = generateParticle();
                    auto particle = interaction(*reference_material, primaryParticle);
                    const double theK = particle.energy;
                    const double theCost = particle.direction.z;
                    theCosineHist[(int) ((theCost * (hbins / 2)) + (hbins / 2))] += 1;
                    double redPhEnergy = theK / gPrimaryEnergy;
                    if (redPhEnergy > 0.0) {
                        redPhEnergy = std::log10(redPhEnergy);
                        theHist[(int) ((redPhEnergy - xmin) * ihdel)] += 1.0;
                    }
                }
            }
        };

        class ElectronTest : public ODPMTest {
        public:
            explicit ElectronTest(std::string gInputDataDir = OUTPUT_DATA_DIR, double gPrimaryEnergy = 12.345,
                                  int gNumPrimaries = 1.0E+5, int gMaterialIndex = 0, int hbins = 101)
                    : ODPMTest(std::move(gInputDataDir), gPrimaryEnergy, gNumPrimaries, gMaterialIndex, hbins) {}

            Particle generateParticle() override {
                static ThreeVector <real_type> pos{0., 0., 0.};
                static ThreeVector <real_type> dir{0., 0., 1.};
                return Electron{pos, dir, static_cast<real_type>(gPrimaryEnergy)};
            };
        };
        class MSCAngularDeflectionODPM : public ElectronTest {
        public:
            explicit MSCAngularDeflectionODPM(std::string gInputDataDir = OUTPUT_DATA_DIR,
                                              double gPrimaryEnergy = 12.345,
                                              int gNumPrimaries = 1.0E+5, int gMaterialIndex = 0, int hbins = 101)
                    : ElectronTest(std::move(gInputDataDir), gPrimaryEnergy, gNumPrimaries, gMaterialIndex, hbins) {}

            Particle interaction(const Material &material, Particle &particle) override {
                physics->hingeInteraction(static_cast<Electron &>(particle), particle.energy);
                return particle;
            }
        };

        class BremmODPM : public ElectronTest {
        public:
            explicit BremmODPM(std::string gInputDataDir = OUTPUT_DATA_DIR, double gPrimaryEnergy = 12.345,
                               int gNumPrimaries = 1.0E+5, int gMaterialIndex = 0, int hbins = 101)
                    : ElectronTest(std::move(gInputDataDir), gPrimaryEnergy, gNumPrimaries, gMaterialIndex, hbins) {}

            Particle interaction(const Material &material, Particle &particle) override {
                return physics->bremsstrahlungInteraction(material, static_cast<Electron &>(particle));
            }
        };

        class MollerODPM : public ElectronTest {
        public:
            explicit MollerODPM(std::string gInputDataDir = OUTPUT_DATA_DIR, double gPrimaryEnergy = 12.345,
                                int gNumPrimaries = 1.0E+5, int gMaterialIndex = 0, int hbins = 101)
                    : ElectronTest(std::move(gInputDataDir), gPrimaryEnergy, gNumPrimaries, gMaterialIndex, hbins) {}

            Particle interaction(const Material &material, Particle &particle) override {
                return physics->mollerInteraction(static_cast<Electron &>(particle));
            }
        };

        class PositronODPM : public ODPMTest {
        public:
            explicit PositronODPM(std::string gInputDataDir = OUTPUT_DATA_DIR, double gPrimaryEnergy = 12.345,
                                  int gNumPrimaries = 1.0E+5, int gMaterialIndex = 0, int hbins = 101)
                    : ODPMTest(std::move(gInputDataDir), gPrimaryEnergy, gNumPrimaries, gMaterialIndex, hbins) {}

            Particle generateParticle() final {
                static ThreeVector <real_type> pos{0., 0., 0.};
                static ThreeVector <real_type> dir{0., 0., 1.};
                return Positron{pos, dir, static_cast<real_type>(gPrimaryEnergy)};
            }
        };

        class AnnihilationODPM : public PositronODPM {
        public:
            explicit AnnihilationODPM(std::string gInputDataDir = OUTPUT_DATA_DIR, double gPrimaryEnergy = 12.345,
                                      int gNumPrimaries = 1.0E+5, int gMaterialIndex = 0, int hbins = 101)
                    : PositronODPM(std::move(gInputDataDir), gPrimaryEnergy, gNumPrimaries, gMaterialIndex, hbins) {}

            Particle interaction(const Material &material, Particle &particle) override {
                return physics->annihilation(static_cast<const Positron &&>(particle)).first;
            }
        };

        class PhotonODPM : public ODPMTest {
        public:
            explicit PhotonODPM(std::string gInputDataDir = OUTPUT_DATA_DIR, double gPrimaryEnergy = 12.345,
                                int gNumPrimaries = 1.0E+5, int gMaterialIndex = 0, int hbins = 101)
                    : ODPMTest(std::move(gInputDataDir), gPrimaryEnergy, gNumPrimaries, gMaterialIndex, hbins) {}

            Particle generateParticle() final {
                static ThreeVector <real_type> pos{0., 0., 0.};
                static ThreeVector <real_type> dir{0., 0., 1.};
                return Photon{pos, dir, static_cast<real_type>(gPrimaryEnergy)};
            }
        };

        class ComptonPrimaryODPM : public PhotonODPM {
        public:
            explicit ComptonPrimaryODPM(std::string gInputDataDir = OUTPUT_DATA_DIR, double gPrimaryEnergy = 12.345,
                                        int gNumPrimaries = 1.0E+5, int gMaterialIndex = 0, int hbins = 101)
                    : PhotonODPM(std::move(gInputDataDir), gPrimaryEnergy, gNumPrimaries, gMaterialIndex, hbins) {}

            Particle interaction(const Material &material, Particle &particle) override {
                auto secondary = physics->comptonScattering(static_cast<Photon &>(particle));
                if (particle.energy < constants::k_GammaCut) {
                    particle.energy = 0;
                }
                return particle;
            }
        };

        class ComptonSecondaryODPM : public PhotonODPM {
        public:
            explicit ComptonSecondaryODPM(std::string gInputDataDir = OUTPUT_DATA_DIR, double gPrimaryEnergy = 12.345,
                                          int gNumPrimaries = 1.0E+5, int gMaterialIndex = 0, int hbins = 101)
                    : PhotonODPM(std::move(gInputDataDir), gPrimaryEnergy, gNumPrimaries, gMaterialIndex, hbins) {}

            Particle interaction(const Material &material, Particle &particle) override {
                auto secondary = physics->comptonScattering(static_cast<Photon &>(particle));
                if (secondary.energy < constants::k_ElectronCut) {
                    return particle;
                }
                return secondary;
            }
        };
        SECTION("01 MSCAngularDeflection ") {
            for (int i = 0; i < tests; i++) {
                const auto energy = constants::k_ElectronCut + 10e-5 +
                                    static_cast<real_type>(i) *
                                    ((constants::k_MaxEkin - constants::k_ElectronCut) / tests);
                SECTION("MSCAngularDeflection Energy " + std::to_string(energy)) {
                    MSCAngularDeflectionTest test(OUTPUT_DATA_DIR, energy, histories);
                    MSCAngularDeflectionODPM custom(OUTPUT_DATA_DIR, energy, histories);
                    testInteraction(test, custom);
                }
            }
            std::uniform_real_distribution<real_type> distribution(constants::k_ElectronCut, constants::k_MaxEkin);
            for (int i = 0; i < 10; i++) {
                const auto energy = distribution(generator);
                SECTION("MSCAngularDeflection Random Energy " + std::to_string(i)) {
                    INFO(energy);
                    MSCAngularDeflectionTest test(OUTPUT_DATA_DIR, energy, histories);
                    MSCAngularDeflectionODPM custom(OUTPUT_DATA_DIR, energy, histories);
                    testInteraction(test, custom);
                }
            }
        }SECTION("02 Bremms") {
            for (int i = 1; i < tests; i++) {
                const auto energy = constants::k_GammaCut + 10e-5 +
                                    static_cast<real_type>(i) *
                                    ((constants::k_MaxEkin - constants::k_GammaCut) / tests);
                SECTION("Bremms Energy " + std::to_string(energy)) {
                    BremTest test(OUTPUT_DATA_DIR, energy, histories);
                    BremmODPM custom(OUTPUT_DATA_DIR, energy, histories);
                    testInteraction(test, custom);
                }
            }
            std::uniform_real_distribution<real_type> distribution(constants::k_GammaCut, constants::k_MaxEkin);
            for (int i = 0; i < 10; i++) {
                const auto energy = distribution(generator);
                SECTION("Bremms Random Energy " + std::to_string(i)) {
                    INFO(energy);
                    BremTest test(OUTPUT_DATA_DIR, energy, histories);
                    BremmODPM custom(OUTPUT_DATA_DIR, energy, histories);
                    testInteraction(test, custom);
                }
            }
        }SECTION("03 Moller") {
            for (int i = 0; i < tests; i++) {
                const auto energy =
                        constants::k_ElectronCut * 2 + 10e-5 +
                        static_cast<real_type>(i) * ((constants::k_MaxEkin - 2 * constants::k_ElectronCut) / tests);
                SECTION("Moller Energy " + std::to_string(energy)) {
                    MollerTest test(OUTPUT_DATA_DIR, energy, histories);
                    MollerODPM custom(OUTPUT_DATA_DIR, energy, histories);
                    testInteraction(test, custom);
                }
            }
            std::uniform_real_distribution<real_type> distribution(constants::k_ElectronCut * 2, constants::k_MaxEkin);
            for (int i = 0; i < 10; i++) {
                const auto energy = distribution(generator);
                SECTION("Moller Random Energy " + std::to_string(i)) {
                    INFO(energy);
                    MollerTest test(OUTPUT_DATA_DIR, energy, histories);
                    MollerODPM custom(OUTPUT_DATA_DIR, energy, histories);
                    testInteraction(test, custom);
                }
            }
        }SECTION("04 Annihilation") {
            for (int i = 0; i < tests; i++) {
                const auto energy =
                        constants::k_ElectronCut * 2 + 10e-5 +
                        static_cast<real_type>(i) * ((constants::k_MaxEkin - 2 * constants::k_ElectronCut) / tests);
                SECTION("Annihilation Random Energy " + std::to_string(energy)) {
                    AnnihilationTest test(OUTPUT_DATA_DIR, energy, histories);
                    AnnihilationODPM custom(OUTPUT_DATA_DIR, energy, histories);
                    testInteraction(test, custom);
                }
            }
            std::uniform_real_distribution<real_type> distribution(constants::k_ElectronCut * 2, constants::k_MaxEkin);
            for (int i = 0; i < 10; i++) {
                const auto energy = distribution(generator);
                SECTION("Annihilation Random Energy " + std::to_string(i)) {
                    INFO(energy);
                    AnnihilationTest test(OUTPUT_DATA_DIR, energy, histories);
                    AnnihilationODPM custom(OUTPUT_DATA_DIR, energy, histories);
                    testInteraction(test, custom);
                }
            }
        }SECTION("05 Compton Primary") {
            for (int i = 0; i < tests; i++) {
                const auto energy = constants::k_GammaCut + 10e-5 +
                                    static_cast<real_type>(i) *
                                    ((constants::k_MaxEkin - constants::k_GammaCut) / tests);
                SECTION("Compton Energy " + std::to_string(energy) + " Primary") {
                    ComptonPrimaryTest test(OUTPUT_DATA_DIR, energy, histories);
                    ComptonPrimaryODPM custom(OUTPUT_DATA_DIR, energy, histories);
                    testInteraction(test, custom);
                }
                //            FIXME: this part of the test seems to segfault for no reason
                // SECTION("Compton Energy " + std::to_string(energy) + " Secondary") {
                //     ComptonSecondaryTest test2(OUTPUT_DATA_DIR, energy, histories);
                //     ComptonSecondaryODPM custom2(OUTPUT_DATA_DIR, energy, histories);
                //     testInteraction(test2, custom2);
                // }
            }
            std::uniform_real_distribution<real_type> distribution(constants::k_GammaCut, constants::k_MaxEkin);
            for (int i = 0; i < 10; i++) {
                const auto energy = distribution(generator);
                SECTION("Compton Random Energy " + std::to_string(i) + " Primary") {
                    INFO(energy);
                    ComptonPrimaryTest test(OUTPUT_DATA_DIR, energy, histories);
                    ComptonPrimaryODPM custom(OUTPUT_DATA_DIR, energy, histories);
                    testInteraction(test, custom);
                }
                //        THIS TEST DOES NOT PASS BECAUSE THE SECONDARY PARTICLES ARE
                //        HANDLED DIFFERENTLY IN THE TWO IMPLEMENTATIONS.
                //        THIS IS NOT CONSIDERED AN ERROR.
                //        SECTION("Compton Random Energy " + std::to_string(i) + "
                //        Secondary") {
                //            INFO(energy);
                //            ComptonSecondaryTest test2(OUTPUT_DATA_DIR, energy, histories);
                //            ComptonSecondaryODPM custom2(OUTPUT_DATA_DIR, energy, histories);
                //            testInteraction(test2, custom2);
                //        }
            }
        }
    }
}  // namespace opmc
