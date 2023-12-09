#define CATCH_CONFIG_MAIN
#include <cuda.h>

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
#include <catch2/catch.hpp>
#include <cstdio>
#include <cuda/api.hpp>
#include <functional>
#include <random>

#include "gpu_unique_ptr.h"
#include "gpu_utils.h"
#include "test_utils.h"

namespace opmc {

static constexpr auto verbose_print = false;
static constexpr auto tests         = 16384;
static constexpr auto block_size    = 256;
static constexpr auto grid_size     = div_rounding_up(tests, block_size);

template <class F>
ODPM_KERNEL static void kernel(F f) {
    const auto i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i < tests) {
        f(i);
    }
}

template <bool enable, typename... Ts>
static constexpr void debugPrintf(Ts... args) {
    if constexpr (enable) {
        printf(args...);
    }
}

#define INVOKE_TEST(arg)                    \
    kernel<<<grid_size, block_size>>>(arg); \
    gpuErrorAssert(cudaDeviceSynchronize())

TEST_CASE("GPU INTERPOLATION TESTS") {
    const auto seed = std::random_device()();
    INFO("SEED " << seed);
    std::default_random_engine rng{seed};
    SimElectronData elData{};
    SimPhotonData phData{};
    SimMaterialData mData{};
    elData.Load(DATA_DIR, 0);
    phData.Load(DATA_DIR, 0);
    mData.Load(DATA_DIR, 0);
    TestMaterial material{DATA_DIR};
    TestPhotonData photonData{DATA_DIR};
    TestElectronData electronData{DATA_DIR};
    // This copies the data needed by the material on the GPU
    auto device = cuda::device::current::get();
    material.initializeGPU(device);
    photonData.initializeGPU(device);
    electronData.initializeGPU(device);
    // this allocates the material on the GPU and copies it
    auto gpu_material      = cuda::make_unique(device, material);
    auto gpu_electron_data = cuda::make_unique(device, electronData);
    auto gpu_photon_data   = cuda::make_unique(device, photonData);

    SECTION("ITr1MFPElastic") {
        std::uniform_real_distribution<real_type> distribution{material.GetITr1MFPElastic().getMinInput(),
                                                               material.GetITr1MFPElastic().getMaxInput()};

        std::vector<real_type> inputs(tests);
        std::generate(inputs.begin(), inputs.end(), [&distribution, &rng]() { return distribution(rng); });
        std::vector<real_type> reference(tests);
        std::transform(inputs.begin(), inputs.end(), reference.begin(), [&elData, &material](const real_type& x) {
            return elData.GetITr1MFPElastic()->GetITr1MFPPerDensity(x, material.id());
        });
        cuda::gpu_array<real_type> gpu_inputs(cuda::device::current::get(), inputs);
        cuda::gpu_array<real_type> gpu_outputs(cuda::device::current::get(), tests);
        auto gpu_function = [in = gpu_inputs.m_data(), out = gpu_outputs.m_data(),
                             material_ptr = gpu_material.get()] __device__(int i) {
            out[i] = material_ptr->iTr1MFPElasticElectron(in[i]);
            debugPrintf<verbose_print>("i = %d, in[i] = %f, out[i] = %f\n", i, in[i], out[i]);
        };

        INVOKE_TEST(gpu_function);

        auto gpu_results = gpu_outputs.to_vector();
        for (int i = 0; i < tests; ++i) {
            INFO("Test #: " << i)
            INFO("input: " << inputs[i])
            REQUIRE(reference[i] == Approx(gpu_results[i]));
        }
    }
    SECTION("MaxScatStrength") {
        std::uniform_real_distribution<real_type> distribution{material.GetMaxScatStrength().getMinInput(),
                                                               material.GetMaxScatStrength().getMaxInput()};

        std::vector<real_type> inputs(tests);
        std::generate(inputs.begin(), inputs.end(), [&distribution, &rng]() { return distribution(rng); });
        std::vector<real_type> reference(tests);
        std::transform(inputs.begin(), inputs.end(), reference.begin(), [&elData, &material](const real_type& x) {
            return elData.GetMaxScatStrength()->GetMaxScatStrength(x);
        });
        cuda::gpu_array<real_type> gpu_inputs(cuda::device::current::get(), inputs);
        cuda::gpu_array<real_type> gpu_outputs(cuda::device::current::get(), tests);
        auto gpu_function = [in = gpu_inputs.m_data(), out = gpu_outputs.m_data(),
                             material_ptr = gpu_material.get()] __device__(int i) {
            out[i] = material_ptr->maxScatteringStrengthElectron(in[i]);
            debugPrintf<verbose_print>("i = %d, in[i] = %f, out[i] = %f\n", i, in[i], out[i]);
        };

        INVOKE_TEST(gpu_function);

        auto gpu_results = gpu_outputs.to_vector();
        for (int i = 0; i < tests; ++i) {
            INFO("Test #: " << i)
            INFO("input: " << inputs[i])
            REQUIRE(reference[i] == Approx(gpu_results[i]));
        }
    }
    SECTION("IMFPMoller") {
        std::uniform_real_distribution<real_type> distribution{material.GetIMFPMoller().getMinInput(),
                                                               material.GetIMFPMoller().getMaxInput()};
        std::vector<real_type> inputs(tests);
        std::generate(inputs.begin(), inputs.end(), [&distribution, &rng]() { return distribution(rng); });
        std::vector<real_type> reference(tests);
        std::transform(inputs.begin(), inputs.end(), reference.begin(), [&elData, &material](const real_type& x) {
            return elData.GetIMFPMoller()->GetIMFPPerDensity(x);
        });
        cuda::gpu_array<real_type> gpu_inputs(cuda::device::current::get(), inputs);
        cuda::gpu_array<real_type> gpu_outputs(cuda::device::current::get(), tests);
        auto gpu_function = [in = gpu_inputs.m_data(), out = gpu_outputs.m_data(),
                             material_ptr = gpu_material.get()] __device__(int i) {
            out[i] = material_ptr->iMFPMollerElectron(in[i]);
            debugPrintf<verbose_print>("i = %d, in[i] = %f, out[i] = %f\n", i, in[i], out[i]);
        };

        INVOKE_TEST(gpu_function);

        auto gpu_results = gpu_outputs.to_vector();
        for (int i = 0; i < tests; ++i) {
            INFO("Test #: " << i)
            INFO("input: " << inputs[i])
            REQUIRE(reference[i] == Approx(gpu_results[i]));
        }
    }
    SECTION("IMFPBrem") {
        std::uniform_real_distribution<real_type> distribution{material.GetIMFPBrem().getMinInput(),
                                                               material.GetIMFPBrem().getMaxInput()};

        std::vector<real_type> inputs(tests);
        std::generate(inputs.begin(), inputs.end(), [&distribution, &rng]() { return distribution(rng); });
        std::vector<real_type> reference(tests);
        std::transform(inputs.begin(), inputs.end(), reference.begin(), [&elData, &material](const real_type& x) {
            return elData.GetIMFPBrem()->GetIMFPPerDensity(x, material.id());
        });
        cuda::gpu_array<real_type> gpu_inputs(cuda::device::current::get(), inputs);
        cuda::gpu_array<real_type> gpu_outputs(cuda::device::current::get(), tests);
        auto gpu_function = [in = gpu_inputs.m_data(), out = gpu_outputs.m_data(),
                             material_ptr = gpu_material.get()] __device__(int i) {
            out[i] = material_ptr->iMFPBremElectron(in[i]);
            debugPrintf<verbose_print>("i = %d, in[i] = %f, out[i] = %f\n", i, in[i], out[i]);
        };

        INVOKE_TEST(gpu_function);

        auto gpu_results = gpu_outputs.to_vector();
        for (int i = 0; i < tests; ++i) {
            INFO("Test #: " << i)
            INFO("input: " << inputs[i])
            REQUIRE(reference[i] == Approx(gpu_results[i]));
        }
    }
    SECTION("StoppingPower") {
        std::uniform_real_distribution<real_type> distribution{material.GetStoppingPower().getMinInput(),
                                                               material.GetStoppingPower().getMaxInput()};

        std::vector<real_type> inputs(tests);
        std::generate(inputs.begin(), inputs.end(), [&distribution, &rng]() { return distribution(rng); });
        std::vector<real_type> reference(tests);
        std::transform(inputs.begin(), inputs.end(), reference.begin(), [&elData, &material](const real_type& x) {
            return elData.GetDEDX()->GetDEDXPerDensity(x, material.id());
        });
        cuda::gpu_array<real_type> gpu_inputs(cuda::device::current::get(), inputs);
        cuda::gpu_array<real_type> gpu_outputs(cuda::device::current::get(), tests);
        auto gpu_function = [in = gpu_inputs.m_data(), out = gpu_outputs.m_data(),
                             material_ptr = gpu_material.get()] __device__(int i) {
            out[i] = material_ptr->stoppingPowerElectron(in[i]);
            debugPrintf<verbose_print>("i = %d, in[i] = %f, out[i] = %f\n", i, in[i], out[i]);
        };

        INVOKE_TEST(gpu_function);

        auto gpu_results = gpu_outputs.to_vector();
        for (int i = 0; i < tests; ++i) {
            INFO("Test #: " << i)
            INFO("input: " << inputs[i])
            REQUIRE(reference[i] == Approx(gpu_results[i]));
        }
    }
    SECTION("SeltzerBerger") {
        std::uniform_real_distribution<real_type> energy_distribution{constants::k_ElectronCut, constants::k_MaxEkin};
        std::uniform_real_distribution<real_type> uniform_distribution{};

        std::vector<real_type> inputs(tests), rng1(tests), rng2(tests), rng3(tests);
        std::generate(inputs.begin(), inputs.end(),
                      [&energy_distribution, &rng]() { return energy_distribution(rng); });
        std::generate(rng1.begin(), rng1.end(), [&uniform_distribution, &rng]() { return uniform_distribution(rng); });
        std::generate(rng2.begin(), rng2.end(), [&uniform_distribution, &rng]() { return uniform_distribution(rng); });
        std::generate(rng3.begin(), rng3.end(), [&uniform_distribution, &rng]() { return uniform_distribution(rng); });

        std::vector<real_type> reference(tests);
        for (int i = 0; i < tests; ++i) {
            reference[i] =
                elData.GetTheSBTables()->SampleEnergyTransfer(inputs[i], material.id(), rng1[i], rng2[i], rng3[i]);
        }
        cuda::gpu_array<real_type> gpu_inputs(cuda::device::current::get(), inputs);
        cuda::gpu_array<real_type> gpu_rng1(cuda::device::current::get(), rng1);
        cuda::gpu_array<real_type> gpu_rng2(cuda::device::current::get(), rng2);
        cuda::gpu_array<real_type> gpu_rng3(cuda::device::current::get(), rng3);
        cuda::gpu_array<real_type> gpu_outputs(cuda::device::current::get(), tests);

        auto gpu_function = [in = gpu_inputs.m_data(), out = gpu_outputs.m_data(), rng1 = gpu_rng1.m_data(),
                             rng2 = gpu_rng2.m_data(), rng3 = gpu_rng3.m_data(),
                             material_ptr = gpu_material.get()] __device__(int i) {
            out[i] = material_ptr->seltzerBergerElectron(in[i], rng1[i], rng2[i], rng3[i]);
            debugPrintf<verbose_print>("i = %d, in[i] = %f, out[i] = %f, rng1[i] = %f, rng2[i] = %f, rng3[i] = %f\n", i,
                                       in[i], out[i], rng1[i], rng2[i], rng3[i]);
        };

        INVOKE_TEST(gpu_function);

        auto gpu_results = gpu_outputs.to_vector();
        for (int i = 0; i < tests; ++i) {
            INFO("Test #: " << i)
            INFO("input: " << inputs[i])
            INFO("rng1 " << rng1[i])
            INFO("rng2 " << rng2[i])
            INFO("rng3 " << rng3[i])
            REQUIRE(reference[i] == Approx(gpu_results[i]));
        }
    }
    SECTION("MollerEnergyTransfer") {
        std::uniform_real_distribution<real_type> energy_distribution{constants::k_ElectronCut * 2,
                                                                      constants::k_MaxEkin};
        std::uniform_real_distribution<real_type> uniform_distribution{};

        std::vector<real_type> inputs(tests), rng1(tests), rng2(tests), rng3(tests);
        std::generate(inputs.begin(), inputs.end(),
                      [&energy_distribution, &rng]() { return energy_distribution(rng); });
        std::generate(rng1.begin(), rng1.end(), [&uniform_distribution, &rng]() { return uniform_distribution(rng); });
        std::generate(rng2.begin(), rng2.end(), [&uniform_distribution, &rng]() { return uniform_distribution(rng); });
        std::generate(rng3.begin(), rng3.end(), [&uniform_distribution, &rng]() { return uniform_distribution(rng); });

        std::vector<real_type> reference(tests);
        for (int i = 0; i < tests; ++i) {
            reference[i] = elData.GetTheMollerTables()->SampleEnergyTransfer(inputs[i], rng1[i], rng2[i], rng3[i]);
        }
        cuda::gpu_array<real_type> gpu_inputs(cuda::device::current::get(), inputs);
        cuda::gpu_array<real_type> gpu_rng1(cuda::device::current::get(), rng1);
        cuda::gpu_array<real_type> gpu_rng2(cuda::device::current::get(), rng2);
        cuda::gpu_array<real_type> gpu_rng3(cuda::device::current::get(), rng3);
        cuda::gpu_array<real_type> gpu_outputs(cuda::device::current::get(), tests);

        auto gpu_function = [in = gpu_inputs.m_data(), out = gpu_outputs.m_data(), rng1 = gpu_rng1.m_data(),
                             rng2 = gpu_rng2.m_data(), rng3 = gpu_rng3.m_data(),
                             electron_data_ptr = gpu_electron_data.get()] __device__(int i) {
            out[i] = electron_data_ptr->mollerEnergyTransfer(in[i], rng1[i], rng2[i], rng3[i]);
            debugPrintf<verbose_print>("i = %d, in[i] = %f, out[i] = %f, rng1[i] = %f, rng2[i] = %f, rng3[i] = %f\n", i,
                                       in[i], out[i], rng1[i], rng2[i], rng3[i]);
        };

        INVOKE_TEST(gpu_function);

        auto gpu_results = gpu_outputs.to_vector();
        for (int i = 0; i < tests; ++i) {
            INFO("Test #: " << i)
            INFO("input: " << inputs[i])
            INFO("rng1 " << rng1[i])
            INFO("rng2 " << rng2[i])
            INFO("rng3 " << rng3[i])
            REQUIRE(reference[i] == Approx(gpu_results[i]));
        }
    }
    SECTION("GoudsmitSaunderson") {
        std::uniform_real_distribution<real_type> energy_distribution{constants::k_ElectronCut, constants::k_MaxEkin};
        std::uniform_real_distribution<real_type> uniform_distribution{};

        std::vector<real_type> inputs(tests), rng1(tests), rng2(tests);
        std::generate(inputs.begin(), inputs.end(),
                      [&energy_distribution, &rng]() { return energy_distribution(rng); });
        std::generate(rng1.begin(), rng1.end(), [&uniform_distribution, &rng]() { return uniform_distribution(rng); });
        std::generate(rng2.begin(), rng2.end(), [&uniform_distribution, &rng]() { return uniform_distribution(rng); });

        std::vector<real_type> reference(tests);
        for (int i = 0; i < tests; ++i) {
            reference[i] = elData.GetTheGSTables()->SampleAngularDeflection(inputs[i], rng1[i], rng2[i]);
        }
        cuda::gpu_array<real_type> gpu_inputs(cuda::device::current::get(), inputs);
        cuda::gpu_array<real_type> gpu_rng1(cuda::device::current::get(), rng1);
        cuda::gpu_array<real_type> gpu_rng2(cuda::device::current::get(), rng2);
        cuda::gpu_array<real_type> gpu_outputs(cuda::device::current::get(), tests);
        auto gpu_function = [in = gpu_inputs.m_data(), out = gpu_outputs.m_data(), rng1 = gpu_rng1.m_data(),
                             rng2 = gpu_rng2.m_data(), material_ptr = gpu_material.get()] __device__(int i) {
            out[i] = material_ptr->angularDeflectionElectron(in[i], rng1[i], rng2[i]);
            debugPrintf<verbose_print>("i = %d, in[i] = %f, out[i] = %f, rng1[i] = %f, rng2[i] = %f\n", i, in[i],
                                       out[i], rng1[i], rng2[i]);
        };

        INVOKE_TEST(gpu_function);

        auto gpu_results = gpu_outputs.to_vector();
        for (int i = 0; i < tests; ++i) {
            INFO("Test #: " << i)
            INFO("input: " << inputs[i])
            INFO("rng1 " << rng1[i])
            INFO("rng2 " << rng2[i])
            REQUIRE(reference[i] == Approx(gpu_results[i]));
        }
    }

    SECTION("IMFPTotalPhoton") {
        std::uniform_real_distribution<real_type> distribution{material.GetIMFPTotalPhoton().GetMinInput(),
                                                               material.GetIMFPTotalPhoton().GetMaxInput()};

        std::vector<real_type> inputs(tests);
        std::generate(inputs.begin(), inputs.end(), [&distribution, &rng]() { return distribution(rng); });
        std::vector<real_type> reference(tests);
        std::transform(inputs.begin(), inputs.end(), reference.begin(), [&phData, &material](const real_type& x) {
            return phData.GetIMFPTotal()->GetIMFPPerDensity(x, material.id());
        });
        cuda::gpu_array<real_type> gpu_inputs(cuda::device::current::get(), inputs);
        cuda::gpu_array<real_type> gpu_outputs(cuda::device::current::get(), tests);
        auto gpu_function = [in = gpu_inputs.m_data(), out = gpu_outputs.m_data(),
                             material_ptr = gpu_material.get()] __device__(int i) {
            out[i] = material_ptr->iMFPTotalPhoton(in[i]);
            debugPrintf<verbose_print>("i = %d, in[i] = %f, out[i] = %f\n", i, in[i], out[i]);
        };

        INVOKE_TEST(gpu_function);

        auto gpu_results = gpu_outputs.to_vector();
        for (int i = 0; i < tests; ++i) {
            INFO("Test #: " << i)
            INFO("input: " << inputs[i])
            REQUIRE(reference[i] == Approx(gpu_results[i]));
        }
    }
    SECTION("IMFPComptonPhoton") {
        std::uniform_real_distribution<real_type> distribution{material.GetIMFPComptonPhoton().GetMinInput(),
                                                               material.GetIMFPComptonPhoton().GetMaxInput()};

        std::vector<real_type> inputs(tests);
        std::generate(inputs.begin(), inputs.end(), [&distribution, &rng]() { return distribution(rng); });
        std::vector<real_type> reference(tests);
        std::transform(inputs.begin(), inputs.end(), reference.begin(), [&phData, &material](const real_type& x) {
            return phData.GetIMFPCompton()->GetIMFPPerDensity(x, material.id());
        });
        cuda::gpu_array<real_type> gpu_inputs(cuda::device::current::get(), inputs);
        cuda::gpu_array<real_type> gpu_outputs(cuda::device::current::get(), tests);
        auto gpu_function = [in = gpu_inputs.m_data(), out = gpu_outputs.m_data(),
                             material_ptr = gpu_material.get()] __device__(int i) {
            out[i] = material_ptr->iMFPComptonPhoton(in[i]);
            debugPrintf<verbose_print>("i = %d, in[i] = %f, out[i] = %f\n", i, in[i], out[i]);
        };

        INVOKE_TEST(gpu_function);

        auto gpu_results = gpu_outputs.to_vector();
        for (int i = 0; i < tests; ++i) {
            INFO("Test #: " << i)
            INFO("input: " << inputs[i])
            REQUIRE(reference[i] == Approx(gpu_results[i]));
        }
    }
    SECTION("IMFPPairProdPhoton") {
        std::uniform_real_distribution<real_type> distribution{material.GetIMFPComptonPhoton().GetMinInput(),
                                                               material.GetIMFPComptonPhoton().GetMaxInput()};

        std::vector<real_type> inputs(tests);
        std::generate(inputs.begin(), inputs.end(), [&distribution, &rng]() { return distribution(rng); });
        std::vector<real_type> reference(tests);
        std::transform(inputs.begin(), inputs.end(), reference.begin(), [&phData, &material](const real_type& x) {
            return phData.GetIMFPCompton()->GetIMFPPerDensity(x, material.id());
        });
        cuda::gpu_array<real_type> gpu_inputs(cuda::device::current::get(), inputs);
        cuda::gpu_array<real_type> gpu_outputs(cuda::device::current::get(), tests);
        auto gpu_function = [in = gpu_inputs.m_data(), out = gpu_outputs.m_data(),
                             material_ptr = gpu_material.get()] __device__(int i) {
            out[i] = material_ptr->iMFPComptonPhoton(in[i]);
            debugPrintf<verbose_print>("i = %d, in[i] = %f, out[i] = %f\n", i, in[i], out[i]);
        };

        INVOKE_TEST(gpu_function);

        auto gpu_results = gpu_outputs.to_vector();
        for (int i = 0; i < tests; ++i) {
            INFO("Test #: " << i)
            INFO("input: " << inputs[i])
            REQUIRE(reference[i] == Approx(gpu_results[i]));
        }
    }
    SECTION("IMFPMaxPhoton") {
        std::uniform_real_distribution<real_type> distribution{photonData.GetIMFPMaxPhoton().GetMinInput(),
                                                               photonData.GetIMFPMaxPhoton().GetMaxInput()};

        std::vector<real_type> inputs(tests);
        std::generate(inputs.begin(), inputs.end(), [&distribution, &rng]() { return distribution(rng); });
        std::vector<real_type> reference(tests);
        std::transform(inputs.begin(), inputs.end(), reference.begin(),
                       [&phData](const real_type& x) { return phData.GetIMFPTotalMax()->GetIMFP(x); });
        cuda::gpu_array<real_type> gpu_inputs(cuda::device::current::get(), inputs);
        cuda::gpu_array<real_type> gpu_outputs(cuda::device::current::get(), tests);
        auto gpu_function = [in = gpu_inputs.m_data(), out = gpu_outputs.m_data(),
                             photon_data_ptr = gpu_photon_data.get()] __device__(int i) {
            out[i] = photon_data_ptr->iMFPMaxPhoton(in[i]);
            debugPrintf<verbose_print>("i = %d, in[i] = %f, out[i] = %f\n", i, in[i], out[i]);
        };

        INVOKE_TEST(gpu_function);

        auto gpu_results = gpu_outputs.to_vector();
        for (int i = 0; i < tests; ++i) {
            INFO("Test #: " << i)
            INFO("input: " << inputs[i])
            REQUIRE(reference[i] == Approx(gpu_results[i]));
        }
    }
    SECTION("KleinNishina") {
        std::uniform_real_distribution<real_type> energy_distribution{constants::k_GammaCut, constants::k_MaxEkin};
        std::uniform_real_distribution<real_type> uniform_distribution{};

        std::vector<real_type> inputs(tests), rng1(tests), rng2(tests), rng3(tests);
        std::generate(inputs.begin(), inputs.end(),
                      [&energy_distribution, &rng]() { return energy_distribution(rng); });
        std::generate(rng1.begin(), rng1.end(), [&uniform_distribution, &rng]() { return uniform_distribution(rng); });
        std::generate(rng2.begin(), rng2.end(), [&uniform_distribution, &rng]() { return uniform_distribution(rng); });
        std::generate(rng3.begin(), rng3.end(), [&uniform_distribution, &rng]() { return uniform_distribution(rng); });

        std::vector<real_type> reference(tests);
        for (int i = 0; i < tests; ++i) {
            reference[i] = phData.GetTheKNTables()->SampleEnergyTransfer(inputs[i], rng1[i], rng2[i], rng3[i]);
        }
        cuda::gpu_array<real_type> gpu_inputs(cuda::device::current::get(), inputs);
        cuda::gpu_array<real_type> gpu_rng1(cuda::device::current::get(), rng1);
        cuda::gpu_array<real_type> gpu_rng2(cuda::device::current::get(), rng2);
        cuda::gpu_array<real_type> gpu_rng3(cuda::device::current::get(), rng3);
        cuda::gpu_array<real_type> gpu_outputs(cuda::device::current::get(), tests);

        auto gpu_function = [in = gpu_inputs.m_data(), out = gpu_outputs.m_data(), rng1 = gpu_rng1.m_data(),
                             rng2 = gpu_rng2.m_data(), rng3 = gpu_rng3.m_data(),
                             photon_data_ptr = gpu_photon_data.get()] __device__(int i) {
            out[i] = photon_data_ptr->kleinNishina(in[i], rng1[i], rng2[i], rng3[i]);
            debugPrintf<verbose_print>("i = %d, in[i] = %f, out[i] = %f, rng1[i] = %f, rng2[i] = %f, rng3[i] = %f\n", i,
                                       in[i], out[i], rng1[i], rng2[i], rng3[i]);
        };

        INVOKE_TEST(gpu_function);

        auto gpu_results = gpu_outputs.to_vector();
        for (int i = 0; i < tests; ++i) {
            INFO("Test #: " << i)
            INFO("input: " << inputs[i])
            INFO("rng1 " << rng1[i])
            INFO("rng2 " << rng2[i])
            INFO("rng3 " << rng3[i])
            REQUIRE(reference[i] == Approx(gpu_results[i]));
        }
    }
}

}  // namespace opmc