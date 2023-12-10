//
// Created by mbarbone on 7/26/21.
//

#ifndef ODPM_CPU_INCLUDE_RUN_H_
#define ODPM_CPU_INCLUDE_RUN_H_

#ifdef RPMALLOC
#include <rpmalloc.h>
#endif

#include <thread>
#include <variant>

#include "beams.h"
#include "geometry.h"
#include "materials.h"
#include "particle_queue.h"
#include "physics.h"
#include "random.h"
#include "utils.h"
#include "initialise_tables.h"

namespace opmc {

    template<typename V>
    class Run {
    private:
        template<class Z>
        class Simulation {
        public:
            ODPM_INLINE constexpr explicit Simulation(Physics<V, Z> &physics) noexcept:
                    physics(physics) {}

            ODPM_INLINE constexpr void operator()(Electron &electron) noexcept { physics.simulate(electron); }

            ODPM_INLINE constexpr void operator()(Photon &photon) noexcept { physics.simulate(photon); };

            ODPM_INLINE constexpr void operator()(Positron &positron) noexcept { physics.simulate(positron); }

        private:
            Physics<V, Z> &physics;
        };

    public:
        Run(const unsigned long seed, const unsigned threads = std::thread::hardware_concurrency())
                : k_threads{threads} {
            for (unsigned long i = 0; i < k_threads; ++i) {
                m_rngs.emplace_back(seed, i);
            }
        }

        template<typename T, typename Z = ParticleQueue>
        void simulate(V &geometry, const T &beam, const Material &referenceMaterial,
                      const ODPMElectronData &electron_data, const ODPMPhotonData &photon_data,
                      std::uint64_t histories) {
            std::cout << "Monte Carlo Start" << std::endl;
            // time the function with chrono
            auto start = std::chrono::steady_clock::now();
            // each thread executes histories/threads
            histories = (histories + k_threads - 1) / k_threads;
            // But it must be at least 1
            m_chunk_size = math::max(math::min(histories, m_chunk_size), 1);
            // chunks to be executed are
            const auto chunks = (histories + m_chunk_size - 1) / m_chunk_size;

#pragma omp parallel for default(none) firstprivate(beam) schedule(static)\
shared(chunks, geometry, referenceMaterial, photon_data, electron_data, std::cout)
            for (unsigned i = 0; i < k_threads; ++i) {
                auto &rng = m_rngs[i];
                Z queue{};
                auto physics = Physics<V, Z>{rng, geometry, queue, referenceMaterial, photon_data, electron_data};
                auto simulation = Simulation<Z>{physics};
                for (unsigned _ = 0; _ < chunks; ++_) {
                    for (unsigned history = 0; history < m_chunk_size; ++history) {
                        queue.push(beam.generateParticle(rng));
                    }
                    while (!queue.empty()) {
                        auto particle = queue.pop();
                        std::visit(simulation, particle);
                    }
                }
            }
            geometry.toGray();
            auto end = std::chrono::steady_clock::now();
            std::chrono::duration<double, std::milli> timeRequired = (end - start);
            std::cout << "Monte Carlo Simulation required " << timeRequired.count() << " milliseconds" << std::endl;
        }

    private:
        const unsigned k_threads;
        std::vector<Random> m_rngs;
        unsigned long m_chunk_size = 16384;

    public:
        Run(Run &&run) noexcept = delete;

        Run(const Run &run) = delete;

        Run &operator=(const Run &run) = delete;

        Run &operator=(Run &&run) noexcept = delete;

        virtual ~Run() = default;
    };  // namespace opmc

}  // namespace opmc
#endif  // ODPM_CPU_INCLUDE_RUN_H_
