//
// Created by mbarbone on 10/26/21.
//

#ifndef ODPM_CPU_INCLUDE_PARTICLE_QUEUE_H_
#define ODPM_CPU_INCLUDE_PARTICLE_QUEUE_H_

#include <absl/container/inlined_vector.h>

#include <deque>
#include <memory_resource>
#include <queue>
#include <utility>
#include <variant>

#include "particles.h"
#include "utils.h"

namespace opmc {

    using SimulationParticles = std::variant<Electron, Photon, Positron>;

    class ParticleQueue {
    public:
        constexpr ODPM_INLINE ParticleQueue() noexcept = default;

        ODPM_INLINE ParticleQueue(const ParticleQueue &) noexcept = default;

        ODPM_INLINE ParticleQueue(ParticleQueue &&) noexcept = default;

        ParticleQueue &operator=(const ParticleQueue &) noexcept = default;

        ParticleQueue &operator=(ParticleQueue &&) noexcept = default;

        template<typename T>
        ODPM_INLINE
        constexpr void push(T &&particle) noexcept {
            m_queue.emplace_back(particle);
        }

        ODPM_INLINE
        SimulationParticles pop() noexcept {
            ODPM_ASSUME(!m_queue.empty());
            auto particle = m_queue.back();
            m_queue.pop_back();
            return particle;
        }

        ODPM_INLINE ODPM_CONST
        bool empty() const noexcept { return m_queue.empty(); }

    private:
        absl::InlinedVector<SimulationParticles, 100000> m_queue{};
    };

    class PriorityParticleQueue {
    public:
        ODPM_INLINE void push(SimulationParticles &&particle) noexcept { m_queue.emplace(particle); }

        ODPM_INLINE SimulationParticles pop() noexcept {
            auto particle = m_queue.top();
            m_queue.pop();
            return particle;
        }

        ODPM_INLINE bool empty() const noexcept { return m_queue.empty(); }

    private:
        struct greater {
            template<typename T, typename V>
            constexpr ODPM_INLINE bool operator()(const T &lhs, const V &rhs) const noexcept {
                return static_cast<const Particle &>(lhs) < static_cast<const Particle &>(rhs);
            }
        };

        constexpr static greater m_greater{};

        struct GreaterParticle {
            constexpr ODPM_INLINE bool operator()(const SimulationParticles &lhs,
                                                  const SimulationParticles &rhs) const noexcept {
                return std::visit(m_greater, lhs, rhs);
            }
        };

        std::priority_queue<SimulationParticles, absl::InlinedVector<SimulationParticles, 1000000>, GreaterParticle>
                m_queue{};
    };

#ifdef GPU
    class GPUParticleQueue {
       public:

        constexpr ODPM_INLINE GPUParticleQueue() : used(0) {}

        virtual ~GPUParticleQueue() {}


        ODPM_INLINE constexpr void push(SimulationParticles&& particle) {
            assert(used < capacity);
            data[used++] = particle;
        }


        ODPM_INLINE constexpr SimulationParticles pop() {
            assert(used > 0);
            return data[--used];
        }


        ODPM_INLINE constexpr bool empty() const { return used <= 0; }

       private:
        unsigned char used;
        static const unsigned char capacity = 64;
        SimulationParticles data[capacity];
    };
#endif

#ifdef TEST

    class DummyQueue {
    public:

        ODPM_INLINE constexpr DummyQueue() = default;


        virtual ~DummyQueue() = default;


        ODPM_INLINE constexpr void push(__attribute__((unused)) SimulationParticles &&particle) {}


        ODPM_INLINE constexpr bool empty() const { return false; }
    };

#endif
}  // namespace opmc
#endif  // ODPM_CPU_INCLUDE_PARTICLE_QUEUE_H_
