//
// Created by mbarbone on 7/28/21.
//

#ifndef ODPM_CPU_INCLUDE_PARTICLES_H_
#define ODPM_CPU_INCLUDE_PARTICLES_H_

#include <ostream>

#include "types.h"
//

namespace opmc {

class Particle {
   public:
    ThreeVector<real_type> position;
    ThreeVector<real_type> direction;
    real_type energy;

    ODPM_INLINE constexpr Particle(ThreeVector<real_type> position, ThreeVector<real_type> direction,
                                   real_type kinetic_energy)
        : position(std::move(position)), direction(std::move(direction)), energy(std::move(kinetic_energy)) {}

    template <typename S1, typename S2>
    ODPM_INLINE constexpr Particle(S1 &&position, S2 &&direction, real_type kinetic_energy)
        : position(std::forward<ThreeVector<real_type>>(position)),
          direction(std::forward<ThreeVector<real_type>>(direction)),
          energy(kinetic_energy) {}

    ODPM_INLINE constexpr Particle() : position(), direction(), energy(0){};

    ODPM_INLINE constexpr Particle(Particle &&) noexcept            = default;

    ODPM_INLINE constexpr Particle(const Particle &)                = default;

    ODPM_INLINE constexpr Particle &operator=(Particle &&) noexcept = default;

    ODPM_INLINE constexpr Particle &operator=(const Particle &)     = default;

    template <typename T, typename V>
    ODPM_INLINE ODPM_PURE constexpr friend bool operator<(const T &lhs, const V &rhs) {
        return lhs.position < rhs.position && lhs.energy < rhs.energy;
    }

    template <typename T, typename V>
    ODPM_INLINE ODPM_PURE constexpr friend bool operator<=(const T &lhs, const V &rhs) {
        return lhs.position <= rhs.position && lhs.energy <= rhs.energy;
    }

    template <typename T, typename V>
    ODPM_INLINE ODPM_PURE constexpr friend bool operator>(const T &lhs, const V &rhs) {
        return lhs.position > rhs.position && lhs.energy > rhs.energy;
    }

    template <typename T, typename V>
    ODPM_INLINE ODPM_PURE constexpr friend bool operator>=(const T &lhs, const V &rhs) {
        return lhs.position >= rhs.position && lhs.energy >= rhs.energy;
    }

    friend std::ostream &operator<<(std::ostream &os, const Particle &particle) {
        os << "position: " << particle.position << " direction: " << particle.direction
           << " energy: " << particle.energy;
        return os;
    }
};

class Electron : public Particle {
   public:
    using Particle::Particle;
};

class Positron : public Electron {
   public:
    using Electron::Electron;
};

class Photon : public Particle {
   public:
    using Particle::Particle;
};

// NOTE: this needs to be trivially allocatable.
struct ParticleState {
    real_type num_tr_1_mfp, num_moller_mfp, inv_moller_mfp, num_brem_mfp, num_tr_1_mfp_0, initial_energy;
    bool is_msc_hinge;

    ODPM_INLINE constexpr ParticleState() noexcept                       = default;

    ODPM_INLINE constexpr ParticleState(ParticleState &&) noexcept       = default;

    ODPM_INLINE constexpr ParticleState(const ParticleState &) noexcept  = default;

    ODPM_INLINE ParticleState &operator=(ParticleState &&) noexcept      = default;

    ODPM_INLINE ParticleState &operator=(const ParticleState &) noexcept = default;

    friend std::ostream &operator<<(std::ostream &os, const ParticleState &state) {
        os << "num_tr_1_mfp: " << state.num_tr_1_mfp << " num_moller_mfp: " << state.num_moller_mfp
           << " inv_moller_mfp: " << state.inv_moller_mfp << " num_brem_mfp: " << state.num_brem_mfp
           << " num_tr_1_mfp_0: " << state.num_tr_1_mfp_0 << " initial_energy: " << state.initial_energy
           << " is_msc_hinge: " << state.is_msc_hinge;
        return os;
    }
};

class StatefulElectron : public Electron, public ParticleState {
   public:
    StatefulElectron(const Electron &electron, const ParticleState &state) : Electron(electron), ParticleState(state) {}
};

class StatefulPositron : public Positron, public ParticleState {
   public:
    StatefulPositron(const Positron &positron, const ParticleState &state) : Positron(positron), ParticleState(state) {}
};

}  // namespace opmc
#endif  // ODPM_CPU_INCLUDE_PARTICLES_H_
