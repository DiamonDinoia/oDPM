//
// Created by mbarbone on 10/26/21.
//

#ifndef ODPM_CPU_INCLUDE_BEAMS_H_
#define ODPM_CPU_INCLUDE_BEAMS_H_

#include <cassert>
#include <fstream>
#include <iostream>

#include "particles.h"
#include "random.h"

#ifdef GPU
#include <cuda/api.hpp>


#endif

namespace opmc {

    template<typename T>
    class PencilBeam {
    public:

        constexpr PencilBeam(PencilBeam &&) = default;

        constexpr PencilBeam(const PencilBeam &) = default;


        PencilBeam &operator=(PencilBeam &&) = default;

        PencilBeam &operator=(const PencilBeam &) = default;

        ODPM_INLINE
        constexpr PencilBeam(const ThreeVector<real_type> &position, const ThreeVector<real_type> &direction,
                             const real_type energy)
                : k_default_particle(position, direction, energy) {}

        ODPM_INLINE ODPM_CONST
        constexpr T generateParticle(__attribute__((unused)) Random &rng) const noexcept { return k_default_particle; }

        friend std::ostream &operator<<(std::ostream &os, const PencilBeam &beam) {
            os << "k_default_particle: " << beam.k_default_particle;
            return os;
        }

    private:
        const T k_default_particle;
    };

    using ElectronPencilBeam = PencilBeam<Electron>;
    using PhotonPencilBeam = PencilBeam<Photon>;

/**
 * Parameters list
 * sad = source to  axis distance in mm
 * ffl = field of view in mm
 * offset = distance from the center of the phantom in mm
 *
 */
    template<typename T>
    class ToyLinac {
    public:

        ODPM_INLINE constexpr explicit ToyLinac(const real_type RZ0, const real_type energy,
                                                const real_type source_axis_distance = 1000,
                                                const real_type fluence_field_length = 100,
                                                const DoubleVector<> offset = {})
                : m_rz0(RZ0),
                  m_source_axis_distance(source_axis_distance),
                  m_fluence_field_length(fluence_field_length),
                  m_offset(offset),
                  m_energy(energy),
                  m_cos_min{source_axis_distance /
                            math::sqrt(
                                    math::pow(source_axis_distance, 2) + 0.5 * math::pow(fluence_field_length, 2))} {};


        constexpr ToyLinac(const ToyLinac &) = default;


        constexpr ToyLinac(ToyLinac &&) = default;


        ToyLinac &operator=(ToyLinac &&) = default;

        ToyLinac &operator=(const ToyLinac &) = default;

        ODPM_CONST
        ODPM_INLINE constexpr T generateParticle(Random &rng) const noexcept {
            ThreeVector<> position;
            ThreeVector<> direction;
            do {
                const real_type cosTheta = m_cos_min + (1 - m_cos_min) * rng.getUniform();
                const real_type sinTheta = math::sqrt(1 - cosTheta * cosTheta);
                const real_type phi = 2 * M_PI * rng.getUniform();
                direction.x = sinTheta * math::cos(phi);
                direction.y = sinTheta * math::sin(phi);
                direction.z = cosTheta;
                position.x = 0;
                position.y = 0;
                position.z = -m_source_axis_distance;
                const real_type s = -position.z / direction.z;
                position.x += s * direction.x;
                position.y += s * direction.y;
                position.z = m_rz0;
            } while (math::abs(position.x) > 0.5 * m_fluence_field_length ||
                     math::abs(position.y) > 0.5 * m_fluence_field_length);
            position.x += m_offset.x;
            position.y += m_offset.y;
            return {position, direction, m_energy};
        }

#ifdef GPU
        void initializeGPU(const cuda::device_t& device) {}
#endif
    private:
        const real_type m_rz0;
        const real_type m_source_axis_distance;
        const real_type m_fluence_field_length;
        const DoubleVector<> m_offset;
        const real_type m_energy;
        const real_type m_cos_min;
    };

    using ToyElectronGun = ToyLinac<Electron>;
    using ToyPhotonGun = ToyLinac<Photon>;

    template<typename T>
    class IRCPhaseSpace {
    public:

        IRCPhaseSpace(const std::string filename, const float distance, const DoubleVector<real_type> &center = {})
                : k_distance(distance), k_center(center) {
            std::ifstream infile(filename);
            if (!infile.is_open()) {
                throw std::runtime_error("Impossible to open phase space file " + filename);
            }
            std::string line;
            std::getline(infile, line);  // skip the first line

            while (std::getline(infile, line)) {
                BeamParticle ray;
                ray.fromString(line, k_delim);
                if (ray.m_fluence == 0) {
                    continue;
                }
                ray.m_particle.position.z = -distance;
                m_particle_list.push_back(ray);
            }
        }

        ODPM_PURE
        T generateParticle(Random &rng) const override {
            static std::uniform_int_distribution<unsigned> distribution{0,
                                                                        static_cast<unsigned>(m_particle_list.size())};
            unsigned trialPoint;
            if (p_previous_ray == nullptr) {
                trialPoint = rng.getSampled(distribution);
                p_previous_ray = &m_particle_list[trialPoint];
            }
            real_type r;
            do {
                trialPoint = rng.getSampled(distribution);
                r = rng.getUniform();
            } while (p_previous_ray->m_fluence * r >= m_particle_list[trialPoint].m_fluence);
            p_previous_ray = &m_particle_list[trialPoint];
            const auto ndotu = k_plane_normal.dot(p_previous_ray->m_particle.direction);
            const auto w = k_plane_point - p_previous_ray->m_particle.position;
            const auto si = k_plane_normal.dot(w) / ndotu;
            auto psi = p_previous_ray->m_particle.position + p_previous_ray->m_particle.direction * si;
            psi.x += k_center.x;
            psi.y += k_center.y;
            psi.z = 0;
            return {psi, p_previous_ray->m_particle.direction, p_previous_ray->m_particle.energy};
        }

    private:
        class BeamParticle {
        public:
            T m_particle;
            float m_fluence;

            BeamParticle(T &particle, float fluence) : m_particle(particle), m_fluence(fluence) {}

            BeamParticle() : m_particle({}, {}, 0), m_fluence(0) {}

            void fromString(const std::string &line, const std::string &delim) {
                auto start = 0U;
                auto end = line.find(delim);

                m_particle.position.x = std::stof(line.substr(start, end - start));
                start = end + delim.length();
                end = line.find(delim, start);

                m_particle.position.y = std::stof(line.substr(start, end - start));
                start = end + delim.length();
                end = line.find(delim, start);

                //            particle.position.z    = std::stof(line.substr(start, end - start));
                start = end + delim.length();
                end = line.find(delim, start);

                m_particle.direction.x = std::stof(line.substr(start, end - start));
                start = end + delim.length();
                end = line.find(delim, start);

                m_particle.direction.y = std::stof(line.substr(start, end - start));
                start = end + delim.length();
                end = line.find(delim, start);

                m_particle.direction.z = std::stof(line.substr(start, end - start));
                start = end + delim.length();
                end = line.find(delim, start);

                m_particle.energy = std::stof(line.substr(start, end - start));
                start = end + delim.length();
                end = line.find(delim, start);

                m_fluence = std::stof(line.substr(start, end));
            }

            friend std::ostream &operator<<(std::ostream &os, const BeamParticle &particle) {
                os << "particle: " << particle.m_particle << " fluence: " << particle.m_fluence;
                return os;
            }
        };

        const ThreeVector<real_type> k_plane_normal{0, 0, 1};
        const ThreeVector<real_type> k_plane_point{0, 0, 0};
        std::vector<BeamParticle> m_particle_list;
        BeamParticle *p_previous_ray = nullptr;
        const float k_distance;
        const DoubleVector<real_type> k_center;
        const std::string k_delim = " ";
    };

}  // namespace opmc
#endif  // ODPM_CPU_INCLUDE_BEAMS_H_
