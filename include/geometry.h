//
// Created by mbarb on 3/23/2020.
//

#ifndef ODPM_VOXEL_MAP_H
#define ODPM_VOXEL_MAP_H

#include <absl/types/optional.h>

#include <experimental/mdspan>

#include "elementary_functions.h"

#include "utils.h"
#include "voxel.h"

#ifdef GPU
#include "gpu_array.h"
#include "gpu_unique_ptr.h"
#endif

namespace opmc {

    template<typename T = unsigned long>
    class DoseDistribution {
    private:
        std::vector<double> m_data;
        using DYNAMIC_3D_EXTENDS =
                std::experimental::extents<std::size_t, std::experimental::dynamic_extent, std::experimental::dynamic_extent,
                        std::experimental::dynamic_extent>;
        using DISTRIBUTION_TYPE = std::experimental::mdspan<double, DYNAMIC_3D_EXTENDS, std::experimental::layout_left>;
        DISTRIBUTION_TYPE m_access_wrapper;
        ThreeVector<> m_voxelSize;
    public:
        constexpr DoseDistribution(T x, T y, T z, ThreeVector<> voxelSize)
                : m_data(x * y * z), m_access_wrapper(m_data.data(), DYNAMIC_3D_EXTENDS(x, y, z)),
                  m_voxelSize(voxelSize) {}

        constexpr DoseDistribution(T x, T y, T z, ThreeVector<> voxelSize, T *data) :
                m_access_wrapper(m_data.data(), DYNAMIC_3D_EXTENDS(x, y, z)),
                m_voxelSize(voxelSize) {
            std::copy(data, data + x * y * z, std::back_inserter(m_data));
        }

        template<typename V>
        constexpr DoseDistribution(T x, T y, T z, ThreeVector<> voxelSize, V &&data)
                : m_data(data), m_access_wrapper(m_data.data(), DYNAMIC_3D_EXTENDS(x, y, z)),
                  m_voxelSize(voxelSize) {}

        constexpr DoseDistribution()
                : m_data(), m_access_wrapper(m_data.data(), DYNAMIC_3D_EXTENDS()), m_voxelSize({}) {}

        constexpr DoseDistribution(const DoseDistribution &other) = default;

        constexpr DoseDistribution(DoseDistribution &&other) noexcept = default;

        constexpr DoseDistribution &operator=(const DoseDistribution &other) {
            std::copy(other.m_data.begin(), other.m_data.end(), std::back_inserter(m_data));
            m_access_wrapper = {m_data.data(), other.m_access_wrapper.extents()};
            return *this;
        };

        constexpr DoseDistribution &operator=(DoseDistribution &&other) noexcept {
            std::swap(m_data, other.m_data);
            m_access_wrapper = {m_data.data(), other.m_access_wrapper.extents()};
            return *this;
        };

        constexpr ODPM_INLINE auto &operator[](T x, T y, T z) { return m_access_wrapper[x, y, z]; }

        constexpr ODPM_INLINE auto &operator[](T x, T y, T z) const { return m_access_wrapper[x, y, z]; }

        constexpr void normalize() noexcept {
            const auto max_value = std::max_element(m_data.begin(), m_data.end());
            std::transform(m_data.begin(), m_data.end(), m_data.begin(),
                           [max_value](const auto &value) { return value / *max_value; });
        }

        constexpr void filter(const float threshold) noexcept {
            // set all elements less than threshold percent of the max to zero
            const auto max_value = *std::max_element(m_data.begin(), m_data.end());
            // calculate the threshold percentage of max
            const auto threshold_value = threshold * max_value;

            std::transform(m_data.begin(), m_data.end(), m_data.begin(),
                           [threshold_value](const auto &value) { return (value > threshold_value) ? value : 0; });
        }

        constexpr auto data() noexcept { return m_data.data(); }

        constexpr ThreeVector<T> dimensions() const noexcept {
            return {m_access_wrapper.extent(0), m_access_wrapper.extent(1), m_access_wrapper.extent(2)};
        }

        constexpr ThreeVector<> voxelSize() const noexcept {
            return m_voxelSize;
        }

    };

    class Geometry {
    protected:
        using DYNAMIC_3D_EXTENDS = std::experimental::dextents<std::size_t, 3>;
        using CONTAINED_T = Voxel;
        template<typename T>
        using MDSPAN_TYPE = std::experimental::mdspan<T, DYNAMIC_3D_EXTENDS, std::experimental::layout_left>;
    };

    class HalfDistanceVoxelCube : public Geometry {
    public:

        HalfDistanceVoxelCube() = default;

        explicit HalfDistanceVoxelCube(const ThreeVector<unsigned long> &dimensions,
                                       const ThreeVector<real_type> &resolution, const std::vector<Material> &materials)
                : m_dimensions(dimensions),
                  m_doses(dimensions.x * dimensions.y * dimensions.z),
                  m_tracks(dimensions.x * dimensions.y * dimensions.z),
                  m_materials(dimensions.x * dimensions.y * dimensions.z),
                  m_doses_span{m_doses.data(), DYNAMIC_3D_EXTENDS(dimensions.x, dimensions.y, dimensions.z)},
                  m_tracks_span{m_tracks.data(), DYNAMIC_3D_EXTENDS(dimensions.x, dimensions.y, dimensions.z)},
                  m_materials_span{m_materials.data(), DYNAMIC_3D_EXTENDS(dimensions.x, dimensions.y, dimensions.z)},
                  m_resolution(resolution),
                  m_inverse_resolution(ThreeVector{1., 1., 1.} / resolution),
                  m_half_resolution(resolution / 2.),
                  m_half_distance(dimensions.x / 2, dimensions.y / 2, 0),
                  m_material_list(materials) {
            for (auto i = 0ULL; i < dimensions.x * dimensions.y * dimensions.z; ++i) {
                m_doses[i] = 0;
                m_tracks[i] = 0;
                m_materials[i] = 0;
            }
        }

        HalfDistanceVoxelCube(HalfDistanceVoxelCube &&) = default;

        HalfDistanceVoxelCube(const HalfDistanceVoxelCube &) = default;

        virtual ~HalfDistanceVoxelCube() = default;

        HalfDistanceVoxelCube &operator=(HalfDistanceVoxelCube &&) = default;

        template<typename T>
        ODPM_INLINE constexpr absl::optional<CONTAINED_T> operator[](const ThreeVector<T> &i) noexcept {
            const auto position = getPosition(i);
            if (ODPM_LIKELY(rangeCheck(position))) {
                return CONTAINED_T{accessSpan(m_doses_span, position), accessSpan(m_tracks_span, position),
                                   getMaterial(position)};
            }
            return absl::nullopt;
        }

        ODPM_INLINE constexpr long numElements() const noexcept { return m_doses_span.size(); }

        ODPM_INLINE constexpr float *origin_dose() const noexcept { return m_doses_span.data_handle(); }

        ODPM_INLINE constexpr std::uint8_t *origin_material() const noexcept { return m_materials_span.data_handle(); }

        ODPM_INLINE constexpr real_type distanceToBoundary(const Particle &particle) const noexcept {
            // Floor the particle position
            // Granularity. Resolution is the voxel size
            auto position = (particle.position + m_half_resolution) * m_inverse_resolution;
            auto relativePosition = particle.position - (position.floor() * m_resolution);
#ifdef VERBOSE
            std::cout << voxel_position << std::endl;
            std::cout << relativePosition << std::endl;
#endif
            const auto &direction = particle.direction;
            ThreeVector<real_type> distance{std::numeric_limits<real_type>::max(),
                                            std::numeric_limits<real_type>::max(),
                                            std::numeric_limits<real_type>::max()};
            if (direction.x > 0) {
                distance.x = (m_half_resolution.x - relativePosition.x) / direction.x;
            }
            if (direction.x < 0) {
                distance.x = -(m_half_resolution.x + relativePosition.x) / direction.x;
            }
            assert(distance.x >= 0);

            if (direction.y > 0) {
                distance.y = (m_half_resolution.y - relativePosition.y) / direction.y;
            }
            if (direction.y < 0) {
                distance.y = -(m_half_resolution.y + relativePosition.y) / direction.y;
            }
            assert(distance.y >= 0);

            if (direction.z > 0) {
                distance.z = (m_half_resolution.z - relativePosition.z) / direction.z;
            }
            if (direction.z < 0) {
                distance.z = -(m_half_resolution.z + relativePosition.z) / direction.z;
            }
            assert(distance.z >= 0);
            // check for small steps and handle them appropriately
            return distance.min();
        }

        std::vector<real_type> doseDepthDistribution() const;

        std::vector<real_type> doseDepthProfile() const;

        void write(const std::string &filename) const;

        const DoseDistribution<unsigned long> doseDistribution() const;

        void setDefaultMaterial(const u_int8_t id) noexcept {
            for (auto &element: m_materials) {
                element = id;
            }
        };

        template<typename V>
        void copyFromGeometry(const V &geometry) {
            static_assert(std::is_same<typeof(*this), V>(), "Geometries must have same type");
            if (numElements() != geometry.numElements()) {
                throw std::runtime_error("Geometries size mismatch");
            }
            for (auto i = 0LL; i < numElements(); ++i) {
                m_materials[i] = geometry.m_materials[i];
                m_doses[i] = 0;
                m_tracks[i] = 0;
            }
        }

        void normalize() noexcept {
            const auto max_value = std::max_element(m_doses.begin(), m_doses.end());
            std::transform(m_doses.begin(), m_doses.end(), m_doses.begin(),
                           [max_value](const auto &value) { return value / *max_value; });
        }

        constexpr auto &dimension() const {
            return m_dimensions;
        };

        constexpr auto &resolution() const {
            return m_resolution;
        }

    protected:
        ThreeVector<unsigned long> m_dimensions;

        std::vector<float> m_doses;
        std::vector<int> m_tracks;
        std::vector<std::uint8_t> m_materials;

        MDSPAN_TYPE<float> m_doses_span;
        MDSPAN_TYPE<int> m_tracks_span;
        MDSPAN_TYPE <std::uint8_t> m_materials_span;

        ThreeVector<real_type> m_resolution;


        ThreeVector<real_type> m_inverse_resolution;
        ThreeVector<real_type> m_half_resolution;
        ThreeVector<unsigned long> m_half_distance;

        std::vector<Material> m_material_list;

        template<typename T>
        ODPM_INLINE constexpr bool rangeCheck(const ThreeVector<T> position) const noexcept {
            if constexpr (std::is_unsigned_v<T>) {
                return position < m_dimensions;
            } else {
                return position >= ThreeVector<T>{0, 0, 0} && position < m_dimensions;
            }
        }

        template<typename T>
        ODPM_INLINE constexpr ThreeVector<long> getPosition(const ThreeVector<T> &i) const noexcept {
            const auto voxel_position = (i + m_half_resolution) * m_inverse_resolution;
            ThreeVector<long> position{static_cast<long>(std::floor(voxel_position.x)),
                                       static_cast<long>(std::floor(voxel_position.y)),
                                       static_cast<long>(std::floor(voxel_position.z))};
            return position + m_half_distance;
        }

        template<typename T, typename V>
        ODPM_INLINE static constexpr T &accessSpan(MDSPAN_TYPE <T> &span, const ThreeVector<V> &position) noexcept {
            return span[position.x, position.y, position.z];
        }

        ODPM_INLINE constexpr const Material &getMaterial(const std::uint8_t id) const noexcept {
            return m_material_list[id];
        }

        template<class V>
        ODPM_INLINE constexpr const Material &getMaterial(const ThreeVector<V> &position) const noexcept {
            return getMaterial(m_materials_span[position.x, position.y, position.z]);
        }
    };

    class VoxelCube : public HalfDistanceVoxelCube {
    public:
        VoxelCube() = default;

        explicit VoxelCube(const ThreeVector<unsigned long> &dimensions, const ThreeVector<real_type> &resolution,
                           const std::vector<Material> &material_list)
                : HalfDistanceVoxelCube(dimensions, resolution, material_list) {};

        VoxelCube(VoxelCube &&) = default;

        VoxelCube(const VoxelCube &) = default;

        virtual ~VoxelCube() = default;

        template<typename T>
        ODPM_INLINE constexpr ThreeVector<T> getPosition(const ThreeVector<T> &i) const noexcept {
            return i * m_inverse_resolution;
        }

        template<typename T>
        ODPM_INLINE constexpr absl::optional<CONTAINED_T> operator[](const ThreeVector<T> &i) noexcept {
            const auto position = static_cast<ThreeVector<int>>(getPosition(i));
            if (ODPM_LIKELY(rangeCheck(position))) {
                return CONTAINED_T{accessSpan(m_doses_span, position), accessSpan(m_tracks_span, position),
                                   getMaterial(position)};
            }
            return absl::nullopt;
        }


        ODPM_INLINE constexpr real_type distanceToBoundary(const Particle &particle) const noexcept {
            // if the machine supports IEEE-754 dividing a non-zero number by ±0.0 gives the correctly-signed infinity and
            // FE_DIVBYZERO is raised. Since I check that IEEE-754 is supported, the branches can be removed.
            const auto minBoundary = m_resolution * getPosition(particle.position).floor();
            const auto maxBoundary = minBoundary + m_resolution;
#if defined(USE_AGGRESSIVE_OPTIMIZATIONS)
            //        // some bitwise trickery to remove branches...
            const real_type boundary_x[2] = {maxBoundary.x, minBoundary.x};
            const real_type boundary_y[2] = {maxBoundary.y, minBoundary.y};
            const real_type boundary_z[2] = {maxBoundary.z, minBoundary.z};
            const ThreeVector boundary{boundary_x[math::signbit(particle.direction.x)],
                                       boundary_y[math::signbit(particle.direction.y)],
                                       boundary_z[math::signbit(particle.direction.z)]};
#else
            const ThreeVector boundary{(particle.direction.x > 0) ? maxBoundary.x : minBoundary.x,
                                       (particle.direction.y > 0) ? maxBoundary.y : minBoundary.y,
                                       (particle.direction.z > 0) ? maxBoundary.z : minBoundary.z};
#endif
            const auto distance = (boundary - particle.position) / particle.direction;

            // note: need to check for micro steps and handle them appropriately
            // If the implementation supports IEEE floating-point arithmetic (IEC 60559)
            // math::min() -> returns:
            // If one of the two arguments is NaN, the value of the other argument is returned.
            // Only if both arguments are NaN, NaN is returned
            // in case one of these is not nan the is correctly returned
            const auto result = distance.min();
            assert(result >= 0 && !std::isnan(result));
            return result;
        }
    };

#ifdef TEST

    class TestHalfDistanceVoxelCube : public HalfDistanceVoxelCube {
       public:
        TestHalfDistanceVoxelCube(const ThreeVector<unsigned long> &dimensions, const ThreeVector<real_type> &resolution,
                                  const std::vector<Material> &material_list)

            : HalfDistanceVoxelCube(dimensions, resolution, material_list) {}

        template <typename T>
        ODPM_INLINE constexpr ThreeVector<long> getPosition(const ThreeVector<T> &i) const {
            return HalfDistanceVoxelCube::getPosition(i);
        }
    };

#endif

}  // namespace opmc

#endif  // ODPM_VOXEL_MAP_H
