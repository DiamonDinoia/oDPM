//
// Created by mbarbone on 7/29/21.
//

#include "geometry.h"

#include <array>
#include <fstream>
#include <iostream>

#include "materials.h"
#include "utils.h"

namespace opmc {
//

    std::vector<real_type> HalfDistanceVoxelCube::doseDepthDistribution() const {
        std::vector<real_type> histogram(m_dimensions.z, 0.);
        for (auto z = 0UL; z < m_dimensions.z; ++z) {
            for (auto y = 0UL; y < m_dimensions.y; ++y) {
                for (auto x = 0UL; x < m_dimensions.x; ++x) {
                    histogram[z] += m_doses_span[x, y, z];
                }
            }
        }
        return histogram;
    }

    std::vector<real_type> HalfDistanceVoxelCube::doseDepthProfile() const {
        std::vector<real_type> histogram(m_dimensions.z, 0.);
        auto x = m_half_distance.x;
        auto y = m_half_distance.y;
        for (auto z = 0UL; z < m_dimensions.z; ++z) {
            histogram[z] = m_doses_span[x, y, z];
        }
//        real_type norm = 1. / *max_element(histogram.begin(), histogram.end());
//        for (auto &item: histogram) {
//            item *= norm;
//        }
        return histogram;
    }

    void HalfDistanceVoxelCube::write(const std::string &filename) const {
        std::cout << "writing dose cube to " << filename << std::endl;
        std::ofstream outStream(filename, std::ios::binary | std::ios::out);
        if (!outStream.good()) {
            throw std::runtime_error("CAN'T OPEN " + filename);
        }
        std::cout << "dimensions " << m_dimensions.x << " " << m_dimensions.y << " " << m_dimensions.z << std::endl;
        auto dimX = static_cast<unsigned_type>(m_dimensions.x);
        outStream.write((char *) &dimX, sizeof(unsigned_type));
        auto dimY = static_cast<unsigned_type>(m_dimensions.y);
        outStream.write((char *) &dimY, sizeof(unsigned_type));
        auto dimZ = static_cast<unsigned_type>(m_dimensions.z);
        outStream.write((char *) &dimZ, sizeof(unsigned_type));
        std::cout << "resolution " << m_resolution << std::endl;
        outStream.write((char *) &m_resolution.x, sizeof(real_type));
        outStream.write((char *) &m_resolution.y, sizeof(real_type));
        outStream.write((char *) &m_resolution.z, sizeof(real_type));
        for (unsigned int z = 0; z < dimZ; ++z) {
            for (unsigned int y = 0; y < dimY; ++y) {
                for (unsigned int x = 0; x < dimX; ++x) {
                    real_type value = m_doses_span[x, y, z];
                    outStream.write((char *) &value, sizeof(real_type));
                }
            }
        }
        outStream.close();
    }

    const DoseDistribution<> HalfDistanceVoxelCube::doseDistribution() const {
        DoseDistribution distribution{m_dimensions.x, m_dimensions.y, m_dimensions.z, m_resolution};
        for (unsigned int x = 0; x < m_dimensions.x; ++x) {
            for (unsigned int y = 0; y < m_dimensions.y; ++y) {
                for (unsigned int z = 0; z < m_dimensions.z; ++z) {
                    const auto dose = m_doses_span[x, y, z];
                    distribution[x, y, z] = dose;
                }
            }
        }
        return distribution;
    }


}  // namespace opmc