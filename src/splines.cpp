//
// Copyright M. Novak: 2021
// Modified by M.Barbone 2021
//

#include <array>
#include <fstream>
#include <sstream>

#include "materials.h"
#include "spline_interpolator.h"

namespace opmc {

void SplineInterpolator::initialize() {
    m_num_data      = m_data.size();

    m_log_min_x     = math::log(m_data[0].x);
    m_inv_log_delta = 1. / (math::log(m_data[1].x) - m_log_min_x);

    int m1          = 1;
    int m2          = m_num_data - 1;
    double s        = 0.0;
    double r        = 0.0;
    for (int m = 0; m < m2; ++m) {
        m_data[m].d = m_data[m + 1].x - m_data[m].x;
        r           = (m_data[m + 1].y - m_data[m].y) / m_data[m].d;
        m_data[m].c = r - s;
        s           = r;
    }
    r            = 0.0;
    s            = 0.0;
    m_data[0].c  = 0.0;
    m_data[m2].c = 0.0;
    for (int m = m1; m < m2; ++m) {
        m_data[m].c = m_data[m].c + r * m_data[m - 1].c;
        m_data[m].b = 2.0 * (m_data[m - 1].x - m_data[m + 1].x) - r * s;
        s           = m_data[m].d;
        r           = s / m_data[m].b;
    }
    int mr = m2 - 1;
    for (int m = m1; m < m2; ++m) {
        m_data[mr].c = (m_data[mr].d * m_data[mr + 1].c - m_data[mr].c) / m_data[mr].b;
        --mr;
    }
    for (int m = 0; m < m2; ++m) {
        s           = m_data[m].d;
        r           = m_data[m + 1].c - m_data[m].c;
        m_data[m].d = r / s;
        m_data[m].c = 3.0 * m_data[m].c;
        m_data[m].b = (m_data[m + 1].y - m_data[m].y) / s - (m_data[m].c + r) * s;
        //    fData[m].fA = fData[m].fY;
    }
}
std::ostream& operator<<(std::ostream& os, const Point5D& point) {
    os << "x: " << point.x << " y: " << point.y << " b: " << point.b << " c: " << point.c << " d: " << point.d;
    return os;
}

}  // namespace opmc