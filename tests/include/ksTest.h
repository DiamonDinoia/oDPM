//
// Created by mbarbone on 7/26/21.
//

#ifndef DPM_CPU_TESTS_INCLUDE_KSTEST_H_
#define DPM_CPU_TESTS_INCLUDE_KSTEST_H_

#include <ostream>
#include <vector>

namespace opmc {
namespace internal {

double KsTestD(const std::vector<double> &vector1, const std::vector<double> &vector2);
}

template <typename T, typename V>
double KsTest(const std::vector<T> &vector1, const std::vector<V> &vector2) {
    std::vector<double> sample1(vector1.begin(), vector1.end());
    std::vector<double> sample2(vector2.begin(), vector2.end());
    return internal::KsTestD(sample1, sample2);
}

}  // namespace opmc
#endif  // DPM_CPU_TESTS_INCLUDE_KSTEST_H_
