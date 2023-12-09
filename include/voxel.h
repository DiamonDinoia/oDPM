//
// Created by mbarb on 3/17/2020.
//
#include <omp.h>

#include "materials.h"

namespace opmc {

// Not the best design for this class, but it must be trivially constructable so that it can be used on the GPU with no
// issues. The geometry should take care of initializing it.

struct Voxel {
    float& dose;
    int& histories;
    const Material& material;


    ODPM_INLINE constexpr Voxel(float& dose, int& histories, const Material& material) noexcept
        : dose(dose), histories(histories), material(material) {}

    template <typename T>
     ODPM_INLINE void addDose(const T new_dose) noexcept {
        const auto casted_dose = static_cast<float>(new_dose);
        assert(new_dose >= 0);
#ifndef __CUDA_ARCH__
#pragma omp atomic update relaxed hint(omp_sync_hint_speculative)
        dose += casted_dose;
#pragma omp atomic update relaxed hint(omp_sync_hint_speculative)
        histories++;
#else
        atomicAdd(&dose, casted_dose);
        atomicAdd(&histories, 1);
#endif
    }
};
}  // namespace opmc

#ifndef ODPM_VOXEL_H
#define ODPM_VOXEL_H

#endif  // ODPM_VOXEL_H
