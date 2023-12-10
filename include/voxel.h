//
// Created by mbarb on 3/17/2020.
//
#include <omp.h>

#include "materials.h"

namespace opmc {

// Not the best design for this class, but it must be trivially constructable so that it can be used on the GPU with no
// issues. The geometry should take care of initializing it.

    struct Voxel {
        double &dose;
        const Material &material;

        template<typename T>
        ODPM_INLINE void addDose(const T new_dose) noexcept {
            assert(new_dose >= 0);
#pragma omp atomic update relaxed hint(omp_sync_hint_speculative)
            dose += new_dose;
        }
    };
}  // namespace opmc

#ifndef ODPM_VOXEL_H
#define ODPM_VOXEL_H

#endif  // ODPM_VOXEL_H
