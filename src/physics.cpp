//
// Created by mbarbone on 7/28/21.
//

#include "physics.h"

namespace opmc {

std::ostream& operator<<(std::ostream& out, const Event value) {
    const char* s = nullptr;
#define PROCESS_VAL(p) \
    case (p):          \
        s = #p;        \
        break;
    switch (value) {
        PROCESS_VAL(NOTHING);
        PROCESS_VAL(BREMSSTRAHLUNG);
        PROCESS_VAL(MOLLER);
        PROCESS_VAL(HINGE);
        PROCESS_VAL(OUT_OF_BOUNDARY);
        PROCESS_VAL(ABSORBED);
        PROCESS_VAL(MICROSTEP);
    }
#undef PROCESS_VAL
    return out << s;
}

}  // namespace opmc