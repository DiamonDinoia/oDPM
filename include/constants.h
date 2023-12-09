//
// Copyright M. Novak: 2021
// Modified by M.Barbone 2021
//

#ifndef ODPM_CPU_INCLUDE_CONSTANTS_H_
#define ODPM_CPU_INCLUDE_CONSTANTS_H_

#include <G4SystemOfUnits.hh>

#include "types.h"
namespace opmc::constants {

constexpr real_type k_PI              = 3.1415926535897932;
constexpr real_type k_EMC2            = 0.510991;
constexpr real_type k_2EMC2           = 2.0 * k_EMC2;
constexpr real_type k_InvEMC2         = 1.0 / k_EMC2;
constexpr real_type k_HalfSqrt2EMC2   = k_EMC2 * 0.7071067812;
constexpr real_type k_ElectronCut     = 0.20;  // 200 [keV]
constexpr real_type k_GammaCut        = 0.05;  //  50 [keV]
constexpr real_type k_MaxEkin         = 21.;   //  21 [MeV]
constexpr real_type k_MinEkin         = 0.001;
constexpr real_type k_lowInvMFPVacuum = 1.0E-20;
constexpr real_type k_EMaxThreshold   = 1.0E-6;
constexpr real_type k_minStepLength   = 1.0E-4;  // [mm]
constexpr real_type k_MscStepSlow     = 5.0;     //   5 [mm]
constexpr real_type k_MscStepShigh    = 10.0;    //   1 [cm]
constexpr real_type k_MscEcross       = 12.0;    //  12 [MeV]
constexpr int k_NumEkinForElectron    = 128;
constexpr int k_NumEkinForPhoton      = 1024;
constexpr real_type k_Gray            = CLHEP::gray;
constexpr real_type k_eV              = CLHEP::eV;

}  // namespace opmc::constants

#endif  // ODPM_CPU_INCLUDE_CONSTANTS_H_
