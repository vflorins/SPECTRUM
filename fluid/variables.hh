/*!
\file variables.hh
\brief A container class for multiple fluid variables
\author Vladimir Florinski

This file is part of the SPECTRUM suite of scientific numerical simulation codes. SPECTRUM stands for Space Plasma and Energetic Charged particle TRansport on Unstructured Meshes. The code simulates plasma or neutral particle flows using MHD equations on a grid, transport of cosmic rays using stochastic or grid based methods. The "unstructured" part refers to the use of a geodesic mesh providing a uniform coverage of the surface of a sphere.
*/

#ifndef SPECTRUM_VARIABLES
#define SPECTRUM_VARIABLES

#include "config.h"
#include "common/physics.hh"

// Possible fluid types:
// GASDYN - compressible gas dynamics (5 vars)
// MHD - ideal MHD (8 vars)
// MHDE - two-fluid ideal MHD (9 vars)
// CGL - anisotropic ideal MHD (9 vars)
// CGLE - anisotropic two-fluid ideal MHD (9 vars)

// TODO: These should go into config.h
#define FLUID_GASDYN 301
#define FLUID_MHD    302
#define FLUID_MHDE   303
#define FLUID_CGL    304
#define FLUID_CGLE   305
#define TURB_ZANK6EQ 321

#define PASSIVE_PRIMARY_01
#define PASSIVE_PRIMARY_02

#define FLUID_PRIMARY_TYPE FLUID_MHD
#define FLUID_PRIMARY_SPECIE Specie::proton

#define FLUID_SECONDARY_TYPE FLUID_GASDYN
#define FLUID_SECONDARY_SPECIE Specie::proton

#define FLUID_TERTIARY_TYPE FLUID_GASDYN
#define FLUID_TERTIARY_SPECIE Specie::proton

#undef FLUID_QUATERNARY_TYPE
#undef FLUID_QUATERNARY_SPECIE

#undef TURBULENCE_TYPE

#if (FLUID_PRIMARY_TYPE == FLUID_GASDYN) || (FLUID_SECONDARY_TYPE == FLUID_GASDYN) || (FLUID_TERTIARY_TYPE == FLUID_GASDYN) || (FLUID_QUATERNARY_TYPE == FLUID_GASDYN)
#include "fluid/conservation_laws_gasdyn.hh"
#endif

#if (FLUID_PRIMARY_TYPE == FLUID_MHD) || (FLUID_SECONDARY_TYPE == FLUID_MHD) || (FLUID_TERTIARY_TYPE == FLUID_MHD) || (FLUID_QUATERNARY_TYPE == FLUID_MHD)
#include "fluid/conservation_laws_mhd.hh"
#endif

#if TURBULENCE_TYPE == TURB_ZANK6EQ
#include "fluid/conservation_laws_zank6eq.hh"
#endif

#if defined(PASSIVE_PRIMARY_01) || defined(PASSIVE_PRIMARY_02) || defined(PASSIVE_SECONDARY_01) || defined (PASSIVE_SECONDARY_02)
#include "fluid/conservation_laws_passive.hh"
#endif

namespace Spectrum {

/*!
\brief A container combining all conserved fluid variables into single contiguous storage
\author Vladimir Florinski
*/
struct ConservedVariables
{

#if FLUID_PRIMARY_TYPE == FLUID_GASDYN
   Gasdyn<FLUID_PRIMARY_SPECIE>::ConservedState var_primary;
#elif FLUID_PRIMARY_TYPE == FLUID_MHD
   MHD<FLUID_PRIMARY_SPECIE>::ConservedState var_primary;
#endif

#ifdef PASSIVE_PRIMARY_01
   Passive::ConservedState passive_primary_01;
#endif

#ifdef PASSIVE_PRIMARY_02
   Passive::ConservedState passive_primary_02;
#endif

#if FLUID_SECONDARY_TYPE == FLUID_GASDYN
   Gasdyn<FLUID_SECONDARY_SPECIE>::ConservedState var_secondary;
#elif FLUID_SECONDARY_TYPE == FLUID_MHD
   MHD<FLUID_SECONDARY_SPECIE>::ConservedState var_secondary;
#endif

#ifdef PASSIVE_SECONDARY_01
   Passive::ConservedState passive_secondary_01;
#endif

#ifdef PASSIVE_SECONDARY_02
   Passive::ConservedState passive_secondary_02;
#endif

#if FLUID_TERTIARY_TYPE == FLUID_GASDYN
   Gasdyn<FLUID_TERTIARY_SPECIE>::ConservedState var_tertiary;
#elif FLUID_TERTIARY_TYPE == FLUID_MHD
   MHD<FLUID_TERTIARY_SPECIE>::ConservedState var_tertiary;
#endif

#ifdef PASSIVE_TERTIARY_01
   Passive::ConservedState passive_tertiary_01;
#endif

#ifdef PASSIVE_TERTIARY_02
   Passive::ConservedState passive_tertiary_02;
#endif

#if FLUID_QUATERNARY_TYPE == FLUID_GASDYN
   Gasdyn<FLUID_QUATERNARY_SPECIE>::ConservedState var_quaternary;
#elif FLUID_QUATERNARY_TYPE == FLUID_MHD
   MHD<FLUID_QUATERNARY_SPECIE>::ConservedState var_quaternary;
#endif

#ifdef PASSIVE_QUATERNARY_01
   Passive::ConservedState passive_quaternary_01;
#endif

#ifdef PASSIVE_QUATERNARY_02
   Passive::ConservedState passive_quaternary_02;
#endif

#if TURBULENCE_TYPE == TURB_ZANK6EQ
   Zank6eq::ConservedState turbulence;
#endif

};

};

#endif
