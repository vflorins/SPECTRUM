/*!
\file compiletime_lists.hh
\brief Defines global compile-time constant options that are available, and some associated constant lists.
\author Vladimir Florinski
\author Juan G Alonso Guzman
\author Lucius Schoenbaum

This file is part of the SPECTRUM suite of scientific numerical simulation codes. SPECTRUM stands for Space Plasma and Energetic Charged particle TRansport on Unstructured Meshes. The code simulates plasma or neutral particle flows using MHD equations on a grid, transport of cosmic rays using stochastic or grid based methods. The "unstructured" part refers to the use of a geodesic mesh providing a uniform coverage of the surface of a sphere.
*/

#ifndef SPECTRUM_COMPILETIME_LISTS_HH
#define SPECTRUM_COMPILETIME_LISTS_HH

#ifdef USE_GSL
#include <gsl/gsl_const_cgsm.h>
#endif

#include "common/compiletime_math.hh"
#include "common/specie.hh"


namespace Spectrum {

//----------------------------------------------------------------------------------------------------------------------------------------------------
// Fields options
//----------------------------------------------------------------------------------------------------------------------------------------------------


//! Coordinate system - defines curvilinear coordinates
enum class CoordinateSystem {
   Cartesian,
   Cylindrical,
   SphericalRTP,
   SphericalRMP,
   Custom,
   Polar,
   PitchAngle, /* magnitude, 'mu', - */
   Anisotropic, /* perpendicular, - , parallel */
};


//----------------------------------------------------------------------------------------------------------------------------------------------------
// Fluid Specie / Convervation Law options
//----------------------------------------------------------------------------------------------------------------------------------------------------


enum class Model {
// GASDYN - compressible gas dynamics (5 vars)
   GasDyn,
// MHD - ideal MHD (8 vars)
   MHD,
// MHDE - two-fluid ideal MHD (9 vars)
   MHDGLM,
   MHDE,
// CGL - anisotropic ideal MHD (9 vars)
   CGL,
// CGLE - anisotropic two-fluid ideal MHD (9 vars)
   CGLE,
};

enum class TurbulenceModel {
   Zank6eq,
};

enum class Form {
   primitive,
   conserved,
   flux,
};

enum class Passivity {
   active,
   passive,
};


//! Indicates that the argument should be std::ratio
#define Ratio typename


//! Converts std::ratio to double
template <typename ratio, typename Float = double>
constexpr Float get_fp() {
   return static_cast<Float>(ratio::num)/static_cast<double>(ratio::den);
}




namespace Impl {

template<bool cond, typename IfTrue, typename IfFalse>
struct Cond_impl;

template<typename IfTrue, typename IfFalse>
struct Cond_impl<true, IfTrue, IfFalse> {
   using Cond = IfTrue;
};

template<typename IfTrue, typename IfFalse>
struct Cond_impl<false, IfTrue, IfFalse> {
   using Cond = IfFalse;
};

}

//! Conditional type
template <bool cond, typename IfTrue, typename IfFalse>
using Cond = typename Impl::Cond_impl<cond, IfTrue, IfFalse>::Cond;


};

#endif
