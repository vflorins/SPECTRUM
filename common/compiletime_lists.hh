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

#include "common/rk_lists.hh"
#include "common/specie.hh"

namespace Spectrum {


//----------------------------------------------------------------------------------------------------------------------------------------------------
// General options
//----------------------------------------------------------------------------------------------------------------------------------------------------


//! Build mode (debug build or production build)
enum class BuildMode {
   debug,
   release
};

//----------------------------------------------------------------------------------------------------------------------------------------------------
// Fields options
//----------------------------------------------------------------------------------------------------------------------------------------------------


enum class CoordinateSystem {
   cartesian,
   polar,
   cylindrical,
   spherical,
   /* In cases where some of the physical coordinates present are not used. */
   none,
   pitchangle, /* magnitude, 'mu', - */
   anisotropic, /* perpendicular, parallel, - */ /* CHECK */
};

//----------------------------------------------------------------------------------------------------------------------------------------------------
// Background options
//----------------------------------------------------------------------------------------------------------------------------------------------------


namespace BackgroundOptions {

//! Derivative method (generally either analytic or numeric)
enum class DerivativeMethod {
   analytic,
   numeric
};


//! Heliospheric current sheet (0: disabled, 1: flat, 2: wavy (Jokipii-Thomas 1981) and static, 3: wavy and time-dependent).
enum class CurrentSheet {
   disabled,
   flat,
   wavy_static,
   wavy_time_dependent
};

//! Magnetic topology region (0: nowhere, 1: same as HCS)
enum class SectoredRegion {
   nowhere,
   HCS
};

//! Correction to Parker Spiral, mainly for polar regions (0: none, 1: Smith-Bieber 1991, 2: Zurbuchen et al. 1997, 3: Schwadron-McComas 2003)
enum class PolarCorrection {
   none,
   Smith_Bieber,
   Zurbuchen_etal,
   Schwadron_McComas
};

//! Latitudinal profile for bulk speed (0: constant, 1: linear step, 2: smooth step)
enum class SpeedLatitudeProfile {
   constant,
   linear_step,
   smooth_step
};

//! Integer exponent of decrease of solar wind speed beyond the termination shock
enum class TermShockSpeedExponent {
   zero,
   one,
   square,
   cube
};

//! What function to use within 'get_ampfactor' (0 = none, 1 = zero, 2 = constant, 3 = scaled)
enum class ModType {
   none,
   zero,
   constant,
   scaled
};

//! Whether to scale relative to s=0 (0) or s=+inf (1)
enum class ModRPos {
   scale_rel_zero,
   scale_rel_inf
};

}


//----------------------------------------------------------------------------------------------------------------------------------------------------
// Trajectory options
//----------------------------------------------------------------------------------------------------------------------------------------------------


enum class TrajectoryId {
   Fieldline,
   Focused,
   Guiding,
   Lorentz,
   Parker,
};


namespace TrajectoryOptions {

//! Direction of trajectory integration (time flow direction)
enum class TimeFlow {
   forward,
   backward
};

//! How many time steps to allow before recording a mirror event
//const int mirror_threshold;

enum class SafetyLevel {
   low,
   medium,
   high
};

// Switch controlling how to calculate mu. "0" means computing it at the end of the step from magnetic moment conservation. "1" means advancing it in time according to the scheme (does not guarantee conservation of MM, but can be used with non-adiabatic terms).
enum class PPerpMethod {
   mag_moment_conservation,
   scheme
};

//! Flag to use gradient and curvature drifts in drift velocity calculation
enum class UseBDrifts {
   none,
   gradient_curvature
};

//! Which stochastic method to use for PA scattering, 0 = Euler, 1 = Milstein, 2 = RK2
enum class StochasticMethod {
   Euler,
   Milstein,
   RK2
};

//! Whether to use constant dmumax or constant dthetamax, 0 = constant dthetamax, 1 = constant dmumax
enum class ConstDmumax {
   constant_dtheta_max,
   constant_dmu_max
};

//! Which method of computation to use for divK: 0 = using direct central FD, 1 = using _spdata.grad quantities
enum class DivkMethod {
   direct,
   gradients
};

}


//----------------------------------------------------------------------------------------------------------------------------------------------------
// Diffusion options
//----------------------------------------------------------------------------------------------------------------------------------------------------


namespace DiffusionOptions {

enum Component {
   perp = 0,
   para = 1,
   mu = 2,
};

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


//! Type indicator to use a default configuration for a Config class
class Default {};






};

#endif
