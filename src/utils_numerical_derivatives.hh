/*!
\file utils_numerical_derivatives.hh
\brief Declares a class to compute the plasma background gradients numerically
\author Vladimir Florinski
\author Juan G Alonso Guzman
\author Lucius Schoenbaum

This file is part of the SPECTRUM suite of scientific numerical simulation codes. SPECTRUM stands for Space Plasma and Energetic Charged particle TRansport on Unstructured Meshes. The code simulates plasma or neutral particle flows using MHD equations on a grid, transport of cosmic rays using stochastic or grid based methods. The "unstructured" part refers to the use of a geodesic mesh providing a uniform coverage of the surface of a sphere.
*/

#ifndef SPECTRUM_UTILS_NUMERICAL_DERIVATIVES_HH
#define SPECTRUM_UTILS_NUMERICAL_DERIVATIVES_HH


// This includes (algorithm, cmath, cstdint, cstring, exception, fstream, vector), data_container, definitions, multi_index, vectors
#include "common/params.hh"
#include "common/physics.hh"
#include "common/matrix.hh"
#include "common/derivativedata.hh"

#include <memory>
#ifdef USE_SILO
#include <silo.h>
#endif

namespace Spectrum {

//----------------------------------------------------------------------------------------------------------------------------------------------------
// NumericalDerivatives class declaration
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*
 * TODO: instantiate trivially if !method==numeric.
 *   In that case, return a trivially constructed ddata, but make sure it has a valid dmax field.
 *
 *
 */

/*!
\brief A stateful helper class for computing derivatives numerically and managing data.
\author Vladimir Florinski
\author Juan G Alonso Guzman
\author Lucius Schoenbaum

Unlike background classes, the numerical derivatives class is stateful.
It must maintain its state in the event of a failure or exception in numerical computations,
an event which must propagate through the entire timestep.
NumericalDerivatives
*/
template <typename Background_>
class NumericalDerivatives {
public:

   using Background = Background_;
   using HConfig = Background::HConfig;
   using Config = HConfig::BackgroundConfig;

private:

//! What fraction of "_dmax" to use to calculate the field increment
   static constexpr double incr_dmax_ratio = Config::incr_dmax_ratio;

//! Rotation angle of local (x,y) plane in numerical derivative evaluation/averaging
   static constexpr double local_rot_ang = M_PI / Config::num_numeric_grad_evals;

//! Sine and cosine of local rotation angle
   static constexpr double sin_lra = sin(local_rot_ang);
   static constexpr double cos_lra = cos(local_rot_ang);

protected:

   // TODO: *******
   //  After refactoring, check whether/when ddata gets
   //  reset after each trajectory advance.
   //  *************

//! Distance reference value, how far the trajectory is allowed to move in one step (persistent)
   static constexpr double dmax0 = Config::dmax0;

//! Derivative data object (transient)
   DerivativeData _ddata;

//! Field-aligned basis (transient)
   GeoMatrix fa_basis;

//! Rotation matrix (transient)
   GeoMatrix rot_mat;

//! Set up the field evaluator
   void Setup(bool construct);

//! Compute the fields at an incremented position or time
   template <typename Coordinates, typename Fields, typename RequestedFields>
   void DirectionalDerivative(int xyz, Coordinates, Fields&, double scale_factor);

//! Compute the field derivatives numerically and statefully
   template <typename Coordinates, typename Fields, typename RequestedFields>
   void EvaluateBackgroundDerivatives(Coordinates&, Fields&);

public:

//! Destructor
   ~NumericalDerivatives() = default;

//! Signal the backend this client no longer needs its service
   void StopServerFront(void);

//! Find "safe" increment in a given direction
   double GetSafeIncr(const GeoVector& dir);

//! Return the derivative data, a diagnostic data structure
   [[nodiscard]]
   DerivativeData GetDerivativeData(void) const;

};

};

// Something like this is needed for templated classes
#include "utils_numerical_derivatives.cc"

#endif
