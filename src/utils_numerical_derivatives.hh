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
#include "background.hh"

#include <memory>
#ifdef USE_SILO
#include <silo.h>
#endif

namespace Spectrum {

//----------------------------------------------------------------------------------------------------------------------------------------------------
// NumericalDerivatives class declaration
//----------------------------------------------------------------------------------------------------------------------------------------------------



template <typename Config>
struct NumericalDerivativesBase {

//! Distance reference value, how far the trajectory is allowed to move in one step (persistent)
   static constexpr double dmax0 = Config::dmax0;

//! What fraction of "_dmax" to use to calculate the field increment
   static constexpr double incr_dmax_ratio = Config::incr_dmax_ratio;

   //! Derivative data object (transient)
   DerivativeData _ddata;

/*!
\author Lucius Schoenbaum
\date 11/25/2025
Reset (re-initialize) the derivative data structure.
 */
   void Reset(double dmax)
   {
      // Initialize "safe" box for derivatives
      _ddata = DerivativeData();
      for (auto xyz = 0; xyz < 3; xyz++) _ddata._dr[xyz] = incr_dmax_ratio * dmax0;
      _ddata._dt = incr_dmax_ratio * dmax0 / c_code;
      _ddata.dmax = dmax;
   }

/*!
\author Juan G Alonso Guzman
\date 10/19/2022
\param[in] dir Direction
\return Safe increment in some direction (potential negative) to stay inside domain
*/
   [[nodiscard]] double GetSafeIncr(const GeoVector& dir)
   {
//FIXME: This is incomplete.
      return incr_dmax_ratio * _ddata.dmax;
   };

/*!
\author Lucius Schoenbaum
\date 08/15/2025
\return DerivativeData, a small struct containing information about the most recent derivative computed, and dmax
\note This information from the background is currently only needed by Diffusion classes
when numerical directional derivatives are computed.
 */
   [[nodiscard]] DerivativeData GetDerivativeData(void) const
   {
      return _ddata;
   }

};



template <typename Background, bool need_numeric>
class NumericalDerivatives;



template <typename Background>
struct NumericalDerivatives<Background, false>: public NumericalDerivativesBase<typename Background::HConfig::BackgroundConfig> {};



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
template <typename Background>
class NumericalDerivatives<Background, true>: public NumericalDerivativesBase<typename Background::HConfig_::BackgroundConfig> {
public:

   using HConfig = Background::HConfig;
   using Config = HConfig::BackgroundConfig;

   using Base = NumericalDerivativesBase<Config>;
   using Base::GetDerivativeData;
   using Base::GetSafeIncr;
   using Base::Reset;
   using Base::dmax0;
   using Base::_ddata;
   using Base::incr_dmax_ratio;

   static constexpr auto specie = HConfig::specie;

private:

//! Rotation angle of local (x,y) plane in numerical derivative evaluation/averaging
   static constexpr double local_rot_ang = M_PI / Config::num_numeric_grad_evals;

//! Sine and cosine of local rotation angle
   static constexpr double sin_lra = sin(local_rot_ang);
   static constexpr double cos_lra = cos(local_rot_ang);

protected:

//! Field-aligned basis (transient)
   GeoMatrix fa_basis;

//! Rotation matrix (transient)
   GeoMatrix rot_mat;

//! Compute the fields at an incremented position or time
   template <typename Coordinates, typename Fields, typename RequestedFields>
   void DirectionalDerivative(int xyz, Coordinates, Fields&, double scale_factor, Background& background);

public:

//! Constructor
   NumericalDerivatives() = default;

//! Destructor
   ~NumericalDerivatives() = default;

//! Compute the field derivatives numerically and statefully
   template <typename Coordinates, typename Fields, typename RequestedFields>
   status_t EvaluateBackgroundDerivatives(Coordinates&, Fields&, Background& background);

};

};

// Something like this is needed for templated classes
#include "utils_numerical_derivatives.cc"

#endif
