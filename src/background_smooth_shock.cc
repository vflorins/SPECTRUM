/*!
\file background_smooth_shock.cc
\brief Implements a smooth shock field background
\author Juan G Alonso Guzman

This file is part of the SPECTRUM suite of scientific numerical simulation codes. SPECTRUM stands for Space Plasma and Energetic Charged particle TRansport on Unstructured Meshes. The code simulates plasma or neutral particle flows using MHD equations on a grid, transport of cosmic rays using stochastic or grid based methods. The "unstructured" part refers to the use of a geodesic mesh providing a uniform coverage of the surface of a sphere.
*/

#include "background_smooth_shock.hh"

namespace Spectrum {

//----------------------------------------------------------------------------------------------------------------------------------------------------
// BackgroundSmoothShock methods
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\author Juan G Alonso Guzman
\date 10/20/2023
*/
template <typename HConfig>
BackgroundSmoothShock<HConfig>::BackgroundSmoothShock(void)
                     : BackgroundShock(bg_name, STATE_NONE)
{
};

/*!
\author Juan G Alonso Guzman
\date 10/20/2023
\param[in] other Object to initialize from

A copy constructor should first first call the Params' version to copy the data container and then check whether the other object has been set up. If yes, it should simply call the virtual method "SetupBackground()" with the argument of "true".
*/
template <typename HConfig>
BackgroundSmoothShock<HConfig>::BackgroundSmoothShock(const BackgroundSmoothShock& other)
                     : BackgroundShock(other)
{
   RAISE_BITS(_status, STATE_NONE);
   if (BITS_RAISED(other._status, STATE_SETUP_COMPLETE)) SetupBackground(true);
};

/*!
\author Juan G Alonso Guzman
\date 05/30/2024
\param [in] x Relative transition region location
\return Relative value of shocked quantity
*/
template <typename HConfig>
double BackgroundSmoothShock<HConfig>::ShockTransition(double x)
{
   double y = x + 0.5;
   if constexpr (HConfig::smooth_discontinuity_order == 0) {
      // not continous
      if (x < -0.5) return 0.0;
      else if (x > 0.5) return 1.0;
      else return y;
   }
   else if constexpr (HConfig::smooth_discontinuity_order == 1){
      // differentiable
      if (x < -0.5) return 0.0;
      else if (x > 0.5) return 1.0;
      else return Sqr(y) * (3.0 - 2.0 * y);
   }
   else if constexpr (HConfig::smooth_discontinuity_order == 2) {
      // twice differentiable
      if (x < -0.5) return 0.0;
      else if (x > 0.5) return 1.0;
      else return Cube(y) * (10.0 - 15.0 * y + 6.0 * Sqr(y));
   }
   else if constexpr (HConfig::smooth_discontinuity_order == 3) {
      // thrice differentiable
      if (x < -0.5) return 0.0;
      else if (x > 0.5) return 1.0;
      else return Sqr(Sqr(y)) * (35.0 - 84.0 * y + 70.0 * Sqr(y) - 20.0 * Cube(y));
   }
   else {
      // smooth
      return 0.5 * (1.0 + tanh(tanh_width_factor * x));
   }
};

/*!
\author Juan G Alonso Guzman
\date 05/30/2024
\param [in] x Relative transition region location 
\return Derivative of relative value of shocked quantity
*/
template <typename HConfig>
double BackgroundSmoothShock<HConfig>::ShockTransitionDerivative(double x)
{
   double y = x + 0.5;
   if constexpr (HConfig::smooth_discontinuity_order == 0) {
      if (x < -0.5) return 0.0;
      if (x > 0.5) return 0.0;
      return 1.0;
   }
   else if constexpr (HConfig::smooth_discontinuity_order == 1){
      if (x < -0.5) return 0.0;
      if (x > 0.5) return 0.0;
      return 6.0 * y * (1.0 - y);
   }
   else if constexpr (HConfig::smooth_discontinuity_order == 2) {
      if (x < -0.5) return 0.0;
      if (x > 0.5) return 0.0;
      return 30.0 * Sqr(y) * (1.0 - 2.0 * y + Sqr(y));
   }
   else if constexpr (HConfig::smooth_discontinuity_order == 3) {
      if (x < -0.5) return 0.0;
      if (x > 0.5) return 0.0;
      return 140.0 * Cube(y) * (1.0 - 3.0 * y + 3.0 * Sqr(y) - 1.0 * Cube(y));
   }
   else {
      return 0.5 * tanh_width_factor * (1.0 - Sqr(tanh(tanh_width_factor * x)));
   }
};

/*!
\author Juan G Alonso Guzman
\date 02/28/2025
\param [in] construct Whether called from a copy constructor or separately

This method's main role is to unpack the data container and set up the class data members and status bits marked as "persistent". The function should assume that the data container is available because the calling function will always ensure this.
*/
template <typename HConfig>
void BackgroundSmoothShock<HConfig>::SetupBackground(bool construct)
{
// The parent version must be called explicitly if not constructing
   if (!construct) BackgroundShock::SetupBackground(false);

// Unpack parameters
   container.Read(width_shock);
   container.Read(dmax_fraction);
};

/*!
\author Juan G Alonso Guzman
\date 05/14/2025
*/
template <typename HConfig>
template <typename Fields>
void BackgroundSmoothShock<HConfig>::EvaluateBackground(BackgroundCoordinates& coords, Fields& fields)
{
   double a1, a2;
   ds_shock = ((coords.Pos() - r0) * n_shock - v_shock * coords.Time()) / width_shock;

   a1 = ShockTransition(ds_shock);
   a2 = 1.0 - a1;

// TODO: Change the way Bvec transitions to depend directly on Uvec to keep some product constant
   if constexpr (Fields::Vel_found()) fields.Vel() = u0 * a1 + u1 * a2;
   if constexpr (Fields::Mag_found()) fields.Mag() = B0 * a1 + B1 * a2;
   if constexpr (Fields::Elc_found()) fields.Elc() = -(fields.Vel() ^ fields.Mag()) / c_code;
   if constexpr (Fields::Iv0_found()) fields.Ev0() = 1.0 * a1 + 2.0 * a2; // same as 2.0 - a1

   LOWER_BITS(_status, STATE_INVALID);
};

/*!
\author Juan G Alonso Guzman
\date 10/20/2023
*/
template <typename HConfig>
template <typename Fields>
void BackgroundSmoothShock<HConfig>::EvaluateBackgroundDerivatives(BackgroundCoordinates& coords, Fields& fields)
{
   if constexpr (HConfig::derivative_method == DerivativeMethod::analytic) {
      if constexpr (Fields::DelVel_found()) {
         fields.gradUvec.Dyadic(n_shock, u0 - u1);
         fields.gradUvec *= ShockTransitionDerivative(ds_shock) / width_shock;
      };
      if constexpr (Fields::DelMag_found()) {
         fields.DelMag().Dyadic(n_shock, B0 - B1);
         fields.DelMag() *= ShockTransitionDerivative(ds_shock) / width_shock;
         GeoVector bhat;
         if constexpr (Fields::HatMag_found()) bhat = fields.HatMag();
         else bhat = fields.Mag() / fields.Mag().Norm();
         fields.DelMag() = fields.DelMag() * bhat;
      };
      if constexpr (Fields::DelElc_found()) {
         fields.gradEvec = -((fields.DelVel() ^ fields.Mag()) + (fields.Vel() ^ fields.DelMag())) / c_code;
      };
      if constexpr (Fields::DotVel_found()) {
         fields.dUvecdt = (ShockTransitionDerivative(ds_shock) * v_shock / width_shock) * (u1 - u0);
      };
      if constexpr (Fields::DotMag_found()) {
         fields.dBvecdt = (ShockTransitionDerivative(ds_shock) * v_shock / width_shock) * (B1 - B0);
      };
      if constexpr (Fields::DotElc_found()) {
         fields.dEvecdt = -((fields.DotVel() ^ fields.Mag()) + (fields.Vel() ^ fields.DotMag())) / c_code;
      };
   }
   else {
      NumericalDerivatives(coords, fields);
   };
};

/*!
\author Juan G Alonso Guzman
\date 02/28/2025
*/
template <typename Fields>
void BackgroundSmoothShock<Fields>::EvaluateDmax(BackgroundCoordinates& coords)
{
   _ddata.dmax = fmin(dmax_fraction * width_shock * fmax(1.0, fabs(ds_shock)), dmax0);
   LOWER_BITS(_status, STATE_INVALID);
};

};
