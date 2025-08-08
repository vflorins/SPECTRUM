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
template <typename Fields>
BackgroundSmoothShock<Fields>::BackgroundSmoothShock(void)
                     : BackgroundShock(bg_name_smooth_shock, 0, STATE_NONE)
{
};

/*!
\author Juan G Alonso Guzman
\date 10/20/2023
\param[in] other Object to initialize from

A copy constructor should first first call the Params' version to copy the data container and then check whether the other object has been set up. If yes, it should simply call the virtual method "SetupBackground()" with the argument of "true".
*/
template <typename Fields>
BackgroundSmoothShock<Fields>::BackgroundSmoothShock(const BackgroundSmoothShock& other)
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
template <typename Fields>
double BackgroundSmoothShock<Fields>::ShockTransition(double x)
{
   double y = x + 0.5;
#if SMOOTH_SHOCK_ORDER == 0 // continous but not differentiable
   if (x < -0.5) return 0.0;
   else if (x > 0.5) return 1.0;
   else return y;
#elif SMOOTH_SHOCK_ORDER == 1 // differentiable
   if (x < -0.5) return 0.0;
   else if (x > 0.5) return 1.0;
   else return Sqr(y) * (3.0 - 2.0 * y);
#elif SMOOTH_SHOCK_ORDER == 2 // twice differentiable
   if (x < -0.5) return 0.0;
   else if (x > 0.5) return 1.0;
   else return Cube(y) * (10.0 - 15.0 * y + 6.0 * Sqr(y));
#elif SMOOTH_SHOCK_ORDER == 3 // thrice differentiable
   if (x < -0.5) return 0.0;
   else if (x > 0.5) return 1.0;
   else return Sqr(Sqr(y)) * (35.0 - 84.0 * y + 70.0 * Sqr(y) - 20.0 * Cube(y));
#else // smooth
   return 0.5 * (1.0 + tanh(tanh_width_factor * x));
#endif
};

/*!
\author Juan G Alonso Guzman
\date 05/30/2024
\param [in] x Relative transition region location 
\return Derivative of relative value of shocked quantity
*/
template <typename Fields>
double BackgroundSmoothShock<Fields>::ShockTransitionDerivative(double x)
{
   double y = x + 0.5;
#if SMOOTH_SHOCK_ORDER == 0
   if (x < -0.5) return 0.0;
   else if (x > 0.5) return 0.0;
   else return 1.0;
#elif SMOOTH_SHOCK_ORDER == 1
   if (x < -0.5) return 0.0;
   else if (x > 0.5) return 0.0;
   else return 6.0 * y * (1.0 - y);
#elif SMOOTH_SHOCK_ORDER == 2
   if (x < -0.5) return 0.0;
   else if (x > 0.5) return 0.0;
   else return 30.0 * Sqr(y) * (1.0 - 2.0 * y + Sqr(y));
#elif SMOOTH_SHOCK_ORDER == 3
   if (x < -0.5) return 0.0;
   else if (x > 0.5) return 0.0;
   else return 140.0 * Cube(y) * (1.0 - 3.0 * y + 3.0 * Sqr(y) - 1.0 * Cube(y));
#else
   return 0.5 * tanh_width_factor * (1.0 - Sqr(tanh(tanh_width_factor * x)));
#endif
};

/*!
\author Juan G Alonso Guzman
\date 02/28/2025
\param [in] construct Whether called from a copy constructor or separately

This method's main role is to unpack the data container and set up the class data members and status bits marked as "persistent". The function should assume that the data container is available because the calling function will always ensure this.
*/
template <typename Fields>
void BackgroundSmoothShock<Fields>::SetupBackground(bool construct)
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
template <typename Fields>
void BackgroundSmoothShock<Fields>::EvaluateBackground(void)
{
   double a1, a2;
   ds_shock = ((_pos - r0) * n_shock - v_shock * _t) / width_shock;

   a1 = ShockTransition(ds_shock);
   a2 = 1.0 - a1;

// TODO: Change the way Bvec transitions to depend directly on Uvec to keep some product constant
   if constexpr (Fields::Vel_found()) _fields.Vel() = u0 * a1 + u1 * a2;
   if constexpr (Fields::Mag_found()) _fields.Mag() = B0 * a1 + B1 * a2;
   if constexpr (Fields::Elc_found()) _fields.Elc() = -(_fields.Vel() ^ _fields.Mag()) / c_code;
   if constexpr (Fields::Iv0_found()) _fields.Ev0() = 1.0 * a1 + 2.0 * a2; // same as 2.0 - a1

   LOWER_BITS(_status, STATE_INVALID);
};

/*!
\author Juan G Alonso Guzman
\date 10/20/2023
*/
template <typename Fields>
void BackgroundSmoothShock<Fields>::EvaluateBackgroundDerivatives(void)
{
#if SMOOTHSHOCK_DERIVATIVE_METHOD == 0
   if constexpr (Fields::DelVel_found()) {
      _fields.gradUvec.Dyadic(n_shock, u0 - u1);
      _fields.gradUvec *= ShockTransitionDerivative(ds_shock) / width_shock;
   };
   if constexpr (Fields::DelMag_found()) {
      _fields.DelMag().Dyadic(n_shock, B0 - B1);
      _fields.DelMag() *= ShockTransitionDerivative(ds_shock) / width_shock;
      GeoVector bhat;
      if constexpr (Fields::HatMag_found()) bhat = _fields.HatMag();
      else bhat = _fields.Mag() / _fields.Mag().Norm();
      _fields.DelMag() = _fields.DelMag() * bhat;
   };
   if constexpr (Fields::DelElc_found()) {
      _fields.gradEvec = -((_fields.DelVel() ^ _fields.Mag()) + (_fields.Vel() ^ _fields.DelMag())) / c_code;
   };
   if constexpr (Fields::DdtVel_found()) {
      _fields.dUvecdt = (ShockTransitionDerivative(ds_shock) * v_shock / width_shock) * (u1 - u0);
   };
   if constexpr (Fields::DdtMag_found()) {
      _fields.dBvecdt = (ShockTransitionDerivative(ds_shock) * v_shock / width_shock) * (B1 - B0);
   };
   if constexpr (Fields::DdtElc_found()) {
      _fields.dEvecdt = -((_fields.DdtVel() ^ _fields.Mag()) + (_fields.Vel() ^ _fields.DdtMag())) / c_code;
   };
#else
   NumericalDerivatives();
#endif
};

/*!
\author Juan G Alonso Guzman
\date 02/28/2025
*/
template <typename Fields>
void BackgroundSmoothShock<Fields>::EvaluateDmax(void)
{
   _ddata.dmax = fmin(dmax_fraction * width_shock * fmax(1.0, fabs(ds_shock)), dmax0);
   LOWER_BITS(_status, STATE_INVALID);
};

};
