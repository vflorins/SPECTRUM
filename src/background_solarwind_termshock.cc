/*!
\file background_solarwind.hh
\brief Implements a plasma background class for the constant speed supersonic wind of a rotating star
\author Vladimir Florinski
\author Juan G Alonso Guzman

This file is part of the SPECTRUM suite of scientific numerical simulation codes. SPECTRUM stands for Space Plasma and Energetic Charged particle TRansport on Unstructured Meshes. The code simulates plasma or neutral particle flows using MHD equations on a grid, transport of cosmic rays using stochastic or grid based methods. The "unstructured" part refers to the use of a geodesic mesh providing a uniform coverage of the surface of a sphere.
*/

#include "background_solarwind_termshock.hh"

namespace Spectrum {

//----------------------------------------------------------------------------------------------------------------------------------------------------
// BackgroundSolarWind methods
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\author Juan G Alonso Guzman
\date 02/22/2023
*/
BackgroundSolarWindTermShock::BackgroundSolarWindTermShock(void)
                            : BackgroundSolarWind(bg_name_solarwind_termshock, 0, MODEL_STATIC)
{
};

/*!
\author Juan G Alonso Guzman
\date 02/22/2023
\param[in] other Object to initialize from

A copy constructor should first first call the Params' version to copy the data container and then check whether the other object has been set up. If yes, it should simply call the virtual method "SetupBackground()" with the argument of "true".
*/
BackgroundSolarWindTermShock::BackgroundSolarWindTermShock(const BackgroundSolarWindTermShock& other)
                            : BackgroundSolarWind(other)
{
   RAISE_BITS(_status, MODEL_STATIC);
   if (BITS_RAISED(other._status, STATE_SETUP_COMPLETE)) SetupBackground(true);
};

/*!
\author Juan G Alonso Guzman
\date 02/22/2023
\param [in] construct Whether called from a copy constructor or separately

This method's main role is to unpack the data container and set up the class data members and status bits marked as "persistent". The function should assume that the data container is available because the calling function will always ensure this.
*/
void BackgroundSolarWindTermShock::SetupBackground(bool construct)
{
// The parent version must be called explicitly if not constructing
   if (!construct) BackgroundSolarWind::SetupBackground(false);
   container.Read(r_TS);
   container.Read(w_TS);
   container.Read(s_TS);

   s_TS_inv = 1.0 / s_TS;
};

/*!
\author Juan G Alonso Guzman
\date 05/19/2025
\param [in] r radial distance
\return Relative value of radial speed
*/
double BackgroundSolarWindTermShock::TermShockTransition(double r)
{
   double y = (r - r_TS) / w_TS;
#if SMOOTH_TERM_SHOCK_ORDER == 0 // continous but not differentiable
   if (y < 0.0) return 0.0;
   else if (y > 1.0) return 1.0;
   else return y;
#elif SMOOTH_TERM_SHOCK_ORDER == 1 // differentiable
   if (y < 0.0) return 0.0;
   else if (y > 1.0) return 1.0;
   else return Sqr(y) * (3.0 - 2.0 * y);
#elif SMOOTH_TERM_SHOCK_ORDER == 2 // twice differentiable
   if (y < 0.0) return 0.0;
   else if (y > 1.0) return 1.0;
   else return Cube(y) * (10.0 - 15.0 * y + 6.0 * Sqr(y));
#elif SMOOTH_TERM_SHOCK_ORDER == 3 // thrice differentiable
   if (y < 0.0) return 0.0;
   else if (y > 1.0) return 1.0;
   else return Sqr(Sqr(y)) * (35.0 - 84.0 * y + 70.0 * Sqr(y) - 20.0 * Cube(y));
#else // smooth
   return 0.5 * (1.0 + tanh(tanh_width_factor * (y - 0.5)));
#endif
};

/*!
\author Juan G Alonso Guzman
\date 05/19/2025
\param [in] r radial distance
\return Derivative of relative value of radial speed
*/
double BackgroundSolarWindTermShock::TermShockTransitionDerivative(double r)
{
   double y = (r - r_TS) / w_TS;
#if SMOOTH_TERM_SHOCK_ORDER == 0
   if (y < 0.0) return 0.0;
   else if (y > 1.0) return 0.0;
   else return 1.0;
#elif SMOOTH_TERM_SHOCK_ORDER == 1
   if (y < 0.0) return 0.0;
   else if (y > 1.0) return 0.0;
   else return 6.0 * y * (1.0 - y);
#elif SMOOTH_TERM_SHOCK_ORDER == 2
   if (y < 0.0) return 0.0;
   else if (y > 1.0) return 0.0;
   else return 30.0 * Sqr(y) * (1.0 - 2.0 * y + Sqr(y));
#elif SMOOTH_TERM_SHOCK_ORDER == 3
   if (y < 0.0) return 0.0;
   else if (y > 1.0) return 0.0;
   else return 140.0 * Cube(y) * (1.0 - 3.0 * y + 3.0 * Sqr(y) - 1.0 * Cube(y));
#else
   return 0.5 * tanh_width_factor * (1.0 - Sqr(tanh(tanh_width_factor * (y - 0.5))));
#endif
};

/*!
\author Juan G Alonso Guzman
\date 05/19/2025
\param[in]  r      radial distance
\param[out] ur_mod modified radial flow
*/
void BackgroundSolarWindTermShock::ModifyUr(const double r, double &ur_mod)
{
   if (r >= r_TS) {
#if SOLARWIND_TERMSHOCK_SPEED_EXPONENT == 1
      if (r > r_TS + w_TS) ur_mod *= s_TS_inv * (r_TS + w_TS) / r;
#elif SOLARWIND_TERMSHOCK_SPEED_EXPONENT == 2
      if (r > r_TS + w_TS) ur_mod *= s_TS_inv * Sqr((r_TS + w_TS) / r);
#else
      if (r > r_TS + w_TS) ur_mod *= s_TS_inv;
#endif
      else ur_mod *= 1.0 + (s_TS_inv - 1.0) * TermShockTransition(r);
   };
};

/*!
\author Juan G Alonso Guzman
\author Swati Sharma
\date 05/19/2025
\param[in]  r      radial distance
*/
double BackgroundSolarWindTermShock::dUrdr(const double r)
{
   if (r >= r_TS) {
#if SOLARWIND_TERMSHOCK_SPEED_EXPONENT == 1
      if (r > r_TS + w_TS) return -_spdata.Uvec.Norm() / r;
#elif SOLARWIND_TERMSHOCK_SPEED_EXPONENT == 2
      if (r > r_TS + w_TS) return -2.0 * _spdata.Uvec.Norm() / r;
#else
      if (r > r_TS + w_TS) return 0.0;
#endif
      else return ur0 * (s_TS_inv - 1.0) * TermShockTransitionDerivative(r) / w_TS;
   }
   else return 0.0;
};

/*!
\author Juan G Alonso Guzman
\date 06/21/2024
\param[in]  r radial distance
\param[out] time lag of propagation from solar surface to current position
*/
double BackgroundSolarWindTermShock::TimeLag(const double r)
{
   if (r < r_TS) return r / ur0;
#if SOLARWIND_TERMSHOCK_SPEED_EXPONENT == 1
   else return (r_TS + s_TS * (Sqr(r) - Sqr(r_TS)) / (2.0 * r_TS)) / ur0;
#elif SOLARWIND_TERMSHOCK_SPEED_EXPONENT == 2
   else return (r_TS + s_TS * (Cube(r) - Cube(r_TS)) / (3.0 * Sqr(r_TS))) / ur0;
#else
   else return (r_TS + s_TS * (r - r_TS)) / ur0;
#endif
};

/*!
\author Juan G Alonso Guzman
\author Swati Sharma
\date 04/30/2025
*/
void BackgroundSolarWindTermShock::EvaluateBackgroundDerivatives(void)
{
#if SOLARWIND_DERIVATIVE_METHOD == 0
   double r;
   GeoMatrix CartesianToSpherical;

   r = (_pos - r0).Norm();
   CartesianToSpherical.Transpose(SphericalToCartesian);

   if (BITS_RAISED(_spdata._mask, BACKGROUND_gradU)) {
// Expression valid only for radial flow
      _spdata.gradUvec = gm_zeros;
      _spdata.gradUvec[0][0] = dUrdr(r);
      _spdata.gradUvec[1][1] = _spdata.Uvec.Norm() / r;
      _spdata.gradUvec[2][2] = _spdata.Uvec.Norm() / r;
// Convert to Cartesian
      _spdata.gradUvec = SphericalToCartesian * _spdata.gradUvec * CartesianToSpherical;
      _spdata.gradUvec.ChangeFromBasis(eprime);
   };
   if (BITS_RAISED(_spdata._mask, BACKGROUND_gradB)) {
// Expressions valid only for radial field ~ 1/r^2
      _spdata.gradBvec = gm_zeros;
      _spdata.gradBvec[0][0] = -2.0 * _spdata.Bvec.Norm() / r;
      _spdata.gradBvec[1][1] = _spdata.Bvec.Norm() / r;
      _spdata.gradBvec[2][2] = _spdata.Bvec.Norm() / r;
      _spdata.gradBmag = gv_zeros;
      _spdata.gradBmag[0] = -2.0 * _spdata.Bvec.Norm() / r;
// Convert to Cartesian
      _spdata.gradBvec = SphericalToCartesian * _spdata.gradBvec * CartesianToSpherical;
      _spdata.gradBvec.ChangeFromBasis(eprime);
      _spdata.gradBmag = SphericalToCartesian * _spdata.gradBmag;
      _spdata.gradBmag.ChangeFromBasis(eprime);
   };
   if (BITS_RAISED(_spdata._mask, BACKGROUND_gradE)) {
      _spdata.gradEvec = -((_spdata.gradUvec ^ _spdata.Bvec) + (_spdata.Uvec ^ _spdata.gradBvec)) / c_code;
   };
   if (BITS_RAISED(_spdata._mask, BACKGROUND_dUdt)) _spdata.dUvecdt = gv_zeros;
   if (BITS_RAISED(_spdata._mask, BACKGROUND_dBdt)) _spdata.dBvecdt = gv_zeros;
   if (BITS_RAISED(_spdata._mask, BACKGROUND_dEdt)) _spdata.dEvecdt = gv_zeros;

#else
   NumericalDerivatives();
#endif
};

/*!
\author Juan G Alonso Guzman
\date 02/23/2024
*/
void BackgroundSolarWindTermShock::EvaluateDmax(void)
{
   BackgroundSolarWind::EvaluateDmax();

// Reduce "dmax" around the shock. This implemenation assumes that "dmax" = "dmax0" near "r_TS" by default.
   double dr_shock = ((_pos - r0).Norm() - r_TS - 0.5 * w_TS) / w_TS;
   _spdata.dmax = fmin(dmax_fraction * w_TS * fmax(1.0, fabs(dr_shock)), _spdata.dmax);
};

};
