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
template <typename Fields>
BackgroundSolarWindTermShock<Fields>::BackgroundSolarWindTermShock(void)
                            : BackgroundSolarWind(bg_name_solarwind_termshock, 0, MODEL_STATIC)
{
};

/*!
\author Juan G Alonso Guzman
\date 02/22/2023
\param[in] other Object to initialize from

A copy constructor should first first call the Params' version to copy the data container and then check whether the other object has been set up. If yes, it should simply call the virtual method "SetupBackground()" with the argument of "true".
*/
template <typename Fields>
BackgroundSolarWindTermShock<Fields>::BackgroundSolarWindTermShock(const BackgroundSolarWindTermShock& other)
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
template <typename Fields>
void BackgroundSolarWindTermShock<Fields>::SetupBackground(bool construct)
{
// The parent version must be called explicitly if not constructing
   if (!construct) BackgroundSolarWind::SetupBackground(false);
   container.Read(r_TS);
   container.Read(w_TS);
   container.Read(s_TS);

   s_TS_inv = 1.0 / s_TS;
   dmax_TS = dmax_fraction * w_TS;
};

/*!
\author Juan G Alonso Guzman
\date 03/14/2024
\param[in]  r      radial distance
\param[out] ur_mod modified radial flow
*/
template <typename Fields>
void BackgroundSolarWindTermShock<Fields>::ModifyUr(const double r, double &ur_mod)
{
   if (r > r_TS) {
#if SOLARWIND_TERMSHOCK_SPEED_EXPONENT == 1
      if (r > r_TS + w_TS) ur_mod *= s_TS_inv * (r_TS + w_TS) / r;
#elif SOLARWIND_TERMSHOCK_SPEED_EXPONENT == 2
      if (r > r_TS + w_TS) ur_mod *= s_TS_inv * Sqr((r_TS + w_TS) / r);
#else
      if (r > r_TS + w_TS) ur_mod *= s_TS_inv;
#endif
      else ur_mod *= 1.0 + (s_TS_inv - 1.0) * (r - r_TS) / w_TS;
   };
};

/*!
\author Juan G Alonso Guzman
\date 05/14/2025
\param[in]  r      radial distance
*/
template <typename Fields>
double BackgroundSolarWindTermShock<Fields>::dUrdr(const double r)
{
   if (r > r_TS) {
#if SOLARWIND_TERMSHOCK_SPEED_EXPONENT == 1
      if (r > r_TS + w_TS) return -_fields.Vel().Norm() / r;
#elif SOLARWIND_TERMSHOCK_SPEED_EXPONENT == 2
      if (r > r_TS + w_TS) return -2.0 * _fields.Vel().Norm() / r;
#else
      if (r > r_TS + w_TS) return 0.0;
#endif
      else return (s_TS_inv - 1.0) * (r_TS / w_TS) * ur0;
   }
   else return 0.0;
};

/*!
\author Juan G Alonso Guzman
\date 06/21/2024
\param[in]  r radial distance
\param[out] time lag of propagation from solar surface to current position
*/
template <typename Fields>
double BackgroundSolarWindTermShock<Fields>::TimeLag(const double r)
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
\date 05/14/2023
*/
template <typename Fields>
void BackgroundSolarWindTermShock<Fields>::EvaluateBackgroundDerivatives(void)
{
#if SOLARWIND_DERIVATIVE_METHOD == 0
   double r;
   GeoVector posprime;
   GeoMatrix rr;

   if constexpr (Fields::DelVel_found()) {
// Expression valid only for radial flow
      posprime = _pos - r0;
      r = posprime.Norm();
      rr.Dyadic(posprime / r);
      _fields.DelVel() = dUrdr(r) * rr + (_fields.Vel().Norm() / r) * (gm_unit - rr);
   };
   if constexpr (Fields::DelMag_found()) {
//TODO: complete
   };
   if constexpr (Fields::DelElc_found()) {
      _fields.DelElc() = -((_fields.DelVel() ^ _fields.Mag()) + (_fields.Vel() ^ _fields.DelMag())) / c_code;
   };
   if constexpr (Fields::DdtVel_found()) _fields.DdtVel() = gv_zeros;
   if constexpr (Fields::DdtMag_found()) _fields.DdtMag() = gv_zeros;
   if constexpr (Fields::DdtElc_found()) _fields.DdtElc() = gv_zeros;

#else
   NumericalDerivatives();
#endif
};

/*!
\author Juan G Alonso Guzman
\date 02/23/2024
*/
template <typename Fields>
void BackgroundSolarWindTermShock<Fields>::EvaluateDmax(void)
{
   BackgroundSolarWind::EvaluateDmax();

// Reduce "dmax" around the shock. This implemenation assumes that "dmax" = "dmax0" near "r_TS" by default.
   double r = (_pos - r0).Norm();
   if (r_TS - dmax0 < r && r < r_TS + w_TS + dmax0) {
      if (r < r_TS) _ddata.dmax += (dmax_TS - dmax0) * (r - r_TS + dmax0) / dmax0;
      else if (r > r_TS + w_TS) _ddata.dmax -= (dmax_TS - dmax0) * (r - r_TS - w_TS - dmax0) / dmax0;
      else _ddata.dmax = dmax_TS;
   };
};

};
