/*!
\file background_dipole.cc
\brief Implements a dipole magnetic field background without a flow
\author Vladimir Florinski
\author Juan G Alonso Guzman

This file is part of the SPECTRUM suite of scientific numerical simulation codes. SPECTRUM stands for Space Plasma and Energetic Charged particle TRansport on Unstructured Meshes. The code simulates plasma or neutral particle flows using MHD equations on a grid, transport of cosmic rays using stochastic or grid based methods. The "unstructured" part refers to the use of a geodesic mesh providing a uniform coverage of the surface of a sphere.
*/

#include "background_dipole.hh"

namespace Spectrum {

//----------------------------------------------------------------------------------------------------------------------------------------------------
// BackgroundDipole methods
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\author Vladimir Florinski
\date 03/25/2022
*/
template <typename HyperParams>
BackgroundDipole<HyperParams>::BackgroundDipole(void)
                : BackgroundBase(bg_name, MODEL_STATIC)
{
};

/*!
\author Vladimir Florinski
\date 03/25/2022
\param[in] other Object to initialize from

A copy constructor should first first call the Params' version to copy the data container and then check whether the other object has been set up. If yes, it should simply call the virtual method "SetupBackground()" with the argument of "true".
*/
template <typename HyperParams>
BackgroundDipole<HyperParams>::BackgroundDipole(const BackgroundDipole& other)
                : BackgroundBase(other)
{
   RAISE_BITS(_status, MODEL_STATIC);
   if (BITS_RAISED(other._status, STATE_SETUP_COMPLETE)) SetupBackground(true);
};

/*!
\author Vladimir Florinski
\date 08/22/2023
\param [in] construct Whether called from a copy constructor or separately

This method's main role is to unpack the data container and set up the class data members and status bits marked as "persistent". The function should assume that the data container is available because the calling function will always ensure this.
*/
template <typename HyperParams>
void BackgroundDipole<HyperParams>::SetupBackground(bool construct)
{
   double r_ref;

// The parent version must be called explicitly if not constructing
   if (!construct) BackgroundBase::SetupBackground(false);
   container.Read(r_ref);
   container.Read(dmax_fraction);
   M = B0 * Cube(r_ref);
};

/*!
\author Vladimir Florinski
\date 03/25/2022
*/
template <typename HyperParams>
void BackgroundDipole<HyperParams>::EvaluateBackground(void)
{
   if constexpr (Fields::Vel_found()) {
      _fields.Vel() = gv_zeros;
   }
   if constexpr (Fields::Mag_found()) {
      GeoVector posprime = _pos - r0;
      double r2 = posprime.Norm();
      double r5 = Cube(r2);
      r2 *= r2;
      r5 *= r2;
      _fields.Mag() = (3.0 * (posprime * M) * posprime - r2 * M) / r5;
   }
   // todo bhat, Bmag?
   if constexpr (Fields::Elc_found()) {
      _fields.Elc() = gv_zeros;
   }
   if constexpr (Fields::Iv0_found()) {
      // todo originally region is 3 variables. Is this one variable or two or three?
      _fields.Iv0() = 1.0;
   }

   LOWER_BITS(_status, STATE_INVALID);
};

/*!
\author Vladimir Florinski
\author Juan G Alonso Guzman
\date 10/14/2022
*/
template <typename HyperParams>
void BackgroundDipole<HyperParams>::EvaluateBackgroundDerivatives(void)
{
   if constexpr (HyperParams::derivative_method == DerivativeMethod::analytic) {
      if constexpr (Fields::DelVel_found())
         _fields.DelVel() = gm_zeros;
      if constexpr (Fields::DelMag_found() || Fields::DelAbsMag_found()) {
         GeoVector posprime = _pos - r0;
         double r2 = posprime.Norm();
         double r5 = Cube(r2);
         r2 *= r2;
         r5 *= r2;
         GeoMatrix mr, rm, rr;
         mr.Dyadic(M,posprime);
         rm.Dyadic(posprime,M);
         rr.Dyadic(posprime);

         auto DelMag = 3.0 * (mr + rm + (M * posprime) * (gm_unit - 5.0 * rr / r2)) / r5;

         if constexpr (Fields::DelMag_found())
            _fields.DelMag() = DelMag;
         if constexpr (Fields::DelAbsMag_found()) {
            // todo how to finesse this?
            auto bhat = _fields.Mag() / _fields.Mag().Norm();
            _fields.DelAbsMag() = DelMag * bhat;
         }
      };
      if constexpr (Fields::DelElc_found())
         _fields.DelElc() = gm_zeros;
      if constexpr (Fields::DdtVel_found())
         _fields.DdtVel() = gv_zeros;
      if constexpr (Fields::DdtMag_found())
         _fields.DdtMag() = gv_zeros;
      if constexpr (Fields::DdtAbsMag_found())
         _fields.DdtAbsMag() = 0.0;
      if constexpr (Fields::DdtElc_found())
         _fields.DdtElc() = gv_zeros;
   }
   else {
      NumericalDerivatives();
   };
};

/*!
\author Vladimir Florinski
\date 03/25/2022
*/
template <typename HyperParams>
void BackgroundDipole<HyperParams>::EvaluateDmax(void)
{
   _ddata.dmax = fmin(dmax_fraction * (_pos - r0).Norm(), dmax0);
   LOWER_BITS(_status, STATE_INVALID);
};

};
