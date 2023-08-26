/*!
\file background_dipole.cc
\brief Implements a dipole magnetic field background without a flow
\author Vladimir Florinski

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
BackgroundDipole::BackgroundDipole(void)
                : BackgroundBase(bg_name_dipole, 0, MODEL_STATIC)
{
};

/*!
\author Vladimir Florinski
\date 03/25/2022
\param[in] other Object to initialize from

A copy constructor should first first call the Params' version to copy the data container and then check whether the other object has been set up. If yes, it should simply call the virtual method "SetupBackground()" with the argument of "true".
*/
BackgroundDipole::BackgroundDipole(const BackgroundDipole& other)
                : BackgroundBase(other)
{
   RAISE_BITS(_status, MODEL_STATIC);
   if(BITS_RAISED(other._status, STATE_SETUP_COMPLETE)) SetupBackground(true);
};

/*!
\author Vladimir Florinski
\date 03/25/2022
\param [in] construct Whether called from a copy constructor or separately

This method's main role is to unpack the data container and set up the class data members and status bits marked as "persistent". The function should assume that the data container is available because the calling function will always ensure this.
*/
void BackgroundDipole::SetupBackground(bool construct)
{
   double r_ref;

// The parent version must be called explicitly if not constructing
   if(!construct) BackgroundBase::SetupBackground(false);
   container.Read(&r_ref);
   container.Read(&dmax_fraction);
   M = B0 * Cube(r_ref) / 2.0;
};

/*!
\author Vladimir Florinski
\date 03/25/2022
*/
void BackgroundDipole::EvaluateBackground(void)
{
   if(BITS_RAISED(_spdata._mask, BACKGROUND_U)) _spdata.Uvec = gv_zeros;
   if(BITS_RAISED(_spdata._mask, BACKGROUND_B)) {
      GeoVector posprime = _pos - r0;
      double r2 = posprime.Norm();
      double r5 = Cube(r2);
      r2 *= r2;
      r5 *= r2;
      _spdata.Bvec = (3.0 * (posprime * M) * posprime - r2 * M) / r5;
   };
   if(BITS_RAISED(_spdata._mask, BACKGROUND_E)) _spdata.Evec = gv_zeros;
   _spdata.region = 1.0;

   LOWER_BITS(_status, STATE_INVALID);
};

/*!
\author Vladimir Florinski
\author Juan G Alonso Guzman
\date 10/14/2022
*/
void BackgroundDipole::EvaluateBackgroundDerivatives(void)
{
#if DIPOLE_DERIVATIVE_METHOD == 0

   if(BITS_RAISED(_spdata._mask, BACKGROUND_gradU)) _spdata.gradUvec = gm_zeros;
   if(BITS_RAISED(_spdata._mask, BACKGROUND_gradB)) {
      GeoVector posprime = _pos - r0;
      double r2 = posprime.Norm();
      double r5 = Cube(r2);
      r2 *= r2;
      r5 *= r2;
      double mdotr = M * posprime;
      GeoMatrix mr, rm, rr;
      mr.Dyadic(M,posprime);
      rm.Dyadic(posprime,M);
      rr.Dyadic(posprime);

// TODO change the second call to Dyadic to Transpose
      _spdata.gradBvec = 3.0 * (mr + rm + mdotr * (gm_unit - 5.0 * rr / r2)) / r5;
   };
   if(BITS_RAISED(_spdata._mask, BACKGROUND_gradE)) _spdata.gradEvec = gm_zeros;
   if(BITS_RAISED(_spdata._mask, BACKGROUND_dUdt)) _spdata.dUvecdt = gv_zeros;
   if(BITS_RAISED(_spdata._mask, BACKGROUND_dBdt)) _spdata.dBvecdt = gv_zeros;
   if(BITS_RAISED(_spdata._mask, BACKGROUND_dEdt)) _spdata.dEvecdt = gv_zeros;

#else
   NumericalDerivatives(); 
#endif
};

/*!
\author Vladimir Florinski
\date 03/25/2022
*/
void BackgroundDipole::EvaluateDmax(void)
{
   _spdata.dmax = fmin(dmax_fraction * (_pos - r0).Norm(), dmax0);
   LOWER_BITS(_status, STATE_INVALID);
};

};
