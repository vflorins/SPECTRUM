/*!
\file background_solarwind.hh
\brief Implements a plasma background class for the constant speed supersonic wind of a rotating star
\author Vladimir Florinski
\author Juan G Alonso Guzman

This file is part of the SPECTRUM suite of scientific numerical simulation codes. SPECTRUM stands for Space Plasma and Energetic Charged particle TRansport on Unstructured Meshes. The code simulates plasma or neutral particle flows using MHD equations on a grid, transport of cosmic rays using stochastic or grid based methods. The "unstructured" part refers to the use of a geodesic mesh providing a uniform coverage of the surface of a sphere.
*/

#include "background_solarwind.hh"

namespace Spectrum {

//----------------------------------------------------------------------------------------------------------------------------------------------------
// BackgroundSolarWind methods
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\author Vladimir Florinski
\author Juan G Alonso Guzman
\date 10/28/2022
*/
BackgroundSolarWind::BackgroundSolarWind(void)
                   : BackgroundBase(bg_name_solarwind, 0, MODEL_STATIC)
{
};

/*!
\author Juan G Alonso Guzman
\date 02/22/2024
*/
BackgroundSolarWind::BackgroundSolarWind(const std::string& name_in, unsigned int specie_in, uint16_t status_in)
                   : BackgroundBase(name_in, specie_in, status_in)
{
};

/*!
\author Vladimir Florinski
\date 01/26/2022
\param[in] other Object to initialize from

A copy constructor should first first call the Params' version to copy the data container and then check whether the other object has been set up. If yes, it should simply call the virtual method "SetupBackground()" with the argument of "true".
*/
BackgroundSolarWind::BackgroundSolarWind(const BackgroundSolarWind& other)
                   : BackgroundBase(other)
{
   RAISE_BITS(_status, MODEL_STATIC);
   if(BITS_RAISED(other._status, STATE_SETUP_COMPLETE)) SetupBackground(true);
};

/*!
\author Vladimir Florinski
\date 01/26/2022
\param [in] construct Whether called from a copy constructor or separately

This method's main role is to unpack the data container and set up the class data members and status bits marked as "persistent". The function should assume that the data container is available because the calling function will always ensure this.
*/
void BackgroundSolarWind::SetupBackground(bool construct)
{
// The parent version must be called explicitly if not constructing
   if(!construct) BackgroundBase::SetupBackground(false);
   container.Read(Omega.Data());
   container.Read(&r_ref);
   container.Read(&dmax_fraction);

// Build the new coordinate system. The z axis is along "omega"
   eprime[2] = UnitVec(Omega);
   eprime[0] = GetSecondUnitVec(eprime[2]);
   eprime[1] = eprime[2] ^ eprime[0];

// Only the first components of velocity is used. The value for B could be negative (for negative cycles).
   ur0 = fabs(u0[0]);
   Br0 = B0[0];
   w0 = Omega.Norm();
};

/*!
\author Vladimir Florinski
\date 01/27/2022
*/
void BackgroundSolarWind::EvaluateBackground(void)
{
// Convert position into flow aligned coordinates and scale to "r_ref"
   posprime = _pos - r0;
   posprime.ChangeToBasis(eprime);
   posprime /= r_ref;

// Position in cylindrical coordinates
   double r, r3;
   r = posprime.Norm();

// Compute the (radial) velocity and convert back to global frame
   if(BITS_RAISED(_spdata._mask, BACKGROUND_U)) {
      _spdata.Uvec = (ur0 / r) * posprime;
      _spdata.Uvec.ChangeFromBasis(eprime);
   };
   
// Compute (Parker spiral) magnetic field and convert back to global frame
   if(BITS_RAISED(_spdata._mask, BACKGROUND_B)) {
      r3 = Cube(r);
      _spdata.Bvec = (Br0 / r3) * (posprime + (w0 * (r - 1.0) * r_ref / ur0) * (posprime[1] * gv_nx - posprime[0] * gv_ny));
      _spdata.Bvec.ChangeFromBasis(eprime);
   };

// Compute electric field, already in global frame. Note that the flags to compute U and B should be enabled in order to compute E.
   if(BITS_RAISED(_spdata._mask, BACKGROUND_E)) _spdata.Evec = -(_spdata.Uvec ^ _spdata.Bvec) / c_code;

   LOWER_BITS(_status, STATE_INVALID);
};

/*!
\author Juan G Alonso Guzman
\author Vladimir Florinski
\date 10/14/2022
*/
void BackgroundSolarWind::EvaluateBackgroundDerivatives(void)
{
#if SOLARWIND_DERIVATIVE_METHOD == 0

   if(BITS_RAISED(_spdata._mask, BACKGROUND_gradU)) {
//TODO: complete
   };
   if(BITS_RAISED(_spdata._mask, BACKGROUND_gradB)) {
//TODO: complete
   };
   if(BITS_RAISED(_spdata._mask, BACKGROUND_gradE)) {
      _spdata.gradEvec = -((_spdata.gradUvec ^ _spdata.Bvec) + (_spdata.Uvec ^ _spdata.gradBvec)) / c_code;
   };
   if(BITS_RAISED(_spdata._mask, BACKGROUND_dUdt)) _spdata.dUvecdt = gv_zeros;
   if(BITS_RAISED(_spdata._mask, BACKGROUND_dBdt)) _spdata.dBvecdt = gv_zeros;
   if(BITS_RAISED(_spdata._mask, BACKGROUND_dEdt)) _spdata.dEvecdt = gv_zeros;

#else
   NumericalDerivatives();
#endif
};

/*!
\author Vladimir Florinski
\date 03/10/2022
*/
void BackgroundSolarWind::EvaluateDmax(void)
{
   _spdata.dmax = fmin(dmax_fraction * (_pos - r0).Norm(), dmax0);
   LOWER_BITS(_status, STATE_INVALID);
};

};
