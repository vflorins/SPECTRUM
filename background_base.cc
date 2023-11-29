/*!
\file background_base.cc
\brief Implements a base class to compute the plasma background
\author Vladimir Florinski

This file is part of the SPECTRUM suite of scientific numerical simulation codes. SPECTRUM stands for Space Plasma and Energetic Charged particle TRansport on Unstructured Meshes. The code simulates plasma or neutral particle flows using MHD equations on a grid, transport of cosmic rays using stochastic or grid based methods. The "unstructured" part refers to the use of a geodesic mesh providing a uniform coverage of the surface of a sphere.
*/

#include "background_base.hh"

namespace Spectrum {

//----------------------------------------------------------------------------------------------------------------------------------------------------
// BackgroundBase methods
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\author Vladimir Florinski
\date 11/24/2020
*/
BackgroundBase::BackgroundBase(void)
              : Params("", 0, STATE_NONE)
{
};

/*!
\author Vladimir Florinski
\date 12/17/2020
\param[in] name_in   Readable name of the class
\param[in] specie_in Particle's specie
\param[in] status_in Initial status
*/
BackgroundBase::BackgroundBase(const std::string& name_in, unsigned int specie_in, uint16_t status_in)
              : Params(name_in, specie_in, status_in)
{
};

/*!
\author Vladimir Florinski
\date 09/26/2021
\param[in] other Object to initialize from

A copy constructor should first first call the Params' version to copy the data container and then check whether the other object has been set up. If yes, it should simply call the virtual method "SetupBackground()" with the argument of "true".
*/
BackgroundBase::BackgroundBase(const BackgroundBase& other)
              : Params(other)
{
// Params' constructor resets all flags
   if(BITS_RAISED(other._status, STATE_SETUP_COMPLETE)) SetupBackground(true);
};

/*!
\author Vladimir Florinski
\date 12/14/2020
\return Maximum distance based on the grid or other properties
*/
double BackgroundBase::GetDmax(void) const
{
   return _spdata.dmax;
};

/*!
\author Vladimir Florinski
\author Juan G Alonso Guzman
\date 07/26/2022
\param[in] int index for which derivative to take (0 = x, 1 = y, 2 = z, else = t)
\note This is a common routine that the derived classes should not change.
*/
void BackgroundBase::DirectionalDerivative(int xyz)
{
   double _t_saved;
   GeoVector _pos_saved;

   if(BITS_LOWERED(_status, STATE_SETUP_COMPLETE)) {
      RAISE_BITS(_status, STATE_INVALID);
      throw ExUninitialized();
   };

// Spatial derivatives
   if(0 <= xyz && xyz <= 2) {

// Save position, compute increment      
      _pos_saved = _pos;
      _spdata_tmp._dr[xyz] = incr_dmax_ratio * _spdata.dmax;
      _pos += _spdata_tmp._dr[xyz] * cart_unit_vec[xyz];
      EvaluateBackground();

// Error in a space increment - try negative increment
      if(BITS_RAISED(_status, STATE_INVALID)) {
         _spdata_tmp._dr[xyz] = -_spdata_tmp._dr[xyz];
         _pos += 2.0 * _spdata_tmp._dr[xyz] * cart_unit_vec[xyz];
         EvaluateBackground();

// Still error in space increment - throw error
         if(BITS_RAISED(_status, STATE_INVALID)) throw ExFieldError();
      };

// Restore position
      _pos = _pos_saved;

      if(BITS_RAISED(_spdata._mask, BACKGROUND_U)) _spdata_tmp.gradUvec[xyz] = (_spdata.Uvec - _spdata_tmp.Uvec) / _spdata_tmp._dr[xyz];
      if(BITS_RAISED(_spdata._mask, BACKGROUND_B)) _spdata_tmp.gradBvec[xyz] = (_spdata.Bvec - _spdata_tmp.Bvec) / _spdata_tmp._dr[xyz];
      if(BITS_RAISED(_spdata._mask, BACKGROUND_E)) _spdata_tmp.gradEvec[xyz] = (_spdata.Evec - _spdata_tmp.Evec) / _spdata_tmp._dr[xyz];
   }

// Time derivatives
   else {

// Save time, compute increment      
      _t_saved = _t;
      _spdata_tmp._dt = incr_dmax_ratio * _spdata.dmax / c_code;
      _t += _spdata_tmp._dt;
      EvaluateBackground();

// Error in a time increment - try negative increment
      if(BITS_RAISED(_status, STATE_INVALID)) {
         _spdata_tmp._dt = -_spdata_tmp._dt;
         _t += 2.0 * _spdata_tmp._dt;
         EvaluateBackground();

// Still error in time increment - throw error
         if(BITS_RAISED(_status, STATE_INVALID)) throw ExFieldError();
      };
// Restore time
      _t = _t_saved;

      if(BITS_RAISED(_spdata._mask, BACKGROUND_U)) _spdata_tmp.dUvecdt = (_spdata.Uvec - _spdata_tmp.Uvec) / _spdata_tmp._dt;
      if(BITS_RAISED(_spdata._mask, BACKGROUND_B)) _spdata_tmp.dBvecdt = (_spdata.Bvec - _spdata_tmp.Bvec) / _spdata_tmp._dt;
      if(BITS_RAISED(_spdata._mask, BACKGROUND_E)) _spdata_tmp.dEvecdt = (_spdata.Evec - _spdata_tmp.Evec) / _spdata_tmp._dt;
   };
};

/*!
\author Vladimir Florinski
\author Juan G Alonso Guzman
\date 10/13/2022
*/
void BackgroundBase::NumericalDerivatives(void)
{
   int xyz;

// Save the mask, u,B,E. This is quicker than using the copy assignment based on _mask
// Note: any background computing a gradXvec or dXvecdt should also have the flag for computing X itself active, otherwise the derivative will be wrong
   _spdata_tmp._mask = _spdata._mask;
   if(BITS_RAISED(_spdata._mask, BACKGROUND_U)) _spdata_tmp.Uvec = _spdata.Uvec;
   if(BITS_RAISED(_spdata._mask, BACKGROUND_B)) _spdata_tmp.Bvec = _spdata.Bvec;
   if(BITS_RAISED(_spdata._mask, BACKGROUND_E)) _spdata_tmp.Evec = _spdata.Evec;
   _spdata_tmp.region = _spdata.region;
   _spdata_tmp.dmax = _spdata.dmax;

// Spatial derivatives. The mask shifting is done to limit the evaluation of the variable to those that require a gradient
   _spdata._mask >>= mask_offset;
   if(BITS_RAISED(_spdata._mask, BACKGROUND_ALL)) {
      for(xyz = 0; xyz < 3; xyz++) DirectionalDerivative(xyz);
   };

// The calculation in DirectionalDerivative() gives gradV[i][j] = dV_j / ds^i. It needs to be transposed to match the standard mathematical notation for the Jacobian.
   if(BITS_RAISED(_spdata._mask, BACKGROUND_U)) _spdata.gradUvec.Transpose();
   if(BITS_RAISED(_spdata._mask, BACKGROUND_B)) _spdata.gradBvec.Transpose();
   if(BITS_RAISED(_spdata._mask, BACKGROUND_E)) _spdata.gradEvec.Transpose();

// Time derivatives. The mask shifting is done to limit the evaluation of the variable to those that require a time derivative.
   _spdata._mask >>= mask_offset;
   if(BITS_RAISED(_spdata._mask, BACKGROUND_ALL)) {
      DirectionalDerivative(3);
   };

// Copy data from _spdata_tmp into _spdata
   _spdata._mask = _spdata_tmp._mask;
   _spdata = _spdata_tmp;
};

/*!
\author Vladimir Florinski
\date 11/25/2020
\param[in] cont_in Container with parameters

This is the default method to set up an object. It should only be defined in the base class (XXXXBase). Derived classes should _not_ modify it! This version always calls the correct virtual "SetupBackground()" method.
*/
void BackgroundBase::SetupObject(const DataContainer& cont_in)
{
   Params::SetContainer(cont_in);
   SetupBackground(false);
};

/*!
\author Vladimir Florinski
\date 09/26/2021
\param [in] construct Whether called from a copy constructor or separately

This method's main role is to unpack the data container and set up the class data members and status bits marked as "persistent". The function should assume that the data container is available because the calling function will always ensure this.
*/
void BackgroundBase::SetupBackground(bool construct)
{
// Only needed in the parent version
   container.Reset();
   RAISE_BITS(_status, STATE_SETUP_COMPLETE);
   LOWER_BITS(_status, STATE_INVALID);

// Unpack parameters
   container.Read(r0.Data());
   container.Read(u0.Data());
   container.Read(B0.Data());
   container.Read(&dmax0);

// Initialize "safe" box for derivatives
   for(int xyz = 0; xyz < 3; xyz++) _spdata._dr[xyz] = incr_dmax_ratio * dmax0;
   _spdata._dt = incr_dmax_ratio * dmax0 / c_code;
   _spdata.dmax = dmax0;
};

/*!
\author Vladimir Florinski
\date 12/08/2021
\note This is only a stub
*/
void BackgroundBase::EvaluateBackground(void)
{
   LOWER_BITS(_status, STATE_INVALID);
};

/*!
\author Vladimir Florinski
\date 10/13/2022
\note This is only a stub
*/
void BackgroundBase::EvaluateBackgroundDerivatives(void)
{
};

/*!
\author Vladimir Florinski
\date 09/15/2021
\note The default method should be good enough for all grid-free backgrounds
*/
void BackgroundBase::EvaluateDmax(void)
{
   _spdata.dmax = dmax0;
   LOWER_BITS(_status, STATE_INVALID);
};

/*!
\author Vladimir Florinski
\date 02/17/2023
\note The default method should be good enough for all grid-free backgrounds
*/
void BackgroundBase::StopServerFront(void)
{
};

/*!
\author Juan G Alonso Guzman
\date 10/19/2022
\param[in] dir Direction
\return Safe increment in some direction (potential negative) to stay inside domain
*/
double BackgroundBase::GetSafeIncr(const GeoVector& dir)
{
   return incr_dmax_ratio * _spdata.dmax;
};

/*!
\author Vladimir Florinski
\author Juan G Alonso Guzman
\date 02/17/2023
\param[in]  t_in   Time
\param[in]  pos_in Position
\param[out] spdata All spatial data
\note This is a common routine that the derived classes should not change.
*/
void BackgroundBase::GetFields(double t_in, const GeoVector& pos_in, SpatialData& spdata)
{
// Check that state setup is complete
   if(BITS_LOWERED(_status, STATE_SETUP_COMPLETE)) {
      RAISE_BITS(_status, STATE_INVALID);
      throw ExUninitialized();
   };

// If "EvaluateDmax()" fails, the state will be set to "STATE_INVALID" and background will not be evaluated
   SetState(t_in, pos_in);
   EvaluateDmax();
   if(BITS_RAISED(_status, STATE_INVALID)) throw ExCoordinates();

// The mask is provided by the caller, we need to copy it into our internal fields structure
   _spdata._mask = spdata._mask;

   EvaluateBackground();
   if(BITS_RAISED(_status, STATE_INVALID)) throw ExFieldError();
   EvaluateBackgroundDerivatives();
   if(BITS_RAISED(_status, STATE_INVALID)) throw ExFieldError();

// Copy the internal fields into arguments and compute derived values for B
   spdata = _spdata;
   if(BITS_RAISED(spdata._mask, BACKGROUND_B)) {
      spdata.Bmag = spdata.Bvec.Norm();
      if(spdata.Bmag < tiny) RAISE_BITS(_status, STATE_INVALID);
      spdata.bhat = spdata.Bvec / spdata.Bmag;
   };
};

};
