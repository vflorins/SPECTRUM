/*!
\file source_other.cc
\brief Implements several classes to compute source terms
\author Juan G Alonso Guzman

This file is part of the SPECTRUM suite of scientific numerical simulation codes. SPECTRUM stands for Space Plasma and Energetic Charged particle TRansport on Unstructured Meshes. The code simulates plasma or neutral particle flows using MHD equations on a grid, transport of cosmic rays using stochastic or grid based methods. The "unstructured" part refers to the use of a geodesic mesh providing a uniform coverage of the surface of a sphere.
*/

#include "source_other.hh"

namespace Spectrum {

//----------------------------------------------------------------------------------------------------------------------------------------------------
// SourceConstant methods
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\author Juan G Alonso Guzman
\date 07/14/2025
*/
SourceConstant::SourceConstant(void)
              : SourceBase(source_name_constant, 0, MODEL_STATIC)
{
};

/*!
\author Juan G Alonso Guzman
\date 07/14/2025
\param[in] other Object to initialize from

A copy constructor should first first call the Params' version to copy the data container and then check whether the other object has been set up. If yes, it should simply call the virtual method "SetupSource()" with the argument of "true".
*/
SourceConstant::SourceConstant(const SourceConstant& other)
              : SourceBase(other)
{
   if (BITS_RAISED(other._status, STATE_SETUP_COMPLETE)) SetupSource(true);
};

/*!
\author Juan G Alonso Guzman
\date 07/14/2025
\param [in] construct Whether called from a copy constructor or separately

This method's main role is to unpack the data container and set up the class data members and status bits marked as "persistent". The function should assume that the data container is available because the calling function will always ensure this.
*/
void SourceConstant::SetupSource(bool construct)
{
// The parent version must be called explicitly if not constructing
   if (!construct) SourceBase::SetupSource(false);
   container.Read(S0);

// Preset source term since it's constant
   SourceTerm = S0;
};

/*!
\author Juan G Alonso Guzman
\date 07/14/2025
*/
void SourceConstant::EvaluateSource(void)
{
// Nothing to do
};

//----------------------------------------------------------------------------------------------------------------------------------------------------
// SourceMomentumInjection methods
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\author Juan G Alonso Guzman
\date 09/09/2025
*/
SourceMomentumInjection::SourceMomentumInjection(void)
                       : SourceBase(source_name_momentum_injection, 0, MODEL_STATIC)
{
};

/*!
\author Juan G Alonso Guzman
\date 09/09/2025
\param[in] other Object to initialize from

A copy constructor should first first call the Params' version to copy the data container and then check whether the other object has been set up. If yes, it should simply call the virtual method "SetupSource()" with the argument of "true".
*/
SourceMomentumInjection::SourceMomentumInjection(const SourceMomentumInjection& other)
                       : SourceBase(other)
{
   if (BITS_RAISED(other._status, STATE_SETUP_COMPLETE)) SetupSource(true);
};

/*!
\author Juan G Alonso Guzman
\date 09/09/2025
\param [in] construct Whether called from a copy constructor or separately

This method's main role is to unpack the data container and set up the class data members and status bits marked as "persistent". The function should assume that the data container is available because the calling function will always ensure this.
*/
void SourceMomentumInjection::SetupSource(bool construct)
{
// The parent version must be called explicitly if not constructing
   if (!construct) SourceBase::SetupSource(false);
   container.Read(p_inj);
   container.Read(rate);

// Initialize momentum deltas to zero
   del_mom = del_mom_old = 0.0;
};

/*!
\author Juan G Alonso Guzman
\date 09/09/2025
*/
void SourceMomentumInjection::EvaluateSource(void)
{
   del_mom = _mom[0] - p_inj;
   if (del_mom * del_mom_old < 0.0) SourceTerm = rate / _dt;
   else SourceTerm = 0.0;
   del_mom_old = del_mom;
};

//----------------------------------------------------------------------------------------------------------------------------------------------------
// SourceSphericalShockInjection methods
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\author Juan G Alonso Guzman
\date 07/14/2025
*/
SourceSphericalShockInjection::SourceSphericalShockInjection(void)
                             : SourceBase(source_name_spherical_shock_injection, 0, MODEL_STATIC)
{
};

/*!
\author Juan G Alonso Guzman
\date 07/14/2025
\param[in] other Object to initialize from

A copy constructor should first first call the Params' version to copy the data container and then check whether the other object has been set up. If yes, it should simply call the virtual method "SetupSource()" with the argument of "true".
*/
SourceSphericalShockInjection::SourceSphericalShockInjection(const SourceSphericalShockInjection& other)
                             : SourceBase(other)
{
   if (BITS_RAISED(other._status, STATE_SETUP_COMPLETE)) SetupSource(true);
};

/*!
\author Juan G Alonso Guzman
\date 07/14/2025
\param [in] construct Whether called from a copy constructor or separately

This method's main role is to unpack the data container and set up the class data members and status bits marked as "persistent". The function should assume that the data container is available because the calling function will always ensure this.
*/
void SourceSphericalShockInjection::SetupSource(bool construct)
{
// The parent version must be called explicitly if not constructing
   if (!construct) SourceBase::SetupSource(false);
   container.Read(p_inj);
   container.Read(r0);
   container.Read(r_sh);
   container.Read(w_sh);
   container.Read(rate);

// Initialize momentum deltas to zero
   del_mom = del_mom_old = 0.0;
};

/*!
\author Juan G Alonso Guzman
\date 07/14/2025
*/
void SourceSphericalShockInjection::EvaluateSource(void)
{
   double r = (_pos - r0).Norm();
   del_mom = _mom[0] - p_inj;
   if (del_mom * del_mom_old < 0.0) {
      if (r_sh <= r && r <= r_sh + w_sh) SourceTerm = rate / _dt;
      else SourceTerm = 0.0;
   }
   else SourceTerm = 0.0;
   del_mom_old = del_mom;
};

};
