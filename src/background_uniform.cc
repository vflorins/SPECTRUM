/*!
\file background_uniform.cc
\brief Implements a simple uniform field background
\author Vladimir Florinski
\author Juan G Alonso Guzman

This file is part of the SPECTRUM suite of scientific numerical simulation codes. SPECTRUM stands for Space Plasma and Energetic Charged particle TRansport on Unstructured Meshes. The code simulates plasma or neutral particle flows using MHD equations on a grid, transport of cosmic rays using stochastic or grid based methods. The "unstructured" part refers to the use of a geodesic mesh providing a uniform coverage of the surface of a sphere.
*/

#include "background_uniform.hh"

namespace Spectrum {

//----------------------------------------------------------------------------------------------------------------------------------------------------
// BackgroundUniform methods
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\author Vladimir Florinski
\date 09/27/2021
*/
template <typename Fields>
BackgroundUniform<Fields>::BackgroundUniform(void)
                 : BackgroundBase(bg_name_uniform, 0, MODEL_STATIC)
{
};

/*!
\author Vladimir Florinski
\date 09/26/2021
\param[in] other Object to initialize from

A copy constructor should first first call the Params' version to copy the data container and then check whether the other object has been set up. If yes, it should simply call the virtual method "SetupBackground()" with the argument of "true".
*/
template <typename Fields>
BackgroundUniform<Fields>::BackgroundUniform(const BackgroundUniform& other)
                 : BackgroundBase(other)
{
   RAISE_BITS(_status, MODEL_STATIC);
   if (BITS_RAISED(other._status, STATE_SETUP_COMPLETE)) SetupBackground(true);
};

/*!
\author Vladimir Florinski
\author Juan G Alonso Guzman
\date 01/04/2024
\param [in] construct Whether called from a copy constructor or separately

This method's main role is to unpack the data container and set up the class data members and status bits marked as "persistent". The function should assume that the data container is available because the calling function will always ensure this.
*/
template <typename Fields>
void BackgroundUniform<Fields>::SetupBackground(bool construct)
{
// The parent version must be called explicitly if not constructing
   if (!construct) BackgroundBase::SetupBackground(false);

// Precompute motional electric field for efficiency
   E0 = -(u0 ^ B0) / c_code;
};

/*!
\author Vladimir Florinski
\author Juan G Alonso Guzman
\date 01/04/2024
*/
template <typename Fields>
void BackgroundUniform<Fields>::EvaluateBackground(void)
{
   if constexpr (Fields::Vel_found()) {
      _fields.Vel() = u0;
   }
   if constexpr (Fields::Mag_found()) {
      _fields.Mag() = B0;
   }
   if constexpr (Fields::Elc_found()) {
      _fields.Elc() = E0;
   }
   if constexpr (Fields::Iv0_found()) {
      _fields.Iv0() = 1.0;
   }
   LOWER_BITS(_status, STATE_INVALID);
};

/*!
\author Juan G Alonso Guzman
\author Vladimir Florinski
\date 10/14/2022
*/
template <typename Fields>
void BackgroundUniform<Fields>::EvaluateBackgroundDerivatives(void)
{
// Spatial derivatives are zero
   if constexpr (Fields::DelVel_found()) _fields.DelVel() = gm_zeros;
   if constexpr (Fields::DelMag_found()) _fields.DelMag() = gm_zeros;
   if constexpr (Fields::DelElc_found()) _fields.DelElc() = gm_zeros;

// Time derivatives are zero
   if constexpr (Fields::DdtVel_found()) _fields.DdtVel() = gv_zeros;
   if constexpr (Fields::DdtMag_found()) _fields.DdtMag() = gv_zeros;
   if constexpr (Fields::DdtElc_found()) _fields.DdtElc() = gv_zeros;
};

};
