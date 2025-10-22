/*!
\file background_uniform.cc
\brief Implements a simple uniform field background
\author Vladimir Florinski
\author Juan G Alonso Guzman

This file is part of the SPECTRUM suite of scientific numerical simulation codes. SPECTRUM stands for Space Plasma and Energetic Charged particle TRansport on Unstructured Meshes. The code simulates plasma or neutral particle flows using MHD equations on a grid, transport of cosmic rays using stochastic or grid based methods. The "unstructured" part refers to the use of a geodesic mesh providing a uniform coverage of the surface of a sphere.
*/

#include "background_uniform.hh"

namespace Spectrum {

using namespace BackgroundOptions;

//----------------------------------------------------------------------------------------------------------------------------------------------------
// BackgroundUniform methods
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\author Vladimir Florinski
\date 09/27/2021
*/
template <typename HConfig>
BackgroundUniform<HConfig>::BackgroundUniform(void)
                 : BackgroundBase(bg_name, MODEL_STATIC)
{
};

/*!
\author Vladimir Florinski
\date 09/26/2021
\param[in] other Object to initialize from

A copy constructor should first first call the Params' version to copy the data container and then check whether the other object has been set up. If yes, it should simply call the virtual method "SetupBackground()" with the argument of "true".
*/
template <typename HConfig>
BackgroundUniform<HConfig>::BackgroundUniform(const BackgroundUniform& other)
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
template <typename HConfig>
void BackgroundUniform<HConfig>::SetupBackground(bool construct)
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
template <typename HConfig>
template <typename Coordinates, typename Fields, typename RequestedFields>
void BackgroundUniform<HConfig>::EvaluateBackground(Coordinates& coords, Fields& fields)
{
   if constexpr (RequestedFields::Fluv_found()) {
      fields.Fluv() = u0;
   }
   if constexpr (RequestedFields::Mag_found()) {
      fields.Mag() = B0;
   }
   if constexpr (RequestedFields::Elc_found()) {
      fields.Elc() = E0;
   }
   if constexpr (RequestedFields::Iv0_found()) {
      fields.Iv0() = 1.0;
   }
   LOWER_BITS(_status, STATE_INVALID);
};

/*!
\author Juan G Alonso Guzman
\author Vladimir Florinski
\date 10/14/2022
*/
template <typename HConfig>
template <typename Coordinates, typename Fields, typename RequestedFields>
void BackgroundUniform<HConfig>::EvaluateBackgroundDerivatives(Coordinates& coords, Fields& fields)
{
// Spatial derivatives are zero
   if constexpr (RequestedFields::DelFluv_found()) fields.DelFluv() = gm_zeros;
   if constexpr (RequestedFields::DelMag_found()) fields.DelMag() = gm_zeros;
   if constexpr (RequestedFields::DelElc_found()) fields.DelElc() = gm_zeros;

// Time derivatives are zero
   if constexpr (RequestedFields::DotFluv_found()) fields.DotFluv() = gv_zeros;
   if constexpr (RequestedFields::DotMag_found()) fields.DotMag() = gv_zeros;
   if constexpr (RequestedFields::DotElc_found()) fields.DotElc() = gv_zeros;
};

};
