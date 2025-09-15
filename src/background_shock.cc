/*!
\file background_shock.cc
\brief Implements a simple planar shock field background
\author Juan G Alonso Guzman

This file is part of the SPECTRUM suite of scientific numerical simulation codes. SPECTRUM stands for Space Plasma and Energetic Charged particle TRansport on Unstructured Meshes. The code simulates plasma or neutral particle flows using MHD equations on a grid, transport of cosmic rays using stochastic or grid based methods. The "unstructured" part refers to the use of a geodesic mesh providing a uniform coverage of the surface of a sphere.
*/

#include "background_shock.hh"
#include "common/print_warn.hh"

namespace Spectrum {

//----------------------------------------------------------------------------------------------------------------------------------------------------
// BackgroundShock methods
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\author Juan G Alonso Guzman
\date 10/20/2023
*/
template <typename HConfig>
BackgroundShock<HConfig>::BackgroundShock(void)
               : BackgroundBase(bg_name, STATE_NONE)
{
};

/*!
\author Juan G Alonso Guzman
\date 10/20/2023
\param[in] other Object to initialize from

A copy constructor should first first call the Params' version to copy the data container and then check whether the other object has been set up. If yes, it should simply call the virtual method "SetupBackground()" with the argument of "true".
*/
template <typename HConfig>
BackgroundShock<HConfig>::BackgroundShock(const BackgroundShock& other)
               : BackgroundBase(other)
{
   RAISE_BITS(_status, STATE_NONE);
   if (BITS_RAISED(other._status, STATE_SETUP_COMPLETE)) SetupBackground(true);
};

/*!
\author Juan G Alonso Guzman
\date 10/20/2023
*/
template <typename HConfig>
BackgroundShock<HConfig>::BackgroundShock(const std::string& name_in, uint16_t status_in)
               : BackgroundBase(name_in, status_in)
{
};

/*!
\author Vladimir Florinski
\author Juan G Alonso Guzman
\date 05/19/2025
\param [in] construct Whether called from a copy constructor or separately

This method's main role is to unpack the data container and set up the class data members and status bits marked as "persistent". The function should assume that the data container is available because the calling function will always ensure this.
*/
template <typename HConfig>
void BackgroundShock<HConfig>::SetupBackground(bool construct)
{
// The parent version must be called explicitly if not constructing
   if (!construct) BackgroundBase::SetupBackground(false);

// Unpack parameters
   container.Read(n_shock);
   container.Read(v_shock);
   container.Read(compression);

// Normalize n_shock
   n_shock.Normalize();

// Normal and tangential components in the shock frame
   GeoVector u0n, u0t, B0n, B0t;
   u0n = (u0 * n_shock - v_shock) * n_shock;
   u0t = u0 - (u0 * n_shock) * n_shock;
   B0n = (B0 * n_shock) * n_shock;
   B0t = B0 - B0n;

   if(u0n * n_shock > 0.0) PrintError(__FILE__, __LINE__, "Upstream normal velocity is in the wrong direction", true);

   u1 = u0n / compression + v_shock * n_shock + u0t;
   B1 = B0n + B0t * compression;
};

/*!
\author Juan G Alonso Guzman
\date 05/14/2025
*/
template <typename HConfig>
template <typename Fields>
void BackgroundShock<HConfig>::EvaluateBackground(Coordinates& coords, Fields& fields)
{
// Upstream
   if ((coords.Pos() - r0) * n_shock - v_shock * coords.Time() > 0) {
      if constexpr (Fields::Vel_found()) fields.Vel() = u0;
      if constexpr (Fields::Mag_found()) fields.Mag() = B0;
      if constexpr (Fields::Iv0_found()) fields.Iv0() = 1.0;
   }
// Downstream
   else {
      if constexpr (Fields::Vel_found()) fields.Vel() = u1;
      if constexpr (Fields::Mag_found()) fields.Mag() = B1;
      if constexpr (Fields::Iv0_found()) fields.Iv0() = 2.0;
   };

   if constexpr (Fields::Elc_found()) fields.Elc() = -(fields.Vel() ^ fields.Mag()) / c_code;

   LOWER_BITS(_status, STATE_INVALID);
};

/*!
\author Juan G Alonso Guzman
\date 10/14/2022
*/
template <typename HConfig>
template <typename Fields>
void BackgroundShock<HConfig>::EvaluateBackgroundDerivatives(Coordinates& coords, Specie& specie, Fields& fields)
{
// Spatial derivatives are zero
   if constexpr (Fields::DelVel_found()) fields.DelVel() = gm_zeros;
   if constexpr (Fields::DelMag_found()) fields.DelMag() = gm_zeros;
   if constexpr (Fields::DelElc_found()) fields.DelElc() = gm_zeros;

// Time derivatives are zero
   if constexpr (Fields::DotVel_found()) fields.DotVel() = gv_zeros;
   if constexpr (Fields::DotMag_found()) fields.DotMag() = gv_zeros;
   if constexpr (Fields::DotElc_found()) fields.DotElc() = gv_zeros;
};

};
