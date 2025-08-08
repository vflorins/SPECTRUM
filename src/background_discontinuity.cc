/*!
\file background_discontinuity.cc
\brief Implements a simple planar MHD discontinuity background
\author Juan G Alonso Guzman

This file is part of the SPECTRUM suite of scientific numerical simulation codes. SPECTRUM stands for Space Plasma and Energetic Charged particle TRansport on Unstructured Meshes. The code simulates plasma or neutral particle flows using MHD equations on a grid, transport of cosmic rays using stochastic or grid based methods. The "unstructured" part refers to the use of a geodesic mesh providing a uniform coverage of the surface of a sphere.
*/

#include "background_discontinuity.hh"

namespace Spectrum {

//----------------------------------------------------------------------------------------------------------------------------------------------------
// BackgroundDiscontinuity methods
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\author Juan G Alonso Guzman
\author Lucius Schoenbaum
\date 08/05/2025
*/
template <typename Fields>
BackgroundDiscontinuity<Fields>::BackgroundDiscontinuity(void)
                       : BackgroundBase(bg_name_discontinuity, 0, STATE_NONE)
{
};

/*!
\author Juan G Alonso Guzman
\date 10/20/2023
\param[in] other Object to initialize from

A copy constructor should first first call the Params' version to copy the data container and then check whether the other object has been set up. If yes, it should simply call the virtual method "SetupBackground()" with the argument of "true".
*/
template <typename Fields>
BackgroundDiscontinuity<Fields>::BackgroundDiscontinuity(const BackgroundDiscontinuity& other)
                       : BackgroundBase(other)
{
   RAISE_BITS(_status, STATE_NONE);
   if (BITS_RAISED(other._status, STATE_SETUP_COMPLETE)) SetupBackground(true);
};

/*!
\author Juan G Alonso Guzman
\date 10/20/2023
*/
template <typename Fields>
BackgroundDiscontinuity<Fields>::BackgroundDiscontinuity(const std::string& name_in, unsigned int specie_in, uint16_t status_in)
               : BackgroundBase(name_in, specie_in, status_in)
{
};

/*!
\author Juan G Alonso Guzman
\date 05/14/2025
\param [in] construct Whether called from a copy constructor or separately

This method's main role is to unpack the data container and set up the class data members and status bits marked as "persistent". The function should assume that the data container is available because the calling function will always ensure this.
*/
template <typename Fields>
void BackgroundDiscontinuity<Fields>::SetupBackground(bool construct)
{
// The parent version must be called explicitly if not constructing
   if (!construct) BackgroundBase::SetupBackground(false);

// Unpack parameters
   container.Read(n_discont);
   container.Read(v_discont);
   container.Read(u1);
   container.Read(B1);

// Normalize n_discont
   n_discont.Normalize();
};

/*!
\author Juan G Alonso Guzman
\date 05/14/2025
*/
template <typename Fields>
void BackgroundDiscontinuity<Fields>::EvaluateBackground(void)
{
// Upstream
   if ((_pos - r0) * n_discont - v_discont * _t > 0) {
      if constexpr (Fields::Vel_found()) _fields.Vel() = u0;
      if constexpr (Fields::Mag_found()) _fields.Mag() = B0;
      if constexpr (Fields::Iv0_found()) _fields.Iv0() = 1.0;
   }
// Downstream
   else {
      if constexpr (Fields::Vel_found()) _fields.Vel() = u1;
      if constexpr (Fields::Mag_found()) _fields.Mag() = B1;
      if constexpr (Fields::Iv0_found()) _fields.Iv0() = 2.0;
   };

   if constexpr (Fields::Elc_found()) _fields.Elc() = -(_fields.Vel() ^ _fields.Mag()) / c_code;

   LOWER_BITS(_status, STATE_INVALID);
};

/*!
\author Juan G Alonso Guzman
\date 10/14/2022
*/
template <typename Fields>
void BackgroundDiscontinuity<Fields>::EvaluateBackgroundDerivatives(void)
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
