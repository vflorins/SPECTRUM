/*!
\file trajectory_fieldline.cc
\brief Defines a class for trajectory following a field line
\author Vladimir Florinski
\author Juan G Alonso Guzman

This file is part of the SPECTRUM suite of scientific numerical simulation codes. SPECTRUM stands for Space Plasma and Energetic Charged particle TRansport on Unstructured Meshes. The code simulates plasma or neutral particle flows using MHD equations on a grid, transport of cosmic rays using stochastic or grid based methods. The "unstructured" part refers to the use of a geodesic mesh providing a uniform coverage of the surface of a sphere.
*/

#include "trajectory_fieldline.hh"

namespace Spectrum {

//----------------------------------------------------------------------------------------------------------------------------------------------------
// TrajectoryFieldline methods
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\author Vladimir Florinski
\author Lucius Schoenbaum
\date 08/16/2025
*/
template <typename HConfig, typename Field_t>
TrajectoryFieldline<HConfig, Field_t>::TrajectoryFieldline(void)
                   : TrajectoryFieldlineBase(traj_name, STATE_NONE)
{
};

/*!
\author Vladimir Florinski
\author Lucius Schoenbaum
\date 08/16/2025
\param[out] slope_pos_istage RK slope for position
\param[out] slope_mom_istage RK slope for momentum
*/
template <typename HConfig, typename Field_t>
void TrajectoryFieldline<HConfig, Field_t>::Slopes(GeoVector& slope_pos_istage, GeoVector& slope_mom_istage)
{
   if constexpr (std::same_as<Field_t, Vel_t>) {
      slope_pos_istage = _coords.Vel()[2] * UnitVec(_fields.Vel());
   }
   else if constexpr (std::same_as<Field_t, Mag_t>) {
      if constexpr(TrajectoryFields::HatMag_found())
         slope_pos_istage = _coords.Vel()[2] * _fields.HatMag();
      else
         slope_pos_istage = _coords.Vel()[2] * UnitVec(_fields.Mag());
   }
   else if constexpr (std::same_as<Field_t, Elc_t>) {
      slope_pos_istage = _coords.Vel()[2] * UnitVec(_fields.Elc());
   }
   else {
      slope_pos_istage = gv_zeros;
   }
   slope_mom_istage = gv_zeros;
};

};
