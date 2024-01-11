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
// TrajectoryLorentz methods
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\author Vladimir Florinski
\date 09/30/2022
*/
TrajectoryFieldline::TrajectoryFieldline(void)
                   : TrajectoryBase(traj_name_fieldline, 0, STATE_NONE, defsize_fieldline)
{
};

/*!
\author Juan G Alonso Guzman
\author Vladimir Florinski
\date 10/13/2022
*/
void TrajectoryFieldline::SetStart(void)
{
// Call the base version of this function.
   TrajectoryBase::SetStart();

// Redefine which field line to follow based on header file parameter
   _spdata._mask = which_field_to_follow;
   spdata0._mask = which_field_to_follow;
};

/*!
\author Vladimir Florinski
\date 09/30/2022
*/
void TrajectoryFieldline::PhysicalStep(void)
{
   dt_physical = cfl_adv_tf * _spdata.dmax / fabs(_vel[2]);
};

/*!
\author Vladimir Florinski
\date 09/30/2022
\param[out] slope_pos_istage RK slope for position
\param[out] slope_mom_istage RK slope for momentum
*/
void TrajectoryFieldline::Slopes(GeoVector& slope_pos_istage, GeoVector& slope_mom_istage)
{
   if(BITS_RAISED(which_field_to_follow, BACKGROUND_U)) slope_pos_istage = _vel[2] * UnitVec(_spdata.Uvec);
   else if(BITS_RAISED(which_field_to_follow, BACKGROUND_B)) slope_pos_istage = _vel[2] * _spdata.bhat;
   else if(BITS_RAISED(which_field_to_follow, BACKGROUND_E)) slope_pos_istage = _vel[2] * UnitVec(_spdata.Evec);
   else slope_pos_istage = gv_zeros;
   slope_mom_istage = gv_zeros;
};

/*!
\author Vladimir Florinski
\date 09/30/2022
\return True if a step was taken

If the state at return contains the TRAJ_TERMINATE flag, the calling program must stop this trajectory. If the state at the end contains the TRAJ_DISCARD flag, the calling program must reject this trajectory (and possibly repeat the trial with a different random number).
*/
bool TrajectoryFieldline::Advance(void)
{
   return RKAdvance();
};

};
