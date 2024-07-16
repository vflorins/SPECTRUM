/*!
\file trajectory_lorentz.cc
\brief Defines a class for trajectory based on Newton-Lorentz equation
\author Vladimir Florinski

This file is part of the SPECTRUM suite of scientific numerical simulation codes. SPECTRUM stands for Space Plasma and Energetic Charged particle TRansport on Unstructured Meshes. The code simulates plasma or neutral particle flows using MHD equations on a grid, transport of cosmic rays using stochastic or grid based methods. The "unstructured" part refers to the use of a geodesic mesh providing a uniform coverage of the surface of a sphere.
*/

#include "trajectory_lorentz.hh"

namespace Spectrum {

//----------------------------------------------------------------------------------------------------------------------------------------------------
// TrajectoryLorentz methods
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\author Vladimir Florinski
\date 12/17/2020
*/
TrajectoryLorentz::TrajectoryLorentz(void)
                 : TrajectoryBase(traj_name_lorentz, 0, STATE_NONE, defsize_lorentz)
{
};

/*!
\author Vladimir Florinski
\date 06/14/2021
*/
void TrajectoryLorentz::SetStart(void)
{
// Call the base version of this function.
   TrajectoryBase::SetStart();

// Lorentz requires E and B fields only
   _spdata._mask = BACKGROUND_E | BACKGROUND_B;
   spdata0._mask = BACKGROUND_E | BACKGROUND_B;

// "p_para" is used to compare p*B to check for mirroring. The check is for <0, so an initial zero will return a false.
   p_para = 0.0;
   t_mirror = 0.0;

// No need to reset "which_fields" because only "basic" fields are needed in Lorentz
};

/*!
\author Vladimir Florinski
\author Juan G Alonso Guzman
\date 12/15/2020
*/
void TrajectoryLorentz::PhysicalStep(void)
{
// Obtain the time step based on the orbit resolution. If B is zero, use a small but finite value of "Omega".
   double Omega = fmax(CyclotronFrequency(_vel.Norm(), _spdata.Bmag, specie), sp_tiny);
   dt_physical = M_2PI / Omega / steps_per_orbit;

// Obtain the grid based time step based on grid max distance.
   dt_physical = fmin(dt_physical, cfl_adv_tl * _spdata.dmax / _vel.Norm());
};

/*!
\author Juan G Alonso Guzman
\date 04/18/2022
\param[out] slope_pos_istage RK slope for position
\param[out] slope_mom_istage RK slope for momentum
*/
void TrajectoryLorentz::Slopes(GeoVector& slope_pos_istage, GeoVector& slope_mom_istage)
{
   slope_pos_istage = _vel;
   slope_mom_istage = q * (_spdata.Evec + (_vel ^ _spdata.Bvec) / c_code);
};

/*!
\author Vladimir Florinski
\date 01/06/2021
\return True if a step was taken

If the state at return contains the TRAJ_TERMINATE flag, the calling program must stop this trajectory. If the state at the end contains the TRAJ_DISCARD flag, the calling program must reject this trajectory (and possibly repeat the trial with a different random number).
*/
bool TrajectoryLorentz::Advance(void)
{
   return RKAdvance();
};

};
