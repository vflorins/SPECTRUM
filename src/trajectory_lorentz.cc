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
template <typename HConfig>
TrajectoryLorentz<HConfig>::TrajectoryLorentz(void)
      : TrajectoryBase(traj_name, STATE_NONE)
{
};

/*!
\author Vladimir Florinski
\date 06/14/2021
*/
template <typename HConfig>
void TrajectoryLorentz<HConfig>::SetStart(void)
{
// Call the base version of this function.
   TrajectoryBase::SetStart();

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
template <typename HConfig>
void TrajectoryLorentz<HConfig>::PhysicalStep(void)
{
   constexpr int steps_per_orbit = HConfig::steps_per_orbit;
   constexpr double cfl_adv = HConfig::cfl_advection;

// Obtain the time step based on the orbit resolution. If B is zero, use a small but finite value of "Omega".
   double Omega = fmax(CyclotronFrequency(_coords.Vel().Norm(), _fields.AbsMag(), specie), sp_tiny);
   dt_physical = M_2PI / Omega / steps_per_orbit;

// Obtain the grid based time step based on grid max distance.
   dt_physical = fmin(dt_physical, cfl_adv * _dmax / _coords.Vel().Norm());
};

/*!
\author Juan G Alonso Guzman
\date 04/18/2022
\param[out] slope_pos_istage RK slope for position
\param[out] slope_mom_istage RK slope for momentum
*/
template <typename HConfig>
void TrajectoryLorentz<HConfig>::Slopes(GeoVector& slope_pos_istage, GeoVector& slope_mom_istage)
{
   slope_pos_istage = _coords.Vel();
   slope_mom_istage = HConfig::q * (_fields.Elc() + (_coords.Vel() ^ _fields.Mag()) / c_code);
};

/*!
\author Vladimir Florinski
\date 01/06/2021
\return True if a step was taken

If the state at return contains the TRAJ_TERMINATE flag, the calling program must stop this trajectory. If the state at the end contains the TRAJ_DISCARD flag, the calling program must reject this trajectory (and possibly repeat the trial with a different random number).
*/
template <typename HConfig>
bool TrajectoryLorentz<HConfig>::Advance(void)
{
   return RKAdvance();
};

};
