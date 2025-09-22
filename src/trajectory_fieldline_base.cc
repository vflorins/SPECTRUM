/*!
\file trajectory_fieldline.cc
\brief Defines a class for trajectory following a field line
\author Vladimir Florinski
\author Juan G Alonso Guzman

This file is part of the SPECTRUM suite of scientific numerical simulation codes. SPECTRUM stands for Space Plasma and Energetic Charged particle TRansport on Unstructured Meshes. The code simulates plasma or neutral particle flows using MHD equations on a grid, transport of cosmic rays using stochastic or grid based methods. The "unstructured" part refers to the use of a geodesic mesh providing a uniform coverage of the surface of a sphere.
*/

#include "trajectory_fieldline_base.hh"

namespace Spectrum {

//----------------------------------------------------------------------------------------------------------------------------------------------------
// TrajectoryFieldlineBase methods
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\author Vladimir Florinski
\author Lucius Schoenbaum
\date 08/16/2025
*/
template <typename Trajectory, typename HConfig>
TrajectoryFieldlineBase<Trajectory, HConfig>::TrajectoryFieldlineBase(void)
                   : TrajectoryBase(traj_name, STATE_NONE)
{
};

/*!
\author Vladimir Florinski
\date 01/28/2022
\param[in] name_in   Readable name of the class
\param[in] status_in Initial status
*/
template <typename Trajectory, typename HConfig>
TrajectoryFieldlineBase<Trajectory, HConfig>::TrajectoryFieldlineBase(const std::string& name_in, uint16_t status_in)
      : TrajectoryBase(name_in, status_in)
{
};


/*!
\author Juan G Alonso Guzman
\author Vladimir Florinski
\author Lucius Schoenbaum
\date 08/16/2025
*/
template <typename Trajectory, typename HConfig>
void TrajectoryFieldlineBase<Trajectory, HConfig>::SetStart(void)
{
// Call the base version of this function.
   TrajectoryBase::SetStart();
};

/*!
\author Lucius Schoenbaum
\date 08/17/2025
This is a stub overriding the method in the lower base class.
 */
template <typename Trajectory, typename HConfig>
void TrajectoryFieldlineBase<Trajectory, HConfig>::Slopes(GeoVector& slope_pos_istage, GeoVector& slope_mom_istage)
{}


/*!
\author Vladimir Florinski
\author Lucius Schoenbaum
\date 08/16/2025
*/
template <typename Trajectory, typename HConfig>
void TrajectoryFieldlineBase<Trajectory, HConfig>::PhysicalStep(void)
{
   dt_physical = cfl_adv_tf * _dmax / fabs(_coords.Vel()[2]);
};

/*!
\author Vladimir Florinski
\author Lucius Schoenbaum
\date 08/16/2025
\return True if a step was taken

If the state at return contains the TRAJ_TERMINATE flag, the calling program must stop this trajectory. If the state at the end contains the TRAJ_DISCARD flag, the calling program must reject this trajectory (and possibly repeat the trial with a different random number).
*/
template <typename Trajectory, typename HConfig>
bool TrajectoryFieldlineBase<Trajectory, HConfig>::Advance(void)
{
   return RKAdvance();
};

};
