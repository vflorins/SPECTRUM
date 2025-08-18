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
template <typename Trajectory, typename Fields>
TrajectoryFieldlineBase<Trajectory, Fields>::TrajectoryFieldlineBase(void)
                   : TrajectoryBase(traj_name_fieldline_base, 0, STATE_NONE, defsize_fieldline)
{
};

/*!
\author Vladimir Florinski
\date 01/28/2022
\param[in] name_in   Readable name of the class
\param[in] specie_in Particle's specie
\param[in] status_in Initial status
\param[in] presize_in Whether to pre-allocate memory for trajectory arrays
*/
template <typename Trajectory, typename Fields>
TrajectoryFieldlineBase<Trajectory, Fields>::TrajectoryFieldlineBase(const std::string& name_in, unsigned int specie_in, uint16_t status_in, bool presize_in)
      : TrajectoryBase(name_in, specie_in, status_in, presize_in)
{
};


/*!
\author Juan G Alonso Guzman
\author Vladimir Florinski
\author Lucius Schoenbaum
\date 08/16/2025
*/
template <typename Trajectory, typename Fields>
void TrajectoryFieldlineBase<Trajectory, Fields>::SetStart(void)
{
// Call the base version of this function.
   TrajectoryBase::SetStart();
};

/*!
\author Lucius Schoenbaum
\date 08/17/2025
This is a stub overriding the method in the lower base class.
 */
template <typename Trajectory, typename Fields>
void TrajectoryFieldlineBase<Trajectory, Fields>::Slopes(GeoVector& slope_pos_istage, GeoVector& slope_mom_istage)
{}


/*!
\author Vladimir Florinski
\author Lucius Schoenbaum
\date 08/16/2025
*/
template <typename Trajectory, typename Fields>
void TrajectoryFieldlineBase<Trajectory, Fields>::PhysicalStep(void)
{
   dt_physical = cfl_adv_tf * _dmax / fabs(_vel[2]);
};

/*!
\author Vladimir Florinski
\author Lucius Schoenbaum
\date 08/16/2025
\return True if a step was taken

If the state at return contains the TRAJ_TERMINATE flag, the calling program must stop this trajectory. If the state at the end contains the TRAJ_DISCARD flag, the calling program must reject this trajectory (and possibly repeat the trial with a different random number).
*/
template <typename Trajectory, typename Fields>
bool TrajectoryFieldlineBase<Trajectory, Fields>::Advance(void)
{
   return RKAdvance();
};

};
