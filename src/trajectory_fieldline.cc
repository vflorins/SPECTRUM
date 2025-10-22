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
// TrajectoryFieldlineBase methods
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\author Vladimir Florinski
\author Lucius Schoenbaum
\date 08/16/2025
*/
template <typename Background, typename Diffusion, typename Field_t>
TrajectoryFieldline<HConfig, Background, Field_t>::TrajectoryFieldline(void)
                   : TrajectoryBase(traj_name, STATE_NONE)
{
};

/*!
\author Vladimir Florinski
\date 01/28/2022
\param[in] name_in   Readable name of the class
\param[in] status_in Initial status
*/
template <typename Background, typename Diffusion, typename Field_t>
TrajectoryFieldline<HConfig, Background, Field_t>::TrajectoryFieldline(const std::string& name_in, uint16_t status_in)
      : TrajectoryBase(name_in, status_in)
{
};


/*!
\author Juan G Alonso Guzman
\author Vladimir Florinski
\author Lucius Schoenbaum
\date 08/16/2025
*/
template <typename Background, typename Diffusion, typename Field_t>
void TrajectoryFieldline<HConfig, Background, Field_t>::SetStart(void)
{
// Call the base version of this function.
   TrajectoryBase::SetStart();
};


/*!
\author Vladimir Florinski
\author Lucius Schoenbaum
\date 08/16/2025
*/
template <typename Background, typename Diffusion, typename Field_t>
void TrajectoryFieldline<HConfig, Background, Field_t>::PhysicalStep(void)
{
   constexpr double cfl_adv = HConfig::cfl_advection;
   // todo review - what is Vel_sys ----> check
   dt_physical = cfl_adv * _dmax / fabs(_coords.VelPara());
};


/*!
\author Vladimir Florinski
\author Lucius Schoenbaum
\date 08/16/2025
\param[out] slope_pos_istage RK slope for position
\param[out] slope_mom_istage RK slope for momentum
*/
template <typename Background, typename Diffusion, typename Field_t>
void TrajectoryFieldline<HConfig, Background, Field_t>::Slopes(GeoVector& slope_pos_istage, GeoVector& slope_mom_istage)
{
   // todo VelPara
   if constexpr (std::same_as<Field_t, Fluv_t>) {
      slope_pos_istage = _coords.VelPara() * _fields.HatFluv();
   }
   else if constexpr (std::same_as<Field_t, Mag_t>) {
      slope_pos_istage = _coords.VelPara() * _fields.HatMag();
   }
   else if constexpr (std::same_as<Field_t, Elc_t>) {
      slope_pos_istage = _coords.VelPara() * _fields.HatElc();
   }
   else {
      slope_pos_istage = gv_zeros;
   }
   slope_mom_istage = gv_zeros;
};

/*!
\author Vladimir Florinski
\author Lucius Schoenbaum
\date 08/16/2025
\return True if a step was taken

If the state at return contains the TRAJ_TERMINATE flag, the calling program must stop this trajectory. If the state at the end contains the TRAJ_DISCARD flag, the calling program must reject this trajectory (and possibly repeat the trial with a different random number).
*/
template <typename Background, typename Diffusion, typename Field_t>
bool TrajectoryFieldline<HConfig, Background, Field_t>::Advance(void)
{
   return RKAdvance();
};

};
