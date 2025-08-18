/*!
\file trajectory_guiding.cc
\brief Defines a class for trajectory based on relativistic guiding center equations
\author Vladimir Florinski
\author Juan G Alonso Guzman

This file is part of the SPECTRUM suite of scientific numerical simulation codes. SPECTRUM stands for Space Plasma and Energetic Charged particle TRansport on Unstructured Meshes. The code simulates plasma or neutral particle flows using MHD equations on a grid, transport of cosmic rays using stochastic or grid based methods. The "unstructured" part refers to the use of a geodesic mesh providing a uniform coverage of the surface of a sphere.
*/

#include "trajectory_guiding.hh"
#include "common/print_warn.hh"

namespace Spectrum {

//----------------------------------------------------------------------------------------------------------------------------------------------------
// TrajectoryGuiding methods
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\author Vladimir Florinski
\date 12/17/2020
*/
template <typename Fields>
TrajectoryGuiding<Fields>::TrajectoryGuiding(void)
                 : TrajectoryGuidingBase(traj_name_guiding, 0, STATE_NONE, defsize_guiding)
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
template <typename Fields>
TrajectoryGuiding<Fields>::TrajectoryGuiding(const std::string& name_in, unsigned int specie_in, uint16_t status_in, bool presize_in)
                 : TrajectoryGuidingBase(name_in, specie_in, status_in, presize_in)
{
};

};
