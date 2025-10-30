/*!
\file trajectory_parker_source.cc
\brief Defines a class for trajectory based on the Parker Transport Equation
\author Juan G Alonso Guzman

This file is part of the SPECTRUM suite of scientific numerical simulation codes. SPECTRUM stands for Space Plasma and Energetic Charged particle TRansport on Unstructured Meshes. The code simulates plasma or neutral particle flows using MHD equations on a grid, transport of cosmic rays using stochastic or grid based methods. The "unstructured" part refers to the use of a geodesic mesh providing a uniform coverage of the surface of a sphere.
*/

#include "trajectory_parker_source.hh"
#include "common/print_warn.hh"

namespace Spectrum {

//----------------------------------------------------------------------------------------------------------------------------------------------------
// TrajectoryParkerSource methods
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\author Juan G Alonso Guzman
\date 07/14/2025
*/
TrajectoryParkerSource::TrajectoryParkerSource(void)
                      : TrajectoryParker(traj_name_parker_source, 0, STATE_NONE, defsize_parker)
{
};

/*!
\author Juan G Alonso Guzman
\date 07/14/2025
\param[in] name_in   Readable name of the class
\param[in] specie_in Particle's specie
\param[in] status_in Initial status
\param[in] presize_in Whether to pre-allocate memory for trajectory arrays
*/
TrajectoryParkerSource::TrajectoryParkerSource(const std::string& name_in, unsigned int specie_in, uint16_t status_in, bool presize_in)
                      : TrajectoryParker(name_in, specie_in, status_in, presize_in)
{
};

/*!
\author Juan G Alonso Guzman
\date 07/14/2025
*/
void TrajectoryParkerSource::SetStart(void)
{
// Call the parent version of this function.
   TrajectoryParker::SetStart();

// Initialize source
   source->GetSource(_t, _pos, _mom, _spdata, 1.0);
};

/*!
\author Juan G Alonso Guzman
\date 07/14/2025
\param[out] slope_pos_istage RK slope for position
\param[out] slope_mom_istage RK slope for momentum
*/
void TrajectoryParkerSource::Slopes(GeoVector& slope_pos_istage, GeoVector& slope_mom_istage, double& slope_amp_istage, double& slope_wgt_istage)
{
// Call the parent version of this function.
   TrajectoryParker::Slopes(slope_pos_istage, slope_mom_istage, slope_amp_istage, slope_wgt_istage);

// Get amplitude
   slope_amp_istage = source->GetSource(_t, _pos, _mom, _spdata, dt);
};

};
