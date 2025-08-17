/*!
\file trajectory_fieldline.hh
\brief Declares a class for trajectory following a field line
\author Vladimir Florinski
\author Juan G Alonso Guzman

This file is part of the SPECTRUM suite of scientific numerical simulation codes. SPECTRUM stands for Space Plasma and Energetic Charged particle TRansport on Unstructured Meshes. The code simulates plasma or neutral particle flows using MHD equations on a grid, transport of cosmic rays using stochastic or grid based methods. The "unstructured" part refers to the use of a geodesic mesh providing a uniform coverage of the surface of a sphere.
*/

#ifndef SPECTRUM_TRAJECTORY_FIELDLINE_HH
#define SPECTRUM_TRAJECTORY_FIELDLINE_HH

#include "trajectory_base.hh"

namespace Spectrum {

//! Readable name of the TrajectoryFieldline class
const std::string traj_name_fieldline = "TrajectoryFieldline";

//! Default initial size
const unsigned int defsize_fieldline = 100000;

//! CFL condition for advection
const double cfl_adv_tf = 0.5;

//! Which field to trace
const uint16_t which_field_to_follow = BACKGROUND_B;

/*!
\brief Field line tracer
\author Juan G Alonso Guzman
\author Vladimir Florinski

Components of "traj_mom" are: unused (x), unused (y), p_para (z)
*/
class TrajectoryFieldline : public TrajectoryBase {

protected:

//! Conversion from (p_x,p_y,p_z) to (p,mu,phi)
   GeoVector ConvertMomentum(void) const override;

//! Compute the RK slopes
   void Slopes(GeoVector& slope_pos_istage, GeoVector& slope_mom_istage) override;

//! Compute the physical time step
   void PhysicalStep(void) override;

//! Take a step
   bool Advance(void) override;

public:

//! Default constructor
   TrajectoryFieldline(void);

//! Copy constructor (class not copyable)
   TrajectoryFieldline(const TrajectoryFieldline& other) = delete;

//! Destructor
   ~TrajectoryFieldline() override = default;

//! Clone function
   CloneFunctionTrajectory(TrajectoryFieldline);

//! Clear the trajectory and start a new one with specified position and momentum
   void SetStart(void) override;
};

//----------------------------------------------------------------------------------------------------------------------------------------------------
// TrajectoryFieldline inline methods
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\author Vladimir Florinski
\date 09/30/2022
\return A vector in the (p,mu,phi) format
\note Not used, but needs to be "overriden" from virtual definition in TrajectoryBase
*/
inline GeoVector TrajectoryFieldline::ConvertMomentum(void) const
{
   return GeoVector(fabs(_mom[2]), 0.0, 0.0);
};

//! Trajectory type
#if TRAJ_TYPE == TRAJ_FIELDLINE
typedef TrajectoryFieldline TrajectoryType;
#endif

};

#endif
