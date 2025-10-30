/*!
\file trajectory_lorentz.hh
\brief Declares a class for trajectory based on Newton-Lorentz equation
\author Vladimir Florinski

This file is part of the SPECTRUM suite of scientific numerical simulation codes. SPECTRUM stands for Space Plasma and Energetic Charged particle TRansport on Unstructured Meshes. The code simulates plasma or neutral particle flows using MHD equations on a grid, transport of cosmic rays using stochastic or grid based methods. The "unstructured" part refers to the use of a geodesic mesh providing a uniform coverage of the surface of a sphere.
*/

#ifndef SPECTRUM_TRAJECTORY_LORENTZ_HH
#define SPECTRUM_TRAJECTORY_LORENTZ_HH

#include "trajectory_base.hh"

namespace Spectrum {

//! Readable name of the TrajectoryLorentz class
const std::string traj_name_lorentz = "TrajectoryLorentz";

//! Default initial size
const unsigned int defsize_lorentz = 100000;

//! CFL condition for advection
const double cfl_adv_tl = 0.1;

//! Number of time steps per one orbit
const unsigned int steps_per_orbit = 100;

//! How many time steps to allow before recording a mirror event
const unsigned int mirror_thresh_lorentz = 300;

/*!
\brief Trajectory tracer for the Newton-Lorentz equation (full orbit)
\author Vladimir Florinski
*/
class TrajectoryLorentz : public TrajectoryBase {

protected:

//! Previous value of parallel momentum (transient)
   double p_para;

//! Time of the last mirror point (transient)
   double t_mirror;

//! Conversion from (p_x,p_y,p_z) to (p,mu,phi)
   GeoVector ConvertMomentum(void) const override;

//! Compute the RK slopes
   void Slopes(GeoVector& slope_pos_istage, GeoVector& slope_mom_istage, double& slope_amp_istage, double& slope_wgt_istage) override;

//! Compute the physical time step
   void PhysicalStep(void) override;

//! Take a step
   bool Advance(void) override;

public:

//! Default constructor
   TrajectoryLorentz(void);

//! Copy constructor (class not copyable)
   TrajectoryLorentz(const TrajectoryLorentz& other) = delete;

//! Destructor
   ~TrajectoryLorentz() override = default;

//! Clone function
   CloneFunctionTrajectory(TrajectoryLorentz);

//! Clear the trajectory and start a new one with specified position and momentum
   void SetStart(void) override;
};

//----------------------------------------------------------------------------------------------------------------------------------------------------
// TrajectoryLorentz inline methods
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\author Vladimir Florinski
\date 09/02/2022
\return A vector in the (p,mu,phi) format
*/
inline GeoVector TrajectoryLorentz::ConvertMomentum(void) const
{
   return GeoVector(_mom.Norm(), (_mom * _spdata.bhat) / _mom.Norm(), 0.0);
};

//! Trajectory type
#if TRAJ_TYPE == TRAJ_LORENTZ
typedef TrajectoryLorentz TrajectoryType;
#endif

};

#endif
