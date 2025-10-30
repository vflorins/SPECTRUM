/*!
\file trajectory_guiding.hh
\brief Declares a class for trajectory based on relativistic guiding center equations
\author Vladimir Florinski
\author Juan G Alonso Guzman

This file is part of the SPECTRUM suite of scientific numerical simulation codes. SPECTRUM stands for Space Plasma and Energetic Charged particle TRansport on Unstructured Meshes. The code simulates plasma or neutral particle flows using MHD equations on a grid, transport of cosmic rays using stochastic or grid based methods. The "unstructured" part refers to the use of a geodesic mesh providing a uniform coverage of the surface of a sphere.
*/

#ifndef SPECTRUM_TRAJECTORY_GUIDING_HH
#define SPECTRUM_TRAJECTORY_GUIDING_HH

#include "trajectory_base.hh"

namespace Spectrum {

// Switch controlling how to calculate p_perp. "0" means computing it at the end of the step from magnetic moment conservation. "1" means advancing it in time according to the scheme (does not guarantee conservation of MM, but can be used with non-adiabatic terms).
#define PPERP_METHOD 1

//! Readable name of the TrajectoryGuiding class
const std::string traj_name_guiding = "TrajectoryGuiding";

//! Default initial size
const unsigned int defsize_guiding = 10000;

//! CFL condition for advection
const double cfl_adv_tg = 0.5;

//! Safety factor for drift-based time step (to modify "drift_vel" with a small fraction of the particle's velocity)
const double drift_safety_tg = 0.5;

//! How many time steps to allow before recording a mirror event
const int mirror_thresh_guiding = 10;

/*!
\brief Trajectory tracer for the relativistic guiding center equations
\author Vladimir Florinski

Components of "traj_mom" are: p_perp (x), unused (y), p_para (z)
*/
class TrajectoryGuiding : public TrajectoryBase {

protected:

//! Magnetic moment (transient)
   double mag_mom;

//! Modified electric field (transient)
   GeoVector Evec_star;

//! Modified magnetic field (transient)
   GeoVector Bvec_star;

//! Drift elocity (transient)
   GeoVector drift_vel;

//! Conversion from (p_perp,0,p_para) to (p,mu,0)
   GeoVector ConvertMomentum(void) const override;

//! Momentum transformation on reflection at a boundary
   void ReverseMomentum(void) override;

//! Compute the modified electromagnetic fields
   void ModifiedFields(void);

//! Compute the drift coefficient
   void DriftCoeff(void);

//! Compute the RK slopes
   void Slopes(GeoVector& slope_pos_istage, GeoVector& slope_mom_istage, double& slope_amp_istage, double& slope_wgt_istage) override;

//! Compute the physical time step
   void PhysicalStep(void) override;

//! Take a step
   bool Advance(void) override;

//! Adjust perp momentum to conserve magnetic moment
   void MomentumCorrection(void) override;

public:

//! Default constructor
   TrajectoryGuiding(void);

//! Constructor with arguments (to speed up construction of derived classes)
   TrajectoryGuiding(const std::string& name_in, unsigned int specie_in, uint16_t status_in, bool presize_in);

//! Copy constructor (class not copyable)
   TrajectoryGuiding(const TrajectoryGuiding& other) = delete;

//! Destructor
   ~TrajectoryGuiding() override = default;

//! Clone function
   CloneFunctionTrajectory(TrajectoryGuiding);

//! Clear the trajectory and start a new one with specified position and momentum
   void SetStart(void) override;
};

//----------------------------------------------------------------------------------------------------------------------------------------------------
// TrajectoryGuiding inline methods
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\author Vladimir Florinski
\date 04/11/2022
\return A vector in the (p,mu,0) format
*/
inline GeoVector TrajectoryGuiding::ConvertMomentum(void) const
{
   return GeoVector(sqrt(Sqr(_mom[0])+Sqr(_mom[2])), _mom[2] / _mom.Norm(), 0.0);
};


/*!
\author Juan G Alonso Guzman
\author Vladimir Florinski
\date 07/07/2023
*/
inline void TrajectoryGuiding::ReverseMomentum(void)
{
   _mom[2] = -_mom[2];
};

//! Trajectory type
#if TRAJ_TYPE == TRAJ_GUIDING
typedef TrajectoryGuiding TrajectoryType;
#endif

};

#endif
